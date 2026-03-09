# ============================================================
# Outer hold-out (70/30) + inner CV tuning (baseline fixed)
# WITH Gompertz–Makeham reference curve:
#   - Fit GM on ages 42–70 (fit_age_range)
#   - Extrapolate to 30–100 (ref_age_range) on a grid (ref_step)
# DO NOT leak: imputer is fit on fit_idx only (within each fold)
# Early stopping uses 10% split from inner-train ONLY
# ============================================================

suppressPackageStartupMessages({
  library(xgboost)
  library(survival)
})

# ----------------------------
# AFT DMatrix helper
# ----------------------------
make_aft_dmatrix <- function(Xm, t, s, base_margin = NULL) {
  d <- xgb.DMatrix(Xm)
  setinfo(d, "label_lower_bound", t)
  setinfo(d, "label_upper_bound", ifelse(s == 1L, t, Inf))
  if (!is.null(base_margin)) setinfo(d, "base_margin", base_margin)
  d
}

# ----------------------------
# Stratified split (by status) for hold-out / early-stop
# p_train + p_valid <= 1; remainder goes to test
# ----------------------------
split_stratified <- function(status, p_train = 0.7, p_valid = 0.1, seed = 2026) {
  set.seed(seed)
  status <- as.integer(status)
  stopifnot(all(status %in% c(0L, 1L)))
  stopifnot(p_train >= 0, p_valid >= 0, p_train + p_valid <= 1)
  
  idx1 <- which(status == 1L)
  idx0 <- which(status == 0L)
  
  split_one <- function(idx) {
    idx <- sample(idx)
    n <- length(idx)
    n_tr <- floor(p_train * n)
    n_va <- floor(p_valid * n)
    list(
      train = idx[1:n_tr],
      valid = idx[(n_tr + 1):(n_tr + n_va)],
      test  = idx[(n_tr + n_va + 1):n]
    )
  }
  
  s1 <- split_one(idx1)
  s0 <- split_one(idx0)
  
  list(
    train = sample(c(s1$train, s0$train)),
    valid = sample(c(s1$valid, s0$valid)),
    test  = sample(c(s1$test,  s0$test))
  )
}

# ----------------------------
# Stratified K-fold split (by status)
# ----------------------------
make_folds_stratified <- function(status, K = 5, seed = 2026) {
  set.seed(seed)
  status <- as.integer(status)
  stopifnot(all(status %in% c(0L, 1L)))
  
  idx1 <- sample(which(status == 1L))
  idx0 <- sample(which(status == 0L))
  
  split_into_K <- function(idx) {
    cut_id <- rep(1:K, length.out = length(idx))
    split(idx, cut_id)
  }
  
  f1 <- split_into_K(idx1)
  f0 <- split_into_K(idx0)
  
  folds <- vector("list", K)
  for (k in 1:K) folds[[k]] <- sort(c(f1[[k]], f0[[k]]))
  folds
}

# ----------------------------
# Median imputer (fit on fit_idx only)
# ----------------------------
fit_imputer_median <- function(X_train) {
  meds <- apply(X_train, 2, function(x) median(x, na.rm = TRUE))
  meds[!is.finite(meds)] <- 0
  meds
}

apply_imputer_median <- function(X, meds) {
  stopifnot(ncol(X) == length(meds))
  for (j in seq_len(ncol(X))) {
    na_idx <- is.na(X[, j])
    if (any(na_idx)) X[na_idx, j] <- meds[j]
  }
  X
}

# ----------------------------
# C-index from risk score
# ----------------------------
cindex_from_risk <- function(time, status, risk) {
  fit <- coxph(Surv(time, status) ~ risk)
  as.numeric(summary(fit)$concordance[1])
}

# ----------------------------
# Enforce monotonic reference curve and invert risk->age (BA)
# ----------------------------
risk_to_BA <- function(risk, age_grid, risk_grid) {
  ok <- is.finite(age_grid) & is.finite(risk_grid)
  age_grid <- age_grid[ok]
  risk_grid <- risk_grid[ok]
  
  ord <- order(age_grid)
  age_grid <- age_grid[ord]
  risk_grid <- risk_grid[ord]
  
  # Monotone non-decreasing risk w.r.t age
  iso <- isoreg(age_grid, risk_grid)
  rg <- iso$yf
  ag <- iso$x
  
  # Invert via interpolation: BA = age such that risk_grid(age) == risk
  approx(x = rg, y = ag, xout = risk, rule = 2, ties = "ordered")$y
}

# ----------------------------
# Fit Gompertz–Makeham: y = A + B * exp(C * age)
# Constrain B>=0, C>=0 to preserve monotonicity
# ----------------------------
fit_gompertz_makeham <- function(age, y) {
  if (!requireNamespace("minpack.lm", quietly = TRUE)) {
    stop("Package 'minpack.lm' is required. Install via install.packages('minpack.lm')")
  }
  
  A0 <- as.numeric(quantile(y, 0.05, na.rm = TRUE)) - 0.01
  y_shift <- pmax(y - A0, 1e-6)
  
  lm0 <- lm(log(y_shift) ~ age)
  C0 <- max(0, as.numeric(coef(lm0)[2]))
  B0 <- exp(as.numeric(coef(lm0)[1]))
  
  fit <- tryCatch(
    minpack.lm::nlsLM(
      y ~ A + B * exp(C * age),
      start = list(A = A0, B = B0, C = C0),
      lower = c(A = -Inf, B = 0, C = 0),
      control = minpack.lm::nls.lm.control(maxiter = 500)
    ),
    error = function(e) NULL
  )
  
  fit
}

# ============================================================
# ONE inner fold training:
# - tr_idx split into fit_idx + es_idx (10% for early stopping)
# - fold test te_idx is ONLY for evaluation
# - GM curve fitted on baseline risk for ages 42-70; extrapolate to 30-100
# ============================================================
train_one_inner_fold_with_GM <- function(dat,
                                         feature_cols,
                                         params_base_fixed,
                                         params_full,
                                         tr_idx,
                                         te_idx,
                                         es_frac = 0.1,
                                         seed = 2026,
                                         censor_check = TRUE,
                                         fit_age_range = c(42, 70),
                                         ref_age_range = c(30, 100),
                                         ref_step = 0.1,
                                         nrounds = 5000,
                                         early_stopping_rounds = 100,
                                         verbose = 0) {
  
  set.seed(seed)
  
  time   <- as.numeric(dat$time_years)
  status <- as.integer(dat$status)
  age0   <- as.numeric(dat$age0)
  
  if (censor_check) {
    stopifnot(all(is.finite(time)), all(time > 0), all(status %in% c(0L, 1L)))
  }
  
  # Matrices
  X_full <- data.matrix(dat[, feature_cols, drop = FALSE])
  X_base <- matrix(age0, ncol = 1)
  colnames(X_base) <- "age0"
  
  # Split inner-train into fit + early-stop valid (stratified)
  sp_es <- split_stratified(status[tr_idx], p_train = 1 - es_frac, p_valid = es_frac, seed = seed + 11)
  fit_idx <- tr_idx[sp_es$train]
  es_idx  <- tr_idx[sp_es$valid]
  
  # Fit imputer on fit only (no leakage)
  meds <- fit_imputer_median(X_full[fit_idx, , drop = FALSE])
  X_full <- apply_imputer_median(X_full, meds)
  
  # ---- 1) Baseline AFT (age only) with early stopping on es_idx ----
  d_fit_base <- make_aft_dmatrix(X_base[fit_idx, , drop = FALSE], time[fit_idx], status[fit_idx])
  d_es_base  <- make_aft_dmatrix(X_base[es_idx,  , drop = FALSE], time[es_idx],  status[es_idx])
  
  bst_base <- xgb.train(
    params = params_base_fixed,
    data = d_fit_base,
    nrounds = nrounds,
    watchlist = list(train = d_fit_base, es = d_es_base),
    early_stopping_rounds = early_stopping_rounds,
    verbose = verbose
  )
  
  base_fit <- predict(bst_base, X_base[fit_idx, , drop = FALSE], outputmargin = TRUE)
  base_es  <- predict(bst_base, X_base[es_idx,  , drop = FALSE], outputmargin = TRUE)
  base_te  <- predict(bst_base, X_base[te_idx,  , drop = FALSE], outputmargin = TRUE)
  
  # ---- 2) Full AFT (biomarkers) with baseline offset; early stop on es_idx ----
  d_fit_full <- make_aft_dmatrix(X_full[fit_idx, , drop = FALSE], time[fit_idx], status[fit_idx], base_margin = base_fit)
  d_es_full  <- make_aft_dmatrix(X_full[es_idx,  , drop = FALSE], time[es_idx],  status[es_idx],  base_margin = base_es)
  d_te_full  <- make_aft_dmatrix(X_full[te_idx,  , drop = FALSE], time[te_idx],  status[te_idx],  base_margin = base_te)
  
  bst_full <- xgb.train(
    params = params_full,
    data = d_fit_full,
    nrounds = nrounds,
    watchlist = list(train = d_fit_full, es = d_es_full),
    early_stopping_rounds = early_stopping_rounds,
    verbose = verbose
  )
  
  # Fold-test predictions (evaluation only)
  margin_te <- predict(bst_full, d_te_full, outputmargin = TRUE)
  risk_te   <- -margin_te
  
  # ---- 3) GM reference curve (baseline) fitted on fit_age_range, extrapolate to ref_age_range ----
  age_grid_ref <- seq(ref_age_range[1], ref_age_range[2], by = ref_step)
  X_age_ref <- matrix(age_grid_ref, ncol = 1)
  colnames(X_age_ref) <- "age0"
  
  base_margin_ref <- predict(bst_base, X_age_ref, outputmargin = TRUE)
  risk_base_ref <- -base_margin_ref
  
  idx_fit_age <- age_grid_ref >= fit_age_range[1] & age_grid_ref <= fit_age_range[2]
  age_fit_gm <- age_grid_ref[idx_fit_age]
  y_fit_gm   <- as.numeric(risk_base_ref[idx_fit_age])
  
  gm_model <- fit_gompertz_makeham(age_fit_gm, y_fit_gm)
  
  if (is.null(gm_model)) {
    # Fallback: linear fit on the fitting age window (reference only)
    lm_ref <- lm(y_fit_gm ~ age_fit_gm)
    risk_grid_ref <- as.numeric(predict(lm_ref, newdata = data.frame(age_fit_gm = age_grid_ref)))
    gm_params <- c(A = NA, B = NA, C = NA)
  } else {
    cf <- coef(gm_model)
    gm_params <- cf
    risk_grid_ref <- as.numeric(cf["A"] + cf["B"] * exp(cf["C"] * age_grid_ref))
  }
  
  # Map risk -> BA via monotone inverse lookup
  BA_te <- risk_to_BA(risk_te, age_grid_ref, risk_grid_ref)
  
  # Evaluation metric for tuning: C-index on fold-test (risk only)
  cidx_te <- cindex_from_risk(time[te_idx], status[te_idx], risk_te)
  
  list(
    cindex = cidx_te,
    best_iter_full = bst_full$best_iteration,
    best_iter_base = bst_base$best_iteration,
    oof = list(margin = margin_te, risk = risk_te, BA = BA_te),
    ref = list(age_grid = age_grid_ref, risk_grid = risk_grid_ref, gm_params = gm_params)
  )
}

# ============================================================
# Inner CV scoring for one params_full (on OUTER_TRAIN only)
# Folds are fixed for fairness across trials
# ============================================================
inner_cv_score_with_GM <- function(dat_outer_train,
                                   feature_cols,
                                   params_base_fixed,
                                   params_full,
                                   folds,
                                   es_frac = 0.1,
                                   seed = 2026,
                                   fit_age_range = c(42, 70),
                                   ref_age_range = c(30, 100),
                                   ref_step = 0.1,
                                   nrounds = 5000,
                                   early_stopping_rounds = 100,
                                   verbose = 0) {
  
  K <- length(folds)
  n <- nrow(dat_outer_train)
  all_idx <- seq_len(n)
  
  fold_c <- rep(NA_real_, K)
  fold_iter_full <- rep(NA_real_, K)
  
  # Optional: store OOF BA/risk on outer_train for diagnostics (not required for tuning)
  oof_risk <- rep(NA_real_, n)
  oof_BA   <- rep(NA_real_, n)
  
  for (k in 1:K) {
    te_idx <- folds[[k]]
    tr_idx <- setdiff(all_idx, te_idx)
    
    res <- train_one_inner_fold_with_GM(
      dat = dat_outer_train,
      feature_cols = feature_cols,
      params_base_fixed = params_base_fixed,
      params_full = params_full,
      tr_idx = tr_idx,
      te_idx = te_idx,
      es_frac = es_frac,
      seed = seed + 1000 + k,
      censor_check = TRUE,
      fit_age_range = fit_age_range,
      ref_age_range = ref_age_range,
      ref_step = ref_step,
      nrounds = nrounds,
      early_stopping_rounds = early_stopping_rounds,
      verbose = verbose
    )
    
    fold_c[k] <- res$cindex
    fold_iter_full[k] <- res$best_iter_full
    
    oof_risk[te_idx] <- res$oof$risk
    oof_BA[te_idx]   <- res$oof$BA
  }
  
  list(
    fold_cindex = fold_c,
    mean_cindex = mean(fold_c, na.rm = TRUE),
    sd_cindex   = sd(fold_c, na.rm = TRUE),
    mean_best_iter_full = mean(fold_iter_full, na.rm = TRUE),
    oof = data.frame(oof_risk = oof_risk, oof_BA = oof_BA)
  )
}

# ============================================================
# Random search space for FULL model only
# ============================================================
sample_params_full <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  eta <- sample(c(0.01, 0.02, 0.03, 0.05), 1)
  depth <- sample(2:6, 1)
  
  # 条件约束：深树必须更保守
  if (depth >= 5 && eta > 0.03) eta <- sample(c(0.01, 0.02, 0.03), 1)
  
  list(
    booster = "gbtree",
    objective = "survival:aft",
    eval_metric = "aft-nloglik",
    
    eta = eta,
    max_depth = depth,
    min_child_weight = sample(c(5, 10, 20, 30, 50), 1),
    
    subsample = sample(c(0.7, 0.8, 0.9, 1.0), 1),
    colsample_bytree = sample(c(0.6, 0.7, 0.8, 0.9, 1.0), 1),
    
    lambda = exp(runif(1, log(1e-2), log(10))),
    alpha  = exp(runif(1, log(1e-3), log(1))),
    gamma  = exp(runif(1, log(1e-6), log(2))),
    
    aft_loss_distribution = "normal",
    aft_loss_distribution_scale = 1
  )
}

# ============================================================
# Tuning on OUTER_TRAIN only (baseline fixed, full tuned)
# ============================================================
tune_full_only_on_outer_train_with_GM <- function(dat_outer_train,
                                                  feature_cols,
                                                  params_base_fixed,
                                                  K = 5,
                                                  n_trials = 50,
                                                  seed = 2026,
                                                  es_frac = 0.1,
                                                  fit_age_range = c(42, 70),
                                                  ref_age_range = c(30, 100),
                                                  ref_step = 0.1,
                                                  nrounds = 5000,
                                                  early_stopping_rounds = 100,
                                                  verbose_inner = 0) {
  
  folds <- make_folds_stratified(as.integer(dat_outer_train$status), K = K, seed = seed)
  
  trials <- vector("list", n_trials)
  
  for (i in 1:n_trials) {
    params_full_i <- sample_params_full(seed = seed + 20000 + i)
    
    sc <- inner_cv_score_with_GM(
      dat_outer_train = dat_outer_train,
      feature_cols = feature_cols,
      params_base_fixed = params_base_fixed,
      params_full = params_full_i,
      folds = folds,
      es_frac = es_frac,
      seed = seed,
      fit_age_range = fit_age_range,
      ref_age_range = ref_age_range,
      ref_step = ref_step,
      nrounds = nrounds,
      early_stopping_rounds = early_stopping_rounds,
      verbose = verbose_inner
    )
    
    trials[[i]] <- list(
      trial = i,
      params_full = params_full_i,
      mean_cindex = sc$mean_cindex,
      sd_cindex = sc$sd_cindex,
      mean_best_iter_full = sc$mean_best_iter_full
    )
    
    message(sprintf(
      "Trial %d/%d | mean_cindex=%.4f (sd=%.4f) | mean_best_iter_full=%.1f",
      i, n_trials, sc$mean_cindex, sc$sd_cindex, sc$mean_best_iter_full
    ))
  }
  
  leaderboard <- do.call(rbind, lapply(trials, function(x) {
    p <- x$params_full
    data.frame(
      trial = x$trial,
      mean_cindex = x$mean_cindex,
      sd_cindex = x$sd_cindex,
      mean_best_iter_full = x$mean_best_iter_full,
      eta = p$eta,
      max_depth = p$max_depth,
      min_child_weight = p$min_child_weight,
      subsample = p$subsample,
      colsample_bytree = p$colsample_bytree,
      lambda = p$lambda,
      alpha = p$alpha,
      gamma = p$gamma,
      aft_scale = p$aft_loss_distribution_scale
    )
  }))
  
  leaderboard <- leaderboard[order(-leaderboard$mean_cindex, leaderboard$sd_cindex), ]
  best_trial <- leaderboard$trial[1]
  best_params_full <- trials[[best_trial]]$params_full
  
  list(
    folds = folds,
    leaderboard = leaderboard,
    best_params_full = best_params_full,
    trials = trials
  )
}

# ============================================================
# Final refit on OUTER_TRAIN with early stopping (10% of outer_train)
# Then evaluate on OUTER_TEST (external 30%)
# WITH GM curve and BA mapping on external test
# ============================================================
fit_final_and_eval_external_with_GM <- function(dat_outer_train,
                                                dat_outer_test,
                                                feature_cols,
                                                params_base_fixed,
                                                params_full_best,
                                                es_frac = 0.1,
                                                seed = 2026,
                                                fit_age_range = c(42, 70),
                                                ref_age_range = c(30, 100),
                                                ref_step = 0.1,
                                                nrounds = 5000,
                                                early_stopping_rounds = 100,
                                                verbose = 0) {
  
  set.seed(seed)
  
  time_tr   <- as.numeric(dat_outer_train$time_years)
  status_tr <- as.integer(dat_outer_train$status)
  age_tr    <- as.numeric(dat_outer_train$age0)
  
  time_te   <- as.numeric(dat_outer_test$time_years)
  status_te <- as.integer(dat_outer_test$status)
  age_te    <- as.numeric(dat_outer_test$age0)
  
  Xtr_full <- data.matrix(dat_outer_train[, feature_cols, drop = FALSE])
  Xte_full <- data.matrix(dat_outer_test[, feature_cols, drop = FALSE])
  
  Xtr_base <- matrix(age_tr, ncol = 1); colnames(Xtr_base) <- "age0"
  Xte_base <- matrix(age_te, ncol = 1); colnames(Xte_base) <- "age0"
  
  # Early-stop split inside outer_train
  sp_es <- split_stratified(status_tr, p_train = 1 - es_frac, p_valid = es_frac, seed = seed + 77)
  fit_idx <- sp_es$train
  es_idx  <- sp_es$valid
  
  # Imputer fit on fit_idx only; apply to outer_train and outer_test
  meds <- fit_imputer_median(Xtr_full[fit_idx, , drop = FALSE])
  Xtr_full <- apply_imputer_median(Xtr_full, meds)
  Xte_full <- apply_imputer_median(Xte_full, meds)
  
  # ---- Baseline ----
  d_fit_base <- make_aft_dmatrix(Xtr_base[fit_idx, , drop = FALSE], time_tr[fit_idx], status_tr[fit_idx])
  d_es_base  <- make_aft_dmatrix(Xtr_base[es_idx,  , drop = FALSE], time_tr[es_idx],  status_tr[es_idx])
  
  bst_base <- xgb.train(
    params = params_base_fixed,
    data = d_fit_base,
    nrounds = nrounds,
    watchlist = list(train = d_fit_base, es = d_es_base),
    early_stopping_rounds = early_stopping_rounds,
    verbose = verbose
  )
  
  base_fit <- predict(bst_base, Xtr_base[fit_idx, , drop = FALSE], outputmargin = TRUE)
  base_es  <- predict(bst_base, Xtr_base[es_idx,  , drop = FALSE], outputmargin = TRUE)
  base_te  <- predict(bst_base, Xte_base, outputmargin = TRUE)
  
  # ---- Full with offset ----
  d_fit_full <- make_aft_dmatrix(Xtr_full[fit_idx, , drop = FALSE], time_tr[fit_idx], status_tr[fit_idx], base_margin = base_fit)
  d_es_full  <- make_aft_dmatrix(Xtr_full[es_idx,  , drop = FALSE], time_tr[es_idx],  status_tr[es_idx],  base_margin = base_es)
  d_ext_full <- make_aft_dmatrix(Xte_full, time_te, status_te, base_margin = base_te)
  
  bst_full <- xgb.train(
    params = params_full_best,
    data = d_fit_full,
    nrounds = nrounds,
    watchlist = list(train = d_fit_full, es = d_es_full),
    early_stopping_rounds = early_stopping_rounds,
    verbose = verbose
  )
  
  # External predictions
  margin_ext <- predict(bst_full, d_ext_full, outputmargin = TRUE)
  risk_ext   <- -margin_ext
  cidx_ext   <- cindex_from_risk(time_te, status_te, risk_ext)
  
  # ---- GM reference curve from baseline; fit 42-70, extrapolate 30-100 ----
  age_grid_ref <- seq(ref_age_range[1], ref_age_range[2], by = ref_step)
  X_age_ref <- matrix(age_grid_ref, ncol = 1); colnames(X_age_ref) <- "age0"
  
  base_margin_ref <- predict(bst_base, X_age_ref, outputmargin = TRUE)
  risk_base_ref <- -base_margin_ref
  
  idx_fit_age <- age_grid_ref >= fit_age_range[1] & age_grid_ref <= fit_age_range[2]
  age_fit_gm <- age_grid_ref[idx_fit_age]
  y_fit_gm   <- as.numeric(risk_base_ref[idx_fit_age])
  
  gm_model <- fit_gompertz_makeham(age_fit_gm, y_fit_gm)
  
  if (is.null(gm_model)) {
    lm_ref <- lm(y_fit_gm ~ age_fit_gm)
    risk_grid_ref <- as.numeric(predict(lm_ref, newdata = data.frame(age_fit_gm = age_grid_ref)))
    gm_params <- c(A = NA, B = NA, C = NA)
  } else {
    cf <- coef(gm_model)
    gm_params <- cf
    risk_grid_ref <- as.numeric(cf["A"] + cf["B"] * exp(cf["C"] * age_grid_ref))
  }
  
  BA_ext <- risk_to_BA(risk_ext, age_grid_ref, risk_grid_ref)
  
  list(
    cindex_external_test = cidx_ext,
    preds_external = data.frame(margin = margin_ext, risk = risk_ext, BA = BA_ext),
    ref = list(age_grid = age_grid_ref, risk_grid = risk_grid_ref, gm_params = gm_params),
    bundle = list(
      bst_base = bst_base,
      bst_full = bst_full,
      feature_cols = feature_cols,
      impute_med = meds
    )
  )
}

# ============================================================
# Master pipeline: outer 70/30 -> tune on 70 -> refit -> eval on 30
# ============================================================
run_outer_holdout_tuning_with_GM <- function(dat_sex,
                                             feature_cols,
                                             params_base_fixed,
                                             K_inner = 5,
                                             n_trials = 50,
                                             seed = 2026,
                                             es_frac = 0.1,
                                             fit_age_range = c(42, 70),
                                             ref_age_range = c(30, 100),
                                             ref_step = 0.1,
                                             nrounds = 5000,
                                             early_stopping_rounds = 100,
                                             verbose_inner = 0,
                                             verbose_final = 0) {
  
  # Outer 70/30 stratified split
  sp_outer <- split_stratified(as.integer(dat_sex$status), p_train = 0.7, p_valid = 0.0, seed = seed)
  idx_outer_train <- sp_outer$train
  idx_outer_test  <- sp_outer$test
  
  dat_outer_train <- dat_sex[idx_outer_train, , drop = FALSE]
  dat_outer_test  <- dat_sex[idx_outer_test,  , drop = FALSE]
  
  # Tune full params on outer_train only
  tune_res <- tune_full_only_on_outer_train_with_GM(
    dat_outer_train = dat_outer_train,
    feature_cols = feature_cols,
    params_base_fixed = params_base_fixed,
    K = K_inner,
    n_trials = n_trials,
    seed = seed,
    es_frac = es_frac,
    fit_age_range = fit_age_range,
    ref_age_range = ref_age_range,
    ref_step = ref_step,
    nrounds = nrounds,
    early_stopping_rounds = early_stopping_rounds,
    verbose_inner = verbose_inner
  )
  
  best_params_full <- tune_res$best_params_full
  
  # Final refit on outer_train, evaluate on external outer_test
  final_res <- fit_final_and_eval_external_with_GM(
    dat_outer_train = dat_outer_train,
    dat_outer_test  = dat_outer_test,
    feature_cols = feature_cols,
    params_base_fixed = params_base_fixed,
    params_full_best  = best_params_full,
    es_frac = es_frac,
    seed = seed + 999,
    fit_age_range = fit_age_range,
    ref_age_range = ref_age_range,
    ref_step = ref_step,
    nrounds = nrounds,
    early_stopping_rounds = early_stopping_rounds,
    verbose = verbose_final
  )
  
  list(
    outer_split = list(train_idx = idx_outer_train, test_idx = idx_outer_test),
    tuning = tune_res,
    best_params_full = best_params_full,
    external = final_res
  )
}

# ============================================================
# Example usage
# ============================================================
# params_base_fixed <- params_base  # baseline is FIXED as you requested
#
# dat_f <- dat %>% dplyr::filter(Sex == 0)
# dat_m <- dat %>% dplyr::filter(Sex == 1)
#
# res_f <- run_outer_holdout_tuning_with_GM(
#   dat_sex = dat_f,
#   feature_cols = feature_cols,
#   params_base_fixed = params_base_fixed,
#   K_inner = 5,
#   n_trials = 80,
#   seed = 2026,
#   es_frac = 0.1,
#   fit_age_range = c(42, 70),
#   ref_age_range = c(30, 100),
#   ref_step = 0.1
# )
# head(res_f$tuning$leaderboard, 10)
# res_f$best_params_full
# res_f$external$cindex_external_test
#
# res_m <- run_outer_holdout_tuning_with_GM(
#   dat_sex = dat_m,
#   feature_cols = feature_cols,
#   params_base_fixed = params_base_fixed,
#   K_inner = 5,
#   n_trials = 80,
#   seed = 2026,
#   es_frac = 0.1,
#   fit_age_range = c(42, 70),
#   ref_age_range = c(30, 100),
#   ref_step = 0.1
# )
# res_m$external$cindex_external_test