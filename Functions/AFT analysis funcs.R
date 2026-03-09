library(minpack.lm)
prepare_aft_table <- function(meas_df,
                              death_df,
                              censor_date = as.Date("2024-12-31"),
                              id_col = "Participant.ID",
                              sex_col = "Sex",
                              inst_col = "Instance",
                              age_base_col = "Age.when.attended.assessment.centre",
                              ac_date_cols = c("Date.of.attending.assessment.centre...Instance.0",
                                               "Date.of.attending.assessment.centre...Instance.1",
                                               "Date.of.attending.assessment.centre...Instance.2",
                                               "Date.of.attending.assessment.centre...Instance.3"),
                              death_age_col  = "Age.at.death...Instance.0",
                              death_date_col = "Date.of.death...Instance.0",
                              yob_col = "Year.of.birth",
                              mob_col = "Month.of.birth",
                              birth_day = 15L,
                              keep_assess_date = TRUE) {
  
  # -----------------------------
  # Helpers
  # -----------------------------
  parse_ymd <- function(x) {
    x_chr <- trimws(as.character(x))
    x_chr[x_chr == "" | x_chr == "NA"] <- NA_character_
    suppressWarnings(as.Date(x_chr, format = "%Y-%m-%d"))
  }
  
  pick_first_non_na <- function(x) {
    x2 <- x[!is.na(x)]
    if (length(x2) == 0) NA else x2[1]
  }
  
  # -----------------------------
  # 0) Column checks
  # -----------------------------
  need_meas <- c(id_col, inst_col, age_base_col)
  miss_meas <- setdiff(need_meas, names(meas_df))
  if (length(miss_meas) > 0) stop("meas_df is missing columns: ", paste(miss_meas, collapse = ", "))
  
  need_death <- c(id_col, death_age_col, death_date_col, yob_col, mob_col)
  miss_death <- setdiff(need_death, names(death_df))
  if (length(miss_death) > 0) stop("death_df is missing columns: ", paste(miss_death, collapse = ", "))
  
  # -----------------------------
  # 1) Baseline measurements: one row per person
  # -----------------------------
  meas_keep <- meas_df %>%
    dplyr::filter(.data[[inst_col]] %in% c(0, 1, 2, 3)) %>%
    dplyr::mutate(inst_tmp = suppressWarnings(as.integer(.data[[inst_col]])))
  
  meas0 <- meas_keep %>%
    dplyr::filter(inst_tmp == 0L) %>%
    dplyr::distinct(.data[[id_col]], .keep_all = TRUE)
  
  missing_ids <- setdiff(unique(meas_keep[[id_col]]), unique(meas0[[id_col]]))
  if (length(missing_ids) > 0) {
    meas_fallback <- meas_keep %>%
      dplyr::filter(.data[[id_col]] %in% missing_ids) %>%
      dplyr::group_by(.data[[id_col]]) %>%
      dplyr::arrange(inst_tmp, .by_group = TRUE) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()
    meas0 <- dplyr::bind_rows(meas0, meas_fallback)
  }
  
  meas0 <- meas0 %>%
    dplyr::mutate(age0 = suppressWarnings(as.numeric(.data[[age_base_col]])))
  
  # -----------------------------
  # 2) Death table: clean and aggregate to one row per person
  #    ALSO bring assessment-centre dates from death_df (they live here in your data)
  # -----------------------------
  # ensure all ac_date_cols exist
  for (cc in ac_date_cols) if (!(cc %in% names(death_df))) death_df[[cc]] <- NA_character_
  
  death_clean <- death_df %>%
    dplyr::transmute(
      !!id_col := .data[[id_col]],
      death_age  = suppressWarnings(as.numeric(.data[[death_age_col]])),
      death_date = parse_ymd(.data[[death_date_col]]),
      yob = suppressWarnings(as.integer(.data[[yob_col]])),
      mob = suppressWarnings(as.integer(.data[[mob_col]])),
      ac0 = parse_ymd(.data[[ac_date_cols[1]]]),
      ac1 = parse_ymd(.data[[ac_date_cols[2]]]),
      ac2 = parse_ymd(.data[[ac_date_cols[3]]]),
      ac3 = parse_ymd(.data[[ac_date_cols[4]]])
    ) %>%
    dplyr::mutate(
      mob = ifelse(is.na(mob) | mob < 1 | mob > 12, 6L, mob),
      yob = ifelse(is.na(yob) | yob < 1800 | yob > 2100, NA_integer_, yob),
      birth_date = dplyr::if_else(
        is.na(yob),
        as.Date(NA),
        as.Date(sprintf("%04d-%02d-%02d", yob, mob, birth_day))
      )
    )
  
  death0 <- death_clean %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      death_age  = pick_first_non_na(death_age),
      death_date = pick_first_non_na(death_date),
      birth_date = pick_first_non_na(birth_date),
      yob        = pick_first_non_na(yob),
      mob        = pick_first_non_na(mob),
      ac0        = pick_first_non_na(ac0),
      ac1        = pick_first_non_na(ac1),
      ac2        = pick_first_non_na(ac2),
      ac3        = pick_first_non_na(ac3),
      .groups = "drop"
    )
  
  # -----------------------------
  # 3) Merge
  # -----------------------------
  dat <- meas0 %>%
    dplyr::left_join(death0, by = id_col)
  
  # baseline assessment date: match by the chosen baseline instance, fallback to ac0
  dat <- dat %>%
    dplyr::mutate(
      assess_date0 = dplyr::case_when(
        inst_tmp == 0L ~ ac0,
        inst_tmp == 1L ~ dplyr::coalesce(ac1, ac0),
        inst_tmp == 2L ~ dplyr::coalesce(ac2, ac0),
        inst_tmp == 3L ~ dplyr::coalesce(ac3, ac0),
        TRUE ~ ac0
      )
    )
  
  # -----------------------------
  # 4) status + age0 refined + time_years
  # -----------------------------
  dat <- dat %>%
    dplyr::mutate(
      status = ifelse(!is.na(death_age) | !is.na(death_date), 1L, 0L),
      
      censor_age = ifelse(
        !is.na(birth_date),
        as.numeric(censor_date - birth_date) / 365.25,
        NA_real_
      ),
      
      death_age_from_date = ifelse(
        !is.na(death_date) & !is.na(birth_date),
        as.numeric(death_date - birth_date) / 365.25,
        NA_real_
      ),
      death_age_final = dplyr::coalesce(death_age, death_age_from_date),
      
      age0_from_date = ifelse(
        !is.na(assess_date0) & !is.na(birth_date),
        as.numeric(assess_date0 - birth_date) / 365.25,
        NA_real_
      ),
      age0 = dplyr::coalesce(age0_from_date, age0),
      
      time_years = dplyr::case_when(
        status == 1L & !is.na(death_date) & !is.na(assess_date0) ~ as.numeric(death_date  - assess_date0) / 365.25,
        status == 0L & !is.na(assess_date0)                     ~ as.numeric(censor_date - assess_date0) / 365.25,
        status == 1L ~ death_age_final - age0,
        status == 0L ~ censor_age     - age0,
        TRUE ~ NA_real_
      )
    )
  
  # -----------------------------
  # 5) Filter invalid rows
  # -----------------------------
  dat <- dat %>%
    dplyr::filter(
      !is.na(age0),
      !is.na(status),
      !is.na(time_years),
      is.finite(time_years),
      time_years > 0,
      !(status == 0L & is.na(assess_date0) & is.na(censor_age)),
      !(status == 1L & is.na(death_date) & is.na(death_age_final))
    ) %>%
    dplyr::select(-dplyr::any_of(c("ac0","ac1","ac2","ac3","age0_from_date")))
  
  dat <- dat %>% dplyr::select(-dplyr::any_of("inst_tmp"))
  
  if (!keep_assess_date) {
    dat <- dat %>% dplyr::select(-dplyr::any_of("assess_date0"))
  }
  
  dat
}





make_aft_dmatrix <- function(Xm, t, s, base_margin = NULL) {
  d <- xgb.DMatrix(Xm)
  setinfo(d, "label_lower_bound", t)
  setinfo(d, "label_upper_bound", ifelse(s == 1L, t, Inf))
  if (!is.null(base_margin)) setinfo(d, "base_margin", base_margin)
  d
}

split_stratified <- function(status, p_train=0.7, p_valid=0.1, seed=2026) {
  set.seed(seed)
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

train_sex_specific_clock <- function(dat_sex,
                                     feature_cols,
                                     params_base,
                                     params_full,
                                     seed=2026,
                                     censor_check=TRUE) {
  
  # 基本向量
  time   <- as.numeric(dat_sex$time_years)
  status <- as.integer(dat_sex$status)
  age0   <- as.numeric(dat_sex$age0)
  
  if (censor_check) {
    stopifnot(all(is.finite(time)), all(time > 0), all(status %in% c(0L,1L)))
  }
  
  # 划分（分层）
  sp <- split_stratified(status, seed=seed)
  tr <- sp$train; va <- sp$valid; te <- sp$test
  
  # 特征矩阵（full：biomarkers）
  X_full <- data.matrix(dat_sex[, feature_cols, drop=FALSE])
  
  # baseline：只用 age0（sex 固定）
  X_base <- matrix(age0, ncol = 1)
  colnames(X_base) <- "age0"
  
  # 缺失填补：按该 sex 子集的中位数
  impute_med <- apply(X_full, 2, function(x) median(x, na.rm = TRUE))
  for (j in seq_len(ncol(X_full))) {
    na_idx <- is.na(X_full[, j])
    if (any(na_idx)) X_full[na_idx, j] <- impute_med[j]
  }
  
  # ---- 1) baseline AFT（age only）----
  dtrain_base <- make_aft_dmatrix(X_base[tr, , drop=FALSE], time[tr], status[tr])
  dvalid_base <- make_aft_dmatrix(X_base[va, , drop=FALSE], time[va], status[va])
  
  bst_base <- xgb.train(
    params = params_base,
    data = dtrain_base,
    nrounds = 5000,
    watchlist = list(train=dtrain_base, valid=dvalid_base),
    early_stopping_rounds = 100,
    verbose = 1
  )
  
  base_tr  <- predict(bst_base, X_base[tr, , drop=FALSE], outputmargin=TRUE)
  base_va  <- predict(bst_base, X_base[va, , drop=FALSE], outputmargin=TRUE)
  base_all <- predict(bst_base, X_base,               outputmargin=TRUE)
  
  # ---- 2) full AFT（biomarkers + base_margin offset）----
  dtrain_full <- make_aft_dmatrix(X_full[tr, , drop=FALSE], time[tr], status[tr], base_margin = base_tr)
  dvalid_full <- make_aft_dmatrix(X_full[va, , drop=FALSE], time[va], status[va], base_margin = base_va)
  
  bst_full <- xgb.train(
    params = params_full,
    data = dtrain_full,
    nrounds = 5000,
    watchlist = list(train=dtrain_full, valid=dvalid_full),
    early_stopping_rounds = 100,
    verbose = 1
  )
  
  # ---- 3) 生成 risk，并做 risk->BA 单调映射（用训练集）----
  d_all <- xgb.DMatrix(X_full)
  setinfo(d_all, "base_margin", base_all)
  resid_all <- predict(bst_full, d_all, outputmargin=TRUE)
  margin_all <-  resid_all
  risk_all <- -margin_all
  
  # 用年龄网格得到 g(age)=E[risk|age]，并强制单调（isoreg）
  age_grid <- seq(floor(min(age0, na.rm=TRUE)), ceiling(max(age0, na.rm=TRUE)), by=0.1)
  
  # 局部平均（窗口0.5岁）
  risk_tr <- risk_all[tr]
  age_tr  <- age0[tr]
  g <- sapply(age_grid, function(a) mean(risk_tr[abs(age_tr - a) < 0.5], na.rm=TRUE))
  ok <- is.finite(g)
  age_grid <- age_grid[ok]; g <- g[ok]
  
  # 单调化（risk 应随年龄上升；若你发现相反，就把 g 取负或改排序）
  iso <- isoreg(age_grid, g)
  risk_grid_mono <- iso$yf  # 单调后的风险曲线
  
  bundle <- list(
    bst_base = bst_base,
    bst_full = bst_full,
    feature_cols = feature_cols,
    impute_med = impute_med,
    age_grid = age_grid,
    risk_grid = risk_grid_mono,
    split = sp   # 可选：保存划分便于复现实验
  )
  
  # 也给你返回 test 上的 C-index（用 risk = -margin）
  library(survival)
  pred_margin_te <- margin_all[te]
  fit <- coxph(Surv(time[te], status[te]) ~ I(-pred_margin_te))
  cindex_te <- as.numeric(summary(fit)$concordance[1])
  
  list(bundle = bundle, cindex_test = cindex_te, 
       tr_set = tr, va_set = va, te_set = te)
}

train_sex_specific_clock_Gompertz <- function(dat_sex,
                                              feature_cols,
                                              params_base,
                                              params_full,
                                              seed = 2026,
                                              censor_check = TRUE,
                                              fit_age_range = c(42, 70),   # 用于拟合 GM 的年龄段
                                              ref_age_range = c(37, 85),   # 外推到的年龄段
                                              ref_step = 0.1) {            # reference 网格步长
  
  # -----------------------
  # 0) 基本向量
  # -----------------------
  time   <- as.numeric(dat_sex$time_years)
  status <- as.integer(dat_sex$status)
  age0   <- as.numeric(dat_sex$age0)
  
  if (censor_check) {
    stopifnot(all(is.finite(time)), all(time > 0), all(status %in% c(0L,1L)))
  }
  
  # 分层划分
  sp <- split_stratified(status, seed = seed)
  tr <- sp$train; va <- sp$valid; te <- sp$test
  
  # 特征矩阵
  X_full <- data.matrix(dat_sex[, feature_cols, drop = FALSE])
  
  # baseline: 只用 age0
  X_base <- matrix(age0, ncol = 1)
  colnames(X_base) <- "age0"
  
  # -----------------------
  # 1) 缺失填补（建议：只用训练集算 median，避免 leakage）
  # -----------------------
  impute_med <- apply(X_full[tr, , drop = FALSE], 2, function(x) median(x, na.rm = TRUE))
  for (j in seq_len(ncol(X_full))) {
    na_idx <- is.na(X_full[, j])
    if (any(na_idx)) X_full[na_idx, j] <- impute_med[j]
  }
  
  # -----------------------
  # 2) baseline AFT（age only）
  # -----------------------
  dtrain_base <- make_aft_dmatrix(X_base[tr, , drop = FALSE], time[tr], status[tr])
  dvalid_base <- make_aft_dmatrix(X_base[va, , drop = FALSE], time[va], status[va])
  
  bst_base <- xgb.train(
    params = params_base,
    data = dtrain_base,
    nrounds = 5000,
    watchlist = list(train = dtrain_base, valid = dvalid_base),
    early_stopping_rounds = 100,
    verbose = 1
  )
  
  base_tr  <- predict(bst_base, X_base[tr, , drop = FALSE], outputmargin = TRUE)
  base_va  <- predict(bst_base, X_base[va, , drop = FALSE], outputmargin = TRUE)
  base_all <- predict(bst_base, X_base,               outputmargin = TRUE)
  
  # -----------------------
  # 3) full AFT（biomarkers + base_margin offset）
  # -----------------------
  dtrain_full <- make_aft_dmatrix(X_full[tr, , drop = FALSE], time[tr], status[tr], base_margin = base_tr)
  dvalid_full <- make_aft_dmatrix(X_full[va, , drop = FALSE], time[va], status[va], base_margin = base_va)
  
  bst_full <- xgb.train(
    params = params_full,
    data = dtrain_full,
    nrounds = 5000,
    watchlist = list(train = dtrain_full, valid = dvalid_full),
    early_stopping_rounds = 100,
    verbose = 1
  )
  
  # -----------------------
  # 4) 生成个体 risk（full 模型输出）
  #    注意：predict(bst_full, d_all) 在设置 base_margin 后通常返回 full margin
  # -----------------------
  d_all <- xgb.DMatrix(X_full)
  setinfo(d_all, "base_margin", base_all)
  margin_all <- predict(bst_full, d_all, outputmargin = TRUE)  # full μ_hat
  risk_all   <- -margin_all
  
  # =========================================================
  # 5) Reference curve: 用 bst_base 的 risk-score vs age
  #    在 42–70 拟合 Gompertz–Makeham（A + B exp(C age)）
  #    外推到 37–85
  # =========================================================
  age_grid_ref <- seq(ref_age_range[1], ref_age_range[2], by = ref_step)
  X_age_ref <- matrix(age_grid_ref, ncol = 1)
  colnames(X_age_ref) <- "age0"
  
  base_margin_ref <- predict(bst_base, X_age_ref, outputmargin = TRUE)
  risk_base_ref <- -base_margin_ref   # 这是 bst_base 的 baseline risk-score
  
  # 只用 42–70 来拟合参数
  idx_fit <- age_grid_ref >= fit_age_range[1] & age_grid_ref <= fit_age_range[2]
  age_fit <- age_grid_ref[idx_fit]
  y_fit   <- as.numeric(risk_base_ref[idx_fit])
  
  # Gompertz–Makeham (用于 risk 标尺的参数化近似): y = A + B * exp(C * age)
  # 约束：B>0, C>=0（保证单调不降）
  fit_gm <- function(age, y) {
    A0 <- as.numeric(quantile(y, 0.05, na.rm = TRUE)) - 0.01
    y_shift <- pmax(y - A0, 1e-6)
    lm0 <- lm(log(y_shift) ~ age)
    C0 <- max(0, as.numeric(coef(lm0)[2]))
    B0 <- exp(as.numeric(coef(lm0)[1]))
    
    fit <- tryCatch(
      nlsLM(
        y ~ A + B * exp(C * age),
        start = list(A = A0, B = B0, C = C0),
        lower = c(A = -Inf, B = 0, C = 0),
        control = nls.lm.control(maxiter = 500)
      ),
      warning = function(w) { message("nlsLM warning: ", conditionMessage(w)); invokeRestart("muffleWarning") },
      error   = function(e) { message("nlsLM error: ", conditionMessage(e)); NULL }
    )
    print(fit)
  }
  gm_model <- fit_gm(age_fit, y_fit)
  
  if (is.null(gm_model)) {
    warning("Gompertz–Makeham 拟合失败：将退化为线性拟合做外推（仅用于 reference curve）。")
    lm_ref <- lm(y_fit ~ age_fit)
    risk_grid_ref <- as.numeric(predict(lm_ref, newdata = data.frame(age_fit = age_grid_ref)))
    gm_params <- c(A = NA, B = NA, C = NA)
  } else {
    cf <- coef(gm_model)
    gm_params <- cf
    risk_grid_ref <- as.numeric(cf["A"] + cf["B"] * exp(cf["C"] * age_grid_ref))
  }
  
  # -----------------------
  # 6) 打包输出
  # -----------------------
  bundle <- list(
    bst_base = bst_base,
    bst_full = bst_full,
    feature_cols = feature_cols,
    impute_med = impute_med,
    
    # reference curve（用于 risk -> BA 反查）
    age_grid = age_grid_ref,
    risk_grid = risk_grid_ref,
    
    # 诊断信息（可选）
    risk_base_ref_raw = as.numeric(risk_base_ref),
    gm_fit_age_range = fit_age_range,
    gm_params = gm_params,
    
    split = sp
  )
  
  # -----------------------
  # 7) test C-index（用 full risk）
  # -----------------------
  library(survival)
  pred_margin_te <- margin_all[te]
  fit <- coxph(Surv(time[te], status[te]) ~ I(-pred_margin_te))
  cindex_te <- as.numeric(summary(fit)$concordance[1])
  
  list(bundle = bundle, cindex_test = cindex_te,
       tr_set = tr, va_set = va, te_set = te)
}








score_BAA_one <- function(new_df, bundle_list) {
  stopifnot(nrow(new_df) >= 1)
  stopifnot(all(c("Sex","age0") %in% names(new_df)))
  
  # 允许多行输入
  out <- vector("list", nrow(new_df))
  
  for (i in seq_len(nrow(new_df))) {
    sex_i <- as.character(as.integer(new_df$Sex[i]))
    b <- bundle_list[[sex_i]]
    if (is.null(b)) stop("Sex=", sex_i, " 没有对应 bundle（检查 Sex 编码是不是 0/1）")
    
    # base: age0
    X_base <- matrix(as.numeric(new_df$age0[i]), ncol=1)
    colnames(X_base) <- "age0"
    base_m <- predict(b$bst_base, X_base, outputmargin=TRUE)
    
    # full: biomarkers（缺列补 NA，再用训练中位数填 NA）
    fc <- b$feature_cols
    row <- new_df[i, , drop=FALSE]
    miss <- setdiff(fc, names(row))
    if (length(miss) > 0) for (m in miss) row[[m]] <- NA_real_
    
    X_full <- data.matrix(row[, fc, drop=FALSE])
    for (j in seq_along(fc)) {
      if (is.na(X_full[1, j])) X_full[1, j] <- b$impute_med[[j]]
    }
    
    dnew <- xgb.DMatrix(X_full)
    setinfo(dnew, "base_margin", base_m)
    resid_m <- predict(b$bst_full, dnew, outputmargin=TRUE)
    
    margin <- base_m + resid_m
    risk <- -margin
    
    BA <- approx(x = b$risk_grid, y = b$age_grid, xout = risk, rule = 2)$y
    BAA <- BA - as.numeric(new_df$age0[i])
    
    out[[i]] <- data.frame(BA=BA, BAA=BAA, risk=risk, pred_margin=margin, Sex=new_df$Sex[i], age0=new_df$age0[i])
  }
  
  bind_rows(out)
}

score_BAA_batch <- function(new_df, bundle_list) {
  stopifnot(nrow(new_df) >= 1)
  stopifnot(all(c("Sex","age0") %in% names(new_df)))
  
  # 兼容 bundle_list 的命名方式
  nms <- names(bundle_list)
  get_bundle <- function(sex01) {
    sex01 <- suppressWarnings(as.integer(sex01))
    if (is.na(sex01) || !(sex01 %in% c(0L, 1L))) stop("Sex 必须是 0/1（0=女,1=男）")
    
    if (!is.null(nms) && all(c("0","1") %in% nms)) {
      return(bundle_list[[as.character(sex01)]])
    }
    if (!is.null(nms) && all(c("female","male") %in% nms)) {
      return(if (sex01 == 1L) bundle_list$male else bundle_list$female)
    }
    stop("bundle_list 的 names 必须是 c('0','1') 或 c('female','male')；当前: ",
         paste(nms, collapse = ", "))
  }
  
  out <- vector("list", nrow(new_df))
  
  for (i in seq_len(nrow(new_df))) {
    b <- get_bundle(new_df$Sex[i])
    if (is.null(b)) stop("Sex=", new_df$Sex[i], " 没有对应 bundle")
    
    # base: age0
    X_base <- matrix(as.numeric(new_df$age0[i]), ncol = 1)
    colnames(X_base) <- "age0"
    base_m <- predict(b$bst_base, X_base, outputmargin = TRUE)
    
    # full: biomarkers（缺列补 NA，再用训练中位数填 NA）
    fc <- b$feature_cols
    row <- new_df[i, , drop = FALSE]
    miss <- setdiff(fc, names(row))
    if (length(miss) > 0) for (m in miss) row[[m]] <- NA_real_
    
    X_full <- data.matrix(row[, fc, drop = FALSE])
    for (j in seq_along(fc)) {
      if (is.na(X_full[1, j])) X_full[1, j] <- b$impute_med[[j]]
    }
    
    dnew <- xgboost::xgb.DMatrix(X_full)
    xgboost::setinfo(dnew, "base_margin", base_m)   # <- 关键：显式 xgboost::
    resid_m <- predict(b$bst_full, dnew, outputmargin = TRUE)
    
    margin <- resid_m
    risk <- -margin
    
    BA <- approx(x = b$risk_grid, y = b$age_grid, xout = risk, rule = 2)$y
    BAA <- BA - as.numeric(new_df$age0[i])
    
    out[[i]] <- data.frame(
      BA = BA, BAA = BAA, risk = risk, pred_margin = margin,
      Sex = new_df$Sex[i], age0 = new_df$age0[i]
    )
  }
  
  dplyr::bind_rows(out)
}


score_BAA_batch_Gompertz <- function(new_df, bundle_list) {
  stopifnot(nrow(new_df) >= 1)
  stopifnot(all(c("Sex","age0") %in% names(new_df)))
  
  nms <- names(bundle_list)
  get_bundle <- function(sex01) {
    sex01 <- suppressWarnings(as.integer(sex01))
    if (is.na(sex01) || !(sex01 %in% c(0L, 1L))) stop("Sex 必须是 0/1（0=女,1=男）")
    
    if (!is.null(nms) && all(c("0","1") %in% nms)) {
      return(bundle_list[[as.character(sex01)]])
    }
    if (!is.null(nms) && all(c("female","male") %in% nms)) {
      return(if (sex01 == 1L) bundle_list$male else bundle_list$female)
    }
    stop("bundle_list 的 names 必须是 c('0','1') 或 c('female','male')；当前: ",
         paste(nms, collapse = ", "))
  }
  
  # GM baseline μ(age) and inverse mapping BA(risk)
  gm_mu <- function(age, A, B, C) {
    # mu = -risk, risk = A + B exp(C age)
    -(A + B * exp(C * age))
  }
  gm_inverse_ba <- function(risk, A, B, C, ba_min, ba_max) {
    if (!is.finite(risk)) return(NA_real_)
    if (!is.finite(A) || !is.finite(B) || !is.finite(C) || B <= 0 || C <= 0) return(NA_real_)
    
    if (risk <= A) return(ba_min)  
    
    ba <- log((risk - A) / B) / C
    
    ba <- max(min(ba, ba_max), ba_min)
    ba
  }
  
  out <- vector("list", nrow(new_df))
  
  for (i in seq_len(nrow(new_df))) {
    b <- get_bundle(new_df$Sex[i])
    if (is.null(b)) stop("Sex=", new_df$Sex[i], " 没有对应 bundle")
    if (is.null(b$gm_params) || any(!c("A","B","C") %in% names(b$gm_params))) {
      stop("bundle 中缺少 gm_params（A,B,C）。请确认你用的是 train_sex_specific_clock_Gompertz() 的输出。")
    }
    
    # reference 年龄范围（自动从 age_grid 推断）
    ba_min <- min(b$age_grid, na.rm = TRUE)
    ba_max <- max(b$age_grid, na.rm = TRUE)
    
    # 取 GM 参数
    A <- as.numeric(b$gm_params["A"])
    B <- as.numeric(b$gm_params["B"])
    C <- as.numeric(b$gm_params["C"])
    
    # --- GM baseline margin (μ_GM) ---
    age_i <- as.numeric(new_df$age0[i])
    base_m <- gm_mu(age_i, A, B, C)   # 关键：不用 bst_base
    
    # --- full biomarkers（缺列补 NA，再用训练中位数填 NA） ---
    fc <- b$feature_cols
    row <- new_df[i, , drop = FALSE]
    miss <- setdiff(fc, names(row))
    if (length(miss) > 0) for (m in miss) row[[m]] <- NA_real_
    
    X_full <- data.matrix(row[, fc, drop = FALSE])
    for (j in seq_along(fc)) {
      if (is.na(X_full[1, j])) X_full[1, j] <- b$impute_med[[j]]
    }
    
    dnew <- xgboost::xgb.DMatrix(X_full)
    xgboost::setinfo(dnew, "base_margin", base_m)
    
    # full μ_hat（在 GM baseline offset 上叠加 biomarker residual）
    margin <- predict(b$bst_full, dnew, outputmargin = TRUE)
    margin <- as.numeric(margin)
    
    risk <- -margin
    
    # --- BA by analytic inverse of GM reference curve ---
    BA  <- gm_inverse_ba(risk, A, B, C, ba_min, ba_max)
    BAA <- BA - age_i
    
    out[[i]] <- data.frame(
      BA = BA, BAA = BAA,
      risk = risk, pred_margin = margin,
      Sex = as.integer(new_df$Sex[i]),
      age0 = age_i
    )
  }
  
  dplyr::bind_rows(out)
}







