score_BAA_batch_Gompertz_fast_hybrid <- function(new_df, bundle_list) {
  require(xgboost)
  require(dplyr)
  
  stopifnot(is.data.frame(new_df))
  stopifnot(nrow(new_df) >= 1)
  stopifnot(all(c("Sex", "age0") %in% names(new_df)))
  
  bundle_names <- names(bundle_list)
  
  # -----------------------------
  # 1) Select the correct sex-specific bundle
  # -----------------------------
  get_bundle <- function(sex01) {
    sex01 <- suppressWarnings(as.integer(sex01))
    
    if (is.na(sex01) || !(sex01 %in% c(0L, 1L))) {
      stop("`Sex` must be coded as 0/1 (0 = female, 1 = male).")
    }
    
    if (!is.null(bundle_names) && all(c("0", "1") %in% bundle_names)) {
      return(bundle_list[[as.character(sex01)]])
    }
    
    if (!is.null(bundle_names) && all(c("female", "male") %in% bundle_names)) {
      return(if (sex01 == 1L) bundle_list$male else bundle_list$female)
    }
    
    stop(
      "`bundle_list` names must be either c('0','1') or c('female','male'). Current names: ",
      paste(bundle_names, collapse = ", ")
    )
  }
  
  # -----------------------------
  # 2) Hybrid baseline margin helper
  # -----------------------------
  gm_margin_from_age <- function(age, A, B, C) {
    -(A + B * exp(C * age))
  }
  
  hybrid_base_margin <- function(age, base_margin_raw, gm_params, fit_age_range) {
    out <- as.numeric(base_margin_raw)
    
    # If Gompertz fitting failed, keep the raw baseline everywhere
    if (any(!is.finite(gm_params))) {
      return(out)
    }
    
    outside_idx <- age < fit_age_range[1] | age > fit_age_range[2]
    if (any(outside_idx)) {
      out[outside_idx] <- gm_margin_from_age(
        age = age[outside_idx],
        A = gm_params["A"],
        B = gm_params["B"],
        C = gm_params["C"]
      )
    }
    
    out
  }
  
  # -----------------------------
  # 3) Stable monotone inverse mapping from risk to BA
  # -----------------------------
  inverse_from_grid <- function(risk, age_grid, risk_grid) {
    out <- rep(NA_real_, length(risk))
    
    if (length(age_grid) == 0 || length(risk_grid) == 0) {
      return(out)
    }
    
    if (!all(is.finite(age_grid)) || !all(is.finite(risk_grid))) {
      return(out)
    }
    
    ok <- is.finite(risk)
    if (!any(ok)) return(out)
    
    iso_fit <- isoreg(age_grid, risk_grid)
    risk_iso <- iso_fit$yf
    keep <- !duplicated(risk_iso)
    
    out[ok] <- approx(
      x = risk_iso[keep],
      y = age_grid[keep],
      xout = risk[ok],
      rule = 2
    )$y
    
    out
  }
  
  # -----------------------------
  # 4) Sanity check Sex values
  # -----------------------------
  sex_values <- sort(unique(as.integer(new_df$Sex)))
  
  if (any(is.na(sex_values)) || !all(sex_values %in% c(0L, 1L))) {
    stop("`Sex` column must contain only 0/1 values.")
  }
  
  result <- vector("list", length(sex_values))
  part_id <- 1L
  
  # -----------------------------
  # 5) Score each sex separately
  # -----------------------------
  for (sex_i in sex_values) {
    idx <- which(as.integer(new_df$Sex) == sex_i)
    df_sub <- new_df[idx, , drop = FALSE]
    b <- get_bundle(sex_i)
    
    if (is.null(b)) stop("No matching bundle found for Sex = ", sex_i)
    if (is.null(b$bst_base)) stop("The bundle is missing `bst_base`.")
    if (is.null(b$bst_full)) stop("The bundle is missing `bst_full`.")
    if (is.null(b$feature_cols) || length(b$feature_cols) == 0) {
      stop("The bundle is missing valid `feature_cols`.")
    }
    if (is.null(b$impute_med) || length(b$impute_med) != length(b$feature_cols)) {
      stop("The bundle is missing valid `impute_med` aligned with `feature_cols`.")
    }
    if (is.null(b$age_grid) || length(b$age_grid) == 0) {
      stop("The bundle is missing a valid `age_grid`.")
    }
    if (is.null(b$risk_grid) || length(b$risk_grid) == 0) {
      stop("The bundle is missing a valid `risk_grid`.")
    }
    if (is.null(b$gm_fit_age_range) || length(b$gm_fit_age_range) != 2) {
      stop("The bundle is missing a valid `gm_fit_age_range`.")
    }
    if (is.null(b$gm_params)) {
      stop("The bundle is missing `gm_params`.")
    }
    if (is.null(b$baa_lm_coef)) {
      stop("The bundle is missing `baa_lm_coef`.")
    }
    
    age_vec <- as.numeric(df_sub$age0)
    
    # -----------------------------
    # 5a) Raw empirical baseline from bst_base
    # -----------------------------
    X_base <- matrix(age_vec, ncol = 1)
    colnames(X_base) <- "age0"
    
    base_margin_raw <- as.numeric(
      predict(b$bst_base, X_base, outputmargin = TRUE)
    )
    
    # -----------------------------
    # 5b) Hybrid offset: raw baseline inside range, Gompertz outside
    # -----------------------------
    base_margin_hybrid <- hybrid_base_margin(
      age = age_vec,
      base_margin_raw = base_margin_raw,
      gm_params = b$gm_params,
      fit_age_range = b$gm_fit_age_range
    )
    
    # -----------------------------
    # 5c) Build full feature matrix
    # -----------------------------
    feature_cols <- b$feature_cols
    missing_cols <- setdiff(feature_cols, names(df_sub))
    if (length(missing_cols) > 0) {
      for (m in missing_cols) {
        df_sub[[m]] <- NA_real_
      }
    }
    
    X_full <- data.matrix(df_sub[, feature_cols, drop = FALSE])
    
    impute_med <- as.numeric(b$impute_med)
    for (j in seq_along(feature_cols)) {
      na_idx <- is.na(X_full[, j])
      if (any(na_idx)) {
        X_full[na_idx, j] <- impute_med[j]
      }
    }
    
    # -----------------------------
    # 5d) Predict full model using the hybrid offset
    # -----------------------------
    dnew <- xgboost::xgb.DMatrix(X_full)
    xgboost::setinfo(dnew, "base_margin", base_margin_hybrid)
    
    pred_margin <- as.numeric(
      predict(b$bst_full, dnew, outputmargin = TRUE)
    )
    
    risk <- -pred_margin
    
    # -----------------------------
    # 5e) Risk -> BA
    # -----------------------------
    BA <- inverse_from_grid(
      risk = risk,
      age_grid = b$age_grid,
      risk_grid = b$risk_grid
    )
    
    raw_BAA <- BA - age_vec
    
    # -----------------------------
    # 5f) Compute PC_BAA using stored residualization coefficients
    # -----------------------------
    beta <- b$baa_lm_coef
    
    if (all(c("beta0", "beta1", "beta2") %in% names(beta))) {
      beta0 <- as.numeric(beta["beta0"])
      beta1 <- as.numeric(beta["beta1"])
      beta2 <- as.numeric(beta["beta2"])
    } else if (all(c("(Intercept)", "age0", "I(age0^2)") %in% names(beta))) {
      beta0 <- as.numeric(beta["(Intercept)"])
      beta1 <- as.numeric(beta["age0"])
      beta2 <- as.numeric(beta["I(age0^2)"])
    } else {
      stop("`baa_lm_coef` does not contain expected coefficient names.")
    }
    
    age_trend <- beta0 + beta1 * age_vec + beta2 * age_vec^2
    PC_BAA <- raw_BAA - age_trend
    
    result[[part_id]] <- data.frame(
      PC_BA = BA,
      raw_BAA = raw_BAA,
      PC_BAA = PC_BAA,
      risk = risk,
      pred_margin = pred_margin,
      base_margin_raw = base_margin_raw,
      base_margin_hybrid = base_margin_hybrid,
      Sex = as.integer(df_sub$Sex),
      age0 = age_vec,
      row_index = idx
    )
    
    part_id <- part_id + 1L
  }
  
  out <- dplyr::bind_rows(result)
  out <- out[order(out$row_index), , drop = FALSE]
  rownames(out) <- NULL
  out$row_index <- NULL
  
  out
}





















