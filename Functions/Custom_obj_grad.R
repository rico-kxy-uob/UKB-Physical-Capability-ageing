custom_age_obj_fast <- function(preds, dtrain, age, window = 1) {
  
  y_true <- xgboost::getinfo(dtrain, "label")
  grad <- preds - y_true
  hess <- rep(1, length(grad))
  
  n <- length(grad)
  
  ord <- order(age)
  age_sorted <- age[ord]
  grad_sorted <- grad[ord]
  
  age_mean_sorted <- numeric(n)
  counts_sorted <- numeric(n)
  
  left <- 1
  right <- 1
  
  for (i in 1:n) {
    while (age_sorted[left] < age_sorted[i] - window && left < n) {
      left <- left + 1
    }
    while (right < n && age_sorted[right + 1] <= age_sorted[i] + window) {
      right <- right + 1
    }
    
    age_mean_sorted[i] <- mean(grad_sorted[left:right])
    counts_sorted[i] <- right - left + 1
  }
  
  age_mean_per_sample <- numeric(n)
  counts_per_sample   <- numeric(n)
  
  age_mean_per_sample[ord] <- age_mean_sorted
  counts_per_sample[ord]   <- counts_sorted
  
  total_mean <- mean(grad)
  if (abs(total_mean) < 1e-6) total_mean <- 1e-6
  
  grad_corrected <-  grad *
    abs((age_mean_per_sample/total_mean))
  
  return(list(grad = grad_corrected, hess = hess))
}

make_age_bin <- function(age, width = 5) {
  floor(age / width)
}

make_age_fair_mae_feval <- function(
    age_valid,
    bin_width = 5,                 
    min_n = 20,                  
    weight_scheme = c("uniform", "inverse_sqrt", "inverse_n")
) {
  weight_scheme <- match.arg(weight_scheme)
  
  if (bin_width == 1) {
    age_bin_valid <- age_valid
  } else {
    age_bin_valid <- floor(age_valid / bin_width) * bin_width
  }
  
  feval_fun <- function(preds, dtrain) {
    y_true <- getinfo(dtrain, "label")
    
    stopifnot(length(y_true) == length(age_bin_valid))
    
    df <- data.frame(
      age_bin = age_bin_valid,
      y_true  = y_true,
      y_pred  = preds
    )
    
    df_sum <- df %>%
      group_by(age_bin) %>%
      summarise(
        n    = n(),
        mae  = mean(abs(y_pred - y_true)),
        .groups = "drop"
      ) %>%
      filter(n >= min_n)
    
    # 如果所有 bin 都被过滤（极端情况），退回整体 MAE
    if (nrow(df_sum) == 0) {
      fair_mae <- mean(abs(preds - y_true))
      return(list(
        metric = "age_fair_mae",
        value  = fair_mae
      ))
    }
    
    if (weight_scheme == "uniform") {
      w <- rep(1 / nrow(df_sum), nrow(df_sum))
    } else if (weight_scheme == "inverse_sqrt") {
      w_raw <- 1 / sqrt(df_sum$n)
      w <- w_raw / sum(w_raw)
    } else if (weight_scheme == "inverse_n") {
      w_raw <- 1 / df_sum$n
      w <- w_raw / sum(w_raw)
    }
    
    fair_mae <- sum(w * df_sum$mae)
    
    return(list(
      metric = "age_fair_mae",
      value  = fair_mae
    ))
  }
  
  return(feval_fun)
}


make_age_fair_rmse_feval <- function(
    age_valid,
    bin_width = 5,                 
    min_n = 20,                  
    weight_scheme = c("uniform", "inverse_sqrt", "inverse_n")
) {
  weight_scheme <- match.arg(weight_scheme)
  
  if (bin_width == 1) {
    age_bin_valid <- age_valid
  } else {
    age_bin_valid <- floor(age_valid / bin_width) * bin_width
  }
  
  feval_fun <- function(preds, dtrain) {
    y_true <- getinfo(dtrain, "label")
    
    stopifnot(length(y_true) == length(age_bin_valid))
    
    df <- data.frame(
      age_bin = age_bin_valid,
      y_true  = y_true,
      y_pred  = preds
    )
    
    df_sum <- df %>%
      group_by(age_bin) %>%
      summarise(
        n    = n(),
        rmse = sqrt(mean((y_pred - y_true)^2)),
        .groups = "drop"
      ) %>%
      filter(n >= min_n)
    
    if (nrow(df_sum) == 0) {
      fair_rmse <- sqrt(mean((preds - y_true)^2))
      return(list(
        metric = "age_fair_rmse",
        value  = fair_rmse
      ))
    }
    
    if (weight_scheme == "uniform") {
      w <- rep(1 / nrow(df_sum), nrow(df_sum))
    } else if (weight_scheme == "inverse_sqrt") {
      w_raw <- 1 / sqrt(df_sum$n)
      w <- w_raw / sum(w_raw)
    } else if (weight_scheme == "inverse_n") {
      w_raw <- 1 / df_sum$n
      w <- w_raw / sum(w_raw)
    }
    
    fair_rmse <- sum(w * df_sum$rmse)
    
    return(list(
      metric = "age_fair_rmse",
      value  = fair_rmse
    ))
  }
  
  return(feval_fun)
}

train_xgb_group_dro <- function(
    X_train,
    y_train,
    group_train,
    X_valid = NULL,
    y_valid = NULL,
    params,
    nrounds_inner = 50,
    T_outer = 10,
    eta_q = 0.1,
    verbose = 1,
    # ===== 新增：outer 选择标准 =====
    eval_fun = c("rmse", "mae", "mse", "custom"),
    custom_eval_fun = NULL,   # function(pred_valid, y_valid, model, dvalid) -> numeric
    minimize = TRUE,          # TRUE 表示越小越好（rmse/mae/mse）；如果你自定义是越大越好可设 FALSE
    # ===== 新增：数值稳定 =====
    eps_q = 1e-3,             # group 权重下限
    center_scale_loss = TRUE  # 是否对 group loss 做 (x-mean)/sd
) {
  eval_fun <- match.arg(eval_fun)
  
  n <- length(y_train)
  if (!is.factor(group_train)) group_train <- factor(group_train)
  if (length(group_train) != n) stop("group_train length must match y_train length")
  
  if (eval_fun == "custom" && is.null(custom_eval_fun)) {
    stop("eval_fun='custom' requires custom_eval_fun.")
  }
  
  # 组信息
  G_levels <- levels(group_train)
  G <- length(G_levels)
  
  # 初始 group 权重：均匀
  q <- rep(1 / G, G)
  names(q) <- G_levels
  
  # 初始样本权重：按 group 权重展开
  w <- q[group_train]
  
  dtrain <- xgb.DMatrix(data = X_train, label = y_train, weight = w)
  dvalid <- NULL
  watchlist <- NULL
  
  if (!is.null(X_valid) && !is.null(y_valid)) {
    dvalid <- xgb.DMatrix(data = X_valid, label = y_valid)
    watchlist <- list(train = dtrain, eval = dvalid)
  } else {
    watchlist <- list(train = dtrain)
  }
  
  # 历史记录
  history <- list(
    q_list = vector("list", T_outer),
    group_loss_list = vector("list", T_outer),
    valid_score = rep(NA_real_, T_outer)
  )
  
  bst <- NULL
  
  # ===== 新增：best outer 追踪 =====
  best_outer <- NA_integer_
  best_score <- if (minimize) Inf else -Inf
  best_model <- NULL
  
  # 计算 valid score 的工具函数
  score_valid <- function(model) {
    if (is.null(dvalid)) return(NA_real_)
    pv <- predict(model, dvalid)
    
    if (eval_fun == "rmse") {
      return(sqrt(mean((pv - y_valid)^2)))
    } else if (eval_fun == "mse") {
      return(mean((pv - y_valid)^2))
    } else if (eval_fun == "mae") {
      return(mean(abs(pv - y_valid)))
    } else {
      # custom: 允许你用模型/矩阵做复杂指标
      return(custom_eval_fun(pv, y_valid, model, dvalid))
    }
  }
  
  for (t in 1:T_outer) {
    if (verbose) cat("\n===== Group DRO outer iteration", t, "=====\n")
    
    # 1) 用当前样本权重训练 / 继续训练 XGBoost
    if (is.null(bst)) {
      bst <- xgb.train(
        params = params,
        data   = dtrain,
        nrounds = nrounds_inner,
        watchlist = watchlist,
        verbose = ifelse(verbose >= 2, 1, 0)
      )
    } else {
      bst <- xgb.train(
        params = params,
        data   = dtrain,
        nrounds = nrounds_inner,
        xgb_model = bst,
        watchlist = watchlist,
        verbose = ifelse(verbose >= 2, 1, 0)
      )
    }
    
    # 1.5) ===== 新增：outer-level valid 评估，并保存 best outer =====
    sc <- score_valid(bst)
    history$valid_score[t] <- sc
    
    if (!is.na(sc)) {
      is_better <- if (minimize) (sc < best_score) else (sc > best_score)
      if (is_better) {
        best_score <- sc
        best_outer <- t
        best_model <- bst
      }
      if (verbose) cat("Outer", t, "valid score =", sc, "\n")
    }
    
    # 2) 在训练集上计算每个 group 的平均 loss（平方误差）
    preds_train <- predict(bst, dtrain)
    residuals <- preds_train - y_train
    loss_i <- residuals^2
    
    df_loss <- data.frame(
      group = group_train,
      loss = loss_i
    )
    
    group_loss <- df_loss %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        n = dplyr::n(),
        mean_loss = mean(loss),
        .groups = "drop"
      )
    
    if (verbose) {
      cat("Group-wise training loss:\n")
      print(group_loss)
    }
    
    history$q_list[[t]] <- q
    history$group_loss_list[[t]] <- group_loss
    
    # 3) exponentiated-gradient 更新 group 权重 q_g
    L_vec <- group_loss$mean_loss
    names(L_vec) <- as.character(group_loss$group)
    L_vec <- L_vec[names(q)]  # 对齐
    
    if (center_scale_loss) {
      L_centered <- L_vec - mean(L_vec)
      s <- sd(L_vec)
      if (is.na(s) || s < 1e-12) {
        L_scaled <- rep(0, length(L_vec))
      } else {
        L_scaled <- L_centered / s
      }
    } else {
      L_scaled <- L_vec
    }
    
    q_raw <- q * exp(eta_q * L_scaled)
    
    # 下限避免塌缩
    q_raw[q_raw < eps_q] <- eps_q
    q <- q_raw / sum(q_raw)
    
    if (verbose) {
      cat("Updated group weights q_g:\n")
      print(q)
    }
    
    # 4) 映射回样本权重并重建 dtrain / watchlist
    w <- q[group_train]
    dtrain <- xgb.DMatrix(data = X_train, label = y_train, weight = w)
    
    if (!is.null(dvalid)) {
      watchlist <- list(train = dtrain, eval = dvalid)
    } else {
      watchlist <- list(train = dtrain)
    }
  }
  
  out <- list(
    model = bst,                 # 最后一轮
    best_model = best_model,     # ✅ 新增：outer 最优模型
    best_outer = best_outer,
    best_score = best_score,
    group_weights = q,
    history = history,
    group_levels = G_levels,
    params = params,
    nrounds_total = nrounds_inner * T_outer,
    eval_fun = eval_fun,
    minimize = minimize
  )
  class(out) <- "xgb_group_dro"
  return(out)
}



make_id_delta_feval <- function(ids_valid, metric = c("rmse", "mae")) 
  {
  metric <- match.arg(metric)
  force(ids_valid)
  function(preds, dtrain) {
    y <- xgboost::getinfo(dtrain, "label")
    stopifnot(length(y) == length(preds), length(y) == length(ids_valid))
    
    df <- data.frame(id = ids_valid, y = y, pred = preds)
    
    df <- df[order(df$id, df$y), ]
    
    dy <- ave(df$y,    df$id, FUN = function(v) c(NA, diff(v)))
    dp <- ave(df$pred, df$id, FUN = function(v) c(NA, diff(v)))
    
    ok <- !is.na(dy) & !is.na(dp)
    err <- dp[ok] - dy[ok]  
    
    value <- if (metric == "rmse") sqrt(mean(err^2)) else mean(abs(err))
    
    list(metric = paste0("id_delta_", metric), value = value)
  }
}

# make_combined_feval <- function(
#     ids_valid,
#     age_valid,
#     bin_width = 5,
#     min_n = 10,
#     weight_scheme = c("uniform", "inverse_sqrt", "inverse_n"),
#     id_metric = c("rmse", "mae"),
#     alpha = 1,
#     beta  = 1,
#     gamma = 1
# ) {
#   weight_scheme <- match.arg(weight_scheme)
#   id_metric <- match.arg(id_metric)
# 
#   feval_id <- make_id_delta_feval(ids_valid = ids_valid, metric = id_metric)
# 
#   feval_agefair <- make_age_fair_rmse_feval(
#     age_valid = age_valid,
#     bin_width = bin_width,
#     min_n = min_n,
#     weight_scheme = weight_scheme
#   )
# 
#   function(preds, dtrain) {
#     y <- xgboost::getinfo(dtrain, "label")
# 
#     v1 <- feval_id(preds, dtrain)$value
#     v2 <- feval_agefair(preds, dtrain)$value
#     v3 <- sqrt(mean((preds - y)^2, na.rm = TRUE))
# 
#     value_c <- abs(alpha * v1 + beta * v2 + gamma * v3)
# 
#     print(paste0('id-rmse: ',v1))
#     print(paste0('age_fair-rmse: ',v2))
#     print(paste0('rmse: ',v3))
# 
#     list(
#       metric = paste0("combined_", id_metric, "_plus_agefair_rmse"),
#       value  = value_c
#     )
#   }
# }

# make_combined_feval <- function(
#     ids_valid,
#     age_valid,
#     bin_width = 5,
#     min_n = 10,
#     weight_scheme = c("uniform", "inverse_sqrt", "inverse_n"),
#     id_metric = c("rmse", "mae"),
#     alpha = 1,
#     beta  = 1,
#     gamma = 0,
#     baseline_rounds = 5,   # 用前多少轮做基准均值；设为1就是“第1轮”
#     eps = 1e-8             # 防止除0
# ) {
#   weight_scheme <- match.arg(weight_scheme)
#   id_metric <- match.arg(id_metric)
# 
#   feval_id <- make_id_delta_feval(ids_valid = ids_valid, metric = id_metric)
# 
#   feval_agefair <- make_age_fair_rmse_feval(
#     age_valid = age_valid,
#     bin_width = bin_width,
#     min_n = min_n,
#     weight_scheme = weight_scheme
#   )
# 
#   # ---- 闭包状态：用于累计前 baseline_rounds 轮的指标，计算基准 ----
#   t <- 0L
#   v1_hist <- numeric(0)
#   v2_hist <- numeric(0)
#   v3_hist <- numeric(0)
# 
#   v1_0 <- NA_real_
#   v2_0 <- NA_real_
#   v3_0 <- NA_real_
# 
#   function(preds, dtrain) {
#     y <- xgboost::getinfo(dtrain, "label")
# 
#     v1 <- feval_id(preds, dtrain)$value
#     v2 <- feval_agefair(preds, dtrain)$value
#     v3 <- sqrt(mean((preds - y)^2, na.rm = TRUE))
# 
#     # 更新轮次与历史
#     t <<- t + 1L
#     if (is.na(v1_0) || is.na(v2_0) || is.na(v3_0)) {
#       v1_hist <<- c(v1_hist, v1)
#       v2_hist <<- c(v2_hist, v2)
#       v3_hist <<- c(v3_hist, v3)
# 
#       if (t >= baseline_rounds) {
#         v1_0 <<- mean(v1_hist, na.rm = TRUE)
#         v2_0 <<- mean(v2_hist, na.rm = TRUE)
#         v3_0 <<- mean(v3_hist, na.rm = TRUE)
#       }
#     }
# 
#     # 如果基准还没建立（前 baseline_rounds 轮），先用当前值当临时基准
#     b1 <- if (!is.na(v1_0)) v1_0 else v1
#     b2 <- if (!is.na(v2_0)) v2_0 else v2
#     b3 <- if (!is.na(v3_0)) v3_0 else v3
# 
#     # 归一化（scale-free）
#     v1n <- v1 / max(b1, eps)
#     v2n <- v2 / max(b2, eps)
#     v3n <- v3 / max(b3, eps)
# 
#     value_c <- abs(alpha * v1n + beta * v2n + gamma * v3n)
# 
#     print(paste0("id-", id_metric, ": ", v1, "  (norm=", round(v1n, 4), ")"))
#     print(paste0("age_fair-rmse: ", v2, "  (norm=", round(v2n, 4), ")"))
#     print(paste0("rmse: ", v3, "  (norm=", round(v3n, 4), ")"))
#     if (!is.na(v1_0)) {
#       print(paste0("baseline (mean first ", baseline_rounds, " rounds): ",
#                    "v1_0=", round(v1_0, 4), ", v2_0=", round(v2_0, 4), ", v3_0=", round(v3_0, 4)))
#     }
# 
#     list(
#       metric = paste0("combined_norm_", id_metric, "_agefair_rmse_rmse_a", alpha, "_b", beta, "_g", gamma),
#       value  = value_c
#     )
#   }
# }


make_id_combo_delta_feval <- function(
    ids_valid,
    inst_valid,
    metric = c("rmse", "mae"),
    weight_scheme = c("uniform", "inverse_sqrt", "inverse_n"),
    allowed_instances = c(0, 2, 3),
    min_n_ids = 5,     # 某组合至少多少个 participant 才纳入（可按需调）
    verbose = FALSE
) {
  metric <- match.arg(metric)
  weight_scheme <- match.arg(weight_scheme)
  
  force(ids_valid)
  force(inst_valid)
  
  function(preds, dtrain) {
    y <- xgboost::getinfo(dtrain, "label")
    stopifnot(length(y) == length(preds),
              length(y) == length(ids_valid),
              length(y) == length(inst_valid))
    
    df <- data.frame(
      id = ids_valid,
      inst = inst_valid,
      y = y,
      pred = preds
    )
    
    # 只保留关心的 instance（0/2/3）
    df <- df[df$inst %in% allowed_instances, , drop = FALSE]
    if (nrow(df) == 0) {
      return(list(metric = paste0("id_combo_delta_", metric), value = NA_real_))
    }
    
    # ---- 1) 每个 (id, inst) 聚合成 1 条，防止重复测量 ----
    df_agg <- df %>%
      dplyr::group_by(id, inst) %>%
      dplyr::summarise(
        y = mean(y),
        pred = mean(pred),
        .groups = "drop"
      )
    
    # ---- 2) 给每个 id 定义组合类型（只考虑 0/2/3 的出现情况）----
    id_combo <- df_agg %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(
        I0 = any(inst == 0),
        I2 = any(inst == 2),
        I3 = any(inst == 3),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        combo = dplyr::case_when(
          I0 & I2 & !I3 ~ "0 & 2",
          I0 & !I2 & I3 ~ "0 & 3",
          !I0 & I2 & I3 ~ "2 & 3",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(combo))
    
    if (nrow(id_combo) == 0) {
      return(list(metric = paste0("id_combo_delta_", metric), value = NA_real_))
    }
    
    df_agg <- df_agg %>%
      dplyr::inner_join(id_combo[, c("id", "combo")], by = "id")
    
    # ---- 3) 在每个 id 内按 instance 排序，计算相邻 delta 误差 ----
    df_delta <- df_agg %>%
      dplyr::arrange(id, inst) %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(
        dy = y - dplyr::lag(y),
        dp = pred - dplyr::lag(pred),
        err = dp - dy
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(err))
    
    # 某些组合理论上每个 id 只有 1 个 delta（0->2 / 0->3 / 2->3）
    # 但若出现异常多点，仍然可以工作（会产生多条 err）。
    
    # ---- 4) 每个 combo 分别计算误差指标（按 id 聚合后再算，避免一个人多条 delta 过度影响）----
    # 这里按“每个 id 的 mean(err^2) / mean(abs(err))”先算，再在 combo 内平均，更符合“按人”比较
    per_id <- df_delta %>%
      dplyr::group_by(combo, id) %>%
      dplyr::summarise(
        id_score = if (metric == "rmse") sqrt(mean(err^2)) else mean(abs(err)),
        .groups = "drop"
      )
    
    combo_sum <- per_id %>%
      dplyr::group_by(combo) %>%
      dplyr::summarise(
        n_ids = dplyr::n(),
        score = mean(id_score),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_ids >= min_n_ids)
    
    if (nrow(combo_sum) == 0) {
      # 退回：不分 combo 的整体（按人）
      overall <- per_id$id_score
      value <- if (length(overall) == 0) NA_real_ else mean(overall)
      return(list(metric = paste0("id_combo_delta_", metric), value = value))
    }
    
    # ---- 5) 对 combo 之间加权汇总 ----
    if (weight_scheme == "uniform") {
      w <- rep(1 / nrow(combo_sum), nrow(combo_sum))
    } else if (weight_scheme == "inverse_sqrt") {
      w_raw <- 1 / sqrt(combo_sum$n_ids)
      w <- w_raw / sum(w_raw)
    } else if (weight_scheme == "inverse_n") {
      w_raw <- 1 / combo_sum$n_ids
      w <- w_raw / sum(w_raw)
    }
    
    value <- sum(w * combo_sum$score)
    
    if (isTRUE(verbose)) {
      print(combo_sum)
      print(data.frame(combo = combo_sum$combo, n_ids = combo_sum$n_ids, w = w))
    }
    
    list(
      metric = paste0("id_combo_delta_", metric, "_w_", weight_scheme),
      value = value
    )
  }
}


make_combined_feval <- function(
    ids_valid,
    inst_valid,                      # ✅ 新增：验证集 instance（0/2/3）
    age_valid,                       # 仍保留：用于 age_fair_rmse
    bin_width = 5,
    min_n = 10,
    weight_scheme = c("uniform", "inverse_sqrt", "inverse_n"),
    id_metric = c("rmse", "mae"),
    alpha = 1,
    beta  = 1,
    gamma = 0,
    baseline_rounds = 5,
    eps = 1e-8,
    n_valid = length(ids_valid),     # ✅ 用于识别 eval
    min_n_ids_combo = 5,             # ✅ combo 至少多少个 id 才纳入
    verbose = FALSE
) {
  weight_scheme <- match.arg(weight_scheme)
  id_metric <- match.arg(id_metric)
  
  # ✅ 换成新的：按组合(0&2 / 0&3 / 2&3)算 id-delta，再按组合加权
  feval_id <- make_id_combo_delta_feval(
    ids_valid  = ids_valid,
    inst_valid = inst_valid,
    metric = id_metric,
    weight_scheme = weight_scheme,
    min_n_ids = min_n_ids_combo
  )
  
  feval_agefair <- make_age_fair_rmse_feval(
    age_valid = age_valid,
    bin_width = bin_width,
    min_n = min_n,
    weight_scheme = weight_scheme
  )
  
  # ---- 闭包状态：仅对 eval 轮次累计 baseline ----
  t <- 0L
  v1_hist <- numeric(0)
  v2_hist <- numeric(0)
  v3_hist <- numeric(0)
  
  v1_0 <- NA_real_
  v2_0 <- NA_real_
  v3_0 <- NA_real_
  
  function(preds, dtrain) {
    
    # ✅ 如果不是验证集（比如 watchlist 里的 train），直接跳过，不影响 early stopping
    if (length(preds) != n_valid) {
      return(list(
        metric = paste0("combined_norm_", id_metric, "_idcombo_agefair_rmse_rmse_a", alpha, "_b", beta, "_g", gamma),
        value  = NA_real_
      ))
    }
    
    y <- xgboost::getinfo(dtrain, "label")
    
    v1 <- feval_id(preds, dtrain)$value
    v2 <- feval_agefair(preds, dtrain)$value
    v3 <- sqrt(mean((preds - y)^2, na.rm = TRUE))
    
    # ✅ 只对 eval 计数，baseline_rounds 不会被 train+eval 双调用搞乱
    t <<- t + 1L
    if (is.na(v1_0) || is.na(v2_0) || is.na(v3_0)) {
      v1_hist <<- c(v1_hist, v1)
      v2_hist <<- c(v2_hist, v2)
      v3_hist <<- c(v3_hist, v3)
      
      if (t >= baseline_rounds) {
        v1_0 <<- mean(v1_hist, na.rm = TRUE)
        v2_0 <<- mean(v2_hist, na.rm = TRUE)
        v3_0 <<- mean(v3_hist, na.rm = TRUE)
      }
    }
    
    b1 <- if (!is.na(v1_0)) v1_0 else v1
    b2 <- if (!is.na(v2_0)) v2_0 else v2
    b3 <- if (!is.na(v3_0)) v3_0 else v3
    
    v1n <- v1 / max(b1, eps)
    v2n <- v2 / max(b2, eps)
    v3n <- v3 / max(b3, eps)
    
    value_c <- abs(alpha * v1n + beta * v2n + gamma * v3n)
    
    if (isTRUE(verbose)) {
      print(paste0("idcombo-", id_metric, ": ", v1, "  (norm=", round(v1n, 4), ")"))
      print(paste0("age_fair-rmse: ", v2, "  (norm=", round(v2n, 4), ")"))
      print(paste0("rmse: ", v3, "  (norm=", round(v3n, 4), ")"))
      if (!is.na(v1_0)) {
        print(paste0(
          "baseline (mean first ", baseline_rounds, " eval rounds): ",
          "v1_0=", round(v1_0, 4), ", v2_0=", round(v2_0, 4), ", v3_0=", round(v3_0, 4)
        ))
      }
    }
    
    list(
      metric = paste0("combined_norm_", id_metric, "_idcombo_agefair_rmse_rmse_a", alpha, "_b", beta, "_g", gamma),
      value  = value_c
    )
  }
}
