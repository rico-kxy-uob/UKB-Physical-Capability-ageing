suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
})

# ============================================================
# ICD10 (Instance 0) -> two-level cause mapping (fixed vectorization)
# ============================================================

# Normalize an ICD10 code to its 3-character block (Letter + 2 digits)
.normalize_to_block3 <- function(code) {
  code <- toupper(str_trim(as.character(code)))
  code <- str_replace_all(code, "\\.", "")
  code <- str_replace_all(code, "[^A-Z0-9]", "")
  ifelse(str_detect(code, "^[A-Z][0-9]{2}"), substr(code, 1, 3), NA_character_)
}

# Vectorized: convert 3-char ICD10 block to integer key
# key = (letter_index * 100) + two_digit_number
.block3_to_key <- function(block3) {
  block3 <- toupper(as.character(block3))
  letter <- substr(block3, 1, 1)
  num2   <- suppressWarnings(as.integer(substr(block3, 2, 3)))
  
  letter_idx <- match(letter, LETTERS) - 1L   # A->0, B->1, ...
  ok <- !is.na(letter_idx) & !is.na(num2)
  
  key <- letter_idx * 100L + num2
  ifelse(ok, key, NA_integer_)
}

# Build a keyed range map: label + range ("C15-C26") -> start/end keys
.build_keyed_map <- function(df, label_col, range_col) {
  df %>%
    mutate(.range = .data[[range_col]]) %>%
    separate(.range, into = c("start", "end"), sep = "-", remove = FALSE) %>%
    mutate(
      start_key = .block3_to_key(start),
      end_key   = .block3_to_key(end),
      .label    = .data[[label_col]]
    ) %>%
    select(.label, start, end, start_key, end_key)
}

# Map codes to a label using a keyed range map (first hit if overlap)
.map_codes_by_range <- function(codes, keyed_map) {
  blk <- .normalize_to_block3(codes)
  k   <- .block3_to_key(blk)
  
  out <- rep(NA_character_, length(codes))
  ok  <- !is.na(k)
  
  if (any(ok)) {
    kk <- k[ok]
    out_ok <- vapply(
      kk,
      function(x) {
        hit <- keyed_map[keyed_map$start_key <= x & x <= keyed_map$end_key, , drop = FALSE]
        if (nrow(hit) == 0) NA_character_ else hit$.label[[1]]
      },
      character(1)
    )
    out[ok] <- out_ok
  }
  out
}

# ============================================================
# Main function
# ============================================================
add_icd10_cause_levels <- function(
    df,
    icd10_col = "death_cause",
    other_label = "Other",
    keep_chapters = c(
      "Chapter II Neoplasms",
      "Chapter V Mental and behavioural disorders",
      "Chapter VI Diseases of the nervous system",
      "Chapter IX Diseases of the circulatory system",
      "Chapter X Diseases of the respiratory system",
      "Chapter XI Diseases of the digestive system"
    ),
    chapter_map = NULL,
    neoplasm_detail_map = NULL
) {
  # Default chapter-level map (broad)
  if (is.null(chapter_map)) {
    chapter_map <- tibble::tribble(
      ~chapter, ~range,
      "Chapter II Neoplasms", "C00-D48",
      "Chapter V Mental and behavioural disorders", "F00-F99",
      "Chapter VI Diseases of the nervous system", "G00-G99",
      "Chapter IX Diseases of the circulatory system", "I00-I99",
      "Chapter X Diseases of the respiratory system", "J00-J99",
      "Chapter XI Diseases of the digestive system", "K00-K93"
    )
  }
  
  # Default neoplasm detailed map (your full list)
  if (is.null(neoplasm_detail_map)) {
    neoplasm_detail_map <- tibble::tribble(
      ~neo_category, ~range,
      "C00-C14 Malignant neoplasms of lip, oral cavity and pharynx", "C00-C14",
      "C15-C26 Malignant neoplasms of digestive organs", "C15-C26",
      "C30-C39 Malignant neoplasms of respiratory and intrathoracic organs", "C30-C39",
      "C40-C41 Malignant neoplasms of bone and articular cartilage", "C40-C41",
      "C43-C44 Melanoma and other malignant neoplasms of skin", "C43-C44",
      "C45-C49 Malignant neoplasms of mesothelial and soft tissue", "C45-C49",
      "C50-C50 Malignant neoplasm of breast", "C50-C50",
      "C51-C58 Malignant neoplasms of female genital organs", "C51-C58",
      "C60-C63 Malignant neoplasms of male genital organs", "C60-C63",
      "C64-C68 Malignant neoplasms of urinary tract", "C64-C68",
      "C69-C72 Malignant neoplasms of eye, brain and other parts of central nervous system", "C69-C72",
      "C73-C75 Malignant neoplasms of thyroid and other endocrine glands", "C73-C75",
      "C76-C80 Malignant neoplasms of ill-defined, secondary and unspecified sites", "C76-C80",
      "C81-C96 Malignant neoplasms of lymphoid, haematopoietic and related tissue", "C81-C96",
      "C97-C97 Malignant neoplasms of independent (primary) multiple sites", "C97-C97",
      "D10-D36 Benign neoplasms", "D10-D36",
      "D37-D48 Neoplasms of uncertain or unknown behaviour", "D37-D48"
    )
  }
  
  # Precompute keyed maps
  chapter_keyed <- .build_keyed_map(chapter_map, "chapter", "range")
  neo_keyed     <- .build_keyed_map(neoplasm_detail_map, "neo_category", "range")
  
  # Pull and clean codes
  codes <- df[[icd10_col]]
  codes <- na_if(str_trim(as.character(codes)), "")
  
  # Broad chapter mapping
  chapter_all <- .map_codes_by_range(codes, chapter_keyed)
  
  # Level 1: only keep your target chapters; otherwise "Other"
  cause_lv1 <- ifelse(!is.na(chapter_all) & chapter_all %in% keep_chapters, chapter_all, other_label)
  
  # Level 2:
  # - Neoplasms -> detailed subtype when available, else fallback text
  # - Non-neoplasms -> chapter
  # - Outside chapter map -> "Other"
  neo_detail <- .map_codes_by_range(codes, neo_keyed)
  
  cause_lv2 <- dplyr::case_when(
    chapter_all == "Chapter II Neoplasms" ~ dplyr::coalesce(
      neo_detail,
      "Chapter II Neoplasms (other/unspecified within C00-D48)"
    ),
    !is.na(chapter_all) ~ chapter_all,
    TRUE ~ other_label
  )
  
  df %>%
    mutate(
      !!icd10_col := codes,
      cause_lv1 = cause_lv1,
      cause_lv2 = cause_lv2
    )
}


# 2. Define comparison calculation function
get_comparison_stats <- function(df, sex_include, cause_code, time_threshold) {
  
  # Landmark filtering
  sub_df <- df %>% filter(follow_up_time > time_threshold,Sex %in% sex_include)
  
  # Define Outcome
  if (cause_code == "All-Cause") {
    sub_df$outcome_event <- sub_df$event
  } else {
    sub_df$outcome_event <- ifelse(sub_df$event == 1 & sub_df$cause_lv1 == cause_code, 1, 0)
  }
  
  n_events <- sum(sub_df$outcome_event)
  if (n_events < 10) return(NULL) 
  
  fit_base <- coxph(Surv(follow_up_time, outcome_event) ~ age_c * Sex, data = sub_df)
  c_base <- summary(fit_base)$concordance[1]
  
  fit_full <- coxph(Surv(follow_up_time, outcome_event) ~ age_c * Sex + age_c * adj_BAA + adj_BAA * Sex, data = sub_df)
  c_full <- summary(fit_full)$concordance[1]
  
  delta_c <- c_full - c_base
  
  lrt_tab <- anova(fit_base, fit_full, test = "LRT")
  p_lrt <- lrt_tab$`Pr(>|Chi|)`[2]
  # p_lrt <- unname(p_lrt)
  print(p_lrt)
  
  return(c(
    Events = n_events,
    C_Base = c_base,
    C_Full = c_full,
    Delta  = delta_c,
    P_LRT  = p_lrt
  ))
}


fmt_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  out <- ifelse(is.na(p), NA_character_,
                ifelse(p < 1e-3, "<0.001",
                       format(round(p, 3), nsmall = 3)))
  out
}


get_comparison_stats_lv2 <- function(df, sex_include, cause_code, time_threshold) {
  
  # Landmark filtering
  sub_df <- df %>% filter(follow_up_time > time_threshold,Sex %in% sex_include)
  
  # Define Outcome
  if (cause_code == "All-Cause") {
    sub_df$outcome_event <- sub_df$event
  } else {
    sub_df$outcome_event <- ifelse(sub_df$event == 1 & sub_df$cause_lv2 == cause_code, 1, 0)
  }
  
  n_events <- sum(sub_df$outcome_event)
  if (n_events < 10) return(NULL) 
  
  fit_base <- coxph(Surv(follow_up_time, outcome_event) ~ age_c * Sex, data = sub_df)
  c_base <- summary(fit_base)$concordance[1]
  
  fit_full <- coxph(Surv(follow_up_time, outcome_event) ~ age_c * Sex + age_c * adj_BAA + adj_BAA * Sex, data = sub_df)
  c_full <- summary(fit_full)$concordance[1]
  
  delta_c <- c_full - c_base
  
  lrt_tab <- anova(fit_base, fit_full, test = "LRT")
  p_lrt <- lrt_tab$`Pr(>|Chi|)`[2]
  # p_lrt <- unname(p_lrt)
  
  return(c(
    Events = n_events,
    C_Base = c_base,
    C_Full = c_full,
    Delta  = delta_c,
    P_LRT  = p_lrt
  ))
}

get_comparison_stats_landmark <- function(df, sex_include, cause_code, time_threshold) {
  
  # Landmark: keep those at risk at tau, then reset time origin to tau
  sub_df <- df %>%
    dplyr::filter(follow_up_time > time_threshold, Sex %in% sex_include) %>%
    dplyr::mutate(t_land = follow_up_time - time_threshold)
  
  # Define outcome (cause-specific: competing deaths treated as censoring at their event time)
  if (cause_code == "All-Cause") {
    sub_df$outcome_event <- as.integer(sub_df$event)
  } else {
    sub_df$outcome_event <- as.integer(sub_df$event == 1 & sub_df$cause_lv1 == cause_code)
  }
  
  n_events <- sum(sub_df$outcome_event, na.rm = TRUE)
  if (n_events < 10) return(NULL)
  
  # Base vs full models (landmark time)
  fit_base <- survival::coxph(survival::Surv(t_land, outcome_event) ~ age_c * Sex, data = sub_df)
  c_base <- summary(fit_base)$concordance[1]
  
  fit_full <- survival::coxph(
    survival::Surv(t_land, outcome_event) ~ age_c * Sex + age_c * adj_BAA + adj_BAA * Sex,
    data = sub_df
  )
  c_full <- summary(fit_full)$concordance[1]
  
  delta_c <- c_full - c_base
  
  lrt_tab <- stats::anova(fit_base, fit_full, test = "LRT")
  p_lrt <- lrt_tab$`Pr(>|Chi|)`[2]
  
  return(c(
    Events = n_events,
    C_Base = c_base,
    C_Full = c_full,
    Delta  = delta_c,
    P_LRT  = p_lrt
  ))
}

