select_columns_by_keyword <- function(df, keywords) {
  pattern <- paste(keywords, collapse = "|")
  selected <- df[, grepl(pattern, names(df)), drop = FALSE]
  
  return(selected)
}


process_longitudinal_right_left_data <- function(df, 
                                         measure_prefix, 
                                         age_prefix = "Age.when.attended.assessment.centre", 
                                         id_var = "Participant.ID", 
                                         sex_var = "Sex",
                                         side = NULL) {
  ## side: right? left? null?
  library(dplyr)
  library(tidyr)
  
  measure_pattern <- paste0(measure_prefix, ".*Instance.(\\d+)")
  age_pattern <- paste0(age_prefix, "...Instance.(\\d+)")
  
  long_data <- df %>%
    select(all_of(c(id_var, sex_var)), starts_with(measure_prefix), starts_with(age_prefix)) %>%
    pivot_longer(
      cols = starts_with(measure_prefix),
      names_to = "instance",
      names_pattern = measure_pattern,
      values_to = "value"
    ) %>%
    pivot_longer(
      cols = starts_with(age_prefix),
      names_to = "age_instance",
      names_pattern = age_pattern,
      values_to = "age"
    ) %>%
    filter(!is.na(value) & !is.na(age) & instance == age_instance) %>%
    mutate(instance = as.numeric(instance),
           age = as.numeric(age),
           Sex = factor(.data[[sex_var]], labels = c("Female","Male"))) %>%
    select(all_of(id_var), Sex, instance, age, value)
  
  if (!is.null(side)) {
    long_data <- long_data %>%
      mutate(side = side)
  }
  
  return(long_data)
}


plot_measure_vs_age <- function(data_long, 
                                measure_name = "Measure", 
                                side_var = "side") {
  library(ggplot2)
  library(patchwork)  
  
  if (side_var %in% names(data_long)) {
    aes_mapping <- aes(x = age, y = value, color = .data[[side_var]])
    group_var <- side_var
  } else {
    aes_mapping <- aes(x = age, y = value)
    group_var <- NULL
  }
  
  p_female <- ggplot(data_long[data_long$Sex == "Female", ], aes_mapping) +
    geom_point(alpha = 0.5) +
    stat_summary(fun = mean, geom = "line", size = 1.2, 
                 aes(group = if (!is.null(group_var)) .data[[group_var]] else 1)) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    labs(title = paste("Female:", measure_name, "vs Age"),
         x = "Age", y = measure_name) +
    theme_minimal()
  
  p_male <- ggplot(data_long[data_long$Sex == "Male", ], aes_mapping) +
    geom_point(alpha = 0.5) +
    stat_summary(fun = mean, geom = "line", size = 1.2, 
                 aes(group = if (!is.null(group_var)) .data[[group_var]] else 1)) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    labs(title = paste("Male:", measure_name, "vs Age"),
         x = "Age", y = measure_name) +
    theme_minimal()
  
  return(p_female + p_male)
}


merge_array_columns <- function(df,
                                array_pattern = "\\.\\.\\.Array\\.[0-9]+$",
                                merge_fun = function(x) {
                                  if (all(is.na(x))) return(NA)
                                  return(mean(x, na.rm = TRUE))
                                }) {
  
  base_names <- unique(sub(array_pattern, "", colnames(df)))
  
  merged_df <- data.frame(matrix(nrow = nrow(df), ncol = length(base_names)))
  colnames(merged_df) <- base_names
  
  for (bn in base_names) {
    target_cols <- grep(paste0("^", bn), colnames(df), value = TRUE)
    
    if (length(target_cols) == 1) {
      merged_df[[bn]] <- df[[target_cols]]
    } else {
      merged_df[[bn]] <- apply(df[, target_cols, drop = FALSE], 1, merge_fun)
    }
  }
  
  return(merged_df)
}


check_array_consistency <- function(df,
                                    array_pattern = "\\.\\.\\.Array\\.[0-9]+$") {
  # Extract the base column names (remove Array suffix)
  base_names <- unique(sub(array_pattern, "", colnames(df)))
  
  # Create result dataframe
  result_df <- data.frame(matrix(nrow = nrow(df), ncol = length(base_names)))
  colnames(result_df) <- paste0(base_names, "_sameMethod")
  
  # Loop through each base name group
  for (bn in base_names) {
    
    # Find all columns belonging to this base name (Array columns)
    target_cols <- grep(paste0("^", bn), colnames(df), value = TRUE)
    
    # If only one Array column exists → Always consistent (unless NA)
    if (length(target_cols) == 1) {
      result_df[[paste0(bn, "_sameMethod")]] <- ifelse(is.na(df[[target_cols]]), NA, 1)
      
    } else {
      # If more than one Array column exists → Compare values row-wise
      result_df[[paste0(bn, "_sameMethod")]] <- apply(df[, target_cols, drop = FALSE], 1, function(x) {
        
        # Case 1: Both values are NA → No measurement → return NA
        if (all(is.na(x))) return(NA)
        
        # Case 2: Only one measurement exists → Consider consistent → return 1
        if (sum(!is.na(x)) == 1) return(1)
        
        # Case 3: Both have values → Check if equal → 1 if equal, 0 if different
        return(ifelse(length(unique(x)) == 1, 1, 0))
      })
    }
  }
  
  return(result_df)
}


prepare_instance_data <- function(data, instance_idx) {
  y_col <- age_cols[instance_idx]
  y <- data[[y_col]]
  
  instance_cols <- colnames(data)[grepl(paste0("\\.Instance\\.", instance_idx, "$"), colnames(data))]
  
  X_cols <- c(instance_cols, general_cols)
  X <- data[, X_cols]
  
  X[] <- lapply(X, function(col) if(is.factor(col) || is.character(col)) as.numeric(as.factor(col)) else col)
  X <- as.matrix(X)
  
  return(list(X=X, y=y))
}





