library(ggplot2)
library(xgboost)
library(Matrix)
library(SHAPforxgboost)
library(dplyr)
library(UpSetR)
library(tidyr)
library(fastDummies)
library(pheatmap)
library(lme4)
library(shapviz)
library(survival)
library(survminer)
library(lubridate)
library(ggsci)
library(readr)
library(survey)

setwd("/rds/projects/s/spillf-systems-mechanobiology-health-disease/RDS-Rico/R/Physical mea ageing clock")
# setwd('C:/Users/Kaiyi/OneDrive - University of Birmingham/Desktop/PhD/PhD/dataset_physical_mea')

bundle_list <- readRDS("aft_clock_sex_specific_bundles_ins023_Gompertz_CV.rds")
names(bundle_list)


################## External TEST ###################
# new_person_df: External new data (at least containing Sex, age0 and as many columns as possible from feature_cols)
# testing_3plusmeas_raw<-read.csv('deep_cleaned_3plusmeas_testing_ins0123_lfm.csv',header=TRUE)
training_2plusmeas_raw<-read.csv('deep_cleaned_2plusmeas_ins023_lfm.csv',header=TRUE)
dim(training_2plusmeas_raw)

test_data_nacols <- training_2plusmeas_raw %>% 
  select(-all_of(colnames(confounders_extra))) %>% 
  filter(Instance %in% c(0,2,3))
colnames(test_data_nacols)

na_percent <- sapply(test_data_nacols, function(x) {
  mean(is.na(x)) * 100
})

# Print NA percentage in descending order
na_percent_sorted <- sort(na_percent, decreasing = TRUE)
print(na_percent_sorted)

# Identify columns with >60% NA
cols_to_drop <- names(na_percent[na_percent > 40])
cat("Columns to drop (>60% NA):\n")
print(cols_to_drop)

test_data_nona_cols<-test_data_nacols[,!colnames(test_data_nacols) %in% cols_to_drop]
dim(test_data_nona_cols)
colnames(test_data_nona_cols)
table(table(test_data_nona_cols$Participant.ID))

test_data <- test_data_nona_cols %>%
  filter(rowMeans(is.na(across(all_of(feature_cols)))) <= 0.20)

careful_ids<- intersect(dat_aft$Participant.ID,test_data$Participant.ID)
test_data <- subset (test_data, !Participant.ID %in%careful_ids)

dim(test_data)
head(test_data)
length(unique(test_data$Participant.ID))
table(table(test_data$Participant.ID))



test_data <- test_data %>%
  mutate(age0 = as.numeric(`Age.when.attended.assessment.centre`))

test_data <- test_data %>%
  dplyr::group_by(Participant.ID) %>%
  dplyr::filter(dplyr::n() >= 2) %>%
  dplyr::ungroup()

dim(test_data)
table(table(test_data$Participant.ID))

test_data_all<-test_data[!test_data$Participant.ID %in% careful_ids,]
dim(test_data_all)

# out_all <- score_BAA_one(test_data_all, bundle_list)
out_all <- score_BAA_batch(test_data_all, bundle_list)
dim(test_data_all)
dim(out_all)

testing_3plusmeas_with_BAandBAA <- cbind(test_data,out_all)

test_data_for_survival<-testing_3plusmeas_with_BAandBAA[,c("Participant.ID" ,
                                                           # "Health.score..England.",
                                                           # "Living.environment.score..England.",
                                                           "Sex","Instance","Age.when.attended.assessment.centre" ,'BAA','BA')]

head(test_data_for_survival)
# test_data_for_survival<-subset(test_data_for_survival,Age.when.attended.assessment.centre<=75)

cor(test_data_for_survival$BAA,test_data_for_survival$Age.when.attended.assessment.centre)

testset_death_data<-read.csv('death_birth_info/test_set_death_cause.csv',header=TRUE)
length(testset_death_data$Date.of.death...Instance.0)

trainset_death_data<-read.csv('death_birth_info/train_set_death_cause.csv',header=TRUE)

length(!is.na(trainset_death_data$Date.of.death...Instance.0)==FALSE)

# testset_death_data <- subset(
#   testset_death_data,
#   !grepl("^[C]", Underlying..primary..cause.of.death..ICD10...Instance.0 %>% as.character())
# )

##################################################
df <- test_data_for_survival
df$BAA_raw <- df$BA - df$Age.when.attended.assessment.centre
df$adj_BAA <- NA_real_

for (s in sort(unique(df$Sex))) {
  idx <- which(df$Sex == s & is.finite(df$BAA_raw) & is.finite(df$Age.when.attended.assessment.centre))
  if (length(idx) < 20) next
  
  fit_q <- lm(BAA_raw ~ Age.when.attended.assessment.centre + I(Age.when.attended.assessment.centre^2), data = df, subset = idx)
  pred <- rep(NA_real_, nrow(df))
  pred[idx] <- predict(fit_q, newdata = df[idx, , drop = FALSE])
  df$adj_BAA[idx] <- df$BAA_raw[idx] - pred[idx]
}

test_data_for_survival <- df

cor.test(test_data_for_survival$Age.when.attended.assessment.centre,test_data_for_survival$adj_BAA)
#################################################

testset_death_data <- subset(
  trainset_death_data,
  !grepl("^[VWXY]", Underlying..primary..cause.of.death..ICD10...Instance.0 %>% as.character())
)
testset_death_data$death_summary_code = substr(testset_death_data$Underlying..primary..cause.of.death..ICD10...Instance.0,1,1)
head(testset_death_data)
dim(testset_death_data)

death_df <- testset_death_data %>%
  select(
    Participant.ID, 
    Year.of.birth, 
    Month.of.birth,
    starts_with("Date.of.death"),
    death_summary_code,
    Underlying..primary..cause.of.death..ICD10...Instance.0 # Added this here
  ) %>%
  pivot_longer(
    cols = starts_with("Date.of.death"),
    names_to = "Instance",
    values_to = "death_date_raw"
  ) %>%
  mutate(
    birth_date = make_date(Year.of.birth, Month.of.birth, 15),
    death_date = as.Date(death_date_raw),
    death_age_calc = as.numeric(difftime(death_date, birth_date, units = "days")) / 365.25
  ) %>%
  group_by(Participant.ID) %>%
  summarise(
    birth_date = first(birth_date),
    death_summary_code = first(death_summary_code), 
    death_age_recorded = min(death_age_calc, na.rm = TRUE),
    death_date_final = min(death_date, na.rm = TRUE),
    # Use first() here to avoid potential length errors during summary
    ICD10_cause = first(Underlying..primary..cause.of.death..ICD10...Instance.0),
    .groups = "drop"
  ) %>%
  mutate(
    death_age_recorded = ifelse(is.infinite(death_age_recorded), NA, death_age_recorded),
    death_date_final = as.Date(ifelse(is.infinite(as.numeric(death_date_final)), NA, death_date_final), origin = "1970-01-01"),
    
    censor_date = as.Date("2025-12-31"),
    censor_age = as.numeric(difftime(censor_date, birth_date, units = "days")) / 365.25,
    
    death_age = ifelse(is.na(death_age_recorded), censor_age, death_age_recorded),
    status = ifelse(is.na(death_age_recorded), 0, 1)
  )
table(death_df$status)

df_long <- test_data_for_survival %>%
  left_join(death_df, by = "Participant.ID")

head(df_long)
length(unique(df_long[df_long$status==1,'Participant.ID']))

source('Functions/ICD10-mapping.R')

death_df_unique <- death_df %>%
  group_by(Participant.ID) %>%
  slice(1) %>% # Takes the first record per person
  ungroup()

# 2. Join and Map
df_long <- test_data_for_survival %>%
  left_join(death_df_unique, by = "Participant.ID") %>%
  mutate(death_cause = as.character(ICD10_cause)) %>%
  add_icd10_cause_levels(icd10_col = "death_cause") %>%
  select(
    Participant.ID, 
    birth_date, 
    death_age, 
    status, 
    death_cause, 
    cause_lv1, 
    cause_lv2
  ) %>%
  distinct()
head(df_long)

table(df_long$status)


date_lookup <- testset_death_data %>%
  select(Participant.ID, starts_with("Date.of.attending.assessment.centre")) %>%
  pivot_longer(
    cols = starts_with("Date.of.attending.assessment.centre"),
    names_to = "Instance",
    values_to = "Date.of.attending"
  ) %>%
  # Extract the instance number (0, 1, 2, 3) from the column name string
  mutate(Instance = as.numeric(gsub(".*\\.", "", Instance))) %>%
  filter(!is.na(Date.of.attending) & Date.of.attending != "")

# 2. Process Time-Varying Data (df_tv) with date-to-decimal conversion
df_tv <- test_data_for_survival %>%
  # Join with the lookup table to get the actual date for each Instance row
  left_join(date_lookup, by = c("Participant.ID", "Instance")) %>%
  # Join with death_df to get birth_date and death_age
  left_join(df_long %>% select(Participant.ID, birth_date, death_age, status,death_cause,cause_lv1,cause_lv2), by = "Participant.ID") %>%
  # Ensure Date.of.attending is in Date format
  mutate(Date.of.attending = as.Date(Date.of.attending)) %>%
  # Calculate exact decimal age at each assessment
  mutate(
    Age.at.assessment = as.numeric(difftime(Date.of.attending, birth_date, units = "days")) / 365.25
  ) %>%
  # If Date.of.attending was missing, fall back to the integer age provided in the data
  mutate(Age.at.assessment = ifelse(is.na(Age.at.assessment), Age.when.attended.assessment.centre, Age.at.assessment)) %>%
  arrange(Participant.ID, Age.at.assessment) %>%
  group_by(Participant.ID) %>%
  mutate(
    BAA_interval = BAA,
    
    # Start of interval: Age at current assessment
    start = Age.at.assessment,
    # Stop of interval: Age at next assessment (lead)
    stop  = lead(Age.at.assessment),
    
    # If no next assessment, interval ends at death_age (death or censoring)
    stop = ifelse(is.na(stop), death_age, stop),
    
    # Ensure stop is always strictly greater than start for Cox model requirements
    stop = ifelse(stop <= start, start + 0.001, stop),
    
    # Event indicator: 1 if the person died and this is their final interval
    event = ifelse(status == 1 & stop == death_age, 1, 0)
  ) %>%
  filter(start < stop) %>% # Filter out any erroneous overlapping intervals
  ungroup()
head(df_tv)
table(df_tv$event)

combo_counts <- df_tv %>%
  group_by(Participant.ID) %>%
  summarise(
    I0 = as.integer(any(Instance == 0, na.rm = TRUE)),
    I2 = as.integer(any(Instance == 2, na.rm = TRUE)),
    I3 = as.integer(any(Instance == 3, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(I0 + I2 + I3 > 0) %>%
  mutate(
    combination = dplyr::case_when(
      I0 == 1 & I2 == 0 & I3 == 0 ~ "Only 0",
      I0 == 0 & I2 == 1 & I3 == 0 ~ "Only 2",
      I0 == 0 & I2 == 0 & I3 == 1 ~ "Only 3",
      I0 == 1 & I2 == 1 & I3 == 0 ~ "0 & 2",
      I0 == 1 & I2 == 0 & I3 == 1 ~ "0 & 3",
      I0 == 0 & I2 == 1 & I3 == 1 ~ "2 & 3",
      I0 == 1 & I2 == 1 & I3 == 1 ~ "0 & 2 & 3",
      TRUE ~ "(No 0/2/3)"
    )
  ) %>%
  count(combination, name = "n_participants") %>%
  arrange(desc(n_participants))

combo_counts

ids_0_2_3 <- df_tv %>%
  group_by(Participant.ID) %>%
  summarise(
    I0 = as.integer(any(Instance == 0, na.rm = TRUE)),
    I2 = as.integer(any(Instance == 2, na.rm = TRUE)),
    I3 = as.integer(any(Instance == 3, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(I0 == 1 & I2 == 0 & I3 == 0) %>%
  arrange(Participant.ID) %>%
  select(Participant.ID)
ids_0_2_3

df_tv_filter <- subset (df_tv, !Participant.ID %in% ids_0_2_3$Participant.ID)

table(df_tv_filter$event)
table(df_tv_filter$cause_lv1)
summary(df_tv_filter$Age.at.assessment)
length(unique(df_tv_filter$Participant.ID))
intersect(df_tv_filter$Participant.ID,dat_aft$Participant.ID)

cox_baa_sex <- coxph(
  Surv(start, stop, event) ~ adj_BAA * Sex 
  # +Sex*Health.score..England. + Sex*Living.environment.score..England.
  # +Health.score..England. + Living.environment.score..England.
  ,
  # data = df_tv[df_tv$Sex==0,],
  data = df_tv_filter,
  cluster = Participant.ID
)
summary(cox_baa_sex)

cox_sex_only <- coxph(
  Surv(start, stop, event) ~ Sex
  # +Sex*Health.score..England. + Sex*Living.environment.score..England.
  # +Health.score..England. + Living.environment.score..England.
  ,
  # data = df_tv[df_tv$Sex==0,],
  data = df_tv_filter,
  cluster = Participant.ID
)

summary(cox_sex_only)

# full model
cox_baa_sex <- coxph(Surv(start, stop, event) ~ Sex + BAA + Sex:BAA,
                     data = df_tv_filter, cluster = Participant.ID)

# joint Wald test for (BAA, BAA:Sex)
b <- coef(cox_baa_sex)[c("BAA", "Sex:BAA")]
V <- vcov(cox_baa_sex)[c("BAA", "Sex:BAA"), c("BAA", "Sex:BAA")]

W <- as.numeric(t(b) %*% solve(V) %*% b)        # chi-square stat
p <- pchisq(W, df = length(b), lower.tail = FALSE)

c(Wald_chisq = W, df = length(b), p_value = p)

df_late <- df_tv_filter %>% filter(stop > (start + 5))  # Adjust this based on your start/stop time units
cox_late <- coxph(Surv(start, stop, event) ~ BAA + Sex, data=df_late, cluster=Participant.ID)
summary(cox_late)

cox_late_sexonly <- coxph(Surv(start, stop, event) ~ Sex, data=df_late, cluster=Participant.ID)
summary(cox_late_sexonly)


valid_ids <- df_tv_filter %>%
  group_by(Participant.ID) %>%
  summarise(
    baseline_age = min(start, na.rm = TRUE),
    max_followup_age = max(stop, na.rm = TRUE),
    total_duration = max_followup_age - baseline_age
  )
dim(valid_ids)

head(valid_ids)

hist(valid_ids$total_duration)
summary(valid_ids$total_duration)

visit_counts <- df_tv_filter %>%
  group_by(Participant.ID) %>%
  summarise(
    n_visits = n_distinct(Instance),
    Sex = first(Sex),
    .groups = "drop"
  )
table(visit_counts$n_visits)

# Attach n_visits back to the long TV dataset
df_tv2 <- df_tv_filter %>%
  left_join(visit_counts[,c('Participant.ID','n_visits')], by = "Participant.ID") %>%
  filter(n_visits %in% c(2, 3))

head(df_tv2)

# write.csv( df_tv2, 'df_tv2.csv',row.names=FALSE)
###############################################################
cor.test(df_tv2$Age.when.attended.assessment.centre,df_tv2$BA)
cor.test(df_tv2$Age.when.attended.assessment.centre,df_tv2$adj_BAA)

library(ggsci)  
library(ggpubr)  

df_plot <- df_tv2 %>%
  mutate(Sex = factor(Sex, levels = c(0, 1), labels = c("Female", "Male"))) %>%
  group_by(Sex) %>%
  sample_n(10000) %>%
  ungroup()

### >70 extrapolation test ###？
longitudinal_BA_CA_corre_plot<-ggplot(df_plot[df_plot$Age.at.assessment>0,], aes(x =Age.at.assessment , y = BA, 
                    color = Sex, fill = Sex)) + 
  geom_point(alpha = 0.4, size = 2) + 
  geom_smooth(method = "lm", formula = y ~ x, 
              se = TRUE,       
              alpha = 0.3,    
              size = 1.2) +    
  stat_cor(method = "pearson", 
           label.x.npc = "left", 
           show.legend = FALSE) + 
  scale_color_nejm() +
  scale_fill_nejm() +
  theme_bw() +
  labs(
    x = "Chronological age",
    y = "Biological age",
    title = "Longitudinal test cohort"
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )
longitudinal_BA_CA_corre_plot

ggsave("Results plots/longitudinal_BA_CA_corre_plot.png", longitudinal_BA_CA_corre_plot, width = 9, height = 6, dpi = 600)

# head(df_tv2,20)
set.seed(1)
sample_ids_death <- sample(unique(df_tv2[df_tv2$n_visits>=3 &df_tv2$event==1,]$Participant.ID), 10)
df_sub_death <- df_tv2[df_tv2$Participant.ID %in% sample_ids_death, ]
df_sub_death$Participant.ID <- as.factor(df_sub_death$Participant.ID)
df_sub_death <- df_sub_death %>%
  mutate(ID_numeric = as.numeric(as.factor(Participant.ID))) %>%
  mutate(Paiticipant = sprintf("P_%02d", ID_numeric))

set.seed(42)
sample_ids <- sample(unique(df_tv2[df_tv2$n_visits>=3 &df_tv2$event==0,]$Participant.ID), 10)
df_sub <- df_tv2[df_tv2$Participant.ID %in% sample_ids, ]
df_sub$Participant.ID <- as.factor(df_sub$Participant.ID)
df_sub <- df_sub %>%
  mutate(ID_numeric = as.numeric(as.factor(Participant.ID))) %>%
  mutate(Paiticipant = sprintf("P_%02d", ID_numeric+10))

longitu_individ_age_traj_death <- ggplot(df_sub_death, aes(x = Age.at.assessment, y = adj_BAA, 
                   group = Participant.ID, 
                   color = Paiticipant)) +

  geom_abline(intercept = 0, slope = 0, color = "red", linetype = "dashed", size = 0.8) +
  geom_line(alpha = 0.8, size = 1) + 
  geom_point(aes(shape = as.factor(Sex)), size = 3) +
  scale_shape_manual(name = "Sex", 
                     values = c(16, 17), 
                     labels = c("Female", "Male")) +
  theme_bw() +
  theme(
    legend.position = "right",
    guides(color = "none")
  ) +
  labs(
    title = "Individual Aging Trajectories (10 Participants with death record)",
    x = "Chronological Age", 
    y = "Biological Age Acceleration(BAA))",
    shape = "Sex"
  ) +
  coord_fixed(ratio = 1) + xlim(c(40, 80))
longitu_individ_age_traj_death 
ggsave("Results plots/longitu_individ_age_traj_death.png", longitu_individ_age_traj_death, width = 9, height = 6, dpi = 600)

longitu_individ_age_traj <- ggplot(df_sub, aes(x = Age.at.assessment, y = adj_BAA, 
                         group = Participant.ID, 
                         color = Paiticipant)) +
  
  geom_abline(intercept = 0, slope = 0, color = "red", linetype = "dashed", size = 0.8) +
  geom_line(alpha = 0.8, size = 1) + 
  geom_point(aes(shape = as.factor(Sex)), size = 3) +
  scale_shape_manual(name = "Sex", 
                     values = c(16, 17), 
                     labels = c("Female", "Male")) +
  theme_bw() +
  theme(
    legend.position = "right",
    guides(color = "none") 
  ) +
  labs(
    title = "Individual Aging Trajectories (10 Participants still alive)",
    x = "Chronological Age", 
    y = "Biological Age Acceleration(BAA))",
    shape = "Sex"
  ) +
  coord_fixed(ratio = 1)
longitu_individ_age_traj
ggsave("Results plots/longitu_individ_age_traj.png", longitu_individ_age_traj, width = 9, height = 6, dpi = 600)

library(survival)
library(survminer)
library(ggsci)

df_tv2$Sex <- factor(df_tv2$Sex, levels = c(0, 1), labels = c("Female", "Male"))
stats_by_sex <- df_tv2 %>%
  group_by(Sex) %>%
  summarise(
    m_baa = mean(adj_BAA, na.rm = TRUE),
    sd_baa = sd(adj_BAA, na.rm = TRUE)
  )

f_mean <- stats_by_sex$m_baa[stats_by_sex$Sex == "Female"]
f_sd   <- stats_by_sex$sd_baa[stats_by_sex$Sex == "Female"]

m_mean <- stats_by_sex$m_baa[stats_by_sex$Sex == "Male"]
m_sd   <- stats_by_sex$sd_baa[stats_by_sex$Sex == "Male"]

new_df_group <- data.frame(
  Sex = factor(c("Female", "Female", "Male", "Male"), levels = c("Female", "Male")),
  adj_BAA = c(f_mean - f_sd, f_mean + f_sd, 
              m_mean - m_sd, m_mean + m_sd)   
)

fit_group <- survfit(cox_tv, newdata = new_df_group)

longitu_surv_curve <- ggsurvplot(
  fit_group, 
  data = new_df_group,
  legend.labs = c("Female: Low BAA", "Female: High BAA", 
                  "Male: Low BAA", "Male: High BAA"),
  palette = c("#FB9A99", "#E41A1C", "#A6CEE3", "#1F78B4"), 
  conf.int = TRUE, 
  xlim = c(45, 82),
  ylim = c(0.80, 1),
  linetype = c("solid", "solid", "dashed", "dashed"), 
  title = "Time-varying Cox based Survival curves by Sex-Specific BAA Levels",
  ggtheme = theme_bw()
)

# longitu_surv_curve$plot <- longitu_surv_curve$plot + 
#   labs(color = "Group", linetype = "Group") 

print(longitu_surv_curve)

ggsave("Results plots/longitu_surv_curve.png", longitu_surv_curve$plot, width = 9, height = 6, dpi = 600)


###################################################################
# 1. Prepare cause of death list
# target_causes <- c("All-Cause", "Non-C")
# cause_names <- c(
#   "All-Cause" = "All-Cause Mortality",
#   "Non-C" = "Non-C: Non-Neoplasms"
# )

table(df_tv2$cause_lv1)
target_causes <- c(
  "All-Cause",
  "Chapter II Neoplasms",
  "Chapter IX Diseases of the circulatory system"
  # "Chapter V Mental and behavioural disorders",
  # "Chapter VI Diseases of the nervous system",
  # "Chapter X Diseases of the respiratory system",
  # "Chapter XI Diseases of the digestive system"
)

cause_labels <- c(
  "All-Cause"                                     = "All-Cause Mortality",
  "Chapter II Neoplasms"                          = "C: Neoplasms (Cancer)",
  "Chapter IX Diseases of the circulatory system" = "I: Circulatory"
  # "Chapter V Mental and behavioural disorders"    = "F: Mental/MBD",
  # "Chapter VI Diseases of the nervous system"     = "G: Nervous System",
  # "Chapter X Diseases of the respiratory system"  = "J: Respiratory",
  # "Chapter XI Diseases of the digestive system"   = "K: Digestive"
)

# 2. Define comparison calculation function
get_comparison_stats <- function(df2, sex_include, cause_code, visit_filter) {
  
  valid_ids <- df2 %>%
    dplyr::filter(Sex %in% sex_include) %>%
    dplyr::group_by(Participant.ID) %>%
    dplyr::summarise(n_visits = dplyr::n_distinct(Instance), .groups = "drop") %>%
    dplyr::filter(n_visits %in% visit_filter) %>%
    dplyr::pull(Participant.ID)
  
  if (length(valid_ids) == 0) return(NULL)
  
  sub_df <- df2 %>%
    dplyr::filter(Participant.ID %in% valid_ids, Sex %in% sex_include) %>%
    dplyr::filter(is.finite(start), is.finite(stop), stop > start)
  
  if (nrow(sub_df) == 0) return(NULL)
  
  # ---- Define outcome_event and force NA -> 0 ----
  if (cause_code == "All-Cause") {
    sub_df$outcome_event <- dplyr::if_else(is.na(sub_df$event), 0, as.integer(sub_df$event))
  } else {
    sub_df$outcome_event <- ifelse(
      !is.na(sub_df$event) & sub_df$event == 1 & sub_df$cause_lv1 == cause_code,
      1L, 0L
    )
  }
  
  # ---- Deduplicated event count (robust to all-NA groups) ----
  n_events <- sub_df %>%
    dplyr::group_by(Participant.ID) %>%
    dplyr::summarise(
      e = as.integer(any(outcome_event == 1, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::summarise(n = sum(e), .groups = "drop") %>%
    dplyr::pull(n)
  
  cat("Cause:", cause_code, "| Visits:", paste(visit_filter, collapse=" "), "| Events:", n_events, "\n")
  if (!is.finite(n_events) || n_events < 5) return(NULL)
  
  tryCatch({
    fit_base <- survival::coxph(
      survival::Surv(start, stop, outcome_event) ~ Sex,
      data = sub_df,
      cluster = Participant.ID
    )
    c_base <- as.numeric(summary(fit_base)$concordance[1])
    
    fit_full <- survival::coxph(
      survival::Surv(start, stop, outcome_event) ~ adj_BAA * Sex,
      data = sub_df,
      cluster = Participant.ID
    )
    c_full <- as.numeric(summary(fit_full)$concordance[1])
    
    delta_c <- c_full - c_base
    
    # 修正 2: 必须显式加载 survey 包或确保环境正确
    # 且 regTermTest 对于 coxph 对象有时需要指定 null 模型
    w_test <- survey::regTermTest(fit_full, ~ adj_BAA + adj_BAA:Sex)
    p_val  <- as.numeric(w_test$p)
    
    # 返回命名向量
    return(c(Events = n_events, 
             C_Base = c_base, 
             C_Full = c_full, 
             Delta  = delta_c, 
             P_Value = p_val))
  }, error = function(e) return(NULL))
}

results_list <- list()

for (cause in target_causes) {
  # Calculate for three time points
  res_23  <- get_comparison_stats(df_tv2, c(0,1), cause, c(2,3))
  print(res_23)
  res_2 <- get_comparison_stats(df_tv2, c(0,1), cause, c(2))
  res_3 <- get_comparison_stats(df_tv2, c(0,1), cause, c(3))
  
  # Check if res is null (for cases with too few events)
  if(is.null(res_23)) next
  
  row_data <- data.frame(
    Cause = cause_labels[cause],
    # Full Follow-up (0)
    Events_23 = res_23["Events"],
    C_Base_23 = round(res_23["C_Base"], 3),
    C_Full_23 = round(res_23["C_Full"], 3),
    Delta_23  = round(res_23["Delta"], 3),
    Wald_p_23 =  res_23["P_Value"],
    
    Events_2 = ifelse(is.null(res_2), NA, res_2["Events"]),
    C_Base_2 = ifelse(is.null(res_2), NA, round(res_2["C_Base"], 3)),
    C_Full_2 = ifelse(is.null(res_2), NA, round(res_2["C_Full"], 3)),
    Delta_2  = ifelse(is.null(res_2), NA, round(res_2["Delta"], 3)),
    Wald_p_2 =  res_2["P_Value"],
    
    Events_3 = ifelse(is.null(res_3), NA, res_3["Events"]),
    C_Base_3 = ifelse(is.null(res_3), NA, round(res_3["C_Base"], 3)),
    C_Full_3 = ifelse(is.null(res_3), NA, round(res_3["C_Full"], 3)),
    Delta_3  = ifelse(is.null(res_3), NA, round(res_3["Delta"], 3)),
    Wald_p_3 =  ifelse(is.null(res_3), NA, (res_3["P_Value"]))
  )
  results_list[[cause]] <- row_data
}

final_comparison_table <- do.call(rbind, results_list)
row.names(final_comparison_table) <- NULL
print(final_comparison_table)

fmt_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  out <- ifelse(is.na(p), NA_character_,
                ifelse(p < 1e-3, "<0.001",
                       format(round(p, 3), nsmall = 3)))
  out
}

final_comparison_table <- final_comparison_table %>%
  mutate(
    Wald_p_23  = fmt_p(Wald_p_23),
    Wald_p_2  = fmt_p(Wald_p_2),
    Wald_p_3 = fmt_p(Wald_p_3)
  )

print(final_comparison_table)

final_comparison_table %>%
  kbl(caption = "Predictive performance of adj_BAA by cause of death and number of assessment visits (Instances 0/2/3)") %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  add_header_above(c(
    " " = 1,
    "Visits >= 2" = 5,
    "Visits = 2"      = 5,
    "Visits = 3"      = 5
  ))

dev.off()

###############################################
############ df_tv2 stats #####################
###############################################
head(df_tv2)

length(unique(df_tv2$Participant.ID))
summary(df_tv2)

unique_people <- df_tv2[!duplicated(df_tv2$Participant.ID), ]
table(Sex = unique_people$Sex, Visits = unique_people$n_visits)

tapply(df_tv2$Age.at.assessment, df_tv2$Sex, summary)

df_tv2 %>%
  group_by(Participant.ID, Sex) %>%
  summarise(n_v = first(n_visits), .groups = "drop") %>%
  { table(.$Sex, .$n_v) }

df_tv2 %>%
  group_by(Participant.ID, Sex) %>%
  summarise(is_dead = any(event == 1, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sex) %>%
  summarise(
    Total_N = n(),
    Deaths = sum(is_dead),
    Death_Rate = paste0(round(mean(is_dead) * 100, 2), "%")
  )

tapply(df_tv2$death_age, df_tv2$Sex, summary)











