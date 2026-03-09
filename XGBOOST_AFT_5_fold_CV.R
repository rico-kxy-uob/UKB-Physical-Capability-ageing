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

setwd("/rds/projects/s/spillf-systems-mechanobiology-health-disease/RDS-Rico/R/Physical mea ageing clock")
# setwd('C:/Users/Kaiyi/OneDrive - University of Birmingham/Desktop/PhD/PhD/dataset_physical_mea')

source('Functions/longitudinal data analysis functions.R')
source('Functions/AFT analysis funcs.R')

extra_1smeas_raw<-read.csv('deep_cleaned_1smeas_ins0123_lfm_for_cox.csv',header=TRUE)
training_2plusmeas_raw<-read.csv('deep_cleaned_2plusmeas_ins023_lfm.csv',header=TRUE)
summary(extra_1smeas_raw)

length(unique(extra_1smeas_raw$Participant.ID))

table(extra_1smeas_raw$Instance)
table(extra_1smeas_raw$Participant.ID %in% extra_1smeas_raw$Participant.ID)

extra_1smeas_raw_ins0<-subset(extra_1smeas_raw,Instance %in% c(0,2,3))
dim(extra_1smeas_raw_ins0)
table(extra_1smeas_raw_ins0$Instance)
# Calculate NA percentage for each column
na_percent <- sapply(extra_1smeas_raw_ins0, function(x) {
  mean(is.na(x)) * 100
})

# Print NA percentage in descending order
na_percent_sorted <- sort(na_percent, decreasing = TRUE)
print(na_percent_sorted)

# Identify columns with >60% NA
cols_to_drop <- names(na_percent[na_percent > 40])
cat("Columns to drop (>60% NA):\n")
print(cols_to_drop)

# Remove these columns
extra_1smeas <- extra_1smeas_raw_ins0[, !(names(extra_1smeas_raw_ins0) %in% cols_to_drop)]

# Print resulting dimension
cat("Original dim:", dim(extra_1smeas_raw_ins0), "\n")
cat("Cleaned dim :", dim(extra_1smeas), "\n")
colnames(extra_1smeas)
confounder_outcome_fields_extra <- c(
  # "Participant.ID", 
  # "Sex", "Instance"，
  # 'Body.mass.index..BMI.',
  # 'Heel_bone_Tscore_', 
  'Hip.circumference',
  # 'Pulse_rate',
  "Trunk.fat.percentage","Arm.fat.percentage",
  "Leg.fat.percentage", 'Body.fat.percentage',
  "Handedness_value",
  # "Pulse.wave.Arterial.Stiffness.index",
  # 'Heel.bone.mineral.density..BMD..T.score..automated',
  # "Trunk.fat.mass","Arm.fat.mass",
  # "Leg.fat.mass"  ,"Whole.body.fat.mass" ,
  # 'Sitting.height',
  'Standing.height',
  # 'Overall.acceleration.average',
  "Health.score..England.", "Living.environment.score..England.",
  "Overall.health.rating","Weight.method",                                             
  "Long.standing.illness..disability.or.infirmity",
  "Falls.in.the.last.year", "Weight.change.compared.with.1.year.ago","Method.of.measuring.blood.pressure")

confounders_extra <- extra_1smeas %>% dplyr::select(all_of(confounder_outcome_fields_extra))

cat("Confounder-outcome dim:", dim(confounders_extra), "\n")

extra_1smeas_death_data<-read.csv('death_birth_info/extra_1meas_death_cause.csv',header=TRUE)

# extra_1smeas_death_data <- subset(
#   extra_1smeas_death_data,
#   !grepl("^[C]", Underlying..primary..cause.of.death..ICD10...Instance.0 %>% as.character())
# )

extra_1smeas_death_data <- subset(
  extra_1smeas_death_data,
  !grepl("^[VWXY]", Underlying..primary..cause.of.death..ICD10...Instance.0 %>% as.character())
)

train_data_extra_1smeas_na <- extra_1smeas%>% 
  select(-all_of(colnames(confounders_extra))) %>% 
  filter(Instance %in% c(0,2,3))
colnames(train_data_extra_1smeas_na)

train_data_extra_1smeas<-train_data_extra_1smeas_na %>%
  filter(if_all(everything(), ~ !is.na(.)))
head(train_data_extra_1smeas)
dim(train_data_extra_1smeas)
summary((train_data_extra_1smeas$Age.when.attended.assessment.centre))

careful_ids <- intersect(train_data_extra_1smeas$Participant.ID,training_2plusmeas_raw$Participant.ID)
careful_ids

colnames(train_data_extra_1smeas)
head(extra_1smeas_death_data)

head(train_data_extra_1smeas)
head(extra_1smeas_death_data)
########################################################
# ===== AFT table =====
dat_aft <- prepare_aft_table(
  meas_df  = train_data_extra_1smeas,
  death_df = extra_1smeas_death_data
)
colnames(dat_aft)
table(is.na(dat_aft$assess_date0))
summary(dat_aft$age0)
dim(dat_aft)
table(dat_aft$Instance)
table(dat_aft$Sex)
table(dat_aft$status)
summary(dat_aft$time_years)
summary(dat_aft$age0)
#####################################################################
library(xgboost)
library(dplyr)

source('Functions/AFT_CV_funcs.R')
dat <- dat_aft

exclude <- c(
  # identifiers / bookkeeping
  "Participant.ID", "Instance",
  
  # target / survival outcome
  "time_years", "status","time_years_age","time_years_date",
  
  # sex and age (you said: ensure not include age columns)
  "Sex","inst_i",
  "Age.when.attended.assessment.centre",
  "age0", "age0_from_date","age_raw",
  
  # dates and derived ages
  "assess_date0", "censor_date",
  "birth_date", "censor_age",
  "death_date", "death_age",
  "death_age_from_date", "death_age_final",
  
  # birth Y/M (you keep these only for time construction)
  "yob", "mob",
  
  # optional intermediate date holders (if present)
  "ac_date0", "ac_date_by_inst"
)

feature_cols <- dat %>%
  select(where(is.numeric)) %>%
  select(-any_of(exclude)) %>%
  names()

# baseline / full parameters
params_base <- list(
  booster="gbtree",
  objective="survival:aft",
  eval_metric="aft-nloglik",
  eta=0.05, max_depth=2,
  min_child_weight=10, subsample=0.9, colsample_bytree=0.7,
  aft_loss_distribution="normal",
  aft_loss_distribution_scale=1.0
)

params_full <- list(
  booster="gbtree",
  objective="survival:aft",
  eval_metric="aft-nloglik",
  eta=0.05, max_depth=5,
  min_child_weight=10, subsample=0.7, colsample_bytree=0.7,
  aft_loss_distribution="normal",
  aft_loss_distribution_scale=1.0
)

# Train separately by sex
# res_f <- train_sex_specific_clock(dat %>% filter(Sex == 0), feature_cols, params_base, params_full, seed=2026)
# res_m <- train_sex_specific_clock(dat %>% filter(Sex == 1), feature_cols, params_base, params_full, seed=2026)



params_base_fixed <- params_base 

dat_f <- dat %>% dplyr::filter(Sex == 0)
dat_m <- dat %>% dplyr::filter(Sex == 1)

res_f_parameter_cv <- run_outer_holdout_tuning_with_GM(
  dat_sex = dat_f,
  feature_cols = feature_cols,
  params_base_fixed = params_base_fixed,
  K_inner = 5,
  n_trials = 100,
  seed = 2026,
  es_frac = 0.1,
  fit_age_range = c(45, 70),
  ref_age_range = c(30, 100),
  ref_step = 0.1
)

head(res_f_parameter_cv$tuning$leaderboard, 10)

res_f_parameter_cv$best_params_full
res_f_parameter_cv$external$cindex_external_test

write.csv(res_f_parameter_cv$tuning$leaderboard,'Functions/female_cv_leaderboard.csv',row.names = FALSE)

res_m_parameter_cv <- run_outer_holdout_tuning_with_GM(
  dat_sex = dat_m,
  feature_cols = feature_cols,
  params_base_fixed = params_base_fixed,
  K_inner = 5,
  n_trials = 100,
  seed = 2026,
  es_frac = 0.1,
  fit_age_range = c(45, 70),
  ref_age_range = c(30, 100),
  ref_step = 0.1
)

head(res_m_parameter_cv$tuning$leaderboard, 10)

res_m_parameter_cv$external$cindex_external_test

write.csv(res_m_parameter_cv$tuning$leaderboard,'Functions/male_cv_leaderboard.csv',row.names = FALSE)










