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
cat("Columns to drop (>40% NA):\n")
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
  eta=0.05, max_depth=2,  #2
  min_child_weight=10, subsample=0.7, colsample_bytree=0.7,
  # eta=0.05, max_depth=5,
  # min_child_weight=10, subsample=0.7, colsample_bytree=0.7,
  aft_loss_distribution="normal",
  aft_loss_distribution_scale=1
)

params_full_male<- list(
  booster="gbtree",
  objective="survival:aft",
  eval_metric="aft-nloglik",
  # eta=0.05, max_depth=5,
  # min_child_weight=30, subsample = 0.7, colsample_bytree=0.6,
  # lambda=0.012, alpha = 0.04, gamma= 0,
  # eta=0.05, max_depth=5,
  # min_child_weight=30, subsample=0.6, colsample_bytree=0.8,
  # eta=0.023, max_depth=3,
  # min_child_weight=20, subsample=0.74, colsample_bytree=0.7,
  # lambda=0.9, alpha = 0.2, gamma= 4.15,
  eta=0.05, max_depth=5,
  min_child_weight=10, subsample=0.7, colsample_bytree=0.7,
  aft_loss_distribution="normal",
  aft_loss_distribution_scale=1
)

params_full_female<- list(
  booster="gbtree",
  objective="survival:aft",
  eval_metric="aft-nloglik",
  # eta=0.05, max_depth=5,
  # min_child_weight=10, subsample=0.7, colsample_bytree=0.7,
  # eta=0.03, max_depth=4,
  # min_child_weight= 5, subsample=0.9, colsample_bytree=0.6,
  # lambda = 0.2, alpha = 0.006 ,gamma= 0.001,
  eta=0.03, max_depth=5,
  min_child_weight= 5, subsample=0.9, colsample_bytree=0.6,
  lambda = 36, alpha = 1.91 ,gamma= 1.21,
  # eta=0.03, max_depth=5,
  # min_child_weight=10, subsample=0.7, colsample_bytree=0.7,
  aft_loss_distribution="normal",
  aft_loss_distribution_scale=1
)


# Train separately by sex
# res_f <- train_sex_specific_clock(dat %>% filter(Sex == 0), feature_cols, params_base, params_full, seed=2026)
# res_m <- train_sex_specific_clock(dat %>% filter(Sex == 1), feature_cols, params_base, params_full, seed=2026)

# gompertz risk curve #
res_f <- train_sex_specific_clock_Gompertz(dat %>% filter(Sex == 0), feature_cols, params_base, params_full_female, seed=2026, 
                                           censor_check = TRUE,
                                           fit_age_range = c(45, 70),   
                                           ref_age_range = c(30, 95),  
                                           ref_step = 0.1)
res_m <- train_sex_specific_clock_Gompertz(dat %>% filter(Sex == 1), feature_cols, params_base, params_full_male, seed=2026, 
                                           censor_check = TRUE,
                                           fit_age_range = c(45, 70),  
                                           ref_age_range = c(30, 95),  
                                           ref_step = 0.1)

res_f$cindex_test
res_m$cindex_test

trainset_f <- res_f$tr_set
trainset_m <- res_m$tr_set
summary(subset(dat,Sex == 0)[trainset_f,'age0'])

summary(subset(dat,Sex == 1)[trainset_m,'age0'])

par(mfrow=c(2,1))
hist(subset(dat,Sex == 0)[trainset_f,'age0'],main='Female age distribution')
hist(subset(dat,Sex == 1)[trainset_f,'age0'],main='Male age distribution')

testset_f <- res_f$te_set
testset_m <- res_m$te_set

bundle_list <- list(`0` = res_f$bundle, `1` = res_m$bundle)
# saveRDS(bundle_list, "aft_clock_sex_specific_bundles_ins023_Gompertz_CV.rds")

table(dat_aft$Sex)
dim(dat_aft)

table(dat_aft$status,dat_aft$Sex)
tapply(dat_aft$age0, dat_aft$Sex, summary)
tapply(dat_aft$death_age_from_date, dat_aft$Sex, summary)
tapply(dat_aft$time_years, dat_aft$Sex, summary)
summary(dat_aft)
############ vis ################
dat_m <- dat_aft %>% filter(Sex == 1) 

not_instant_death_index <- which(dat_m$Sex == 1 & dat_m$time_years >= 5)

# b <- res_m$bundle   # male bundle (including bst_full, feature_cols, impute_med)
b <- bundle_list$`1`
bst_full <- b$bst_full
fc <- b$feature_cols

# Prepare X (impute missing values with median from training set)
X_m <- data.matrix(dat_m[, fc, drop = FALSE])
for (j in seq_along(fc)) {
  na_idx <- is.na(X_m[, j])
  if (any(na_idx)) X_m[na_idx, j] <- b$impute_med[[j]]
}

# If you use sex-specific internal split, b$split$test is more convenient
train_m <- b$split$train  # This is the index within dat_m (recommended)

set.seed(1)
n_show <- 10000
idx_show <- intersect(sample(train_m, min(n_show, length(train_m))),not_instant_death_index)

X_show <- X_m[idx_show, , drop = FALSE]
dim(X_show)
# Calculate TreeSHAP (for AFT: contributions on log-time / margin scale)
shap_contrib <- predict(bst_full, X_show, predcontrib = TRUE)
shap_mat <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]  # Remove BIAS

# risk = -margin， so we negate the SHAP values here
sv_risk <- shapviz(-shap_mat, X = X_show)
colnames(sv_risk)[2]<-'FEV1.'
male_shap_FEV1.<-sv_dependence(sv_risk, "FEV1.", color_var = "auto")+
  theme_classic()+
  labs(title = "Male: Dependence plot of FEV1..")  
male_shap_FEV1.

male_shap_all<-sv_importance(sv_risk, kind = "beeswarm", show_numbers = TRUE,max_display = 30)+
  theme_classic()+
  labs(title = "Male: SHAP importance of features")  
male_shap_all
ggsave("Results plots/male_sv_importance_beeswarm_top30.png", male_shap_all, width = 12, height = 8, dpi = 600)

male_shap_BMI<-sv_dependence(sv_risk, "Body.mass.index..BMI.", color_var = "auto")+
  theme_classic()+
  labs(title = "Male: Dependence plot of Body.mass.index..BMI.")  
male_shap_BMI
ggsave("Results plots/male_sv_importance_BMI_auto.png", male_shap_BMI, width = 9, height = 6, dpi = 600)

male_shap_impedance_leg<-sv_dependence(sv_risk, "Impedance.of.leg", color_var = "auto")+
  theme_classic()+ xlim(100,450)+
  labs(title = "Male: Dependence plot of Impedance.of.leg")  
male_shap_impedance_leg
ggsave("Results plots/male_sv_importance_impedance_leg_auto.png", male_shap_impedance_leg, width = 9, height = 6, dpi = 600)

male_shap_BP<-sv_dependence(sv_risk, "SystolicBP", color_var = 'auto')+
  theme_classic()+
  labs(title = "Male: Dependence plot of Sys/Dia Blood Pressure")  
male_shap_BP
# male_shap_BP<-sv_dependence(sv_risk, "DiastolicBP", color_var = 'SystolicBP')+
#   theme_classic()+
#   labs(title = "Male: Dependence plot of Sys/Dia Blood Pressure")  

ggsave("Results plots/male_sv_importance_blood_pressure.png", male_shap_BP, width = 9, height = 6, dpi = 600)

male_shap_BP_diaauto<-sv_dependence(sv_risk, "DiastolicBP", color_var = 'auto')+
  theme_classic()+
  labs(title = "Male: Dependence plot of Dia BP and waist circumference")  
male_shap_BP_diaauto

ggsave("Results plots/male_sv_importance_BP_diaauto.png", male_shap_BP_diaauto, width = 9, height = 6, dpi = 600)


sv_risk$X$Hip.circumference
Leg.fat.free.mass
age_grid <- seq(40, 85, by = 1)
X_age <- matrix(age_grid, ncol = 1); colnames(X_age) <- "age0"

base_margin <- predict(res_m$bundle$bst_base, X_age, outputmargin = TRUE)
plot(age_grid, -base_margin, type = "l", xlab = "Age", ylab = "Baseline risk (-margin)")

##########################################
dat_f <- dat_aft %>% filter(Sex == 0)

not_instant_death_index <- which(dat_f$Sex == 0 & dat_f$time_years >5)

b <- bundle_list$`0`
bst_full <- b$bst_full
fc <- b$feature_cols

X_f <- data.matrix(dat_f[, fc, drop = FALSE])
for (j in seq_along(fc)) {
  na_idx <- is.na(X_f[, j])
  if (any(na_idx)) X_f[na_idx, j] <- b$impute_med[[j]]
}

train_f <- b$split$train
set.seed(1)
idx_show <- intersect(sample(train_f, min(10000, length(train_f))),not_instant_death_index)

X_show <- X_f[idx_show, , drop = FALSE]
dim(X_show)

shap_contrib <- predict(bst_full, X_show, predcontrib = TRUE)
shap_mat <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]

sv_risk <- shapviz(-shap_mat, X = X_show)
colnames(sv_risk)[2]<-'FEV1.'
female_shap_all<-sv_importance(sv_risk, kind = "beeswarm", show_numbers = TRUE,max_display = 30)+
  theme_classic()+d
  labs(title = "Female: SHAP importance of features")  
female_shap_all
ggsave("Results plots/female_sv_importance_beeswarm_top30.png", female_shap_all, width = 12, height = 8, dpi = 600)

female_shap_BMI<-sv_dependence(sv_risk, "Body.mass.index..BMI.", color_var = "auto")+
  theme_classic()+
  labs(title = "Female: Dependence plot of Body.mass.index..BMI.")  
female_shap_BMI
ggsave("Results plots/female_sv_importance_BMI_auto.png", female_shap_BMI, width = 9, height = 6, dpi = 600)

female_shap_impedance_leg<-sv_dependence(sv_risk, "Impedance.of.leg", color_var = "auto")+
  theme_classic()+xlim(100,450)+
  labs(title = "Female: Dependence plot of Impedance.of.leg")  
female_shap_impedance_leg
ggsave("Results plots/female_sv_importance_impedance_leg_auto.png", female_shap_impedance_leg, width = 9, height = 6, dpi = 600)


# female_shap_BP<-sv_dependence(sv_risk, "DiastolicBP", color_var = 'SystolicBP')+
#   theme_classic()+
#   labs(title = "Female: Dependence plot of Dia/Sys Blood Pressure.")  

female_shap_BP<-sv_dependence(sv_risk, "SystolicBP", color_var = 'DiastolicBP')+
  theme_classic()+
  labs(title = "Female: Dependence plot of Sys/Dia Blood Pressure.")  
female_shap_BP
ggsave("Results plots/female_sv_importance_Blood_pressure.png", female_shap_BP, width = 9, height = 6, dpi = 600)

female_shap_BP_sysauto<-sv_dependence(sv_risk, "SystolicBP", color_var = 'auto')+
  theme_classic()+
  labs(title = "Female: Dependence plot of Sys Bp and pulse rate.")  
female_shap_BP_sysauto
ggsave("Results plots/female_sv_importance_BP_sysauto.png", female_shap_BP_sysauto, width = 9, height = 6, dpi = 600)

female_shap_BP_diaauto<-sv_dependence(sv_risk, "DiastolicBP", color_var = 'auto')+
  theme_classic()+
  labs(title = "Female: Dependence plot of Dia Bp and waist circumference.")  
female_shap_BP_diaauto
# ggsave("Results plots/female_sv_importance_BP_diaauto.png", female_shap_BP_diaauto, width = 9, height = 6, dpi = 600)

sv_risk$X$Leg.fat.free.mass
sv_risk$X$Body.mass.index..BMI.



## male ##
age_grid <- seq(40, 85, by = 0.2)
X_age <- matrix(age_grid, ncol = 1); colnames(X_age) <- "age0"
base_margin <- predict(bundle_list$`1`$bst_base, X_age, outputmargin = TRUE)

par(mfrow=c(2,1))
plot(age_grid, -base_margin, type = "l", xlab = "Age", ylab = "Baseline risk (-margin)",main='Male')

plot(bundle_list$`1`$age_grid,bundle_list$`1`$risk_grid,xlab = "Age", ylab = "Risk (-margin)",main='Male')
dev.off()

p <- bundle_list$`1`$gm_params
A <- as.numeric(p["A"])
B <- as.numeric(p["B"])
C <- as.numeric(p["C"])
A
B
C
# Age range for plotting
library(ggsci)

age_grid <- seq(37, 85, by = 0.1)

# Gompertz–Makeham reference risk curve: r(age) = A + B * exp(C * age)
risk_gm <- A + B * exp(C * age_grid)

X_age <- matrix(age_grid, ncol = 1)
colnames(X_age) <- "age0"
risk_base_raw <- -predict(bundle_list$`1`$bst_base, X_age, outputmargin = TRUE)

# Build a tidy plotting table
df_curve <- tibble(
  age = age_grid,
  `GM reference` = as.numeric(risk_gm),
  `bst_base raw (-margin)` = as.numeric(risk_base_raw)
) %>%
  pivot_longer(-age, names_to = "curve", values_to = "risk") %>%
  mutate(
    curve = factor(curve, levels = c("GM reference", "bst_base raw (-margin)"))
  )

gm_label <- sprintf(
  "GM reference: r(age) = %.3f + %.3f * exp(%.3f * age)",
  A, B, C
)

p <- ggplot(df_curve, aes(x = age, y = risk, color = curve)) +
  geom_line(linewidth = 1.05, alpha = 0.95) +
  geom_vline(xintercept = c(40, 70), linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
  scale_color_npg(name = NULL) +
  labs(
    title = "Male: Gompertz–like extrapolated baseline risk curve",
    subtitle = gm_label,
    x = "Age",
    y = "Baseline reference risk"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

p
ggsave("Results plots/Male_mortality_curve.png", p, width = 9, height = 6, dpi = 600)

# bundle_list$`1`$bst_base


## female ##
age_grid <- seq(40, 85, by = 0.2)
X_age <- matrix(age_grid, ncol = 1); colnames(X_age) <- "age0"
base_margin <- predict(bundle_list$`0`$bst_base, X_age, outputmargin = TRUE)

par(mfrow=c(2,1))
plot(age_grid, -base_margin, type = "l", xlab = "Age", ylab = "Baseline risk (-margin)",main='Female')

plot(bundle_list$`0`$age_grid,bundle_list$`0`$risk_grid,xlab = "Age", ylab = "Risk (-margin)",main='Female')
dev.off()

p <- bundle_list$`0`$gm_params
A <- as.numeric(p["A"])
B <- as.numeric(p["B"])
C <- as.numeric(p["C"])
A
B
C

age_grid <- seq(37, 85, by = 0.1)

# Gompertz–Makeham reference risk curve: r(age) = A + B * exp(C * age)
risk_gm <- A + B * exp(C * age_grid)

X_age <- matrix(age_grid, ncol = 1)
colnames(X_age) <- "age0"
risk_base_raw <- -predict(bundle_list$`0`$bst_base, X_age, outputmargin = TRUE)

# Build a tidy plotting table
df_curve <- tibble(
  age = age_grid,
  `GM reference` = as.numeric(risk_gm),
  `bst_base raw (-margin)` = as.numeric(risk_base_raw)
) %>%
  pivot_longer(-age, names_to = "curve", values_to = "risk") %>%
  mutate(
    curve = factor(curve, levels = c("GM reference", "bst_base raw (-margin)"))
  )

gm_label <- sprintf(
  "GM reference: r(age) = %.3f + %.3f * exp(%.3f * age)",
  A, B, C
)

p_f <- ggplot(df_curve, aes(x = age, y = risk, color = curve)) +
  geom_line(linewidth = 1.05, alpha = 0.95) +
  geom_vline(xintercept = c(40, 70), linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
  scale_color_npg(name = NULL) +
  labs(
    title = "Female: Gompertz–like extrapolated baseline risk curve",
    subtitle = gm_label,
    x = "Age",
    y = "Baseline reference risk"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

p_f

ggsave("Results plots/Female_mortality_curve.png", p_f, width = 9, height = 6, dpi = 600)

############### Male vs. female GM ###############
age_grid <- seq(37, 85, by = 0.1)

p_m <- bundle_list$`1`$gm_params
risk_gm_m <- as.numeric(p_m["A"]) + as.numeric(p_m["B"]) * exp(as.numeric(p_m["C"]) * age_grid)

df_male <- tibble(
  age = age_grid,
  risk = risk_gm_m,
  gender = "Male"
)

p_f_params <- bundle_list$`0`$gm_params
risk_gm_f <- as.numeric(p_f_params["A"]) + as.numeric(p_f_params["B"]) * exp(as.numeric(p_f_params["C"]) * age_grid)

df_female <- tibble(
  age = age_grid,
  risk = risk_gm_f,
  gender = "Female"
)

df_combined <- bind_rows(df_male, df_female)

p_combined <- ggplot(df_combined, aes(x = age, y = risk, color = gender)) +
  geom_line(linewidth = 1.2, alpha = 0.9) +
  geom_vline(xintercept = c(40, 70), linetype = "dashed", linewidth = 0.5, color = "grey50") +
  scale_color_npg() +
  labs(
    title = "Comparison of Gompertz–like Baseline Risk",
    subtitle = "Extrapolated curves for Male vs Female",
    x = "Age",
    y = "GM Reference Risk",
    color = "Gender"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_combined)

ggsave("Results plots/Combined_gender_mortality_curve.png", p_combined, width = 9, height = 6, dpi = 600)

####### sample skew check >70 ##########
df <- dat_aft   # Your actual training data frame
# If you have train index tr, use df_tr <- df[tr,]
df_tr <- df  

plot(jitter(df_tr$age0, amount=0.1), df_tr$status,
     pch=16, cex=0.5,
     xlab="age0", ylab="status (1=death)",
     main="Deaths along age0 (train)")
abline(v=70, lty=2)

br <- seq(floor(min(df_tr$age0)), ceiling(max(df_tr$age0)), by=1)
df_tr$age_bin <- cut(df_tr$age0, breaks=br, right=FALSE)

tabN <- as.integer(table(df_tr$age_bin))
tabD <- tapply(df_tr$status, df_tr$age_bin, sum)
tabP <- tabD / tabN

mid <- br[-length(br)] + 0.5
ok <- is.finite(tabP)

# library(Cairo)
# 
# CairoPNG("Results plots/Sample_skew_check.png",
#          width = 1800, height = 2200, res = 300)
par(mfrow = c(3,1), mar = c(4,4,2,1))
plot(mid[ok], tabN[ok], type="h", lwd=2, xlab="Age (bin)", ylab="N",
     main="Sample size per 1-year age bin"); abline(v=70, lty=2)
plot(mid[ok], tabD[ok], type="h", lwd=2, xlab="Age (bin)", ylab="Deaths",
     main="Deaths per 1-year age bin"); abline(v=70, lty=2)
plot(mid[ok], tabP[ok], type="b", pch=16, xlab="Age (bin)", ylab="Death proportion",
     main="Death proportion per 1-year age bin"); abline(v=70, lty=2)

dev.off()

dev.off()
plot(df_tr$age0, df_tr$time_years,
     pch=16, cex=0.4,
     col=ifelse(df_tr$status==1, "red", "grey50"),
     xlab="age0", ylab="time_years",
     main="Follow-up time vs CA (red=death)")
abline(v=70, lty=2)


###############################################################
# bundle_list <- readRDS("aft_clock_sex_specific_bundles.rds")
# 
# bundle_list <- readRDS("aft_clock_sex_specific_bundles_ins023.rds")
# bundle_list <- readRDS("aft_clock_sex_specific_bundles_ins023_Gompertz.rds")
bundle_list <- readRDS("aft_clock_sex_specific_bundles_ins023_Gompertz_CV.rds")
names(bundle_list)

bundle_list$`0`$split$train

trainset_f <- bundle_list$`0`$split$train
trainset_m <- bundle_list$`1`$split$train
summary(subset(dat,Sex == 0)[trainset_f,'age0'])

summary(subset(dat,Sex == 1)[trainset_m,'age0'])

par(mfrow=c(2,1))
hist(subset(dat,Sex == 0)[trainset_f,'age0'],main='Female age distribution')
hist(subset(dat,Sex == 1)[trainset_f,'age0'],main='Male age distribution')

testset_f <- bundle_list$`0`$split$test
testset_m <- bundle_list$`1`$split$test

###################################################################
ICD10_coding <- read_tsv('coding ICD10.tsv')

head(ICD10_coding,20)

##############################################################
################### Test set Performance #######################
################################################################
# new_person_df: External new data (at least containing Sex, age0 and as many columns as possible from feature_cols)
dat_f <- dat %>% dplyr::filter(Sex == 0)
testset_own_f <- dat_f[testset_f, , drop = FALSE]
dim(testset_own_f)

dat_m <- dat %>% dplyr::filter(Sex == 1)
testset_own_m <- dat_m[testset_m, , drop = FALSE]

test_data_own <- rbind(testset_own_f,testset_own_m)

test_data <- test_data_own %>% 
  select(-all_of(colnames(confounder_outcome_fields_extra))) %>% 
  filter(Instance %in% c(0,2,3))
colnames(test_data)

dim(test_data)
colnames(test_data)
head(test_data)
length(unique(test_data$Participant.ID))

table(test_data$Sex)
summary(test_data$age0)

test_data_all <- test_data 

dim(test_data_all)
# out_all <- score_BAA_one(test_data_all, bundle_list)
# out_all <- score_BAA_batch(test_data_all, bundle_list)
out_all <- score_BAA_batch_Gompertz(test_data_all, bundle_list)
dim(test_data_all)
dim(out_all)

testset_own_with_BAandBAA <- cbind(test_data,out_all)

testset_own_for_survival<-testset_own_with_BAandBAA[,c("Participant.ID" ,"censor_age",
                                                           # "Health.score..England.",
                                                           # "Living.environment.score..England.",
                                                           "Sex","Instance",'status','death_age_from_date',
                                                           "age0",'time_years','BAA','BA')]

head(testset_own_for_survival)
summary(testset_own_for_survival$age0)

plot(testset_own_for_survival$age0,testset_own_for_survival$BA)
cor.test(testset_own_for_survival$age0,testset_own_for_survival$BA)

plot(testset_own_for_survival$age0,testset_own_for_survival$BAA)
cor.test(testset_own_for_survival$age0,testset_own_for_survival$BAA)

testset_own_death_data<-read.csv('death_birth_info/extra_1meas_death_cause.csv',header=TRUE)
length(!is.na(testset_own_death_data$Date.of.death...Instance.0)==FALSE)

# testset_death_data <- subset(
#   testset_death_data,
#   !grepl("^[C]", Underlying..primary..cause.of.death..ICD10...Instance.0 %>% as.character())
# )

testset_own_death_data <- subset(
  testset_own_death_data,
  !grepl("^[VWXY]", Underlying..primary..cause.of.death..ICD10...Instance.0 %>% as.character())
)

head(testset_own_death_data)
testset_own_death_data$Underlying..primary..cause.of.death..ICD10...Instance.0

source('Functions/ICD10-mapping.R')

death_df <- testset_own_death_data %>%
  select(Participant.ID, Age.at.death...Instance.0, 
         Year.of.birth,Underlying..primary..cause.of.death..ICD10...Instance.0) %>%
  mutate(
    death_cause = Underlying..primary..cause.of.death..ICD10...Instance.0,
  ) %>%
  select(Participant.ID, Year.of.birth,death_cause)

# ---------- 2. Merge baseline covariates ----------
analysis_df_death_cause <- testset_own_for_survival %>%
  select(Participant.ID, event=status, censor_age, baseline_age = age0, 
         death_age_from_date, follow_up_time = time_years,
         Sex, BAA, BA 
         # , Health.score..England., Living.environment.score..England.
         ) %>%
  left_join(death_df, by = "Participant.ID")

analysis_df_death_cause <- add_icd10_cause_levels(analysis_df_death_cause, icd10_col = "death_cause")
analysis_df_death_cause %>% count(cause_lv1, sort = TRUE)
analysis_df_death_cause %>% count(cause_lv2, sort = TRUE)

head(analysis_df_death_cause,20)
table(analysis_df_death_cause$event)

# linear regression BA ~ baseline_age ---
df <- analysis_df_death_cause
df$BAA_raw <- df$BA - df$baseline_age
df$adj_BAA <- NA_real_

for (s in sort(unique(df$Sex))) {
  idx <- which(df$Sex == s & is.finite(df$BAA_raw) & is.finite(df$baseline_age))
  if (length(idx) < 20) next
  
  fit_q <- lm(BAA_raw ~ baseline_age + I(baseline_age^2), data = df, subset = idx)
  pred <- rep(NA_real_, nrow(df))
  pred[idx] <- predict(fit_q, newdata = df[idx, , drop = FALSE])
  df$adj_BAA[idx] <- df$BAA_raw[idx] - pred[idx]
}

analysis_df_death_cause <- df

cor.test(analysis_df_death_cause$baseline_age,analysis_df_death_cause$adj_BAA)

# follow-up years filter
analysis_df <- analysis_df_death_cause %>% filter(!is.na(follow_up_time) & follow_up_time >0)
dim(analysis_df)

plot_data <- analysis_df %>%
  filter(!is.na(Sex)) %>%
  mutate(Sex_label = factor(Sex, levels = c(0, 1), labels = c("Female", "Male")))

testset_BAA_density_plot <- ggplot(plot_data, aes(x = adj_BAA, fill = Sex_label, color = Sex_label)) +
  geom_density(alpha = 0.4, size = 0.8) +
  geom_vline(data = plot_data %>% group_by(Sex_label) %>% summarise(m = mean(adj_BAA, na.rm=T)),
             aes(xintercept = m, color = Sex_label), linetype = "dashed", size = 0.6) +
  
  scale_fill_nejm() +
  scale_color_nejm() +
  
  labs(
    x = "Biological Age Acceleration (BAA)",
    y = "Density",
    fill = "Sex",
    color = "Sex",
    title = "Distribution of BAA by Sex"
  ) +
  
  theme_bw() + 
  theme(
    panel.grid = element_blank(),           
    legend.position = c(0.85, 0.85),        
    legend.background = element_blank()
    # axis.text = element_text(size = 10, color = "black"),
    # axis.title = element_text(size = 12, face = "bold"),
    # plot.title = element_text(hjust = 0.5, size = 14) 
  )

print(testset_BAA_density_plot)

ggsave("Results plots/testset_BAA_density_plot.png", testset_BAA_density_plot, width = 9, height = 6, dpi = 600)

with(analysis_df, tapply(adj_BAA, Sex, mean, na.rm=TRUE))
with(analysis_df, tapply(adj_BAA, Sex, median, na.rm=TRUE))

cor(analysis_df$BAA, analysis_df$baseline_age, use="complete.obs") 

# ---------- 3. Kaplan-Meier & Cox baseline model ----------
surv_obj <- Surv(time = analysis_df$follow_up_time, event = analysis_df$event)

# KM overall
km_fit <- survfit(surv_obj ~ 1, data = analysis_df)
plot(km_fit, xlab = "Years since baseline", ylab = "Survival probability", main = "Kaplan-Meier Overall")

dev.off()
# KM by Sex
km_sex <- survfit(surv_obj ~ Sex, data = analysis_df)
plot(km_sex, xlab = "Years since baseline", ylab = "Survival", main = "KM by Sex")
legend("topright", legend = levels(as.factor(analysis_df$Sex)), lty = 1:2, col = 1:2)

analysis_df$age_c    <- analysis_df$baseline_age -57

cox_full_withBAA <- coxph(Surv(follow_up_time, event) ~ age_c*Sex + age_c*adj_BAA + adj_BAA*Sex
                          # + Sex*health_c + Sex*living_c
                          # + BAA*health_c + BAA*living_c
                          # + health_c + living_c
                          ,
                          data=analysis_df)
summary(cox_full_withBAA)


cox_full_withoutBAA <- coxph(Surv(follow_up_time, event) ~ age_c*Sex + age_c + Sex
                             # + Sex*health_c + Sex*living_c
                             # + health_c + living_c
                             ,
                             data=analysis_df)
summary(cox_full_withoutBAA)

anova(cox_full_withBAA, cox_full_withoutBAA, test="LRT")


################################################################
library(survival)
library(dplyr)
library(kableExtra)
table(analysis_df$cause_lv1)
# table(analysis_df$cause_lv2)
target_causes <- c(
  "All-Cause",
  "Chapter II Neoplasms",
  "Chapter IX Diseases of the circulatory system",
  "Chapter V Mental and behavioural disorders",
  "Chapter VI Diseases of the nervous system",
  "Chapter X Diseases of the respiratory system",
  "Chapter XI Diseases of the digestive system"
)

cause_labels <- c(
  "All-Cause"                                     = "All-Cause Mortality",
  "Chapter II Neoplasms"                          = "C: Neoplasms (Cancer)",
  "Chapter IX Diseases of the circulatory system" = "I: Circulatory",
  "Chapter V Mental and behavioural disorders"    = "F: Mental/MBD",
  "Chapter VI Diseases of the nervous system"     = "G: Nervous System",
  "Chapter X Diseases of the respiratory system"  = "J: Respiratory",
  "Chapter XI Diseases of the digestive system"   = "K: Digestive"
)

results_list <- list()

for (cause in target_causes) {
  # res_0  <- get_comparison_stats(analysis_df,c(0,1) ,cause, 0)
  # res_5  <- get_comparison_stats(analysis_df,c(0,1) ,cause, 5)
  # res_10 <- get_comparison_stats(analysis_df,c(0,1) ,cause, 10)
  
  res_0  <- get_comparison_stats_landmark(analysis_df,c(0,1) ,cause, 0)
  res_5  <- get_comparison_stats_landmark(analysis_df,c(0,1) ,cause, 5)
  res_10 <- get_comparison_stats_landmark(analysis_df,c(0,1) ,cause, 10)
  
  print(res_0)
  
  row_data <- data.frame(
    Cause     = as.character(cause_labels[cause]),
    
    # Full Follow-up
    Events_0  = res_0[["Events"]],
    C_Base_0  = round(as.numeric(res_0[["C_Base.C"]]), 3),
    C_Full_0  = round(as.numeric(res_0[["C_Full.C"]]), 3),
    Delta_0   = round(as.numeric(res_0[["Delta.C"]]), 3),
    P_LRT_0  =  as.numeric(res_0[["P_LRT"]]),
    
    # > 5 Years
    Events_5  = res_5[["Events"]],
    C_Base_5  = round(as.numeric(res_5[["C_Base.C"]]), 3),
    C_Full_5  = round(as.numeric(res_5[["C_Full.C"]]), 3),
    Delta_5   = round(as.numeric(res_5[["Delta.C"]]), 3),
    P_LRT_5  =  as.numeric(res_5[["P_LRT"]]),
    
    # > 10 Years
    Events_10 = res_10[["Events"]],
    C_Base_10 = round(as.numeric(res_10[["C_Base.C"]]), 3),
    C_Full_10 = round(as.numeric(res_10[["C_Full.C"]]), 3),
    Delta_10  = round(as.numeric(res_10[["Delta.C"]]), 3),
    P_LRT_10  =  as.numeric(res_10[["P_LRT"]]),
    
    stringsAsFactors = FALSE
  )
  results_list[[cause]] <- row_data
}

final_comparison_table <- do.call(rbind, results_list)
row.names(final_comparison_table) <- NULL

final_comparison_table <- final_comparison_table %>%
  mutate(
    P_LRT_0  = fmt_p(P_LRT_0),
    P_LRT_5  = fmt_p(P_LRT_5),
    P_LRT_10 = fmt_p(P_LRT_10)
  )

print(final_comparison_table)

final_comparison_table %>%
  kbl(caption = "Performance of BAA across Specific Causes of Death and Follow-up Windows") %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  add_header_above(c(" " = 1,
                     "Full Follow-up" = 5,
                     "> 5 Years" = 5,
                     "> 10 Years" = 5))



library(tidycmprsk) 

head(analysis_df,20)

sub_df_FG <- analysis_df %>%
  mutate(status_cr = case_when(
    is.na(event) | event == 0 ~ "censor",
    event == 1 & cause_lv1 == 'Chapter IX Diseases of the circulatory system' ~ "target",
    event == 1 & cause_lv1 != 'Chapter IX Diseases of the circulatory system' ~ "competing",
    TRUE ~ "censor"
  )) %>%
  mutate(status_cr = factor(status_cr, levels = c("censor", "target", "competing")))

fit_base_FG <- tidycmprsk::crr(Surv(follow_up_time, status_cr) ~ Sex, data = sub_df_FG)
fit_full_FG <- tidycmprsk::crr(Surv(follow_up_time, status_cr) ~ adj_BAA * Sex, data = sub_df_FG)

summary(fit_base_FG)
summary(fit_full_FG)


get_comparison_stats_FG <- function(df, sex_include, cause_code, time_threshold) {
  
  sub_df <- df %>% 
    filter(Sex %in% sex_include, follow_up_time > time_threshold)
  
  sub_df <- sub_df %>%
    mutate(status_cr = case_when(
      is.na(event) | event == 0 ~ "censor",
      event == 1 & cause_lv1 == cause_code ~ "target",
      event == 1 & cause_lv1 != cause_code ~ "competing",
      TRUE ~ "censor"
    )) %>%
    mutate(status_cr = factor(status_cr, levels = c("censor", "target", "competing")))
  
  n_events <- sum(sub_df$status_cr == "target")
  cat("Checking Cause:", cause_code, "| Events found:", n_events, "\n")
  
  if (n_events < 10) return(NULL) 
  
  tryCatch({
    fit_base <- tidycmprsk::crr(Surv(follow_up_time, status_cr) ~ Sex, data = sub_df)
    fit_full <- tidycmprsk::crr(Surv(follow_up_time, status_cr) ~ adj_BAA * Sex, data = sub_df)
    
    c_base <- survival::concordance(fit_base)$concordance
    c_full <- survival::concordance(fit_full)$concordance
    
    p_val <- tidycmprsk::tidy(fit_full) %>% 
      filter(term %in% c("adj_BAA", "adj_BAA:Sex")) %>% 
      pull(p.value) %>% 
      min(na.rm = TRUE)
    
    return(c(Events = n_events, C_Base = c_base, C_Full = c_full, 
             Delta = c_full - c_base, P_Value = p_val))
  }, error = function(e) {
    cat("Model Error for", cause_code, ":", e$message, "\n")
    return(NULL)
  })
}

results_list <- list()

for (cause in target_causes) {
  res_0  <- get_comparison_stats_FG(analysis_df, c(0,1), cause, 0)
  # res_5  <- get_comparison_stats_FG(analysis_df, c(0,1), cause, 5)
  # res_10 <- get_comparison_stats_FG(analysis_df, c(0,1), cause, 10)
  
  safe_extract <- function(res, field, digits = 3) {
    if (is.null(res) || is.na(res[field])) return(NA)
    val <- as.numeric(res[field])
    if (field == "P_Value") return(val) 
    return(round(val, digits))
  }
  print(res_0)

  if (!is.null(res_0)) {
    row_data <- data.frame(
      Cause     = as.character(cause_labels[cause]),
      
      # Full Follow-up
      Events_0  = safe_extract(res_0, "Events", 0),
      C_Base_0  = safe_extract(res_0, "C_Base"),
      C_Full_0  = safe_extract(res_0, "C_Full"),
      Delta_0   = safe_extract(res_0, "Delta"),
      P_LRT_0   = safe_extract(res_0, "P_Value"),
      
      # # > 5 Years
      # Events_5  = safe_extract(res_5, "Events", 0),
      # C_Base_5  = safe_extract(res_5, "C_Base"),
      # C_Full_5  = safe_extract(res_5, "C_Full"),
      # Delta_5   = safe_extract(res_5, "Delta"),
      # P_LRT_5   = safe_extract(res_5, "P_Value"),
      # 
      # # > 10 Years
      # Events_10 = safe_extract(res_10, "Events", 0),
      # C_Base_10 = safe_extract(res_10, "C_Base"),
      # C_Full_10 = safe_extract(res_10, "C_Full"),
      # Delta_10  = safe_extract(res_10, "Delta"),
      # P_LRT_10  = safe_extract(res_10, "P_Value"),
      
      stringsAsFactors = FALSE
    )
    results_list[[cause]] <- row_data
  }
}

final_comparison_table <- do.call(rbind, results_list)
row.names(final_comparison_table) <- NULL

final_comparison_table <- final_comparison_table %>%
  mutate(across(starts_with("P_LRT"), ~ ifelse(is.na(.), "--", fmt_p(.))))

final_comparison_table %>%
  kbl(format = "latex", 
      booktabs = TRUE, 
      caption = "Performance of BAA across Specific Causes of Death and Follow-up Windows (Fine-Gray Sensitivity Analysis)",
      col.names = c("Cause", rep(c("Events", "C-Base", "C-Full", "Delta C", "P-value"), 3))) %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  add_header_above(c(" " = 1,
                     "Full Follow-up" = 5,
                     "> 5 Years" = 5,
                     "> 10 Years" = 5))



######### sex specific ################
results_list_male <- list()

for (cause in target_causes) {
  res_0  <- get_comparison_stats(analysis_df,c(0) ,cause, 0)
  res_5  <- get_comparison_stats(analysis_df,c(0) ,cause, 5)
  res_10 <- get_comparison_stats(analysis_df,c(0),cause, 10)
  
  row_data <- data.frame(
    Cause     = as.character(cause_labels[cause]),
    
    # Full Follow-up
    Events_0  = res_0[["Events"]],
    C_Base_0  = round(as.numeric(res_0[["C_Base.C"]]), 3),
    C_Full_0  = round(as.numeric(res_0[["C_Full.C"]]), 3),
    Delta_0   = round(as.numeric(res_0[["Delta.C"]]), 3),
    P_LRT_0  =  as.numeric(res_0[["P_LRT"]]),
    
    # > 5 Years
    Events_5  = res_5[["Events"]],
    C_Base_5  = round(as.numeric(res_5[["C_Base.C"]]), 3),
    C_Full_5  = round(as.numeric(res_5[["C_Full.C"]]), 3),
    Delta_5   = round(as.numeric(res_5[["Delta.C"]]), 3),
    P_LRT_5  =  as.numeric(res_5[["P_LRT"]]),
    
    # > 10 Years
    Events_10 = res_10[["Events"]],
    C_Base_10 = round(as.numeric(res_10[["C_Base.C"]]), 3),
    C_Full_10 = round(as.numeric(res_10[["C_Full.C"]]), 3),
    Delta_10  = round(as.numeric(res_10[["Delta.C"]]), 3),
    P_LRT_10  =  as.numeric(res_10[["P_LRT"]]),
    
    stringsAsFactors = FALSE
  )
  results_list_male[[cause]] <- row_data
}

final_comparison_table <- do.call(rbind, results_list_male)
row.names(final_comparison_table) <- NULL

fmt_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  out <- ifelse(is.na(p), NA_character_,
                ifelse(p < 1e-3, "<0.001",
                       format(round(p, 3), nsmall = 3)))
  out
}

final_comparison_table <- final_comparison_table %>%
  mutate(
    P_LRT_0  = fmt_p(P_LRT_0),
    P_LRT_5  = fmt_p(P_LRT_5),
    P_LRT_10 = fmt_p(P_LRT_10)
  )

print(final_comparison_table)

final_comparison_table %>%
  kbl(caption = "Performance of BAA across Specific Causes of Death and Follow-up Windows") %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  add_header_above(c(" " = 1,
                     "Full Follow-up" = 5,
                     "> 5 Years" = 5,
                     "> 10 Years" = 5))

#########################################################
# table(analysis_df$cause_lv2)
target_causes_cancer <- c(
  "C15-C26 Malignant neoplasms of digestive organs",
  "C30-C39 Malignant neoplasms of respiratory and intrathoracic organs",
  "C50-C50 Malignant neoplasm of breast",
  "C51-C58 Malignant neoplasms of female genital organs",
  "C60-C63 Malignant neoplasms of male genital organs",
  "C64-C68 Malignant neoplasms of urinary tract",
  "C81-C96 Malignant neoplasms of lymphoid, haematopoietic and related tissue"
)

cause_names_cancer <- c(
  "C15-C26 Malignant neoplasms of digestive organs" =
    "C15–C26: Digestive cancers",
  "C30-C39 Malignant neoplasms of respiratory and intrathoracic organs" =
    "C30–C39: Respiratory/intrathoracic cancers",
  "C50-C50 Malignant neoplasm of breast" =
    "C50: Breast cancer",
  "C51-C58 Malignant neoplasms of female genital organs" =
    "C51–C58: Female genital cancers",
  "C60-C63 Malignant neoplasms of male genital organs" =
    "C60–C63: Male genital cancers",
  "C64-C68 Malignant neoplasms of urinary tract" =
    "C64–C68: Urinary tract cancers",
  "C81-C96 Malignant neoplasms of lymphoid, haematopoietic and related tissue" =
    "C81–C96: Haematological malignancies"
)

# 2. Define comparison calculation function
results_list <- list()

for (cause in target_causes_cancer) {
  res_0  <- get_comparison_stats_lv2(analysis_df,c(0,1) ,cause, 0)
  res_5  <- get_comparison_stats_lv2(analysis_df,c(0,1) ,cause, 5)
  res_10 <- get_comparison_stats_lv2(analysis_df,c(0,1),cause, 10)
  
  row_data <- data.frame(
    Cause     = as.character(cause_names_cancer[cause]),
    
    # Full Follow-up
    Events_0  = res_0[["Events"]],
    C_Base_0  = round(as.numeric(res_0[["C_Base.C"]]), 3),
    C_Full_0  = round(as.numeric(res_0[["C_Full.C"]]), 3),
    Delta_0   = round(as.numeric(res_0[["Delta.C"]]), 3),
    P_LRT_0  =  as.numeric(res_0[["P_LRT"]]),
    
    # > 5 Years
    Events_5  = res_5[["Events"]],
    C_Base_5  = round(as.numeric(res_5[["C_Base.C"]]), 3),
    C_Full_5  = round(as.numeric(res_5[["C_Full.C"]]), 3),
    Delta_5   = round(as.numeric(res_5[["Delta.C"]]), 3),
    P_LRT_5  =  as.numeric(res_5[["P_LRT"]]),
    
    # > 10 Years
    Events_10 = res_10[["Events"]],
    C_Base_10 = round(as.numeric(res_10[["C_Base.C"]]), 3),
    C_Full_10 = round(as.numeric(res_10[["C_Full.C"]]), 3),
    Delta_10  = round(as.numeric(res_10[["Delta.C"]]), 3),
    P_LRT_10  =  as.numeric(res_10[["P_LRT"]]),
    
    stringsAsFactors = FALSE
  )
  results_list[[cause]] <- row_data
}

final_comparison_table <- do.call(rbind, results_list)
row.names(final_comparison_table) <- NULL

final_comparison_table <- final_comparison_table %>%
  mutate(
    P_LRT_0  = fmt_p(P_LRT_0),
    P_LRT_5  = fmt_p(P_LRT_5),
    P_LRT_10 = fmt_p(P_LRT_10)
  )

print(final_comparison_table)

final_comparison_table %>%
  kbl(caption = "Performance of BAA across Specific Causes of Death and Follow-up Windows") %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  add_header_above(c(" " = 1,
                     "Full Follow-up" = 5,
                     "> 5 Years" = 5,
                     "> 10 Years" = 5))

##################################################

df_plot <- analysis_df %>%
  mutate(Sex = factor(Sex, levels = c(0, 1), labels = c("Female", "Male"))) %>%
  group_by(Sex) %>%
  sample_n(3000) %>%
  ungroup()
hist(df_plot$baseline_age)
Testset_BA_CA_corre_plot<-ggplot(df_plot, aes(x =baseline_age , y = BA, 
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
    title = "Testset cohort (2000 random selected samples)"
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )
Testset_BA_CA_corre_plot

ggsave("Results plots/Testset_BA_CA_corre_plot.png", Testset_BA_CA_corre_plot, width = 9, height = 6, dpi = 600)

dev.off()


##############################################
##############################################
###############################################
###############################################

################## External TEST ###################
# new_person_df: External new data (at least containing Sex, age0 and as many columns as possible from feature_cols)
# testing_3plusmeas_raw<-read.csv('deep_cleaned_3plusmeas_testing_ins0123_lfm.csv',header=TRUE)
training_2plusmeas_raw<-read.csv('deep_cleaned_2plusmeas_ins023_lfm.csv',header=TRUE)

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

dim(test_data)
head(test_data)
length(unique(test_data$Participant.ID))

table(table(test_data$Participant.ID))

intersect(dat_aft$Participant.ID,test_data$Participant.ID)

test_data <- test_data %>%
  mutate(age0 = as.numeric(`Age.when.attended.assessment.centre`))

test_data <- test_data %>%
  dplyr::group_by(Participant.ID) %>%
  dplyr::filter(dplyr::n() >= 2) %>%
  dplyr::ungroup()

dim(test_data)
table(table(test_data$Participant.ID))

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

testset_death_data<-read.csv('death_birth_info/test_set_death_cause.csv',header=TRUE)
trainset_death_data<-read.csv('death_birth_info/train_set_death_cause.csv',header=TRUE)

length(!is.na(trainset_death_data$Date.of.death...Instance.0)==FALSE)

# testset_death_data <- subset(
#   testset_death_data,
#   !grepl("^[C]", Underlying..primary..cause.of.death..ICD10...Instance.0 %>% as.character())
# )

testset_death_data <- subset(
  trainset_death_data,
  !grepl("^[VWXY]", Underlying..primary..cause.of.death..ICD10...Instance.0 %>% as.character())
)

# testset_death_data$death_summary_code = substr(testset_death_data$Underlying..primary..cause.of.death..ICD10...Instance.0,1,1)

head(testset_death_data)
dim(testset_death_data)
death_df <- testset_death_data %>%
  select(Participant.ID, 
         Year.of.birth, 
         Month.of.birth,
         starts_with("Date.of.death"),
         death_summary_code) %>%
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
    # Key point here: retain death_summary_code
    death_summary_code = first(death_summary_code), 
    death_age_recorded = min(death_age_calc, na.rm = TRUE),
    death_date_final = min(death_date, na.rm = TRUE),
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

df_long <- test_data_for_survival %>%
  left_join(death_df, by = "Participant.ID")

head(test_data_for_survival)

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
  left_join(death_df %>% select(Participant.ID, birth_date, death_age, status,death_summary_code), by = "Participant.ID") %>%
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
table(df_tv_filter$death_summary_code)
summary(df_tv_filter$Age.at.assessment)
length(unique(df_tv_filter$Participant.ID))
intersect(df_tv_filter$Participant.ID,dat_aft$Participant.ID)

cox_baa_sex <- coxph(
  Surv(start, stop, event) ~ BAA + Sex 
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
  data = df_tv,
  cluster = Participant.ID
)

summary(cox_sex_only)

# full model
cox_baa_sex <- coxph(Surv(start, stop, event) ~ Sex + BAA + Sex:BAA,
                     data = df_tv, cluster = Participant.ID)

# joint Wald test for (BAA, BAA:Sex)
b <- coef(cox_baa_sex)[c("BAA", "Sex:BAA")]
V <- vcov(cox_baa_sex)[c("BAA", "Sex:BAA"), c("BAA", "Sex:BAA")]

W <- as.numeric(t(b) %*% solve(V) %*% b)        # chi-square stat
p <- pchisq(W, df = length(b), lower.tail = FALSE)

c(Wald_chisq = W, df = length(b), p_value = p)

df_late <- df_tv %>% filter(stop > (start + 5))  # Adjust this based on your start/stop time units
cox_late <- coxph(Surv(start, stop, event) ~ BAA + Sex, data=df_late, cluster=Participant.ID)
summary(cox_late)

cox_late_sexonly <- coxph(Surv(start, stop, event) ~ Sex, data=df_late, cluster=Participant.ID)
summary(cox_late_sexonly)

hist(df_tv$stop - df_tv$start)
summary(df_tv$stop - df_tv$start)
# 
# ###################################################################
# # 1. Prepare cause of death list
# target_causes <- c("All-Cause", "Non-C")
# cause_names <- c(
#   "All-Cause" = "All-Cause Mortality",
#   "Non-C" = "Non-C: Non-Neoplasms"
# )
# 
# # 2. Define comparison calculation function
# get_comparison_stats <- function(df, sex_include, cause_code, time_threshold) {
#   
#   # A. Select subset based on Participant ID
#   # Logic: Find individuals whose total follow-up time (max stop - min start) exceeds the threshold
#   valid_ids <- df %>%
#     group_by(Participant.ID) %>%
#     summarise(
#       baseline_age = min(start, na.rm = TRUE),
#       max_followup_age = max(stop, na.rm = TRUE),
#       total_duration = max_followup_age - baseline_age
#     ) %>%
#     filter(total_duration > time_threshold) %>%
#     pull(Participant.ID)
#   
#   # B. Extract subsequent records for these individuals
#   sub_df <- df %>%
#     filter(Participant.ID %in% valid_ids, Sex %in% sex_include) %>%
#     group_by(Participant.ID) %>%
#     # 只保留基线 + threshold 之后的时间段
#     mutate(baseline_age = min(start, na.rm = TRUE)) %>%
#     filter(stop > (baseline_age + time_threshold)) %>%
#     # Recalibrate start time: if start is still before threshold, truncate it to threshold
#     mutate(start = pmax(start, baseline_age + time_threshold)) %>%
#     ungroup()
#   
#   if (nrow(sub_df) == 0) return(NULL)
#   
#   # C. 定义 Outcome (Cause-specific)
#   if (cause_code == "All-Cause") {
#     sub_df$outcome_event <- sub_df$event
#   } else {
#     # 只有当 event=1 且 死因代码不是 'C' 时才计为事件
#     sub_df$outcome_event <- ifelse(sub_df$event == 1 & !(sub_df$death_summary_code %in% c("C")), 1, 0)
#   }
#   
#   # Count events (deduplicated count; each ID contributes at most 1 event in this time period)
#   n_events <- sum(sub_df$outcome_event, na.rm = TRUE)
#   cat("Cause:", cause_code, "| Threshold:", time_threshold, "| Events:", n_events, "\n")
#   
#   if (is.na(n_events) || n_events < 10) return(NULL) 
#   
#   # D. Run models
#   tryCatch({
#     fit_base <- coxph(Surv(start, stop, outcome_event) ~ Sex, data = sub_df, cluster = Participant.ID)
#     c_base <- summary(fit_base)$concordance[1]
#     
#     fit_full <- coxph(Surv(start, stop, outcome_event) ~ BAA + Sex, data = sub_df, cluster = Participant.ID)
#     c_full <- summary(fit_full)$concordance[1]
#     
#     delta_c <- c_full - c_base
#     
#     return(c(Events = n_events, C_Base = c_base, C_Full = c_full, Delta = delta_c))
#   }, error = function(e) return(NULL))
# }
# 
# # Subsequent results_list loop logic remains unchanged
# 
# results_list <- list()
# 
# for (cause in target_causes) {
#   # Calculate for three time points
#   res_0  <- get_comparison_stats(df_tv_filter, c(0,1), cause, 0)
#   print(res_0)
#   res_5  <- get_comparison_stats(df_tv_filter, c(0,1), cause, 5)
#   res_10 <- get_comparison_stats(df_tv_filter, c(0,1), cause, 10)
#   
#   # Check if res is null (for cases with too few events)
#   if(is.null(res_0)) next
#   
#   row_data <- data.frame(
#     Cause = cause_names[cause],
#     # Full Follow-up (0)
#     Events_0 = res_0["Events"],
#     C_Base_0 = round(res_0["C_Base.C"], 3),
#     C_Full_0 = round(res_0["C_Full.C"], 3),
#     Delta_0  = round(res_0["Delta.C"], 3),
#     
#     # > 10 Years
#     Events_10 = ifelse(is.null(res_5), NA, res_5["Events"]),
#     C_Base_10 = ifelse(is.null(res_5), NA, round(res_5["C_Base.C"], 3)),
#     C_Full_10 = ifelse(is.null(res_5), NA, round(res_5["C_Full.C"], 3)),
#     Delta_10  = ifelse(is.null(res_5), NA, round(res_5["Delta.C"], 3)),
#     
#     # > 15 Years
#     Events_15 = ifelse(is.null(res_10), NA, res_10["Events"]),
#     C_Base_15 = ifelse(is.null(res_10), NA, round(res_10["C_Base.C"], 3)),
#     C_Full_15 = ifelse(is.null(res_10), NA, round(res_10["C_Full.C"], 3)),
#     Delta_15  = ifelse(is.null(res_10), NA, round(res_10["Delta.C"], 3))
#   )
#   results_list[[cause]] <- row_data
# }
# 
# final_comparison_table <- do.call(rbind, results_list)
# row.names(final_comparison_table) <- NULL
# 
# print(final_comparison_table)
# 
# final_comparison_table %>%
#   kbl(caption = "Performance of BAA across Specific Causes of Death and Follow-up Windows") %>%
#   kable_classic(full_width = F, html_font = "Arial") %>%
#   add_header_above(c(" " = 1, "Full Follow-up" = 4, "> 10 Years" = 4, "> 15 Years" = 4))
# 
# dev.off()










































































########################################################
########################################################


######### Train ###################################
dat <- dat_aft

exclude <- c(
  "Participant.ID", "Instance",
  "Age.when.attended.assessment.centre", 
  "Age.at.death...Instance.0", "Date.of.death...Instance.0",
  "Year.of.birth", "Month.of.birth",
  "birth_date", "censor_age", "death_date", "death_age",
  "death_age_from_date", "death_age_final", "age0", "yob",'mob',
  "time_years", "status"
)

feature_cols <- dat %>%
  select(where(is.numeric)) %>%
  select(-any_of(exclude)) %>%
  names()

# Feature matrix and labels
X <- data.matrix(dat[, feature_cols, drop = FALSE])
time   <- as.numeric(dat$time_years)
hist(time)
status <- as.integer(dat$status)

# sanity check (recommended to run and check)
stopifnot(nrow(X) == length(time), nrow(X) == length(status))

dall <- xgb.DMatrix(X)
xgboost::setinfo(dall, "label_lower_bound", time)
xgboost::setinfo(dall, "label_upper_bound", ifelse(status == 1L, time, Inf))

##############################
set.seed(2026)

idx1 <- which(status == 1L)  # deaths
idx0 <- which(status == 0L)  # censored

split_one <- function(idx, p_train = 0.6, p_valid = 0.3) {
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

train_idx <- c(s1$train, s0$train)
valid_idx <- c(s1$valid, s0$valid)
test_idx  <- c(s1$test,  s0$test)

# Shuffle
train_idx <- sample(train_idx)
valid_idx <- sample(valid_idx)
test_idx  <- sample(test_idx)

# Check if event proportions in each set are close to the overall proportion
prop.table(table(status))
prop.table(table(status[train_idx]))
prop.table(table(status[valid_idx]))
prop.table(table(status[test_idx]))

make_aft_dmatrix <- function(Xm, t, s) {
  d <- xgb.DMatrix(Xm)
  xgboost::setinfo(d, "label_lower_bound", t)
  xgboost::setinfo(d, "label_upper_bound", ifelse(s == 1L, t, Inf))
  d
}

dtrain <- make_aft_dmatrix(X[train_idx, , drop = FALSE], time[train_idx], status[train_idx])
dvalid <- make_aft_dmatrix(X[valid_idx, , drop = FALSE], time[valid_idx], status[valid_idx])
dtest  <- make_aft_dmatrix(X[test_idx,  , drop = FALSE], time[test_idx],  status[test_idx])

params <- list(
  booster = "gbtree",
  objective = "survival:aft",
  eval_metric = "aft-nloglik",
  eta = 0.05,
  max_depth = 4,
  min_child_weight = 10,
  subsample = 0.7,
  colsample_bytree = 0.7,
  aft_loss_distribution = "normal",
  aft_loss_distribution_scale = 1.0
)

bst <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 5000,
  watchlist = list(train = dtrain, valid = dvalid),
  early_stopping_rounds = 100,
  verbose = 1
)

pred_test <- predict(bst, X[test_idx, , drop = FALSE])
rmse_event_only <- function(pred, t, s) {
  ii <- which(s == 1L)
  sqrt(mean((pred[ii] - t[ii])^2))
}

rmse_event_only(pred_test, time[test_idx], status[test_idx])
summary(time[test_idx])
summary(pred_test)

summary(dat_aft$time_years[dat_aft$status == 1])
summary(dat_aft$time_years[dat_aft$status == 0])


pred_time   <- predict(bst, X[test_idx, , drop = FALSE])
pred_margin <- predict(bst, X[test_idx, , drop = FALSE], outputmargin = TRUE)

summary(pred_time)
summary(pred_margin)
summary(exp(pred_margin))

cor(pred_time, exp(pred_margin), use="complete.obs")

library(survival)

cindex <- function(t, s, pred_margin) {
  fit <- coxph(Surv(t, s) ~ I(-pred_margin))
  as.numeric(summary(fit)$concordance[1])
}

cindex(time[test_idx], status[test_idx], pred_margin)
"Age.when.attended.assessment.centre" %in% feature_cols
base_fit <- coxph(Surv(time[test_idx], status[test_idx]) ~ dat$age0[test_idx] + dat$Sex[test_idx])
summary(base_fit)$concordance[1]
##########################


set.seed(1)

n_show <- 10000
idx <- sample.int(nrow(X), min(n_show, nrow(X)))

X_show <- X[idx, , drop = FALSE]

shap_contrib <- predict(bst, X_show, predcontrib = TRUE, approxcontrib = FALSE)

shap_mat <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]

sv <- shapviz(shap_mat, X = X_show)
sv_importance(sv, kind = "beeswarm",max_display = 30,show_numbers = TRUE)
sv_dependence(sv, "Body.mass.index..BMI.", color_var = 'auto')
sv_dependence(sv, "Body.mass.index..BMI.", color_var = 'Age.when.attended.assessment.centre')

sv$X$Age.when.attended.assessment.centre





##############











