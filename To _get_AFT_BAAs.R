library(xgboost)
library(dplyr)
library(tidyr)
library(readr)


### prepare predict function ###
source('Functions/AFT analysis funcs.R')

### input xgboost-aft trained bundles ###
bundle_list_24pm <- readRDS("UKB_XGBOOST_AFT_clock_of_24_physical_measurements.rds")
bundle_list_25cb <- readRDS("UKB_XGBOOST_AFT_clock_of_25_blood_biomarkers.rds")
bundle_list_10pm <- readRDS("UKB_XGBOOST_AFT_clock_of_10_physical_measurements.rds")

### sex specifc models: 0 for female and 1 for male
names(bundle_list_24pm)
names(bundle_list_24pm$'0')

### check the features each model needs 
bundle_list_10pm$'0'$feature_cols
bundle_list_24pm$'0'$feature_cols
bundle_list_25cb$'0'$feature_cols



######################################################################
# So new dataframe with the participants need not only the features but also 
# one column 'Sex'for gender (0 or 1, for female and male respectively) 
# and one column 'age0' for age when the measurements were collected..
#######################################################################

# Recommend age from 37 to 82 as the model was trained on participants in this
# age range, with minumum age of 27 for male and 30 for female
# and maximum age of 92 for both, as this is the extrapolation limit of the
# model. In other word, the model predicted BA would be between 27 and 92 for
# male and between 30 and 92 for female.


#######################################################################################
physical_demo<- read_csv("Demo dataset/Physical_measurements_demo.csv")

physical_base_BAAs <- score_BAA_batch_Gompertz_fast_hybrid(physical_demo, bundle_list_10pm)
# PC_BA is the predicted biological age, and PC_BAA is the predicted biological age acceleration (BAA)
head(physical_base_BAAs)

Physical_enhanced_BAAs <- score_BAA_batch_Gompertz_fast_hybrid(physical_demo, bundle_list_24pm)
head(Physical_enhanced_BAAs)

#########################################################################################
blood_demo<- read_csv("Demo dataset/Blood_demo.csv")

AFT_blood_BAAs <- score_BAA_batch_Gompertz_fast_hybrid(blood_demo, bundle_list_25cb)
head(AFT_blood_BAAs,10)




