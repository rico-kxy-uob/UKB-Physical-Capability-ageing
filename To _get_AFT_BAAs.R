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
bundle_list_24pm$'0'$feature_cols
bundle_list_25cb$'0'$feature_cols
bundle_list_25cb$'0'$feature_cols




