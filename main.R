library(dplyr)
library(glmnet)
library(caret)

########################################################################
## initial exploratory lasso analysis of raw data

data <- read.csv("data/TablaLPS_IL18_20180217.csv")
