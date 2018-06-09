## Machine learning on the CD4 T cells reconstitution in a cohort of VIH patients
library(dplyr)
library(reshape2)
library(ggplot2)
library(randomForest)
library(glmnet)        ## lasso and ridge regression
library(caret)
library(party)

## function for calculating r^2  as 1 - rss/tss
r_sq <- function(observed,predicted) {
        rss <- sum((observed-predicted)^2)
        tss <- sum((observed-mean(observed))^2)
        return( 1 - (rss/tss) )
}

#################################################################################
## Data Preprocessing                                                          ##
################################################################################

## reading data
data <- read.csv("../data/TablaLPS_IL18_20180217.csv")
head(data)


## preprocessing raw data
## drop categorical data
data <- data[,sapply(data, is.numeric)]

## list of removed variables
drop <- c("Expediente",
          "DeltaCD4W52atleast150",
          "DeltaCD4_W52",
          "CD4TCabove200_W052",
          "CD4TCabove350_W052",
          "Numero_consecutivo")

## removing list of variables
input <- data[,!(colnames(data)%in%drop)]

## input definition
#input <- as.matrix(input)

## output definition
input$DeltaCD4_W52 <- data$DeltaCD4_W52
colnames(input)

###############################################################################
initialTime <- Sys.time()
forest <- cforest(DeltaCD4_W52~.,data=input)
forest.pred <- predict(forest,newdata=input)
finalTime <- Sys.time() - initialTime
finalTime


