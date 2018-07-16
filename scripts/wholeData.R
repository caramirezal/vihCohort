library(dplyr)
library(ggplot2)

## reading data as character variables for preprocessing
data <- read.csv("../data/TablaCompletaCohorteINER_C08-04_20170630.csv",
                 header = TRUE,
                 colClasses = "character",
                 skipNul = TRUE,
                 check.names = FALSE)

data$"eventosClinicos" <- data[,163] 
data <- data[,-163]

## assign numeric class if there is at least one entry that
## can be coerced to numeric to the dataframe column
classes <- sapply(data, function(x) sum(!is.na(as.numeric(x))))
classes[classes>10] <- "numeric"
for (i in 1:ncol(data)) {
        if ( classes[i] == "numeric" ) {
                data[,i] <- as.numeric(data[,i])
        } 
}


###################################################################################
## reading top 20 coeficients from lasso
coeficients <- read.csv("../data/lassoCoefFormtatted.csv")

## selecting top 20 coefficients from lasso 
sharedTop20 <- coeficients$X[1:20][sapply(coeficients$X[1:20],function(x) x %in% names(data) )]
data.filtered <- data[,c(as.character(sharedTop20),"Delta_CD4_year1")]

## selecting complete cases 
data.complete <- data.filtered[complete.cases(data.filtered),]
dim(data.complete) 

## vector to store prediction values for validation
results <- numeric(nrow(data.complete))

#################################################################################
## Linear regression Model

## running Leave-One-Out (LOOC) Validation
for (i in 1:nrow(data.complete)) {
        model <- lm(Delta_CD4_year1~.,data=data.complete[-i,])
        results[i] <- predict(model,newdata = data.complete[i,])
}

mss <- mean((data.complete$Delta_CD4_year1-results)^2)
mse <- sqrt(mss)
mse

## plot data
theme_set(theme_light())
linReg <- data.frame("observed"=data.complete$Delta_CD4_year1,
                     "predicted"=results)
g <- ggplot(data = linReg, aes(x=observed,y=predicted))
g <- g + geom_point(colour="steelblue")
g <- g + geom_smooth(method = "lm",colour="red")
g <- g + theme(text=element_text(size = 22,face="bold"))
g <- g + labs(x="Observed",y="Predicted")
plot(g)