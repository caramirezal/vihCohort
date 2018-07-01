library(dplyr)
library(randomForest)

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

## drop observations with no output
data.s <- data[!is.na(data$Delta_CD4_year1),] 
dim(data.s)

## Inspect the number of variables with at least some threshold number 
## of complete cases
numKeepVar <- function(df) {
        int <- 0:nrow(df)
        result <- rep(0,length(int))
        for (i in int) {
                result[i] <- sum(sapply(df,function(x) i <= sum(!is.na(x))))
        }
        result
}
keepVar <- numKeepVar(data.s)
plot(0:nrow(data.s),
     keepVar,
     xlab = "Complete cases",
     ylab = "Variables", 
     type = "l",
     col="steelblue",
     lwd=2)



## keep variables with more than 2/3 complete cases of the total of observations
th <- (2/3)*nrow(data.s)
keep <- sapply(data.s,function(x) th <= sum(!is.na(x))) 
data.p <- data.s[,keep]
dim(data.p)

## drop near zero variance columns
data.p <- data.p[,-nearZeroVar(data.p)]
dim(data.p)

output <- data.p$Delta_CD4_year1

## list of removed variables
drop <- c("Expediente",
          "DeltaCD4W52atleast150",
          "DeltaCD4_W52",
          "CD4TCabove200_W052",
          "CD4TCabove350_W052",
          "CD4TCabove200_W052",
          "\xff\xfeNumero_consecutivo",
          "CD4_S52",
          "Delta_CD4_year1")

## removing list of variables
input <- data.p[,!(colnames(data.p)%in%drop)]




