library(dplyr)

data <- read.csv("data/TablaLPS_IL18_20180217.csv")

data <- data[,sapply(data, is.numeric)]
data <- select(data,-Expediente)
input <- as.matrix(select(data,-DeltaCD4_W52))
response <- data$DeltaCD4_W52

lasso <- cv.glmnet(input,response)
plot(lasso)
pred.vals <- predict(lasso,
                     newx = input,
                     s = lasso$lambda.min)

plot(response,pred.vals)
abline(0,1,col="red")

coef(lasso)

