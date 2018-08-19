## lasso implementation on VIH patients cohort data
library(glmnet)
library(dplyr)

vih_data <- read.csv("../data/cleandata.csv", stringsAsFactors = FALSE)
str(vih_data)
## No NA values are presented in data
dim(vih_data[!complete.cases(vih_data),])

input <- select(vih_data, -Delta_CD4_year1)
input <- select(input, -CD4_S0)
input <- select(input, -CD4_S52)
input <- model.matrix(~., data = input, 
                      contrasts.arg = sapply(input, is.factor))
head(input)
output <-vih_data$Delta_CD4_year1

res <- numeric(nrow(input))
for (i in 1:nrow(input)) {
        lambda.cv <- cv.glmnet(x=input[-i,], y = output[-i])$lambda.1se
        lasso <- glmnet(x=input[-i,], y = output[-i], lambda = lambda.cv)
        prediction <- predict(lasso, newx = input, type = "response", s = lambda.cv)
        res[i] <- prediction[i]
}
coef(lasso)
plot(output,res, pch=20, col="steelblue", cex=2)
abline(0,1, col="red", lwd=2.5)
