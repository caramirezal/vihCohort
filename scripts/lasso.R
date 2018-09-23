## lasso implementation on VIH patients cohort data
library(glmnet)
library(dplyr)
library(ggplot2)

vih_data <- read.csv("../data/cleandata.csv", stringsAsFactors = FALSE)
str(vih_data)
## No NA values are presented in data
dim(vih_data[!complete.cases(vih_data),])

## processing data for lasso
input <- vih_data[, ! sapply(vih_data, function(x) class(x)=="character") ]
input <- select(input, -Delta_CD4_year1)
input <- select(input, -CD4_S0)
input <- select(input, -CD4_S52)
input <- model.matrix(~., data = input, 
                      contrasts.arg = sapply(input, is.factor))
write.table(input, "../data/model_matrix.tsv", sep = "\t")
str(input)
output <- vih_data$Delta_CD4_year1

######################################################################################
## leave-one-out validation

## vector to store predictions
res <- numeric(nrow(input))

## matrix to store lasso coefficients
lasso_coefs <- matrix(0, nrow(input), ncol(input)+1)

## perform leave-one-out validation
for (i in 1:nrow(input)) {
        lambda.cv <- cv.glmnet(x=input[-i,], y = output[-i])$lambda.1se
        lasso <- glmnet(x=input[-i,], y = output[-i], lambda = lambda.cv)
        prediction <- predict(lasso, newx = input, type = "response", s = lambda.cv)
        res[i] <- prediction[i]
        lasso_coefs[i,] <- as.vector(coef(lasso)) 
}

## plot predicted vs target values
validation <- data.frame("lasso_prediction"=res,
                         "values"=output)
theme_set(theme_light())
p <- ggplot(validation, aes(x=values, y=lasso_prediction)) + 
        geom_point(colour="steelblue", size= 2.5) + 
        geom_abline(slope = 1,colour="red",size=1) +
        labs(x="Delta TCD4 values", y="LASSO_prediction") +
        theme(text = element_text(face="bold", size = 18))
plot(p)
jpeg("../figures/lasso.jpeg")
plot(p)
dev.off()

## calculate mean coefficient values
colnames(lasso_coefs) <- rownames(coef(lasso))
mean_coef <- apply(lasso_coefs, 2, mean)
sd_coef <- apply(lasso_coefs, 2, sd)
summary_coefs <- data.frame("coefficient"=colnames(lasso_coefs),
                            "mean"=mean_coef,
                            "sd"=sd_coef) %>%
                        arrange(desc(abs(mean)))
write.csv(summary_coefs, "../data/lasso_only_numeric.csv", row.names=FALSE)
