## linear regression model on top 20 variables with highest coeficiente in LASSO 
## (absolute values)
library(dplyr)
library(ggplot2)

## load LASSO coeficients
coefs <- read.csv("../data/lasso_only_numeric.csv", 
                  stringsAsFactors = FALSE)
## drop intercept coeficient
coefs <- filter(coefs, coefficient != "(Intercept)")

top20vars <- coefs$coefficient[1:20]


## load cleaned data
vih_data <- read.csv("../data/cleandata.csv", stringsAsFactors = FALSE)

## selecting CD4 increase variable as output for lm
vih_data$"output" <- vih_data$Delta_CD4_year1

## processing data for lasso
input <- vih_data[, names(vih_data) %in% c(top20vars, "output")]

## perform linear model
lModel <- lm(output~., data = input)
preds <- predict(lModel, newdata = input)

## process data
res <- data.frame("value"= vih_data$output, "prediction"= preds)

theme_set(theme_light())
g <- ggplot(res, aes(x=value, y=prediction)) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="True Value", y="Prediction") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "lm")
plot(g)

