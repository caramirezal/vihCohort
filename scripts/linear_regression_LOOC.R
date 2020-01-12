## Linear regression cross validation
library(dplyr)
library(ggplot2)

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

y <- input$output
input <- select(input, -output)

preds <- numeric(nrow(input))
## perform LOOC with linear model
for (i in 1:nrow(input)){
        ## Learning over n-1 samples
        lModel <- lm(y[-i]~., data = input[-i, ])
        preds[i] <- predict(lModel, newdata = input)[i]
}

## plot without coloured IRIS patients
res <- data.frame("value"= vih_data$output, "prediction"= preds)
theme_set(theme_light())
g <- ggplot(res, aes(x=prediction, y=value)) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.9) + 
        labs(x="Predicted value", y="Del TCD4 values") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "lm")
plot(g)

## plot with coloured IRIS patients
with_without_IRIS <- sapply(vih_data$with_without_IRIS, 
                            function(x) ifelse(x=="with", "IRIS", "no symptoms"))
res <- data.frame("value"= vih_data$output, 
                  "prediction"= preds, 
                  "iris"=with_without_IRIS)
theme_set(theme_light())
g <- ggplot(res, aes(x=prediction, y=value)) +     ## adjust point color manually
        geom_point() + 
        labs(x="Predicted value", y="Delta TCD4 values")   +     ## remove legend title
        geom_smooth(method = 'lm')

plot(g)

lin_mod <- lm(value~prediction, res)


#####################################################################################
## CD4 T cells are not correlated to Delta CD4
fit <- lm(CocCD4_CD8_S52 ~ CD4_S52, data = vih_data)
summary(fit)
## Adding CD4 prediction
fit.df <- mutate(vih_data, 
                 prediction = predict(fit, newdata = vih_data)
)

theme_set(theme_light())
ggplot(data = vih_data,
       aes(x = CocCD4_CD8_S52, y = CD4_S52)) + 
        geom_point() + geom_smooth(method = 'lm')


