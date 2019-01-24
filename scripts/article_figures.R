## Construction of the figures for the article
library(glmnet)
library(dplyr)
library(ggplot2)
library(reshape2)


##################################################################################
## data loading and processing

vih_data <- read.csv("../data/cleandata.csv", stringsAsFactors = FALSE)
str(vih_data)
names(vih_data)


## reading table description meta data
description <- read.csv("../data/Desripción de columnas.csv",
                        stringsAsFactors = FALSE,header = FALSE,
                        nrows = 10)
str(description)

## loading data
description <- read.csv("../data/Desripción de columnas.csv",
                        stringsAsFactors = FALSE,
                        header = FALSE)
description <- select(description,V1:V2)
## ""
description$V2 <- gsub(" ","",description$V2)
str(description)
dim(description)
tail(description)
## second column in this table contains metadata
## for the selection of meaningful variables
keepVars <- logical(nrow(description))
## drop column specifed as No
keepVars <- ifelse(description$V2=="Sí" | description$V2=="?"
                   ,TRUE,FALSE)
description$V2 <- keepVars
str(description)
## keeping meaningful var names
meaningfulVars <- description$V1[keepVars]

## prereading data
tabla <- read.csv("../data/TablaCompletaCohorteINER_C08-04_20180810.csv",
                  stringsAsFactors = FALSE)
str(tabla)
dim(tabla)

## drop log cv redundant with cv values
tabla.p <- tabla[,-grep("LogCV_",colnames(tabla))]
## drop cv values with values < 50 == NA at week S052
tabla.p <- tabla.p[!is.na(tabla$CV_S052),]
## Drop non meaningful data specified in metadata
tabla.p <- tabla.p[, colnames(tabla.p) %in% meaningfulVars]
dim(tabla.p)
str(tabla.p)

## drop samples without output delta CD4
tabla.ov <- tabla.p[!is.na(tabla.p$Delta_CD4_year1),]
dim(tabla.ov)

###################################################################################
## 1. Viral load plot

## selecting viral loads
cv <- tabla.ov[, grep("^CV_S", names(tabla.ov), value = TRUE)]
names(cv) <- gsub("CV_S", "", names(cv))
cv.m <- melt(cv)
cv.m <- cv.m[complete.cases(cv.m),]
cv.m <- mutate(cv.m, variable=as.integer(as.character(variable)))

theme_set(theme_light())
p <- ggplot(cv.m, aes(x=variable, y=log10(value))) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="Week", y=" Log_10(HIV RNA copies/mL)") +
        theme(text = element_text(face="bold", size = 18)) +
        scale_x_continuous(breaks=c(0,4,8,12,24,39,52,104)) +     ## change ticks interval labels
        geom_smooth()
plot(p)

##################################################################################
## 2. Test binarize IRIS

vih_data <- read.csv("../data/cleandata.csv", stringsAsFactors = FALSE)
str(vih_data)
## No NA values are presented in data
dim(vih_data[!complete.cases(vih_data),])

## processing data for lasso
vih_data$with_without_IRIS <- sapply(vih_data$with_without_IRIS, 
                                     function(x) ifelse(x=="with", 1, 0))
input <- vih_data[, ! sapply(vih_data, function(x) class(x)=="character") ]
input <- select(input, -Delta_CD4_year1)
input <- select(input, -CD4_S0)
input <- select(input, -CD4_S52)
input <- as.matrix(input)
str(input)
head(input)
write.table(input, "../data/model_matrix_plus_iris.tsv", sep = "\t")
str(input)
output <- vih_data$Delta_CD4_year1

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
        labs(x="Delta TCD4 values", y="LASSO") +
        theme(text = element_text(face="bold", size = 18))
plot(p)
jpeg("../figures/lasso_plus_iris.jpeg")
plot(p)
dev.off()

## calculate mean coefficient values
colnames(lasso_coefs) <- rownames(coef(lasso))
mean_coef <- apply(lasso_coefs, 2, mean)
sd_coef <- apply(lasso_coefs, 2, sd)
summary_coefs <- data.frame("coefficient"=colnames(lasso_coefs),
                            "mean"=mean_coef,
                            "sd"=sd_coef)     %>%
                                      arrange(desc(abs(mean)))
write.csv(summary_coefs, "../data/lasso_only_numeric_plus_iris.csv", row.names=FALSE)

#####################################################################################
## 3. Linear regression model on top 20 variables with highest coeficiente in LASSO 
## (absolute values)

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

## plot without coloured IRIS patients
res <- data.frame("value"= vih_data$output, "prediction"= preds)
theme_set(theme_light())
g <- ggplot(res, aes(x=value, y=prediction)) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.9) + 
        labs(x="True Value", y="Prediction") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "lm")
plot(g)

## plot with coloured IRIS patients
with_without_IRIS <- sapply(vih_data$with_without_IRIS, 
                            function(x) ifelse(x=="with", "IRIS", "no symptoms"))
res <- data.frame("value"= vih_data$output, "prediction"= preds, "iris"=with_without_IRIS)
theme_set(theme_light())
g <- ggplot(res, aes(x=value, y=prediction, color=as.factor(iris))) + 
        scale_color_manual(values = c("green", "black")) +    ## adjust point color manually
        geom_point(size= 2.5) + 
        labs(x="True Value", y="Prediction") +
        theme(text = element_text(face="bold", size = 18)) +
        theme(legend.position = c(0.16, 0.85)) +            ## change legend position
        theme(legend.title = element_blank())  +           ## remove legend title
        geom_abline(slope = 1, size=1, colour="red")

plot(g)



summary(lModel)

##################################################################################
## 4. CD4 T dynamic

## selecting TCD4 columns
cd4_names <- grep("CD4_S", names(vih_data), value = TRUE)
cd4_names

## subsetting data and setting data
cd4 <- select(vih_data, cd4_names)
names(cd4) <- c("0", "8", "12", "24","39", "52") 
head(cd4)
cd4.m <- melt(cd4)
cd4.m <- mutate(cd4.m, variable=as.integer(as.character(variable)))

## plotting the data
theme_set(theme_light())
p <- ggplot(cd4.m, aes(x=variable, y=log10(value))) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="Week", y="Log_10(CD4 T cells cells/mm^3 blood)") +
        theme(text = element_text(face="bold", size = 18)) +
        scale_x_continuous(breaks=c(0,8,12,24,39,52)) +     ## change ticks interval labels
        geom_smooth(method = "loess") 
plot(p)

##################################################################################
## 4.1 CD4/CD8 ratio dynamic
## selecting CD4/CD8 ratio columns
cociente_names <- grep("^CocCD4_CD8_S", names(vih_data), value = TRUE)
cociente_names

## subsetting data and setting data
cociente <- select(vih_data, cociente_names)
names(cociente) <- c("0", "8", "12", "24", "39", "52") 
head(cociente)
cociente.m <- melt(cociente)
cociente.m <- mutate(cociente.m, 
                     variable=as.integer(as.character(variable)))

## plotting the log data
theme_set(theme_light())
p <- ggplot(cociente.m, aes(x=variable, y=log10(value))) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="Week", y="Log_10(CD4/CD8 ratio)") +
        theme(text = element_text(face="bold", size = 18)) +
        scale_x_continuous(breaks=c(0,8,12,24,39,52)) +     ## change ticks interval labels
        geom_smooth(method = "loess")
plot(p)

## plotting the data
theme_set(theme_light())
p <- ggplot(cociente.m, aes(x=variable, y=value)) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="Week", y="CD4/CD8 ratio") +
        theme(text = element_text(face="bold", size = 18)) +
        scale_x_continuous(breaks=c(0,8,12,24,39,52)) +     ## change ticks interval labels
        geom_smooth(method = "loess")
plot(p)


##################################################################################
## 6. % CD4 T dynamics
## selecting Percentage of CD4 T columns
cd4_perc_names <- grep("CD4porcentajeS", names(vih_data), value = TRUE)
cd4_perc_names

## subsetting data and setting data
cd4_perc <- select(vih_data, cd4_perc_names)
names(cd4_perc) <- c("0", "8", "12", "24","39", "52") 
head(cd4_perc)
cd4_perc.m <- melt(cd4_perc)
cd4_perc.m <- mutate(cd4_perc.m, variable=as.integer(as.character(variable)))

## plotting the data
theme_set(theme_light())
p <- ggplot(cd4_perc.m, aes(x=variable, y=value)) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="Week", y="% CD4 T cells") +
        theme(text = element_text(face="bold", size = 18)) +
        scale_x_continuous(breaks=c(0,8,12,24,39,52)) +         ## change ticks interval labels
        geom_smooth(method = "loess")
plot(p)


###################################################################################
## figure 5.1

##  Linear regression of CD4porcentajeS52 vs Delta_CD4_year1
theme_set(theme_light())
g <- ggplot(vih_data, aes(x=CD4porcentajeS52, y=Delta_CD4_year1)) +
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21) + 
        labs(x="CD4 percentage week 52", y=" Delta CD4 (Year 1)") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "lm")
plot(g)


#Y= DELTA_CD4_S52, X=Creatinina S24
theme_set(theme_light())
g <- ggplot(vih_data, aes(y=Delta_CD4_year1, x=CreatininaS24)) +
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21) + 
        labs(x="Creatinine (Week 24)", y=" Delta CD4 week 52") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "lm")
plot(g)

## CreatininaS0
theme_set(theme_light())
g <- ggplot(vih_data, 
            aes(y=Delta_CD4_year1, x=CreatininaS0)) +
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21) + 
        labs(x="Creatinine (Week 0)", y=" % CD4 (Year 1)") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "lm")
plot(g)


##############################################################################
## Creatinina

creatinina <- select(vih_data, CreatininaS0, CreatininaS24)
names(creatinina) <- c("0", "24")
creatinina.m <- melt(creatinina)
creatinina.m

g <- ggplot(creatinina.m, aes(x=variable, y=value)) + 
        geom_boxplot() + 
        geom_jitter(position = position_jitter(0.3)) +
        labs(x="Week", y="Creatinine levels") +
        theme(text = element_text(face="bold", size = 18))
plot(g)

##########################################################################
## leave-one-out validation LASSO

## binarizing IRIS
wwiris <- sapply(vih_data$with_without_IRIS, 
                 function(x) ifelse(x=="with","IRIS","No symptoms")) 

## processing data for lasso
input <- vih_data[, ! sapply(vih_data, function(x) class(x)=="character") ]
input <- select(input, -Delta_CD4_year1)
input <- select(input, -CD4_S0)
input <- select(input, -CD4_S52)
input <- as.matrix(input)
str(input)
head(input)
write.table(input, "../data/model_matrix.tsv", sep = "\t")
str(input)
output <- vih_data$Delta_CD4_year1

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
                         "values"=output,
                         "iris"=as.factor(wwiris))
theme_set(theme_light())
p <- ggplot(validation, aes(x=values, y=lasso_prediction, color= iris)) + 
        geom_point(size= 2.5) +
        scale_colour_manual(values=c("green", "black")) +
        geom_abline(slope = 1, size=1, colour="red") +
        labs(x="Delta TCD4 values", y="LASSO") +
        theme(legend.position = c(0.15, 0.85)) +            ## change legend position
        theme(legend.title = element_blank())  +           ## remove legend title
        theme(text = element_text(face="bold", size = 18))
plot(p)



