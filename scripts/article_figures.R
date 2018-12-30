## Construction of the figures for the article
library(glmnet)
library(dplyr)
library(ggplot2)
library(reshape2)

vih_data <- read.csv("../data/cleandata.csv", stringsAsFactors = FALSE)
str(vih_data)
names(vih_data)

##################################################################################
## figure 2 dynamic of viral load, T CD4, % CD4 and cocCD4_CD8

## figure 2.a. Viral load dynamic
## selecting viral load columns
library(dplyr)
library(reshape2)
library(ggplot2)

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
        labs(x="Week", y=" Log (HIV RNA copies/mL)") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth()
plot(p)

##################################################################################
## figure 2.a. CD4 T dynamic
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
p <- ggplot(cd4.m, aes(x=variable, y=value)) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="Week", y="Number of CD4 T cells", title = "CD4 T cells dynamic") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "loess")
plot(p)

##################################################################################
## figure 2.c. Percentage of CD4 T dynamics
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
        labs(x="Week", y="Percentage of CD4 T cells", title = "% CD4 T cells dynamic") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "loess")
plot(p)

##################################################################################
## figure 2 dynamic of viral load, T CD4, % CD4 and cocCD4_CD8
## figure 2.d. CD4/CD8 ratio dynamic
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

## plotting the data
theme_set(theme_light())
p <- ggplot(cociente.m, aes(x=variable, y=log(value))) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="Week", y="CD4/CD8 ratio", title = "CD4/CD8 ratio dynamic") +
        theme(text = element_text(face="bold", size = 18)) +
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
g <- ggplot(vih_data, aes(y=CD4porcentajeS52, x=CreatininaS24)) +
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21) + 
        labs(x="Creatinine (Week 24)", y=" % CD4 (Year 1)") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "lm")
plot(g)

## CreatininaS0
theme_set(theme_light())
g <- ggplot(vih_data, 
            aes(y=CD4porcentajeS52, x=CreatininaS0)) +
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

##################################################################################
## Test binarize IRIS
library(glmnet)
library(dplyr)
library(ggplot2)

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
                            "sd"=sd_coef) %>%
        arrange(desc(abs(mean)))
write.csv(summary_coefs, "../data/lasso_only_numeric_plus_iris.csv", row.names=FALSE)
