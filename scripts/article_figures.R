## Construction of the figures for the article
library(glmnet)
library(dplyr)
library(ggplot2)
library(reshape2)
library(corrplot)
library(googlesheets)
library(forcats)
library(gplots)
library(plyr)
library(pheatmap)
library(heatmap3)
library(partykit)


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
cv.m <- filter(cv.m, variable!=4)
cv.m <- mutate(cv.m, variable=as.vector(variable))
cv.m <- filter(cv.m, variable != 104)

theme_set(theme_light())
p <- ggplot(cv.m, aes(x=variable, y=log10(value))) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="Week", y=" Log_10(HIV RNA copies/mL)") +
        theme(text = element_text(face="bold", size = 18)) +
        scale_x_continuous(breaks=c(0,8,12,24,39,52)) +     ## change ticks interval labels
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
p <- ggplot(validation, aes(x=lasso_prediction, y=values)) + 
        geom_point(colour="steelblue", size= 2.5) + 
        geom_abline(slope = 1,colour="red",size=1) +
        labs(x="LASSO", y="Delta TCD4 values") +
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
## Bar graph of the top 20 variables ranked by LASSO
sheet <- gs_title("Para heat map 20190416")

labels.df <- gs_read(sheet, ws="Labels", col_names=FALSE)
labels <- labels.df$X4
labels

## load LASSO coeficients
coefs <- read.csv("../data/lasso_only_numeric.csv", 
                  stringsAsFactors = FALSE)
## drop intercept coeficient
coefs <- filter(coefs, coefficient != "(Intercept)")
coefs <- mutate(coefs, coefficient=plyr::mapvalues(coefs$coefficient,
                                                   from=labels.df$X1,
                                                   to=labels.df$X4))

theme_set(theme_light())
coefs[1:20,] %>% 
        mutate(coefficient= fct_reorder(coefficient, log(abs(mean)+1))) %>%
        ggplot(aes(x=coefficient, y=log(abs(mean)+1))) +
           geom_bar(stat="identity") +
           coord_flip() +
           xlab("") +
           theme(legend.position = "None",
                 axis.title.x = element_blank(),
                 axis.text.y = element_text(face="bold", size=10))

coefs.ord <- arrange(coefs[1:35,], desc(abs(mean))) 
ggplot(data=coefs.ord, aes(x=coefficient, y=log(abs(mean)))) + geom_bar(stat="identity")

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
g <- ggplot(res, aes(x=prediction, y=value, color=as.factor(iris))) + 
        scale_color_manual(values = c("green", "black")) +    ## adjust point color manually
        geom_point(size= 2.5) + 
        labs(x="Predicted value", y="Delta TCD4 values") +
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
p <- ggplot(cd4.m, aes(x=variable, y=value)) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="Week", y="CD4 T cells cells/mm^3 blood") +
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
        labs(x="Creatinine Week 24", y=" Delta CD4 Week 52") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "lm")
plot(g)

## CreatininaS0
theme_set(theme_light())
g <- ggplot(vih_data, 
            aes(y=Delta_CD4_year1, x=CreatininaS0)) +
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21) + 
        labs(x="Creatinine Week 0", y="Delta CD4 Week 52") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "lm")
plot(g)


##############################################################################
## Creatinina

creatinina <- select(vih_data, CreatininaS0, CreatininaS24)
names(creatinina) <- c("0", "24")
creatinina.m <- melt(creatinina)
creatinina.m$"PairNumber" <- rep(seq(1,nrow(creatinina.m)/2),2)

g <- ggplot(creatinina.m, aes(x=variable, y=value)) + 
        geom_point() +
        geom_line(aes(group = PairNumber)) +
        labs(x="Week", y="Creatinine levels") +
        theme(text = element_text(face="bold", size = 22))
plot(g)

##############################################################################
## CD4 % week 0 and 24 paired sample plot

cd4 <- select(vih_data, CD4porcentajeS0, CD4porcentajeS24)
names(cd4) <- c("0", "24")
cd4.m <- melt(cd4)
cd4.m$"PairNumber" <- rep(seq(1,nrow(cd4.m)/2),2)

g <- ggplot(cd4.m, aes(x=variable, y=value)) + 
        geom_point() +
        geom_line(aes(group = PairNumber)) +
        labs(x="Week", y="% CD4 T cells") +
        theme(text = element_text(face="bold", size = 22))
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
        labs(y="Delta TCD4 values", x="LASSO") +
        theme(legend.position = c(0.15, 0.85)) +            ## change legend position
        theme(legend.title = element_blank())  +           ## remove legend title
        theme(text = element_text(face="bold", size = 18))
plot(p)

################################################################################################
Sys.setlocale('LC_ALL','C') ## fix non english characters error

ssheet <- gs_title("Para heat map")
basal <- gs_read(ssheet, ws = "Data")
labels <- gs_read(ssheet, ws="Labels", col_names=FALSE)
names(labels) <- c("var_names","labels")


## computes p-values for pair-wise variable correlations
cor.mtest <- function(mat, ...) {
        mat <- as.matrix(mat)
        n <- ncol(mat)
        p.mat<- matrix(NA, n, n)
        diag(p.mat) <- 0
        for (i in 1:(n - 1)) {
                for (j in (i + 1):n) {
                        tmp <- cor.test(mat[, i], mat[, j], ...)
                        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
                }
        }
        colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
        p.mat
}

basal.m <- as.matrix(basal)
names(basal.m) <- labels$labels

## plot correlation between variables
cor_vars <- cor(basal.m, use = "complete.obs", method = "spearman")
heatmap(cor_vars, cexRow = 0.5, labCol = FALSE)
corrplot(cor_vars, 
         method = "circle",
         type="upper",
         diag = FALSE,
         order="hclust",
         tl.col="black",    ## label color
         tl.srt=90,
         tl.cex = 0.45,
         pch.cex = 10,
         tl.pos = "td")

## saving correlations table
write.table(cor_vars, 
            "../data/correlations.tsv", 
            sep = "\t")

## calculate p-values
pvals <- cor.mtest(basal, method="spearman", alternative="two.sided")


## saving pvals
rownames(pvals) <- colnames(pvals)
write.table(pvals, 
            "../data/pvals.tsv",
            sep = "\t")


corrplot(cor_vars, 
         p.mat = pvals,
         sig.level = c(.001, .01, .05), 
         pch.cex = 1.25,
         cl.cex = 1.25,
         tl.cex = 1.25,
         insig = "label_sig", 
         pch.col = "white",
         method = "circle",
         type="upper",
         diag = FALSE,
         order="hclust",
         tl.col="black",    ## label color
         tl.srt=45)         ## column label angle

########################################################################################################################
## Loading data
sheet <- gs_title("Para heat map 20190416")
basales <- gs_read(sheet, ws = "Data")
str(basales)

## performing correlation calculation 
basales_cor <- cor(basales, use = "complete.obs")
#basales_pvals <- cor.mtest(basales)

## Load labels and colors
categories <- gs_read(sheet, ws="colors")
str(categories)
unique(categories$Color)
## check de color-name
sum(categories$Abreviated_name == colnames(basales_cor)) ## column variables and colors are not in the same order
categories_sort <- arrange(categories, Abreviated_name) ## sort categories
basales_cor_sort <- basales_cor[sort(colnames(basales_cor)), sort(colnames(basales_cor))]
categories_sort <- mutate(categories_sort, cat_int= paste(Category, Intensity, sep = " "))
color_indexes <- order(unique(categories_sort$cat_int))
color_ord_legend <- unique(categories_sort$Color)[color_indexes]

pheatmap(basales_cor_sort, 
         cutree_cols = 3,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         color = colorpanel(30,"yellow","red"))

library(dendextend)
library(ComplexHeatmap)

## Heatmap implementation with complexHeatmap Library
row_dend = hclust(dist(basales_cor_sort)) # row clustering
col_dend = hclust(dist(t(basales_cor_sort))) # column clustering
Heatmap(basales_cor_sort,
        cluster_rows = color_branches(row_dend, k = 3),
        cluster_columns = color_branches(col_dend, k = 3),
        show_column_names = FALSE)

## plot heatmap
par(mar=c(7, 1, 0, 0))
tags <- unique(categories_sort$cat_int)[color_indexes]
tags[1] <- "Age"      ## change Age normal
jpeg("../figures/heatmap_v2.jpg")
heatmap.2(basales_cor_sort, trace = "none",
          RowSideColors = categories_sort$Color, ColSideColors = categories_sort$Color, 
          labRow=NA, density.info = "none", key = FALSE, margins=c(12,8), cexCol = 1)
#coords <- locator(1)     ## locates coordinates in the heatmap using cursor
#coords

legend(.01, 0.99, legend = tags, 
       col= color_ord_legend, 
       lty = 5, lwd = 5, cex = 0.58, box.lty = 0, text.font = 2)
dev.off()

## heatmap with column and row labels (difficult to see)
heatmap.2(basales_cor_sort, trace = "none",
          RowSideColors = categories_sort$Color, ColSideColors = categories_sort$Color, 
          density.info = "none")


############################################################################

coefs <- read.table("../data/lasso_only_numeric.csv",
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    sep = ",")
str(coefs)

##################################################################################
## Adaptative lasso implementation

## Performing adaptative lasso
alasso <- glmnet(
        x = input, 
        y = output,
        alpha = 1,
        penalty.factor = 1 / (1 + abs(lasso_coefs)) 
)

## extracting and sorting lasso coefs 
values <- coef(alasso) %>% rowSums 
alasso.df <-  data.frame(
        'variables' = names(values),
        'value' = values
)

## saving adapative lasso results
alasso.df %>% 
        arrange(desc(abs(value))) %>%
        write.table(file = '../data/alasso_coefs.tsv',
                    row.names = FALSE,
                    sep = '\t')

###################################################################################
## random forest implementation
input <- data.frame(input)
forest <- cforest(output~., data = input, )

var_imp <- varimp(forest)
cforest.df <- data.frame(
        'variables' = names(var_imp),
        'value' = var_imp, 
        stringsAsFactors = FALSE
)

## saving random forest results
cforest.df %>%
        arrange(desc(value)) %>%
        write.table(file = '../data/forest_var_imp.tsv',
                    row.names = FALSE,
                    sep = '\t')

