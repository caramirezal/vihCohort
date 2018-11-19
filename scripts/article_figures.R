## Construction of the figures for the article
library(glmnet)
library(dplyr)
library(ggplot2)
library(reshape2)

vih_data <- read.csv("../data/cleandata.csv", stringsAsFactors = FALSE)
str(vih_data)
names(vih_data)

##################################################################################
## figure 1 dynamic of viral load, T CD4, % CD4 and cocCD4_CD8
## figure 1.a. Viral load dynamic
## selecting viral load columns
cv_names <- grep("^CV_S0", names(vih_data), value = TRUE)
cv_names

## subsetting data and setting data
cv <- select(vih_data, cv_names)
names(cv) <- c("0", "8", "12", "39", "52") 
head(cv)
cv.m <- melt(cv)
cv.m <- mutate(cv.m, variable=as.integer(as.character(variable)))

## plotting the data
theme_set(theme_light())
p <- ggplot(cv.m, aes(x=variable, y=log(value))) + 
        geom_point(fill="black", colour="white", 
                   size= 2.5, pch=21, alpha=0.3) + 
        labs(x="Week", y="log(viral_load)", title = "Viral load dynamic") +
        theme(text = element_text(face="bold", size = 18)) +
        geom_smooth(method = "loess")
plot(p)

##################################################################################
## figure 1.a. CD4 T dynamic
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
## figure 1.c. Percentage of CD4 T dynamics
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
## figure 1 dynamic of viral load, T CD4, % CD4 and cocCD4_CD8
## figure 1.d. CD4/CD8 ratio dynamic
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

