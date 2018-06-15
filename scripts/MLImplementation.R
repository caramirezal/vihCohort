## This script contains the R functions to implement Machine Learning, parameter selection,
## and plotting results in technicalReportML.Rmd 

## Dependencies
library(randomForest)
library(caret)
library(party)
library(gridExtra)
library(reshape2)
library(glmnet)        

## function for calculating r^2  as 1 - rss/tss
r_sq <- function(observed,predicted) {
        rss <- sum((observed-predicted)^2)
        tss <- sum((observed-mean(observed))^2)
        return( 1 - (rss/tss) )
}

#########################################################################################
## data preprocessing

## reading data
data <- read.csv("../data/TablaLPS_IL18_20180217.csv")

## drop categorical data
data <- data[,sapply(data, is.numeric)]

## list of removed variables
drop <- c("Expediente",
          "DeltaCD4W52atleast150",
          "DeltaCD4_W52",
          "CD4TCabove200_W052",
          "CD4TCabove350_W052",
          "Numero_consecutivo")

## removing list of variables
input <- data[,!(colnames(data)%in%drop)]

## input definition
input <- as.matrix(input)

## input to dataframe
input.df <- as.data.frame(input)

## output definition
output <- data$DeltaCD4_W52

###########################################################################################
## mtry variation

## perform cforest with mtry variation
forestmtry <- function(output,
                       input,
                       nsims=30,
                       method="randomForest",
                       mtryInterval=seq(2,10,1)) {
        
        if ( class(input) == "matrix" ) {
                input <- as.data.frame(input)
        } 
        if ( class(input) != "data.frame") {
                stop("Input must be a data.frame or matrix")
        }
        
        if ( ! method %in% c("cforest","randomForest") ) {
                stop("Method must be 'cforest' or 'randomForest'")
        }
        
        results <- matrix(0,nsims*length(mtryInterval),3)
        results <- as.data.frame(results)
        colnames(results) <- c("Accuracy","mtry","performance")
        
        for (i in 1:length(mtryInterval)){
                for (j in 1:nsims) {
                        
                        if ( method == "cforest" ) {
                                ## simulation parameters
                                intrain <- createDataPartition(output,p = 0.8,list = FALSE)
                                cforestControl <- cforest_control(teststat = "quad",
                                                                  testtype = "Univ", 
                                                                  mincriterion = 0, 
                                                                  ntree = 50, 
                                                                  mtry = mtryInterval[i], ## variation
                                                                  replace = FALSE)
                                
                                
                                ## perform simulation
                                initialTime <- Sys.time()
                                forest <- cforest(output[intrain]~.,
                                                  data=as.data.frame(input[intrain,]),
                                                  control = cforestControl)
                                finalTime <- Sys.time() - initialTime
                                
                                ## perform prediction on test sample
                                forest.pred <- predict(forest,
                                                       newdata=input[-intrain,])
                                mss <- mean((forest.pred-output[-intrain])^2)
                                mse <- sqrt(mss)
                        }
                        
                        if ( method == "randomForest" ) {
                                intrain <- createDataPartition(output,p = 0.8,list = FALSE)
                                
                                ## perform randomForest
                                initialTime <- Sys.time()
                                rforest <- randomForest(output[intrain]~.,
                                                        data=input[intrain,])
                                finalTime <- Sys.time() - initialTime
                                
                                ## perform predictions
                                rforest.pred <- predict(rforest,
                                                        newdata = input[-intrain,])
                                
                                ## calcualte mse values
                                mss <- mean((rforest.pred-output[-intrain])^2)
                                mse <- sqrt(mss)
                                mse
                                
                        }
                        
                        ## index to store simulation data
                        index <- (nsims*(i-1)) + j
                        
                        ## storing results
                        results[index,1] <- mse
                        results[index,2] <- mtryInterval[i]
                        results[index,3] <- finalTime
                } 
        }
        
        ## plot accuracy test
        theme_set(theme_bw())
        g <- ggplot(results, aes(x=mtry,y=Accuracy)) + geom_point() 
        g <- g + geom_boxplot(aes(group=mtry,fill="steelblue"))
        g <- g + theme(legend.position = "none") 
        g <- g + labs(title =paste(method,"- Accuracy"),y="MSE")
        g <- g + theme(title = element_text(face="bold"),
                       axis.title.x = element_text(face="bold"),
                       axis.title.y = element_text(face="bold"))
        
        ## plot performance time
        theme_set(theme_bw())
        g1 <- ggplot(results, aes(x=mtry,y=performance)) + geom_point() 
        g1 <- g1 + geom_boxplot(aes(group=mtry,       ## x axis cathegorical values  
                                    fill=132))        ## 132 -> steelblue color
        g1 <- g1 + theme(legend.position = "none")    ## removing legend
        g1 <- g1 + labs(title =paste(method,"- Waiting time"),
                        y="Time (secs)")
        g1 <- g1 + theme(title = element_text(face="bold"),  
                         axis.title.x = element_text(face="bold"),
                         axis.title.y = element_text(face="bold"))
        
        grid.arrange(g,g1,nrow=1)
        
} 

## performs mtry parameter variation in randomForest
jpeg("../figures/randforest_mtry_nsims=100_method=RandomForest.jpg")
forestmtry(output = output,
           input = input,
           nsims = 100,
           method = "randomForest",
           mtryInterval = seq(2,10,1))      ## parameters values 
dev.off()

## performs ntree parameter variation in cforest
jpeg("../figures/randforest_mtry_nsims=100_method=cforest.jpg")
forestmtry(output = output,
           input = input,
           nsims = 100,
           method = "cforest",
           mtryInterval = seq(2,10,1))      ## parameters values 
dev.off()

############################################################################################
## ntree paramter variation

## perform cforest with ntree variation
forestntree <- function(output,
                        input,
                        method="randomForest",
                        nsims=30,
                        mtry=3,
                        ntreeInterval=seq(20,100,10)) {
        
        if ( class(input) == "matrix" ) {
                input <- as.data.frame(input)
        } 
        if ( class(input) != "data.frame") {
                warning("Input must be a data.frame or matrix")
        }
        if ( ! method %in% c("cforest","randomForest") ) {
                stop("Method must be 'cforest' or 'randomForest'")
        }
        
        results <- matrix(0,nsims*length(ntreeInterval),3)
        results <- as.data.frame(results)
        colnames(results) <- c("Accuracy","ntree","performance")
        
        for (i in 1:length(ntreeInterval)){
                for (j in 1:nsims) {
                        
                        ## Use cforest algorithm
                        if ( method == "cforest" ) {
                                ## initial settings
                                intrain <- createDataPartition(output,p = 0.8,list = FALSE)
                                cforestControl <- cforest_control(teststat = "quad",
                                                                  testtype = "Univ", 
                                                                  mincriterion = 0, 
                                                                  ntree = ntreeInterval[i], 
                                                                  mtry = 3,
                                                                  replace = FALSE)
                                
                                ## perform cforest
                                initialTime <- Sys.time()
                                forest <- cforest(output[intrain]~.,data=input[intrain,],
                                                  control = cforestControl)
                                finalTime <- Sys.time() - initialTime
                                
                                ## calculate predictions
                                forest.pred <- predict(forest,newdata=input[-intrain,])
                                
                                ## calculate mse
                                mss <- mean((forest.pred-output[-intrain])^2)
                                mse <- sqrt(mss)
                        }
                        
                        ## Use randomForest method
                        if ( method == "randomForest" ) {
                                intrain <- createDataPartition(output,p = 0.8,list = FALSE)
                                initialTime <- Sys.time()
                                forest <- randomForest(output[intrain]~.,
                                                       data=input[intrain,],
                                                       ntree = ntreeInterval[i])
                                finalTime <- Sys.time() - initialTime
                        }
                        
                        ## perform prediction on test sample
                        forest.pred <- predict(forest,
                                               newdata=input[-intrain,])
                        mss <- mean((forest.pred-output[-intrain])^2)
                        mse <- sqrt(mss)
                        
                        ## index to store simulation data
                        index <- (nsims*(i-1)) + j
                        
                        ## storing results
                        results[index,1] <- mse
                        results[index,2] <- ntreeInterval[i]
                        results[index,3] <- finalTime
                } 
        }
        

        ## plot Accuracy MSE values
        theme_set(theme_bw())
        g <- ggplot(results, aes(x=ntree,y=Accuracy)) + geom_point() 
        g <- g + geom_boxplot(aes(group=ntree,fill="steelblue"))
        g <- g + theme(legend.position = "none") 
        g <- g + labs(title =paste(method,"- Accuracy"),y="MSE")
        g <- g + theme(title = element_text(face="bold"),
                       axis.title.x = element_text(face="bold"),
                       axis.title.y = element_text(face="bold"))
        
        ## plot waiting times
        theme_set(theme_bw())
        g1 <- ggplot(results, aes(x=ntree,y=performance)) + geom_point() 
        g1 <- g1 + geom_boxplot(aes(group=ntree,       ## x axis cathegorical values  
                                    fill=132))        ## 132 -> steelblue color
        g1 <- g1 + theme(legend.position = "none")    ## removing legend
        g1 <- g1 + labs(title =paste(method,"- Waiting time"),
                        y="Time (secs)")
        g1 <- g1 + theme(title = element_text(face="bold"),  
                         axis.title.x = element_text(face="bold"),
                         axis.title.y = element_text(face="bold"))
        
        grid.arrange(g,g1,nrow=1)
        
}

## performs ntree variation in randomForest method
jpeg("../figures/randforest_ntree_nsims=100_method=RandomForest.jpg")
forestntree(output = output,input = input,
           nsims = 100,method = "randomForest",
           ntreeInterval = seq(20,100,10))   ## parameter values
dev.off()

## performs ntree variation in cforest method
jpeg("../figures/randforest_ntree_nsims=100_method=cforest.jpg")
forestntree(output = output,input = input,
           nsims = 100,method = "cforest",
           ntreeInterval = seq(20,100,10))   ## parameter values
dev.off()


######################################################################################
## Compare cforest versus randomForest functions for different train
## sizes

nsims <- 100
p = seq(0.5,0.9,0.05)
input.df <- as.data.frame(input)
results <- matrix(0,nsims*length(p),7)
results <- as.data.frame(results)
colnames(results) <- c("Size",
                       "cforest_accuracy",
                       "cforest_rsq",
                       "cforest_time",
                       "rforest_accuracy",
                       "rforest_rsq",
                       "rforest_time")

for ( i in 1:length(p) ) {
        for ( j in 1:nsims ) {
                
                ###########################################################
                ## cforest method
                
                ## simulation parameters
                intrain <- createDataPartition(output,p = p[i],list = FALSE)
                cforestControl <- cforest_control(teststat = "quad",
                                                  testtype = "Univ", 
                                                  mincriterion = 0, 
                                                  ntree = 45, 
                                                  mtry = 3,
                                                  replace = FALSE)
                
                
                ## perform simulation
                initialTime <- Sys.time()
                forest <- cforest(output[intrain]~.,
                                  data=input.df[intrain,],
                                  control = cforestControl)
                finalTime <- Sys.time() - initialTime
                
                ## perform prediction on test sample
                forest.pred <- predict(forest,
                                       newdata=input.df[-intrain,])
                mss <- mean((forest.pred-output[-intrain])^2)
                msecforest <- sqrt(mss)
                
                ## calculation of R^2
                rsqcforest <- r_sq(output[-intrain],forest.pred)
                
                index <- (i - 1)*nsims + j
                results$Size[index] <- p[i]
                results$cforest_accuracy[index] <- msecforest
                results$cforest_rsq[index] <- rsqcforest
                results$cforest_time[index] <- finalTime
                
                ########################################################
                ## randomForest method
                initialTime <- Sys.time()
                #intrain <- createDataPartition(output,p = p[i],list = FALSE)
                forest <- randomForest(output[intrain]~.,
                                       data=input[intrain,],
                                       ntree = 50,
                                       mtry=4)
                finalTime <- Sys.time() - initialTime
                
                ## perform prediction on test sample
                forest.pred <- predict(forest,
                                       newdata=input.df[-intrain,])
                
                ## calculate accuracy
                mss <- mean((forest.pred-output[-intrain])^2)
                mse <- sqrt(mss)
                
                ## calculate rsq
                rsq <- r_sq(output[-intrain],forest.pred)
                
                ## store results
                results$rforest_accuracy[index] <- mse
                results$rforest_rsq[index] <- rsq
                results$rforest_time[index] <- finalTime
                
        }
}

## stroing results
write.csv(results,"../data/forestTrainSizes.csv",row.names = FALSE)


## preprocessing previuos results
results <- read.csv("../data/forestTrainSizes.csv")
results <- melt(results,id="Size")

## calculating mean of mse and R^2 values
means <- with(results,aggregate(value,list(Size,variable),mean) )
## standard deviations of mse and R^2 values
sd <- with(results, aggregate(value,list(Size,variable),sd) )
results.summary <- merge(means,sd,by=c("Group.1","Group.2"))
colnames(results.summary) <- c("size","technique","mean","sd")


accuracy <- results.summary[grepl("accuracy",results.summary$technique),]

theme_set(theme_light())

## plotting accuracy
g <- ggplot(accuracy, aes(x=size,y=mean,colour=technique)) 
g <- g +  geom_line()
g <- g + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01) 
g <- g + geom_point(size=3.5)
g <- g + labs(x="Size", y="MSE")
g <- g + theme(text = element_text(size=16,face="bold"),
               axis.line = element_line(colour = 'black', size = 0.7))
g <- g + labs(title="Accuracy")
g <- g + theme(legend.position = "none")


time <- results.summary[grepl("time",results.summary$technique),]

## plotting performance time
g1 <- ggplot(time, aes(x=size,y=mean,colour=technique)) 
g1 <- g1 +  geom_line()
g1 <- g1 + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01) 
g1 <- g1 + geom_point(size=3.5)
g1 <- g1 + labs(x="Size", y="MSE")
g1 <- g1 + theme(text = element_text(size=16,face="bold"),
                 axis.line = element_line(colour = 'black', size = 0.7))
g1 <- g1 + labs(title="Performance time")
g1 <- g1 + theme(legend.position = "none")

## saving plot
jpeg("../figures/forestTrainSizes.jpg")
grid.arrange(g,g1,nrow=1)
dev.off()


###################################################################################
## Lasso simulations

## number of simulations
nsims <- 100

## train testing sizes
sizes <- seq(0.5,0.90,0.05)

## matrix to store results of simulations
results <- matrix(0,length(sizes)*nsims,4)
results <- as.data.frame(results)
colnames(results) <- c("mse_lasso","r_sq_lasso",
                       "size","time_lasso")


## iteration of lasso, ridge, and random forest regression simulations
for ( i in 1:length(sizes) ) {
        ## train test definition
        intrain <- createDataPartition(output,p=sizes[i],list = FALSE)
        
        for ( j in 1:nsims ) {
                ##########################################################
                ## lasso 
                
                ## lasso Lambda optimal calculation
                initialTime <- Sys.time()
                lasso <- cv.glmnet(x=input[intrain,],y = output[intrain])
                totalTime <- Sys.time() - initialTime 
                
                y_predicted <- predict(lasso, newx=input[-intrain,],
                                               c=lasso$lambda.min,
                                               type="response")
                
                mss_lasso <- mean((y_predicted-output[-intrain])^2)
                mse_lasso <- sqrt(mss_lasso)
                
                r_sq_lasso <- r_sq(observed=output[-intrain],
                                   predicted=y_predicted)
                
                index <- (i-1)*nsims + j
                results[index, "mse_lasso"] <- mse_lasso
                results[index,"r_sq_lasso"] <- r_sq_lasso
                results[index,"size"] <- sizes[i]
                results[index,"time_lasso"] <- totalTime
        }
}

head(results)
write.csv(results,"../data/lassoTrainSizes.csv",
          row.names = FALSE)

results.m <- melt(results,id="size")
head(results.m)

## calculating mean of mse and R^2 values
means <- with(results.m,aggregate(value,list(size,variable),mean) )
## standard deviations of mse and R^2 values
sd <- with(results.m, aggregate(value,list(size,variable),sd) )
results.summary <- merge(means,sd,by=c("Group.1","Group.2"))
colnames(results.summary) <- c("size","technique","mean","sd")
results.summary

## select accuracy values for plot
mse.summary <- results.summary[grepl("mse",results.summary$technique),] 

## plot mse values
theme_set(theme_light())
g <- ggplot(mse.summary, aes(x=size,y=mean,colour=technique)) 
g <- g +  geom_line()
g <- g + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01) 
g <- g + geom_point(size=3.5)
g <- g + labs(x="Size", y="MSE")
g <- g + theme(text = element_text(size=26,face="bold"),
               axis.line = element_line(colour = 'black', size = 0.7))
g <- g + labs(title="Lasso accuracy")
g <- g + theme(legend.position = "none")
jpeg("../figures/lassoMSEtrainSize.jpg")
plot(g)
dev.off()

## select accuracy values for plot
rsq.summary <- results.summary[grepl("r_sq",results.summary$technique),] 

## plot rsq values
theme_set(theme_light())
g <- ggplot(rsq.summary, aes(x=size,y=mean,colour=technique)) 
g <- g +  geom_line()
g <- g + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01) 
g <- g + geom_point(size=3.5)
g <- g + labs(x="Size", y="R^2")
g <- g + theme(text = element_text(size=26,face="bold"),
               axis.line = element_line(colour = 'black', size = 0.7))
g <- g + labs(title="Lasso Variance explained")
g <- g + theme(legend.position = "none")
jpeg("../figures/lassoRsqtrainSize.jpg")
plot(g)
dev.off()

time.summary <- results.summary[grepl("time",results.summary$technique),]
## plot rsq values
theme_set(theme_light())
g <- ggplot(time.summary, aes(x=size,y=mean,colour=technique)) 
g <- g +  geom_line()
g <- g + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd), width=.01) 
g <- g + geom_point(size=3.5)
g <- g + labs(x="Size", y="Time (secs)")
g <- g + theme(text = element_text(size=26,face="bold"),
               axis.line = element_line(colour = 'black', size = 0.7))
g <- g + labs(title="Lasso Performance time")
g <- g + theme(legend.position = "none")
jpeg("../figures/lassoTimeTrainSize.jpg")
plot(g)
dev.off()