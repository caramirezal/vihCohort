## Random forest implementation on CD4 T cells reconstitution of a cohort of VIH patients

library(randomForest)
library(caret)
library(party)
library(gridExtra)

## function for calculating r^2  as 1 - rss/tss
r_sq <- function(observed,predicted) {
        rss <- sum((observed-predicted)^2)
        tss <- sum((observed-mean(observed))^2)
        return( 1 - (rss/tss) )
}

##########################################################################
## data preprocessing

## reading data
data <- read.csv("~/scripts/vihCohort/data/TablaLPS_IL18_20180217.csv")
head(data)


## preprocessing raw data
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

## output definition
output <- data$DeltaCD4_W52
colnames(input)


########################################################################
## simple random forest usinf cforest package
set.seed(3333)

train <- createDataPartition(output,p = 0.8,list = FALSE)

initialTime <- Sys.time()
cforestControl <- cforest_control(teststat = "quad",
                                      testtype = "Univ", mincriterion = 0, ntree = 50, mtry = 3,
                                      replace = FALSE)
forest <- cforest(output[train]~.,data=as.data.frame(input[train,]),
                  control = cforestControl)
forest.pred <- predict(forest,newdata=as.data.frame(input[-train,]))
mss <- mean(sum(forest.pred-output[-train])^2)
mse <- sqrt(mss)
mse
finalTime <- Sys.time() - initialTime
finalTime

plot(forest.pred,output[-train],
     pch=20,cex=1.5,col="steelblue")
abline(0,1, col="red",lwd=2)

## parameters: mtry, ntree, teststat, testtype 

##############################################################################
## mtry variation

## perform cforest with mtry variation
cforestmtry <- function(output,
                        input,
                        nsims=30,
                        mtryInterval=seq(2,10,1)) {
        
        if ( class(input) == "matrix" ) {
                input <- as.data.frame(input)
        } 
        if ( class(input) != "data.frame") {
                warning("Input must be a data.frame or matrix")
        }
        
        results <- matrix(0,nsims*length(mtryInterval),3)
        results <- as.data.frame(results)
        colnames(results) <- c("Accuracy","mtry","performance")
        
        for (i in 1:length(mtryInterval)){
                for (j in 1:nsims) {
                        
                        ## simulation parameters
                        train <- createDataPartition(output,p = 0.8,list = FALSE)
                        cforestControl <- cforest_control(teststat = "quad",
                                                          testtype = "Univ", 
                                                          mincriterion = 0, 
                                                          ntree = 50, 
                                                          mtry = mtryInterval[i], ## variation
                                                          replace = FALSE)
                        
                        
                        ## perform simulation
                        initialTime <- Sys.time()
                        forest <- cforest(output[train]~.,
                                          data=input[train,],
                                          control = cforestControl)
                        finalTime <- Sys.time() - initialTime
                        
                        ## perform prediction on test sample
                        forest.pred <- predict(forest,
                                               newdata=input[-train,])
                        mss <- mean(sum(forest.pred-output[-train])^2)
                        mse <- sqrt(mss)
                        
                        ## index to store simulation data
                        index <- (nsims*(i-1)) + j
                        
                        ## storing results
                        results[index,1] <- mse
                        results[index,2] <- mtryInterval[i]
                        results[index,3] <- finalTime
                } 
        }
        
        # return(results)
        
        theme_set(theme_bw())
        g <- ggplot(results, aes(x=mtry,y=Accuracy)) + geom_point() 
        g <- g + geom_boxplot(aes(group=mtry,fill="steelblue"))
        g <- g + theme(legend.position = "none") 
        g <- g + labs(title ="\t\t\t Accuracy",y="MSE")
        g <- g + theme(title = element_text(face="bold"),
                  axis.title.x = element_text(face="bold"),
                  axis.title.y = element_text(face="bold"))
        
        theme_set(theme_bw())
        g1 <- ggplot(results, aes(x=mtry,y=performance)) + geom_point() 
        g1 <- g1 + geom_boxplot(aes(group=mtry,       ## x axis cathegorical values  
                                    fill=132))        ## 132 -> steelblue color
        g1 <- g1 + theme(legend.position = "none")    ## removing legend
        g1 <- g1 + labs(title ="\t\t Waiting time",y="Time (secs)")
        g1 <- g1 + theme(title = element_text(face="bold"),  
                  axis.title.x = element_text(face="bold"),
                  axis.title.y = element_text(face="bold"))
        
        grid.arrange(g,g1,nrow=1)
        
} 

## cforest accuracy does not depends on mtry
cforestmtry(output,as.data.frame(input),nsims = 30)


##############################################################################
## ntree variation

## perform cforest with ntree variation
cforestntree <- function(output,
                        input,
                        nsims=30,
                        mtry=3,
                        ntreeInterval=seq(20,100,10)) {
        
        if ( class(input) == "matrix" ) {
                input <- as.data.frame(input)
        } 
        if ( class(input) != "data.frame") {
                warning("Input must be a data.frame or matrix")
        }
        
        results <- matrix(0,nsims*length(ntreeInterval),3)
        results <- as.data.frame(results)
        colnames(results) <- c("Accuracy","ntree","performance")
        
        for (i in 1:length(ntreeInterval)){
                for (j in 1:nsims) {
                        
                        ## simulation parameters
                        train <- createDataPartition(output,p = 0.8,list = FALSE)
                        initialTime <- Sys.time()
                        cforestControl <- cforest_control(teststat = "quad",
                                                          testtype = "Univ", 
                                                          mincriterion = 0, 
                                                          ntree = ntreeInterval[i],  ## variation
                                                          replace = FALSE)
                        finalTime <- Sys.time() - initialTime
                        
                        ## perform simulation
                        initialTime <- Sys.time()
                        forest <- cforest(output[train]~.,
                                          data=input[train,],
                                          control = cforestControl)
                        finalTime <- Sys.time() - initialTime
                        
                        ## perform prediction on test sample
                        forest.pred <- predict(forest,
                                               newdata=input[-train,])
                        mss <- mean(sum(forest.pred-output[-train])^2)
                        mse <- sqrt(mss)
                        
                        ## index to store simulation data
                        index <- (nsims*(i-1)) + j
                        
                        ## storing results
                        results[index,1] <- mse
                        results[index,2] <- ntreeInterval[i]
                        results[index,3] <- finalTime
                } 
        }
        
        # return(results)
        
        theme_set(theme_bw())
        g <- ggplot(results, aes(x=ntree,y=Accuracy)) + geom_point() 
        g <- g + geom_boxplot(aes(group=ntree,fill="steelblue"))
        g <- g + theme(legend.position = "none") 
        g <- g + labs(title ="\t\t\t Accuracy",y="MSE")
        g <- g + theme(title = element_text(face="bold"),
                       axis.title.x = element_text(face="bold"),
                       axis.title.y = element_text(face="bold"))
        
        theme_set(theme_bw())
        g1 <- ggplot(results, aes(x=ntree,y=performance)) + geom_point() 
        g1 <- g1 + geom_boxplot(aes(group=ntree,       ## x axis cathegorical values  
                                    fill=132))        ## 132 -> steelblue color
        g1 <- g1 + theme(legend.position = "none")    ## removing legend
        g1 <- g1 + labs(title ="\t\t Waiting time",y="Time (secs)")
        g1 <- g1 + theme(title = element_text(face="bold"),  
                         axis.title.x = element_text(face="bold"),
                         axis.title.y = element_text(face="bold"))
        
        grid.arrange(g,g1,nrow=1)
        
} 

## cforest does not depends on ntree
cforestntree(output,input,nsims = 50)

#######################################################################
## simple random forest using RandomForest package
set.seed(3333)

initialTime <- Sys.time()
rforest <- randomForest(output[train]~.,
                        data=as.data.frame(input[train,]))
rforest.pred <- predict(rforest,newdata = as.data.frame(input[-train,]))
mss <- mean(sum(rforest.pred-output[-train])^2)
mse <- sqrt(mss)
mse
finalTime <- Sys.time() - initialTime
finalTime

plot(rforest.pred,output[-train],
     pch=20,cex=1.5,col="steelblue")
abline(0,1, col="red",lwd=2)
