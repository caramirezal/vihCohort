## This script contains the R functions to implement Machine Learning, parameter selection,
## and plotting results in technicalReportML.Rmd 

## function for calculating r^2  as 1 - rss/tss
r_sq <- function(observed,predicted) {
        rss <- sum((observed-predicted)^2)
        tss <- sum((observed-mean(observed))^2)
        return( 1 - (rss/tss) )
}

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
                                mss <- mean(sum(forest.pred-output[-intrain])^2)
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
                                mss <- mean(sum(rforest.pred-output[-intrain])^2)
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
                                mss <- mean(sum(forest.pred-output[-intrain])^2)
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
                        mss <- mean(sum(forest.pred-output[-intrain])^2)
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
                mss <- mean(sum(forest.pred-output[-intrain])^2)
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
                mss <- mean(sum(forest.pred-output[-intrain])^2)
                mse <- sqrt(mss)
                
                ## calculate rsq
                rsq <- r_sq(output[-intrain],forest.pred)
                
                ## store results
                results$rforest_accuracy[index] <- mse
                results$rforest_rsq[index] <- rsq
                results$rforest_time[index] <- finalTime
                
        }
}

## writing results
#write.csv(results,"../data/forestTrainSizes.csv",row.names = FALSE)

