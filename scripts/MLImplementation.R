## This script contains the R functions to implement Machine Learning, parameter selection,
## and plotting results in technicalReportML.Rmd 

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
        
        # return(results)
        
        theme_set(theme_bw())
        g <- ggplot(results, aes(x=mtry,y=Accuracy)) + geom_point() 
        g <- g + geom_boxplot(aes(group=mtry,fill="steelblue"))
        g <- g + theme(legend.position = "none") 
        g <- g + labs(title =paste(method,"- Accuracy"),y="MSE")
        g <- g + theme(title = element_text(face="bold"),
                       axis.title.x = element_text(face="bold"),
                       axis.title.y = element_text(face="bold"))
        
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