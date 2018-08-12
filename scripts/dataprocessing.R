library(dplyr)


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


################################################################################
## prereading data
tabla <- read.csv("../data/TablaCompletaCohorteINER_C08-04_20180810.csv",
                  stringsAsFactors = FALSE)
str(tabla)
dim(tabla)

## drop log cv redundant with cv values
tabla.p <- tabla[,-grep("LogCV_",colnames(tabla))]
## drop cv values with values < 50 == NA
tabla.p <- tabla.p[!is.na(tabla$CV_S052),]
## Drop non meaningful data specified in metadata
tabla.p <- tabla.p[, colnames(tabla.p) %in% meaningfulVars]
dim(tabla.p)
str(tabla.p)

## drop samples without output delta CD4
tabla.ov <- tabla.p[!is.na(tabla.p$Delta_CD4_year1),]
dim(tabla.ov)

## selecting columns without NAs for the model
NA_nb <- sapply(tabla.ov, function(x) sum(is.na(x))) 
varNames <- names(NA_nb)[NA_nb==0]
length(varNames)
varNames
## data for the model
data.m <- tabla.ov[,varNames]
dim(data.m[complete.cases(data.m),])
str(data.m)

## Addtional cleaning
## dropping date fecha_S0
data.m <- select(data.m, -Fecha_S0)
data.m$IRIS_manifestation_2 <- sapply(data.m$IRIS_manifestation_2,
                                      function(x) ifelse(x=="","Ninguno",x))
## change non recorded to ninguno
table(data.m$IRIS_manifestation_2)
str(data.m)
dim(data.m)

write.csv(data.m,"../data/cleandata.csv",row.names = FALSE)
