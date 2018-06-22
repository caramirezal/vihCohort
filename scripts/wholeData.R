
data <- read.csv("../data/TablaCompletaCohorteINER_C08-04_20170630.csv",
                   header = TRUE,check.names = FALSE,
                   nrows = 3,
                   colClasses = "character",skipNul = TRUE)

data.r <- readLines("../data/TablaCompletaCohorteINER_C08-04_20170630.csv",n = 3)

print("hola")
