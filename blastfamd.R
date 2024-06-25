rm(list=ls())

#Run Factorial Analysis of Mixed Data (FAMD) on blast data to effect a dimensionality reduction.

require(FactoMineR)
require(vcd)
require(factoextra)

data.blast <- read.csv("BLAST_ePAL_SeqID_NoGPS.csv")

data <- data.blast[,-c(1:3)]

for(col in 1:dim(data)[2]) if(typeof(data[,col])=='character'){data[,col]<-as.factor(data[,col])}

# Identify columns with NA values
columns_with_na <- colnames(data)[apply(data, 2, anyNA)]

# Remove columns with NA values
data_clean <- data[, !colnames(data) %in% columns_with_na]

data <- data_clean
rm(data_clean)

# Identify columns with only one level
single_level_factors <- sapply(data, function(x) is.factor(x) && length(levels(x)) < 2)

# Remove columns with only one level
data <- data[, !single_level_factors]

famd <- FAMD(data)

#Multiple Correspondence Analysis:
mca <- MCA(data)
# ^^ probably need to run just on the categorical variables


