rm(list=ls())

#Run PCA on blast data to effect a dimensionality reduction.

library(ggplot2)  	#Visualization
library(factoextra)  	#PCA visualization
library(caret)		#one-hot encoding for categorical variables

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


# Define the formula for one-hot encoding
dummies_model <- dummyVars(~ ., data = data)
data_encoded <- predict(dummies_model, newdata = data)
data_encoded <- as.data.frame(data_encoded)

scaled_data <- scale(data_encoded)

pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)

summary(pca_result)


# Scree plot
fviz_eig(pca_result)

# Biplot
fviz_pca_biplot(pca_result, repel = TRUE)

# Extracting principal components
principal_components <- pca_result$x
head(principal_components)




