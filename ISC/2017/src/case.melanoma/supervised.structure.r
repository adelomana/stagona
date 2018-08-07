#################################################################
# @Course: Systems Biology of Disease 2017                                                                            
# @Rscript: supervised.r                                                                                                                  
# @Author: Adrian de Lomana and Chris Plaisier                                                            
# This source code is distributed under GNU GPL v3.0    
#################################################################

###
### These comments define the steps to perform supervised learning from melanoma non-malignant single cell transcriptomes
###

# 0.1. loading external packages. if library is missing in your computer, use the R command: install.packages('name of the library')
library(Rtsne) # implementation of t-SNE algorithm
library(RColorBrewer) # library to access easily multiple colors for plotting
library(tictoc) # library to profile execution time
library(scatterplot3d) # library for static 3D plotting
library(randomForest) # a library that includes an implementation of a random forest algorithm

# 0.2. user-specific definition of paths.
# please navigate into "case.melanoma" through the "Files" tab on the bottom right, and set that directory as your working directory
dataDirectory='../../data/case.melanoma/formatted/' 

# 1. reading the data and metadata for malignant cells
print('reading and treating data...')s
dataFilePath=paste(dataDirectory,'nonMalignant.2kgenes.data.testing.csv',sep='')
metadataFilePath=paste(dataDirectory,'nonMalignant.2kgenes.tumorMetadata.testing.csv',sep='')
originalData=read.csv(dataFilePath,header=TRUE,row.names=1)
# names(expression)=make.names(names(expression)) <--- formatting required for smooth running of random forests

# 2. dimensionality reduction of original data

# 3. supervised learning on cell type classify data
# 3.1. read prior information about cell types, both profiles and labels from other tumor samples
immuneProfilesDataFile=paste(dataDirectory,'nonMalignant.2kgenes.data.training.csv',sep='')
rawData=read.csv(immuneProfilesDataFile,row.names=1,header=T)
names(immuneProfiles)=make.names(names(immuneProfiles))

knownImmuneLabelsFile=paste(dataDirectory,'nonMalignant.2kgenes.immuneMetadata.training.csv',sep='')
immuneLabels=read.csv(knownImmuneLabelsFile,row.names=1,header=T)
knownLabels=immuneLabels$immune.label

# 3.2. perform learning
# consider random forest, a good compromise between accuracy and computational speed. But feel free to use any other method!

# 3.3. plotting the results

# 3.4. confusion matrix and classification features

# confusion matrix

# list best genes for classification

# best classifiers distributions












