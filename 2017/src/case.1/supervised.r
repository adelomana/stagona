#################################################################
# @Course: Systems Biology of Disease 2017                                                                            
# @Rscript: supervised.r                                                                                                                  
# @Author: Adrian de Lomana and Chris Plaisier                                                            
# This source code is distributed under GNU GPL v3.0    
#################################################################

###
### This script performs data supervised learning from melanoma non-malignant single cell transcriptomes
###

# 0.1. loading external packages. if library is missing in your computer, use the R command: install.packages('name of the library')

# 0.2. user-specific definition of paths.
# these lines should be edited accordingly to your working directory. Type "getwd()" to know where you're at in the tree of directories
sourceDirectory='/Users/alomana/github/stagona/2017/src/case.1/' # this line should be edited to match your hierarchy of dirs
dataDirectory='/Users/alomana/github/stagona/2017/data/case.1/formatted/' # this line should be edited to match your hierarchy of dirs
setwd(sourceDirectory) # no need to edit this line :-)

# 1. reading the data and metadata for malignant cells
print('reading and treating data...')
dataFilePath=paste(dataDirectory,'nonMalignant.8kgenes.data.csv',sep='')
metadataFilePath=paste(dataDirectory,'nonMalignant.8kgenes.tumorMetadata.csv',sep='')
originalData=read.csv(dataFilePath,header=TRUE,row.names=1)
expression=as.data.frame(t(originalData)) # transposing the original data into the appropriate form: 2,678 observations in a 8,000 dimensional space
tumorMetadata=read.csv(metadataFilePath,header=TRUE,row.names=1)

# explore data, run a t-SNE and color based on tumor

# classify data