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
library(Rtsne) # implementation of t-SNE algorithm
library(RColorBrewer) # library to access easily multiple colors for plotting
library(tictoc) # library to profile execution time
library(scatterplot3d) # library for static 3D plotting
library(caret) # a nice library to consider for supervised learning

# 0.2. user-specific definition of paths.
# please navigate into "case.melanoma" through the "Files" tab on the bottom right, and set that directory as your working directory
dataDirectory='../../data/case.melanoma/formatted/' 

# 1. reading the data and metadata for malignant cells
print('reading and treating data...')
dataFilePath=paste(dataDirectory,'nonMalignant.2kgenes.data.prediction.csv',sep='')
metadataFilePath=paste(dataDirectory,'nonMalignant.2kgenes.tumorMetadata.prediction.csv',sep='')
originalData=read.csv(dataFilePath,header=TRUE,row.names=1)
expression=as.data.frame(t(originalData)) # transposing the original data into the appropriate form: 2,249 observations in a 8,000 dimensional space
tumorMetadata=read.csv(metadataFilePath,header=TRUE,row.names=1)

# 2. dimensionality reduction of original data
# 2.1. setting some variables for plotting
tumorLabels=as.character(tumorMetadata$tumor.label)
plottingColors=brewer.pal(length(unique(tumorLabels)),'Dark2')
names(plottingColors)=unique(tumorLabels)
# 2.2. high resolution of t-SNE
print('running 2D t-SNE with large perplexity...')
tic()
results2D=Rtsne(expression,dims=2,perplexity=50,verbose=TRUE,theta=0) 
toc()
# plotting
pdf('figure.non.malignantCells.tsne.2d.pdf')
plot(results2D$Y,main='tSNE of non malignant cells, p=50',col=plottingColors[tumorLabels],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, otherwise the PDF file of the figure won't be ready to be opened and you'll get an error

######## to be removed
dataFilePath=paste(dataDirectory,'nonMalignant.2kgenes.immuneMetadata.prediction.csv',sep='')
immuneMetadata=read.csv(dataFilePath,header=TRUE,row.names=1)
labels=as.character(immuneMetadata$immune.label)
plottingColors=brewer.pal(length(unique(labels)),'Dark2')
names(plottingColors)=unique(labels)
pdf('figure.non.malignantCells.immuneLabels.pdf')
plot(results2D$Y,main='tSNE of non malignant cells, p=50',col=plottingColors[labels],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(labels),fill=plottingColors[unique(labels)])
dev.off() # don't forget this command, otherwise the PDF file of the figure won't be ready
#############

# 3. supervised learning on cell type classify data
# 3.1. read prior information about cell types, both profiles and labels
immuneProfilesDataFile=paste(dataDirectory,'nonMalignant.2kgenes.data.learning.csv',sep='')
rawData=read.csv(immuneProfilesDataFile,row.names=1,header=T)
immuneProfiles=as.data.frame(t(rawData))

knownImmuneLabelsFile=paste(dataDirectory,'nonMalignant.2kgenes.immuneMetadata.learning.csv',sep='')
immuneLabels=read.csv(knownImmuneLabelsFile,row.names=1,header=T)
knownLabels=immuneLabels$immune.label

# 3.2. perform learning
fitControl=trainControl(method="cv",number=10,verboseIter=TRUE,classProbs=TRUE)
learningMethod='rFerns'

tic()
model=train(immuneProfiles,knownLabels,method=learningMethod,trControl=fitControl,verbose=TRUE) # tuneLength = 4,metric = 'ROC'
predictions=predict.train(object=model,expression,type="raw")
toc()
table(predictions)

# 3.3. plotting the results
figureFileName=paste('figure.non.malignantCells.classification.',learningMethod,'.pdf',sep='')
pdf(figureFileName)
predictedLabels=as.character(predictions)
plottingColors=brewer.pal(length(unique(predictedLabels)),'Dark2')
names(plottingColors)=unique(predictedLabels)
plot(results2D$Y,main='learned classification',col=plottingColors[predictedLabels],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('bottomleft',legend=unique(predictedLabels),fill=plottingColors[unique(predictedLabels)])
dev.off() # don't forget this command, otherwise the PDF file of the figure won't be ready

# 3.4. describe learning: error and features (gene markers, possibly a heatmap of markers)
# for students, consider pca and learning biomarkers in 200D, then define biomarkers through t-test
vimp <- varImp(model, scale = TRUE)
vimp <- vimp[order(vimp$Overall,decreasing = TRUE),drop = FALSE]







