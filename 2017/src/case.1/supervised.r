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
# these lines should be edited accordingly to your working directory. Type "getwd()" to know where you're at in the tree of directories
sourceDirectory='/Users/alomana/github/stagona/2017/src/case.1/' # this line should be edited to match your hierarchy of dirs
dataDirectory='/Users/alomana/github/stagona/2017/data/case.1/formatted/' # this line should be edited to match your hierarchy of dirs
setwd(sourceDirectory) # no need to edit this line :-)

# 1. reading the data and metadata for malignant cells
print('reading and treating data...')
dataFilePath=paste(dataDirectory,'nonMalignant.8kgenes.data.prediction.csv',sep='')
metadataFilePath=paste(dataDirectory,'nonMalignant.8kgenes.tumorMetadata.prediction.csv',sep='')
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
results2D=Rtsne(expression,dims=2,perplexity=50,verbose=TRUE,theta=0) # this step takes a bit more than 5 minutes in my laptop...
toc()
# plotting
pdf('figure.non.malignantCells.tsne.2d.pdf')
plot(results2D$Y,main='tSNE of non malignant cells, p=50',col=plottingColors[tumorLabels],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, otherwise the PDF file of the figure won't be ready to be opened and you'll get an error

######## to be removed
dataFilePath=paste(dataDirectory,'nonMalignant.8kgenes.immuneMetadata.prediction.csv',sep='')
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
immuneProfilesDataFile=paste(dataDirectory,'nonMalignant.8kgenes.data.learning.csv',sep='')
rawData=read.csv(immuneProfilesDataFile,row.names=1,header=T)
immuneProfiles=as.data.frame(t(rawData))

knownImmuneLabelsFile=paste(dataDirectory,'nonMalignant.8kgenes.immuneMetadata.learning.csv',sep='')
immuneLabels=read.csv(knownImmuneLabelsFile,row.names=1,header=T)

# 3.2. prepare data for learning
miniProfiles=immuneProfiles[1:400,]
miniLabels=as.factor(immuneLabels$immune.label[1:400])

# 3.2. perform learning
library(parallel)
library(doParallel)
workingCores <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(workingCores)
fitControl <- trainControl(method = "cv",
                           number = 10,
                           allowParallel = TRUE)
learningMethod='rf'

tic()
model=train(miniProfiles,miniLabels,method=learningMethod,trControl = fitControl)
predictions=predict.train(object=model,expression,type="raw")
toc()
table(predictions)

# 3.3. plotting the results
figureFileName=paste('classification.',learningMethod,'.pdf',sep='')
pdf(figureFileName)
predictedLabels=as.character(predictions)
plottingColors=brewer.pal(length(unique(predictedLabels)),'Dark2')
names(plottingColors)=unique(predictedLabels)
plot(results2D$Y,main='learned classification',col=plottingColors[predictedLabels],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(predictedLabels),fill=plottingColors[unique(predictedLabels)])
dev.off() # don't forget this command, otherwise the PDF file of the figure won't be ready


# check random forest and support vector machines, also bayes and self organizing maps
# knn, 398.64 sec elapsed
# lda, 1281.324 sec elapsed
# random forest, rf
# neural networks, nnet, 
# parallel random forest, parRF
# support vector machines
# bayes, nb, 
# self organizing maps

# consider using 2k genes if time is above 5 min














# Classify using Single Cell Predictor
scp1 = cor(immuneProfiles,expression, method='spearman')
rownames(scp1)
scp1_calls = sapply(colnames(scp1), function(x) { rownames(scp1)[which(scp1[,x]==max(scp1[,x]))] })
scp2_calls = sapply(names(scp1_calls), function(x) { ifelse(scp1[scp1_calls[x],x]>=0.1, scp1_calls[x], NA) } )
table(scp2_calls) # Tabulation of cells to each type

######################################
### Significance of predcitions by ###
### resampling of genes, and calls ###
### based on significance.         ###
###                                ###
### Cell type calls made based on: ###
###  1. p-value <= 0.05            ###
###  2. R >= 0.1                   ###
######################################
permutations = 100
m1 = matrix(nrow=ncol(scp1), ncol=4)
rownames(m1) = colnames(scp1)
colnames(m1) = c('Call','Cor','p_value','Final Call')
# For each single cell, and bulk tumor
tic() # Takes about 7 minutes
for(i in colnames(scp1)) {
  call = NA
  sub = sapply(1:permutations, function(j) { cor(centroids.immune[,scp1_calls[i]], originalData[sample(rownames(originalData),nrow(centroids.immune)), i], method='spearman') })
  m1[i,1] = scp1_calls[i]
  m1[i,2] = scp1[scp1_calls[i],i]
  m1[i,3] = length(which(sub >= scp1[scp1_calls[i],i]))/permutations
  if(m1[i,3]<=0.05 && scp1[scp1_calls[i],i]>=0.1) {
    call = scp1_calls[i]
  }
  m1[i,4] = call
}     
toc()

table(m1[,4])


# plot results in 3d interactive tsne
###########################################
### Colorize using immune cell analysis ###
###########################################
# Prepare colors for plotting
colorInterpolation = colorRampPalette(brewer.pal(10,'Paired'))
col2 = colorInterpolation(length(sort(unique(m1[,4]))))
names(col2) = sort(unique(m1[,4]))
cols2 = as.character(col2[m1[,4]])

# Static 3D scatter plot of tSNE
pdf('tSNE_3D_tumorColorized.pdf')
library(scatterplot3d)
scatterplot3d(rtsne3d$Y[,1], rtsne3d$Y[,2], rtsne3d$Y[,3], color=cols1, pch=19, type='p', xlab='tSNE1', ylab='tSNE2', zlab='tSNE3')
legend('topleft', legend=unique(tumorMetadata[,1]), fill=col1[unique(tumorMetadata[,1])], border=F)
dev.off()

# Interactive 3D scatter plot
par3d(windowRect=c(50,50,700,700))
plot3d(rtsne3d$Y[,1], rtsne3d$Y[,2], rtsne3d$Y[,3], col=cols2, type='p', xlab='tSNE1',ylab='tSNE2',zlab='tSNE3')
legend3d('topright', legend=sort(unique(m1[,4])), fill=col2[sort(unique(m1[,4]))])
#play3d(spin3d(), duration=10)
# Need to close window by hand









