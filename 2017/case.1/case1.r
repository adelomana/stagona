#################################################################
# @Course: Systems Biology of Disease                                                                            
# @Rscript: case1.r                                                                                                                  
# @Author: Adrian de Lomana and Chris Plaisier                                                            
# This source code is distributed under GNU GPL v3.0                                                  
# the text of which is available at:                            
# http://choosealicense.com/licenses/gpl-3.0/                   
#################################################################

###
### This script performs data exploration and classification of  single cell transcriptomes from melanoma malignant cells.
###

# 0.1. loading external packages. if library is missing in your computer, use the R command: install.packages('name of the library')
library(Rtsne) # implementation of t-SNE algorithm
library(RColorBrewer) # library to access easily multiple colors for plotting
library(scatterplot3d) # library for static 3D plotting
library(rgl) # 3D Visualization Using OpenGL. If it fails in your Mac OS X, make sure you have X11 installed (https://www.xquartz.org/)
library(tictoc) # a library to profile execution time

# 0.2. defining working directory
setwd('/Users/alomana/github/stagona/2017/case.1/') # this line should be edited accordingly to your working directory. Type "getwd()" to know where you're at

# 1. reading the data and metadata for malignant cells
print('reading and treating data...')
expression=read.csv('data/malignant.8k.genes.data.csv',header=TRUE,row.names=1)
transposedExpression=as.data.frame(t(expression)) # transposing the original data in the appropriate form
tumorMetadata=read.csv('data/malignant.8k.genes.tumorMetadata.csv',header=TRUE,row.names=1)

# 2. dimensionality reduction analysis
# 2.0. setting some variables for plotting
tumorLabels=as.character(tumorMetadata$tumor.label)
plottingColors=brewer.pal(length(unique(tumorLabels)),'Dark2')
names(plottingColors)=unique(tumorLabels)
# 2.1. PCA
print('running PCA...')
tic()
pcaResults=prcomp(transposedExpression)
toc()
# plotting
pdf('figure.malignantCells.pca.pdf')
plot(pcaResults$x[,'PC1'],pcaResults$x[,'PC2'],main='PCA of malignant cells',col=plottingColors[tumorLabels],pch=19,xlab='PC 1',ylab='PC 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, it closes the PDF file, so you can cleanly open it!
                                        
# 2.2. t-SNE
print('running t-SNE...')
# 2.2.1. low resolution of t-SNE
print('running 2D t-SNE with small perplexity...')
tic()
results=Rtsne(transposedExpression,dims=2,perplexity=5,verbose=TRUE)
toc()
# plotting
pdf('figure.malignantCells.tsne.2d.lowres.pdf')
plot(results$Y,main='tSNE of malignant cells, p=5',col=plottingColors[tumorLabels],pch=19,xlab='tSNE 1',ylab='tSNE 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, it closes the PDF file, so you can cleanly open it!
                                        
# 2.2.2. high resolution of t-SNE
print('running 2D t-SNE with large perplexity...')
tic()
results=Rtsne(transposedExpression,dims=2,perplexity=50,verbose=TRUE)
toc()
# plotting
pdf('figure.malignantCells.tsne.2d.hires.pdf')
plot(results$Y,main='tSNE of malignant cells, p=50',col=plottingColors[tumorLabels],pch=19,xlab='tSNE 1',ylab='tSNE 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, it closes the PDF file, so you can cleanly open it!
# 2.2.3. high resolution of t-SNE showing 3D
print('running 3D t-SNE with large perplexity...')
tic()
results=Rtsne(transposedExpression,dims=3,perplexity=50,verbose=TRUE)
toc()
# static 3D scatter plot of tSNE
pdf('figure.malignantCells.tsne.3d.hires.pdf')
scatterplot3d(results$Y[,1],results$Y[,2],results$Y[,3],color=plottingColors[tumorLabels],xlab='tSNE 1',ylab='tSNE 2',zlab='tSNE 3',pch=19) 
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, it closes the PDF file, so you can cleanly open it!
                                        
# interactive 3D scatter plot. You will need to close the interactive plot window when you're done with it
par3d(windowRect=c(50,50,700,700))
plot3d(results$Y[,1],results$Y[,2],results$Y[,3],col=plottingColors[tumorLabels],xlab='tSNE 1',ylab='tSNE 2',zlab='tSNE 3',pch=19) 
legend3d('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])

# let's make a play from the data! Again, you will need to close the interactive plot window when you're done with it
par3d(windowRect=c(50,50,700,700))
plot3d(results$Y[,1],results$Y[,2],results$Y[,3],col=plottingColors[tumorLabels],xlab='tSNE 1',ylab='tSNE 2',zlab='tSNE 3',pch=19) 
legend3d('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
play3d(spin3d(),duration=20)

