#################################################################
# @Course: Systems Biology of Disease                           #
# @Rscript: dimensionalityReduction.R                           #
# @Author: Adrian de Lomana                                     #
#                                                               #
# This source code is distributed under GNU GPL v3.0            #
# the text of which is available at:                            #
# http://choosealicense.com/licenses/gpl-3.0/                   #
#################################################################

###
### this script performs the analysis of single cell tumor cells
###

# 0.1. loading external packages. if library is missing in your computer, use the R command: install.packages('name of the library')
library(Rtsne) # implementation of t-SNE algorithm
library(RColorBrewer) # library to access easily multiple colors for plotting
library(scatterplot3d) # library for static 3D plotting
library(rgl) # 3D Visualization Using OpenGL

# 0.2. user defined variables
setwd('/Users/alomana/gDrive2/education/isb.summer.course.2016/src/pmid27124452') # this line should be edited accordingly to your working directory. Type "getwd()" to know where you're at.
inputFileName='data/formatted/tumorCells.data.csv'
metadataTumor='data/formatted/tumorCells.metadata.csv'

# 1. reading the data and metadata
print('reading and treating data...')
originalData=read.csv(inputFileName,header=TRUE,row.names=1)
metadata=read.csv(metadataTumor,header=TRUE,row.names=1)
# 1.1. transposing the original data in file into appropriate expression data frame
expression=as.data.frame(t(originalData))
                                        
# 2. dimensionality reduction analysis
tumorLabels = as.character(metadata$tumor.label)
plottingColors=brewer.pal(length(unique(tumorLabels)),'Dark2')
names(plottingColors)=unique(tumorLabels)
                                        
# 2.1. PCA
print('running PCA...')
pcaResults=prcomp(expression)
# plotting
pdf('figure.pca.pdf')
plot(pcaResults$x[,'PC1'],pcaResults$x[,'PC2'],main='PCA',col=plottingColors[tumorLabels],pch=19,xlab='PC 1',ylab='PC 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, it closes the PDF file, so you can cleanly open it!
                                        
# 2.2. t-SNE
print('running t-SNE...')
                                        
# 2.2.1. low resolution of t-SNE
print('running 2D t-SNE with small perplexity...')
results=Rtsne(expression,dims=2,theta=0.,perplexity=5,verbose=TRUE)
# plotting
pdf('figure.tsne.2d.lowres.pdf')
plot(results$Y,main='tSNE, p=5',col=plottingColors[tumorLabels],pch=19,xlab='tSNE 1',ylab='tSNE 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, it closes the PDF file, so you can cleanly open it!
                                        
# 2.2.2. high resolution of t-SNE
print('running 2D t-SNE with large perplexity...')
results=Rtsne(expression,dims=2,theta=0.,perplexity=50,verbose=TRUE)
# plotting
pdf('figure.tsne.2d.hires.pdf')
plot(results$Y,main='tSNE, p=50',col=plottingColors[tumorLabels],pch=19,xlab='tSNE 1',ylab='tSNE 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, it closes the PDF file, so you can cleanly open it!
                                        
# 2.2.3. high resolution of t-SNE showing 3D
print('running 3D t-SNE with large perplexity...')
results=Rtsne(expression,dims=3,theta=0.,perplexity=50,verbose=TRUE)
# static 3D scatter plot of tSNE
pdf('figure.tsne.3d.hires.pdf')
scatterplot3d(results$Y[,1],results$Y[,2],results$Y[,3],color=plottingColors[tumorLabels],xlab='tSNE 1',ylab='tSNE 2',zlab='tSNE 3',pch=19) 
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, it closes the PDF file, so you can cleanly open it!
                                        
# interactive 3D scatter plot (you will need to close the interactive plot window when you're done with it)
par3d(windowRect=c(50,50,700,700))
plot3d(results$Y[,1],results$Y[,2],results$Y[,3],col=plottingColors[tumorLabels],xlab='tSNE 1',ylab='tSNE 2',zlab='tSNE 3',pch=19) 
legend3d('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
legend3d('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)]) # do legend twice because of bug in RGL

# let's make a play from the data!
par3d(windowRect=c(50,50,700,700))
plot3d(results$Y[,1],results$Y[,2],results$Y[,3],col=plottingColors[tumorLabels],xlab='tSNE 1',ylab='tSNE 2',zlab='tSNE 3',pch=19) 
legend3d('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
legend3d('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)]) # do legend twice because of bug in RGL
play3d(spin3d(),duration=20)
