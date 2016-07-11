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
setwd('/Users/alomana/dropbox.alomana/Dropbox/ISB Summer Course 2016/Day 2_Tuesday/Day 2_Afternoon/finalCode/melanoma')
inputFileName='data/tumorCells.data.csv'
metadataTumor='data/tumorCells.metadata.csv'
                                        # 1. reading the data and metadata
print('reading and treating data...')
originalData=read.csv(inputFileName,header=TRUE,row.names=1)
metadata=read.csv(metadataTumor,header=TRUE,row.names=1)
                                        # 1.1. transposing the original data in file into appropriate expression data frame
expression=as.data.frame(t(originalData))
                                        # 2. dimensionality reduction analysis
pdf('resultsFigures.pdf')
tumorLabels = as.character(metadata$tumor.label)
plottingColors=brewer.pal(length(unique(tumorLabels)),'Dark2')
names(plottingColors)=unique(tumorLabels)
                                        # 2.1. PCA
print('running PCA...')
maxN=min(dim(expression))
trimmedExpression = expression[1:maxN,1:maxN]
pcaResults = princomp(trimmedExpression)$scores[,1:2]

plot(pcaResults,main='PCA',col=plottingColors[tumorLabels],pch=19,xlab='PC 1',ylab='PC 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
                                        # 2.2. t-SNE
print('running t-SNE...')
                                        # 2.2.1. low resolution of t-SNE
results=Rtsne(expression,dims=2,theta=0.,perplexity=5,verbose=TRUE) 
plot(results$Y,main='tSNE, p=5',col=plottingColors[tumorLabels],pch=19,xlab='tSNE 1',ylab='tSNE 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
                                        # 2.2.2. high resolution of t-SNE
results=Rtsne(expression,dims=2,theta=0.,perplexity=50,verbose=TRUE) 
plot(results$Y,main='tSNE, p=50',col=plottingColors[tumorLabels],pch=19,xlab='tSNE 1',ylab='tSNE 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
                                        # 2.2.3. high resolution of t-SNE showing 3D
results=Rtsne(expression,dims=3,theta=0.,perplexity=50,verbose=TRUE)
                                        # static 3D scatter plot of tSNE
scatterplot3d(results$Y[,1],results$Y[,2],results$Y[,3],color=plottingColors[tumorLabels],xlab='tSNE 1',ylab='tSNE 2',zlab='tSNE 3',pch=19) 
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, it closes the PDF file, so you can cleanly open it!
                                        # interactive 3D scatter plot
par3d(windowRect=c(50,50,700,700))
plot3d(results$Y[,1],results$Y[,2],results$Y[,3],col=plottingColors[tumorLabels],xlab='tSNE 1',ylab='tSNE 2',zlab='tSNE 3',pch=19) 
legend3d('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
legend3d('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)]) # do legend twice because of bug in RGL
movie3d(spin3d(),dev = rgl.cur(),duration=10,movie='movie.tSNE.tumorCells', dir=getwd()) # it will generate a file named movie.tSNE.tumorCells.gif at your current working dir to be opened in your web browser

print('... analysis completed!')
