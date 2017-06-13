#################################################################
# @Course: Systems Biology of Disease 2017                                                                            
# @Rscript: case1.r                                                                                                                  
# @Author: Adrian de Lomana and Chris Plaisier                                                            
# This source code is distributed under GNU GPL v3.0    
#################################################################

###
### This script performs data exploration and classification of  single cell transcriptomes from melanoma malignant cells.
###

# 0.1. loading external packages. if library is missing in your computer, use the R command: install.packages('name of the library')
library(Rtsne) # implementation of t-SNE algorithm
library(RColorBrewer) # library to access easily multiple colors for plotting
library(scatterplot3d) # library for static 3D plotting
library(rgl) # 3D Visualization Using OpenGL. If it fails in your Mac OS X, make sure you have X11 installed (https://www.xquartz.org/)
library(tictoc) # library to profile execution time
library(cluster) # library with functions to perform clustering

# 0.2. defining working directory
setwd('/Users/alomana/github/stagona/2017/case.1/') # this line should be edited accordingly to your working directory. Type "getwd()" to know where you're at in the tree of directories
#setwd('/Users/adriandelomana/git/stagona/2017/case.1/')

# 1. reading the data and metadata for malignant cells
print('reading and treating data...')
originalData=read.csv('data/malignant.8k.genes.data.csv',header=TRUE,row.names=1)
expression=as.data.frame(t(originalData)) # transposing the original data into the appropriate form: 977 observations in a ~8k dimensional space
tumorMetadata=read.csv('data/malignant.8k.genes.tumorMetadata.csv',header=TRUE,row.names=1)

# 2. dimensionality reduction analysis
# 2.0. setting some variables for plotting
tumorLabels=as.character(tumorMetadata$tumor.label)
plottingColors=brewer.pal(length(unique(tumorLabels)),'Dark2')
names(plottingColors)=unique(tumorLabels)
# 2.1. PCA
print('running PCA...')
tic()
results=prcomp(expression) # this step takes about a minute in my laptop...
toc()
# plotting
pdf('figure.malignantCells.pca.pdf')
plot(results$x[,'PC1'],results$x[,'PC2'],main='PCA of malignant cells',col=plottingColors[tumorLabels],pch=19,xlab='PC 1',ylab='PC 2')
legend('bottomright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, otherwise the PDF file of the figure won't be ready to be opened and you'll get an error
                                        
# 2.2. t-SNE
print('running t-SNE...')
# 2.2.1. low resolution of t-SNE
print('running 2D t-SNE with small perplexity...')
tic()
results=Rtsne(expression,dims=2,perplexity=5,verbose=TRUE) # this step takes about a minute in my laptop...
toc()
# plotting
pdf('figure.malignantCells.tsne.2d.low.perplexity.pdf')
plot(results$Y,main='tSNE of malignant cells, p=5',col=plottingColors[tumorLabels],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, otherwise the PDF file of the figure won't be ready to be opened and you'll get an error
                                        
# 2.2.2. high resolution of t-SNE
print('running 2D t-SNE with large perplexity...')
tic()
results2D=Rtsne(expression,dims=2,perplexity=50,verbose=TRUE,theta=0) # this step takes a bit more than a minute in my laptop...
toc()
# plotting
pdf('figure.malignantCells.tsne.2d.hires.pdf')
plot(results2D$Y,main='tSNE of malignant cells, p=50',col=plottingColors[tumorLabels],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, otherwise the PDF file of the figure won't be ready to be opened and you'll get an error
# 2.2.3. high resolution of t-SNE showing 3D
print('running 3D t-SNE with large perplexity...')
tic()
results3D=Rtsne(expression,dims=3,perplexity=50,verbose=TRUE,theta=0) # this step takes a bit more than a minute in my laptop...
toc()
# static 3D scatter plot of tSNE
pdf('figure.malignantCells.tsne.3d.hires.pdf')
scatterplot3d(results3D$Y[,1],results3D$Y[,2],results3D$Y[,3],color=plottingColors[tumorLabels],xlab='tSNE Component 1',ylab='tSNE 2',zlab='tSNE 3',pch=19) 
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
dev.off() # don't forget this command, otherwise the PDF file of the figure won't be ready to be opened and you'll get an error
                                        
# interactive 3D scatter plot. You will need to close the interactive plot window when you're done with it to get back to the code
par3d(windowRect=c(50,50,700,700))
plot3d(results3D$Y[,1],results3D$Y[,2],results3D$Y[,3],col=plottingColors[tumorLabels],xlab='tSNE Component 1',ylab='tSNE Component2',zlab='tSNE Component 3',pch=19) 
legend3d('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])

# let's make a movie from the data! Again, you will need to close the interactive window when you're done with it to get back here
par3d(windowRect=c(50,50,700,700))
plot3d(results3D$Y[,1],results3D$Y[,2],results3D$Y[,3],col=plottingColors[tumorLabels],xlab='tSNE Component 1',ylab='tSNE Component 2',zlab='tSNE Component 3',pch=19) 
legend3d('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])
play3d(spin3d(),duration=20)

# 3. clustering in low dimentional space

# 3.1. definining optimal number of clusters
bestPartition=NbClust(results3D$Y,distance="euclidean",min.nc=2,max.nc=10,method="centroid",index="all")

possibilities=c(4:6)
results=c()
for(k in possibilities) {
  groups=kmeans(results2D$Y,centers=k)
  results=c(results,groups)
}
print(results)
silhouette()
apply(possibilities,2,function(x) print(x))
#transposedObject=as.data.frame(t(results3D$Y)) # surprisingly pvclust clusters columns, not rows
#fit=pvclust(transposedObject,method.hclust="ward.D2",method.dist="euclidean",nboot=10)
#plot(fit) # dendogram with p values
