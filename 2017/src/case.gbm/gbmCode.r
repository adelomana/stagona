#################################################################
# @Course: Systems Biology of Disease                           #
# @Rscript: gbmCode.R                                             #
# @Version: 1                                                   #
# @Author: Chris Plaisier and Adrian de Lomana                  #
# @Sponsored by:                                                #
# Institute for Systems Biology                                 #
# 1441 North 34th Street                                        #
# Seattle, Washington  98103-8904                               #
# (216) 732-2139                                                #
#                                                               #
# Copyright (C) 2016 by Institute for Systems Biology,          #
# Seattle, Washington, USA.  All rights reserved.               #
#                                                               #
# This source code is distributed under the GNU Lesser          #
# General Public License, the text of which is available at:    #
#   http://www.gnu.org/copyleft/lesser.html                     #
#################################################################

#####################################
### Libraries needed for analysis ###
#####################################
library(RColorBrewer)
library(scatterplot3d)
library(rgl)
library(Rtsne)

##############################
### Load data for analysis ###
##############################
# Set working directory

# please navigate into "case.gbm" through the "Files" tab on the bottom right, and set that directory as your working directory
# Read in expression data from GEO GSE57872 for glioblastoma single cells
d1 = read.csv('../../data/case.gbm/GSE57872_GBM_data_matrix.csv', header=1, row.names=1)

##################################
### Exploring single cell data ###
##################################
# Prepare to colorize based on single cell tumor origin
colNames1 = sapply(colnames(d1)[1:430], function(x) { strsplit(x,'_')[[1]][1] } )
colNames1[which(colNames1=='MGH264')] = 'MGH26' # MGH264 is MGH26
col1 = brewer.pal(5, "Dark2")
names(col1) = unique(colNames1)
cols1 = c(as.character(paste(col1[colNames1],'66',sep='')), col1)
cex1 = c(rep(1,430),rep(2,5))

#################
## Compute PCA ##
#################
pca1 = prcomp(t(d1))
 
# Static 2D plot of PCA
pdf('pca_2D_tumorColorized.pdf')
plot(pca1$x[,'PC1'], pca1$x[,'PC2'], main="GBM single cell PCA", col=cols1, pch=19, cex=cex1)
legend('topright', legend=unique(colNames1), fill=col1[unique(colNames1)], border=F, title='Tumor')
dev.off()

# Static 3D scatter plot of PCA
pdf('pca_3D_tumorColorized.pdf')
scatterplot3d(pca1$x[,1], pca1$x[,2], pca1$x[,3], color=cols1, pch=19, type='p', cex.symbols=cex1, xlab='PCA1', ylab='PCA2', zlab='PCA3')
legend('topleft', legend=unique(colNames1), fill=col1[unique(colNames1)], border=F, title='Tumor')
dev.off()

# Interactive 3D scatter plot of PCA
par3d(windowRect=c(50,50,700,700))
plot3d(pca1$x[,1], pca1$x[,2], pca1$x[,3], col=cols1, type='s', size=cex1, xlab='PCA1', ylab='PCA2', zlab='PCA3')
legend3d('topright', legend=unique(colNames1), fill=col1[unique(colNames1)])
legend3d('topright', legend=unique(colNames1), fill=col1[unique(colNames1)]) # Do legend twice because of bug in RGL
play3d(spin3d(), duration=20)
# Need to close window by hand

##################
## Compute tSNE ##
##################
# 2D tSNE
rtsne2d = Rtsne(t(as.matrix(d1)))
# Static 2D plot of tSNE
pdf('tSNE_2D_tumorColorized.pdf')
plot(rtsne2d$Y, main="GBM single cell tSNE", col=cols1, pch=19, cex=cex1)
legend('topright', legend=unique(colNames1), fill=col1[unique(colNames1)], border=F, title='Tumor')
dev.off()

# 3D tSNE
rtsne3d = Rtsne(t(as.matrix(d1)),dim=3)
# Static 3D scatter plot of tSNE
pdf('tSNE_3D_tumorColorized.pdf')
library(scatterplot3d)
scatterplot3d(rtsne3d$Y[,1], rtsne3d$Y[,2], rtsne3d$Y[,3], color=cols1, pch=19, type='p', cex.symbols=cex1, xlab='tSNE1', ylab='tSNE2', zlab='tSNE3')
legend('topleft', legend=unique(colNames1), fill=col1[unique(colNames1)], border=F, title='Tumor')
dev.off()

# Interactive 3D scatter plot
par3d(windowRect=c(50,50,700,700))
plot3d(rtsne3d$Y[,1], rtsne3d$Y[,2], rtsne3d$Y[,3], col=cols1, type='s', size=cex1,xlab='tSNE1',ylab='tSNE2',zlab='tSNE3')
legend3d('topright', legend=unique(colNames1), fill=col1[unique(colNames1)])
legend3d('topright', legend=unique(colNames1), fill=col1[unique(colNames1)]) # Do legend twice because of bug in RGL
play3d(spin3d(), duration=20)
# Need to close window by hand


#############################
### supervised learning ###
#############################

# separate training and testing data sets

#  perform learning

# plotting the results

# confusion matrix and classification features

# confusion matrix

# list best genes for classification

# best classifiers distributions


