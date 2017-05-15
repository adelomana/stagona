#################################################################
# @Course: Systems Biology of Disease                           #
# @Rscript: case2.R                                             #
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
setwd('/Users/alomana/github/stagona/2017/case.2/')

# Read in expression data from GEO GSE57872 for glioblastoma single cells
d1 = read.csv('data/GSE57872_GBM_data_matrix.csv', header=1, row.names=1)

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
play3d(spin3d(), duration=5)
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
play3d(spin3d(), duration=5)
# Need to close window by hand


#############################
### Single cell predictor ###
#############################
# Load up TCGA GBM sub-type centroids
centroids1 = read.csv('data/gbmTCGACentroids.csv',header=1,row.names=1)

# Comupte 
scp1 = cor(centroids1,d1[rownames(centroids1),],method='spearman')
scp1_calls = sapply(colnames(scp1), function(x) { rownames(scp1)[which(scp1[,x]==max(scp1[,x]))] })
table(scp1_calls)


######################################
### Significance of predcitions by ###
### resampling of genes, and calls ###
### based on significance.         ###
######################################
permutations = 100
m1 = matrix(nrow=ncol(scp1), ncol=9)
rownames(m1) = colnames(scp1)
colnames(m1) = c(paste(rep(rownames(scp1),each=2),c('cor','p_value'),sep='.'),'Calls')
# For each single cell, and bulk tumor
for(i in colnames(scp1)) {
    calls = c()
    sub = sapply(1:permutations, function(j) { cor(centroids1, d1[sample(rownames(d1),nrow(centroids1)), i], method='spearman') })
    for(j in 1:4) {
        m1[i,(j*2-1)] = scp1[j,i]
        m1[i,(j*2)] = length(which(sub[j,] >= scp1[j,i]))/permutations
        if(m1[i,(j*2)]<=0.05) {
            calls = c(calls,rownames(scp1)[j])
        }
    }
    if(length(calls)==0) {
        calls = c('NA')
    }
    m1[i,9] = paste(calls,collapse=' ')
}     

table(m1[,9])
write.csv(m1,'samplingResults.csv')


#####################################
### Making a table of percentages ###
#####################################
d2 = d1[,1:430]
colNames1 = sapply(colnames(d2), function(x) { strsplit(x,'_')[[1]][1] } )
colNames1[which(colNames1=='MGH264')] = 'MGH26'
m2 = matrix(ncol=length(unique(m1[,9])),nrow=length(unique(colNames1)))
rownames(m2) = unique(colNames1)
colnames(m2) = unique(m1[,9])
for(i in rownames(m2)) {
    for(j in colnames(m2)) {
        m2[i,j] = sum(m1[which(colNames1==i),9]==j)/length(m1[which(colNames1==i),9])
    }
}

par(mar=c(3, 4, 4, 2))
col2 = brewer.pal(length(unique(colnames(m2))),"Paired")
names(col2) = unique(colnames(m2))
pdf('subtypesSingleCells.pdf')
barplot(t(m2), beside=T, legend.text=rownames(t(m2)), col=col2, ylim=c(0,1), ylab='Percent of Cells')
dev.off()


#################################
### Explore cells using calls ###
#################################
cols2 = as.character(col2[m1[1:430,9]])

##################
## Compute tSNE ##
##################
# 2D tSNE
rtsne2d.sub = Rtsne(t(as.matrix(d1[rownames(centroids1),1:430])))

# Static 2D plot of tSNE
pdf('tSNE_2D_subtypeColorized.pdf', width=10, height=6)
par(mfrow=c(1,2))
plot(rtsne2d.sub$Y[1:430,], main="GBM single cell tSNE", col=cols2, pch=19,cex=cex1)
plot(rtsne2d.sub$Y[1:430,], main="GBM single cell tSNE", col=cols1[1:430], pch=19, cex=cex1[1:430])
dev.off()

# 3D tSNE
rtsne3d.sub = Rtsne(t(as.matrix(d1[rownames(centroids1),1:430])),dim=3)

# Static 3D scatter plot of tSNE
pdf('tSNE_3D_subtypeColorized.pdf')
scatterplot3d(rtsne3d.sub$Y[1:430,1], rtsne3d.sub$Y[1:430,2], rtsne3d.sub$Y[1:430,3], color=cols2[1:430], pch=19, type='p', cex.symbols=cex1[1:430], xlab='tSNE1', ylab='tSNE2', zlab='tSNE3')
legend('topleft', legend=names(col2), fill=col2, border=F, title='Subtype')
dev.off()

# Interactive 3D scatter plot
par3d(windowRect=c(50,50,700,700))
plot3d(rtsne3d.sub$Y[1:430,1], rtsne3d.sub$Y[1:430,2], rtsne3d.sub$Y[1:430,3], col=cols2[1:430], type='s', size=cex1[1:430], xlab='tSNE1', ylab='tSNE2', zlab='tSNE3')
legend3d('topright', legend=names(col2), fill=col2)
legend3d('topright', legend=names(col2), fill=col2) # Do legend twice because of bug in RGL
play3d(spin3d(), duration=5)
# Need to close window by hand

