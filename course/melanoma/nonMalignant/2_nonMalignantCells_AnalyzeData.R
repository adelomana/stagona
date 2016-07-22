#################################################################
# @Course: Systems Biology of Disease                           #
# @Rscript: dimensionalityReduction.R                           #
# @Author: Adrian de Lomana and Chris Plaisier                  #
#                                                               #
# This source code is distributed under GNU GPL v3.0            #
# the text of which is available at:                            #
# http://choosealicense.com/licenses/gpl-3.0/                   #
#################################################################

################################################################################
### This script explores the expression of single cells from melanoma tumors ###
### and performs classification of immune cells from single cell             ###
### transcriptomes of melanoma tumors.                                       ###
################################################################################

########################
### Loading packages ###
######################## 
# Loading external packages. if library is missing in your computer, use the R command: install.packages('name of the library')
library(Rtsne) # implementation of t-SNE algorithm
library(RColorBrewer) # library to access easily multiple colors for plotting
library(scatterplot3d) # library for static 3D plotting
library(rgl) # 3D Visualization Using OpenGL
library(tictoc)
library(scatterplot3d)

##############################
### Loadata for analysis ###
##############################
# Set your working directory
setwd("/Users/adriandelomana/gDrive/education/isb.summer.course.2016/src/melanoma/nonMalignant")

# Reading the data and metadata
print('reading and treating data...')
originalData = read.csv('data/formatted/immuneCells.8k.genes.data.csv', header=TRUE, row.names=1)
tumorMetadata = read.csv('data/formatted/immuneCells.8k.genes.tumorMetadata.csv', header=TRUE, row.names=1)

# Load centroids for immune, CAF, and Endothelial cells
centroids.immune = read.csv('data/formatted/immuneEtcCentroids.csv',row.names=1,header=T)

# Select genes that match to centroid genes
expression = originalData[rownames(centroids.immune),]

##################################
### Exploring single cell data ###
##################################
# Prepare to colorize based on single cell tumor origin
colorInterpolation = colorRampPalette(brewer.pal(9,'Set1'))
col1 = colorInterpolation(length(sort(unique(tumorMetadata[,1]))))
names(col1) = sort(unique(tumorMetadata[,1]))
cols1 = as.character(col1[tumorMetadata[,1]])

##################
## Compute tSNE ##
##################

# 3D tSNE
tic()
rtsne3d = Rtsne(t(as.matrix(expression)),dim=3,perplexity=50)
toc()

# Interactive 3D scatter plot
par3d(windowRect=c(50,50,700,700))
plot3d(rtsne3d$Y[,1], rtsne3d$Y[,2], rtsne3d$Y[,3], col=cols1, type='s', size=0.5, xlab='tSNE1',ylab='tSNE2',zlab='tSNE3')
legend3d('topright', legend=unique(tumorMetadata[,1]), fill=col1[unique(tumorMetadata[,1])])
#play3d(spin3d(), duration=5)
# Need to close window by hand

################################
### Classify to immune cells ###
################################

# Classify using Single Cell Predictor
scp1 = cor(centroids.immune, expression, method='spearman')
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

###########################################
### Colorize using immune cell analysis ###
###########################################
# Prepare colors for plotting
colorInterpolation = colorRampPalette(brewer.pal(10,'Paired'))
col2 = colorInterpolation(length(sort(unique(m1[,4]))))
names(col2) = sort(unique(m1[,4]))
cols2 = as.character(col2[m1[,4]])

# Interactive 3D scatter plot
par3d(windowRect=c(50,50,700,700))
plot3d(rtsne3d$Y[,1], rtsne3d$Y[,2], rtsne3d$Y[,3], col=cols2, type='s', size=0.5, xlab='tSNE1',ylab='tSNE2',zlab='tSNE3')
legend3d('topright', legend=sort(unique(m1[,4])), fill=col2[sort(unique(m1[,4]))])
#play3d(spin3d(), duration=10)
# Need to close window by hand

#################################
### Show only one immune type ###
#################################
# Redo coloring with only one 
colorInterpolation = colorRampPalette(brewer.pal(10,'Paired'))
col2 = colorInterpolation(length(sort(unique(m1[,4]))))
names(col2) = sort(unique(m1[,4]))
cols2 = as.character(col2[m1[,4]])
# Replace DENDA2 with any label to see the plot colored
# with only that immune cell type
cols2[which(m1[,4]!='BASO1')] = NA 

par3d(windowRect=c(50,50,700,700))
plot3d(rtsne3d$Y[,1], rtsne3d$Y[,2], rtsne3d$Y[,3], col=cols2, xlab='tSNE 1', ylab='tSNE 2', zlab='tSNE 3',pch=19, type='s', size=0.5) 
legend3d('topright', legend=sort(unique(m1[,4])), fill=col2)
#play3d(spin3d(),duration=10)

