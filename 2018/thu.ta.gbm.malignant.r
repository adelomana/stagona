# ..............................................................#
# @Course: Systems Biology of Disease                           #
# @Rscript: thu.ta.gbm.malignant.r                              #
# @Version: 3.2                                                 #
# @Authors: Adrian Lopez Garcia de Lomana and Chris Plaisier    #
#                                                               #
# @Sponsored by:                                                #
# Institute for Systems Biology                                 #
# 401 Terry Avenue North                                        #
# Seattle, WA 98109-5263                                        #
#                                                               #
# This source code is distributed under the                     #
# GNU General Public License v3.0,                              #
# the text of which is available at:                            #
# https://www.gnu.org/licenses/gpl-3.0.en.html                  #
# ..............................................................#

rm(list=ls())

# libraries
library(RColorBrewer)
library(Rtsne)
library(Seurat)
library(magrittr) 
library(dplyr)

# set working directory and load data
setwd("~/scratch/")
d1 = read.csv('GSE57872_GBM_data_matrix.ALO.txt', header=1, row.names=1,sep="\t")
bulk = read.csv('GSE57872_GBM_data_matrix.ALO.bulk.txt', header=1, row.names=1,sep="\t")

# transposing expression data
library(data.table)
gbm.expression <- transpose(d1)
colnames(gbm.expression) <- rownames(d1)
rownames(gbm.expression) <- colnames(d1)

gbm.bulk.expression <- transpose(bulk)
colnames(gbm.bulk.expression) <- rownames(bulk)
rownames(gbm.bulk.expression) <- colnames(bulk)

# brew colors for students
patientMetadata = sapply(rownames(gbm.expression), function(x) { strsplit(x,'\\.')[[1]][1] } )
plottingColors=brewer.pal(length(unique(patientMetadata)),'Dark2')
names(plottingColors)=unique(patientMetadata)

save(gbm.expression, patientMetadata, file = "/Users/adriandelomana/scratch/gbm.RData")
load("data/gbm.RData")

# PCA
pca1 = prcomp(gbm.expression)
plot(pca1$x[,'PC1'], pca1$x[,'PC2'], main="GBM single cell PCA",pch=19,col=plottingColors[patientMetadata],xlab='PC 1',ylab='PC 2')
legend('topright',legend=unique(patientMetadata),fill=plottingColors)

# tSNE
rtsne2d = Rtsne(gbm.expression,perplexity=50,verbose=TRUE,theta=0)
plot(rtsne2d$Y, main="GBM single-cell tSNE",pch=19,col=plottingColors[patientMetadata],xlab='tSNE 1',ylab='tSNE 2')
legend('topright',legend=unique(patientMetadata),fill=plottingColors)

# BIOMARKERS

# create object
# careful! Seurat requires count, TPM or FPKM. Don't feed log2, or the negative values will be removed!
fpkm=2**d1
gbm <- CreateSeuratObject(raw.data = fpkm, min.cells = 3, min.genes = 200,is.expr=0.01,do.scale=FALSE,do.center = FALSE)

# data normalization
gbm <- NormalizeData(object = gbm, normalization.method = "LogNormalize")

# find high variable genes
gbm <- FindVariableGenes(object = gbm, mean.function = ExpMean, dispersion.function = LogVMR)
length(x = gbm@var.genes)

# scale data
gbm <- ScaleData(object = gbm)

# pca. consider using 50 many PCs as possible, this data set is challenging!
gbm <- RunPCA(object = gbm, pc.genes = gbm@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute=20)

PrintPCA(object = gbm, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = gbm, pcs.use = 1:2)
PCAPlot(object = gbm, dim.1 = 1, dim.2 = 2)

PCHeatmap(object = gbm, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCElbowPlot(object = gbm)

# clustering cells
gbm <- FindClusters(object = gbm, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
gbm <- RunTSNE(object = gbm, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = gbm)

# find markers for every cluster compared to all remaining cells, report
gbm.markers <- FindAllMarkers(object = gbm, only.pos = TRUE)
View(gbm.markers %>% group_by(cluster) %>% top_n(5,avg_logFC))
gbm.markers %>% group_by(cluster) %>% top_n(1,avg_logFC)

cluster1.markers <- FindMarkers(object = gbm, ident.1 = 1, test.use = "roc", only.pos = TRUE)
VlnPlot(object = gbm, features.plot = c("NKG7", "CCL4", "CD8A", "CD6"))


