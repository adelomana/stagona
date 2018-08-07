# ..............................................................#
# @Course: Systems Biology of Disease                           #
# @Rscript: thu.ta.mela.non.malignant.r                         #
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

library(Rtsne)
library(tictoc)
library(Seurat)
library(RColorBrewer)
library(magrittr) 
library(dplyr)

# DIMENSIONALITY REDUCTION AND CELL STATES IN NON-MALIGNANT CELLS
#setwd("~/Google Drive/education/ISC/2018/sandbox")

#dataFilePath=paste('~/scratch/','nonMalignant.2kgenes.data.prediction.csv',sep='')
#originalData=read.csv(dataFilePath,header=TRUE,row.names=1)
#nonMalignantCellsExpression=as.data.frame(t(originalData)) # transposing the original data into the appropriate form: 2,250 observations in a 8,000 dimensional space

#metadataFilePath=paste('~/scratch/','nonMalignant.2kgenes.tumorMetadata.prediction.csv',sep='')
#tumorMetadata=read.csv(metadataFilePath,header=TRUE,row.names=1)
#tumorLabels=as.character(tumorMetadata$tumor.label)
#plottingColors=brewer.pal(length(unique(tumorLabels)),'Set3')
#names(plottingColors)=unique(tumorLabels)

#immuneFilePath=paste("~/scratch/","nonMalignant.2kgenes.immuneMetadata.prediction.csv",sep="")
#immuneMetadata=read.csv(immuneFilePath,header=TRUE,row.names=1)
#immuneLabels=as.character(immuneMetadata$immune.label)
#plottingColors=brewer.pal(length(unique(immuneLabels)),'Set3')
#names(plottingColors)=unique(immuneLabels)

#save(nonMalignantCellsExpression, tumorLabels, immuneLabels, file = "/Users/adriandelomana/scratch/nonMalignantCells.2k.RData")
load("data/nonMalignantCells.2k.RData")

tic()
results2D=Rtsne(nonMalignantCellsExpression,dims=2,perplexity=50,verbose=TRUE,theta=0) # this step takes a bit more than a minute in my laptop...
toc()

plot(results2D$Y,main='tSNE of non malignant cells, p=50',pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')

plottingColors=brewer.pal(length(unique(tumorLabels)),'Set3')
names(plottingColors)=unique(tumorLabels)
plot(results2D$Y,main='tSNE of non malignant cells, p=50',pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2',col=plottingColors[tumorLabels])
legend('topright',legend=unique(tumorLabels),fill=plottingColors[unique(tumorLabels)])

plottingColors=brewer.pal(length(unique(immuneLabels)),'Set3')
names(plottingColors)=unique(immuneLabels)
plot(results2D$Y,main='tSNE of non malignant cells, p=50',pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2',col=plottingColors[immuneLabels])
legend('topright',legend=unique(immuneLabels),fill=plottingColors[unique(immuneLabels)])

# CHARACTHERIZATION OF NON-MALIGNANT CELL STATES

# prepare data for Seurat
library(data.table)
transpose.non.malignant.expression <- transpose(nonMalignantCellsExpression)
colnames(transpose.non.malignant.expression) <- rownames(nonMalignantCellsExpression)
rownames(transpose.non.malignant.expression) <- colnames(nonMalignantCellsExpression)

# create Seurat object
scso <- CreateSeuratObject(raw.data = transpose.non.malignant.expression, min.cells = 3, min.genes = 200,is.expr=0.01,do.scale=FALSE,do.center = FALSE)

# data normalization
scso <- NormalizeData(object = scso, normalization.method = "LogNormalize")

# find high variable genes
scso <- FindVariableGenes(object = scso, mean.function = ExpMean, dispersion.function = LogVMR)
length(x = scso@var.genes)

# scale data
scso <- ScaleData(object = scso)

# pca. consider using 50 many PCs as possible, this data set is challenging!
scso <- RunPCA(object = scso, pc.genes = scso@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute=50)

PrintPCA(object = scso, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = scso, pcs.use = 1:2)
PCAPlot(object = scso, dim.1 = 1, dim.2 = 2)

PCHeatmap(object = scso, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCElbowPlot(object = scso)

# clustering cells
scso <- FindClusters(object = scso, reduction.type = "pca", dims.use = 1:50, resolution = 0.6, print.output = 0, save.SNN = TRUE)
scso <- RunTSNE(object = scso, dims.use = 1:50, do.fast = TRUE)
TSNEPlot(object = scso)

# find markers for every cluster compared to all remaining cells, report
scso.markers <- FindAllMarkers(object = scso, only.pos = TRUE)
View(scso.markers %>% group_by(cluster) %>% top_n(5,avg_logFC))
scso.markers %>% group_by(cluster) %>% top_n(1,avg_logFC)

cluster0.markers <- FindMarkers(object = scso, ident.1 = 0, test.use = "roc", only.pos = TRUE)
VlnPlot(object = scso, features.plot = c("NKG7", "CCL4", "CD8A", "CD6"))

# nice tsne plot for biomarkers
FeaturePlot(object = scso, features.plot = c("NKG7", "IL7R", "IRF8", "CD68", "IFITM3"), cols.use = c("grey", "blue"), reduction.use = "tsne")

# nice heatmap of biomarkers 
top10 <- scso.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = scso, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

# DIFFERENTIAL CELLS ====
library(stringr)
uniqueTumorLabels=unique(tumorLabels)

patientCellTypes <- function(x){
  tableOfCells=table(immuneLabels[str_detect(tumorLabels,x)])
  plot(tableOfCells)
  if (tableOfCells["Bcells"][[1]] > 10) {
    ratio=tableOfCells["Tcells"][[1]]/tableOfCells["Bcells"][[1]]
  }else{
    ratio=NaN
  }
  print(c(tableOfCells,ratio,x))
  info=c(x,ratio)
  return(info)
}
response=apply(matrix(uniqueTumorLabels),1,patientCellTypes)
x=response[1,]
y=as.numeric(response[2,])
df=data.frame(x,y)
plot(df,xlab='Patient',ylab='T/B-cell ratio')
