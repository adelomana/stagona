# ..............................................................#
# @Course: Systems Biology of Disease                           #
# @Rscript: thu.walk.r                                          #
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

# 0. Preliminaries ====

# 0.1. Load required packages for today
library(RColorBrewer)   # 1. library to access easily multiple colors for plotting
library(Rtsne)          # 2. R implementation of the t-SNE algorithm
library(tictoc)         # 3. library to profile execution time
library(Seurat)         # 4. library to single-cell analysis
library(magrittr)       # 5. library for introducing pipe syntax: https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html
library(dplyr)          # 6. useful to manipulate data frames

# 0.1. Set working directory
# use the appropriately edited line below or RStudio graphical interphace: Lower right quadrant "Files" --> navigate to Thursday materials --> "More" --> "Set as Working Directory"
#setwd("/Users/myname/Documents/Day 4...") # edit this line to match your local folder

# 1. DIMENSIONALITY REDUCTION ====

# 1.0. load expression data and patient metadata ----
print('loading data...')
load("data/malignantCells.2k.RData")
View(malignantCellsExpression[1:5,1:5])

# 1.1. PCA ----
print('running PCA...')
tic()
results=prcomp(malignantCellsExpression) # this step takes about 1 minute in my laptop...
toc()
plot(results$x[,'PC1'],results$x[,'PC2'],main='PCA of malignant cells',pch=19,xlab='PC 1',ylab='PC 2')

# associate patient metadata to plotting colors
plottingColors=brewer.pal(length(unique(malignantCellsPatientMetadata)),'Dark2')
names(plottingColors)=unique(malignantCellsPatientMetadata)

View(malignantCellsPatientMetadata)
View(unique(malignantCellsPatientMetadata))
View(plottingColors)
View(plottingColors[malignantCellsPatientMetadata])

# plot again with patient metadata being a different color 
# put legend out
plot(results$x[,'PC1'],results$x[,'PC2'],main='PCA of malignant cells',col=plottingColors[malignantCellsPatientMetadata],pch=19,xlab='PC 1',ylab='PC 2')
legend('bottomright',legend=unique(malignantCellsPatientMetadata),fill=plottingColors)

# 1.2. tSNE ----
print('running 2D t-SNE with large perplexity...')
tic()
results2D=Rtsne(malignantCellsExpression,dims=2,perplexity=50,verbose=TRUE,theta=0) # this step takes a bit more than a minute in my laptop...
toc()
plot(results2D$Y,main='tSNE of malignant cells, p=50',col=plottingColors[malignantCellsPatientMetadata],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(malignantCellsPatientMetadata),fill=plottingColors)

# what if we use other tSNE parameters?
tic()
results2D=Rtsne(malignantCellsExpression,dims=2,perplexity=5,verbose=TRUE,theta=0) # this step takes a bit more than a minute in my laptop...
toc()
plot(results2D$Y,main='tSNE of malignant cells, p=5',col=plottingColors[malignantCellsPatientMetadata],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(malignantCellsPatientMetadata),fill=plottingColors)

# what if we would've used less genes, would we be able to separate patients? Say we get 20 genes...
reducedNumberOfGenes=20
reducedSetExpression=malignantCellsExpression[,1:reducedNumberOfGenes]
tic()
results2Dreduced=Rtsne(reducedSetExpression,dims=2,perplexity=50,verbose=TRUE,theta=0) # this step takes a bit more than a minute in my laptop...
toc()
plot(results2Dreduced$Y,main='tSNE of malignant cells, p=50',col=plottingColors[malignantCellsPatientMetadata],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(malignantCellsPatientMetadata),fill=plottingColors)

# what if we would've used less genes, would we be able to separate patients? OK, let's get 200 genes...
reducedNumberOfGenes=200
reducedSetExpression=malignantCellsExpression[,1:reducedNumberOfGenes]
tic()
results2Dreduced=Rtsne(reducedSetExpression,dims=2,perplexity=50,verbose=TRUE,theta=0) # this step takes a bit more than a minute in my laptop...
toc()
plot(results2Dreduced$Y,main='tSNE of malignant cells, p=50',col=plottingColors[malignantCellsPatientMetadata],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(malignantCellsPatientMetadata),fill=plottingColors)

# what if we would've used less genes, would we be able to separate patients? Hmmm, what about the 20 most variable genes...
reducedNumberOfGenes=200
reducedSetExpression=malignantCellsExpression[,1:reducedNumberOfGenes]
tic()
results2Dreduced=Rtsne(reducedSetExpression,dims=2,perplexity=50,verbose=TRUE,theta=0) # this step takes a bit more than a minute in my laptop...
toc()
plot(results2Dreduced$Y,main='tSNE of malignant cells, p=50',col=plottingColors[malignantCellsPatientMetadata],pch=19,xlab='tSNE Component 1',ylab='tSNE Component 2')
legend('topright',legend=unique(malignantCellsPatientMetadata),fill=plottingColors)

# 2. BIOMARKERS ====

# 2.0. load data ----
load("data/pbmc.RData") 
# Originally retrieved from: https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
# using the code:
#pbmc.data=Read10X(data.dir="~/scratch/filtered_gene_bc_matrices/hg19/")
#save(pbmc.data, file = "~/scratch/pbmc.RData")

# 2.1. initialize a Seurat object ----
# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_PBMC")

# 2.2. QC ----
# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

# 2.3. data normalization ----
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 2.4. find high variable genes ----
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)

# 2.5. regress out confounding variables ----
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

# 2.6. PCA ----
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)

PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

# 2.7. select significant PCs ----
#pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = TRUE) # very slow. From 36 to xx
#JackStrawPlot(object = pbmc, PCs = 1:12)
PCElbowPlot(object = pbmc)

# 2.8. cluster cells ----
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)

# 2.9. tSNE ----
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = pbmc)

# 2.10. find biomarkers
# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3),min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

# plot expression of biomarkers
#cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(object = pbmc, features.plot = c("LYZ", "S100A9"))

# plot distribution of biomarkers in cell across cell-types
FeaturePlot(object = pbmc, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14","FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"),reduction.use = "tsne")

# show a heatmap of biomarkers per cell type
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

# 2.11. recapitulate biology ----
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells","FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)

# 3. IMPUTATION ====

