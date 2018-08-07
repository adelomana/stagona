# ..............................................................#
# @Course: Systems Biology of Disease                           #
# @Rscript: thu.ta.mela.malignant.r                             #
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


# 1. MELANOMA MALIGNANT BIOMARKERS ====

library(Seurat)         
library(magrittr)
library(dplyr)

# 1.0. load expression data and patient metadata ----
print('loading data...')
setwd("~/Google Drive/education/ISC/2018/sandbox")
load("data/malignantCells.2k.RData")

# prepare data for Seurat
library(data.table)
transpose.melanoma.expression <- transpose(malignantCellsExpression)
colnames(transpose.melanoma.expression) <- rownames(malignantCellsExpression)
rownames(transpose.melanoma.expression) <- colnames(malignantCellsExpression)

# create object
mela <- CreateSeuratObject(raw.data = transpose.melanoma.expression, min.cells = 3, min.genes = 200,is.expr=0.01,do.scale=FALSE,do.center = FALSE)

# 3. DATA NORMALIZATION
mela <- NormalizeData(object = mela, normalization.method = "LogNormalize")

# 4. HIGH VARIABLE GENES
mela <- FindVariableGenes(object = mela, mean.function = ExpMean, dispersion.function = LogVMR)
length(x = mela@var.genes)

# 5. scale data
mela <- ScaleData(object = mela)

# 6. PCA
mela <- RunPCA(object = mela, pc.genes = mela@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

PrintPCA(object = mela, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = mela, pcs.use = 1:2)
PCAPlot(object = mela, dim.1 = 1, dim.2 = 2)

PCHeatmap(object = mela, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCElbowPlot(object = mela)

# 7. clustering cells
mela <- FindClusters(object = mela, reduction.type = "pca", dims.use = 1:20, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

# 8. tSNE
mela <- RunTSNE(object = mela, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = mela)

# 9. finding biomarkers
#  find all markers of cluster 1
cluster1.markers <- FindMarkers(object = mela, ident.1 = 1)
print(x = head(x = cluster1.markers, n = 5))

# # find all markers distinguishing cluster 0 from clusters 1, 2 and 3
cluster0.markers <- FindMarkers(object = mela, ident.1 = 0, ident.2 = c(1,2,3))
print(x = head(x = cluster0.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
mela.markers <- FindAllMarkers(object = mela, only.pos = TRUE)

mela.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

cluster0.markers <- FindMarkers(object = mela, ident.1 = 0, test.use = "roc", only.pos = TRUE)
VlnPlot(object = mela, features.plot = c("APOC2", "APOC1", "FN1", "CREB3L2"))

# nice tsne plot for biomarkers
FeaturePlot(object = mela, features.plot = c("APOC2", "FN1", "LEPRE1", "DCT", "NENF", "RGS1"), cols.use = c("grey", "blue"), reduction.use = "tsne")

top10 <- mela.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = mela, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
