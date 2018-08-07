# ..............................................................#
# @Course: Systems Biology of Disease                           #
# @Rscript: thu.apply.mela.non.malignant.r                      #
# @Version: 3.2                                                 #
# @Authors: Adrian Lopez Garcia de Lomana                       #
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

# 0. Syntax hints ----
# find here some really great syntax hints about general R commands:
# http://github.com/rstudio/cheatsheets/raw/master/base-r.pdf

# and very advanced coding tips:
# https://www.rstudio.com/resources/cheatsheets/

# 1.0. load data
load("data/nonMalignantCells.2k.RData")
# consider using labels to brew colors for plotting --> get the syntax from the walkthrough script

# 1.1. MELANOMA NON-MALIGNANT DIMENSIONALITY REDUCTION
# Is there heterogeneity/cell states in non-malignant cells?
# consider brewing colors using this syntax
# plottingColors=brewer.pal(length(unique(labels to use)),'Set3')
# names(plottingColors)=unique(labels to use)

# 1.2. MELANOMA NON-MALIGNANT BIOMARKERS ----

# Which genes are differentially expressed between cell clusters?
# Do cell clusters match immune cell types? Which ones?
  
# Seurat uses a matrix of genes as rows, we need to transpose expression and create a Seurat object
library(data.table)
transpose.nonMalignantCellsExpression <- transpose(nonMalignantCellsExpression)
colnames(transpose.nonMalignantCellsExpression) <- rownames(nonMalignantCellsExpression)
rownames(transpose.nonMalignantCellsExpression) <- colnames(nonMalignantCellsExpression)
# create Seurat object
scso <- CreateSeuratObject(raw.data = transpose.nonMalignantCellsExpression, min.cells = 3, min.genes = 200,is.expr=0.01,do.scale=FALSE,do.center = FALSE)
# data is expression, so now you are ready for normalization step and beyond!