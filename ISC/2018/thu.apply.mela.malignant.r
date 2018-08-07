# ..............................................................#
# @Course: Systems Biology of Disease                           #
# @Rscript: thu.apply.mela.malignant.r                          #
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

# 0. Syntax hints ----
# find here some really great syntax hints about general R commands:
# http://github.com/rstudio/cheatsheets/raw/master/base-r.pdf

# and very advanced coding tips:
# https://www.rstudio.com/resources/cheatsheets/

# 1.0. load data and metadata ----
load("data/malignantCells.2k.RData")

# 1.1. MELANOMA MALIGNANT DIMENSIONALITY REDUCTION ----
# which is an optimal value for perplexity? Why? What does perplexity control?

# consider brewing colors using this syntax
# plottingColors=brewer.pal(length(unique(labels to use)),'Set3')
# names(plottingColors)=unique(labels to use)

# 1.2. MELANOMA MALIGNANT BIOMARKERS ----
# Seurat uses a matrix of genes as rows, we need to transpose expression and create a Seurat object
library(data.table)
transpose.melanoma.expression <- transpose(malignantCellsExpression)
colnames(transpose.melanoma.expression) <- rownames(malignantCellsExpression)
rownames(transpose.melanoma.expression) <- colnames(malignantCellsExpression)
# create Seurat object
mela <- CreateSeuratObject(raw.data = transpose.melanoma.expression, min.cells = 3, min.genes = 200,is.expr=0.01,do.scale=FALSE,do.center = FALSE)
# data is expression, so now you are ready for normalization step and beyond!