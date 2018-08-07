# ..............................................................#
# @Course: Systems Biology of Disease                           #
# @Rscript: thu.ta.gbm.malignant.r                              #
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

# libraries
library(RColorBrewer)
library(Rtsne)
library(Seurat)
library(magrittr) 
library(dplyr)

# 0. Syntax hints ----
# find here some really great syntax hints about general R commands:
# http://github.com/rstudio/cheatsheets/raw/master/base-r.pdf

# and very advanced coding tips:
# https://www.rstudio.com/resources/cheatsheets/

# set working directory and load data. Make sure you inspect the new variables loaded!
load("data/gbm.RData")

# Is there heterogeneity/cell states in GBM cells?
# consider brewing colors using this syntax
# plottingColors=brewer.pal(length(unique(labels to use)),'Set3')
# names(plottingColors)=unique(labels to use)

# PCA

# tSNE

# BIOMARKERS

# create Seurat object. careful! Seurat requires count, TPM or FPKM. Don't feed log2, or the negative values will be removed!
# Seurat uses a matrix of genes as rows, do we need to transpose expression or not to create a Seurat object?

# data normalization

# find high variable genes

# scale data

# pca. consider using 20 PCs as default

# cluster cells

# find markers for every cluster compared to all remaining cells, report
