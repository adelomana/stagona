# ..............................................................#
# @Course: Systems Biology of Disease                           #
# @Rscript: thu.walk.r                                          #
# @Version: 3.1                                                 #
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

# Check package installation ----
required_bioconductor_packages <- c("GEOquery")
required_cran_packages <- c("RColorBrewer", "Rtsne", "tictoc", "magrittr", "dplyr", "Seurat","data.table")
required_packages <- c(required_bioconductor_packages, required_cran_packages)

# Load all required packages
# sapply -- apply a function to all the elements of a vector like object
# this checks that all required packages are installed i.e. require("ConsensusClusterPlus")
sapply(X = required_packages, FUN = require, character.only = TRUE, quietly=TRUE)

# Install required R packages ----
# Bioconductor installation
source("https://bioconductor.org/biocLite.R") 
biocLite(required_bioconductor_packages) 

# CRAN installation
install.packages(required_cran_packages)