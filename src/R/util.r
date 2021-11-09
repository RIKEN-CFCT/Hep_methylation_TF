#==============================================================================================#

# load functions

#==============================================================================================#
source("./src/functions.r")

#==============================================================================================#

# load libraries

#==============================================================================================#
devtools::install_github("takahirosuzuki1980/InfiniumDiffMetMotR")
library("InfiniumDiffMetMotR")

library("ggplot2")
library("ggfortify")
library("cowplot")
library("reshape2")
library("circlize")
library("RColorBrewer")
library("tidyverse")
library("dplyr")
library("tidyr")
library("ComplexHeatmap")
library("gridExtra")
library("rGREAT")
library("rrvgo")
library("UpSetR")
library("genefilter")
library("rtracklayer")