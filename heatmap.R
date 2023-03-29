library(dplyr)
library(textshape)
library(tidyverse)
library(ggplot2)
library(data.table)
library(plyr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

#give directions to relevant folders and convert data to matrix
datafile="C:/directions/to/file.csv"

dataframe <- read.csv(datafile, header=TRUE, sep=",", stringsAsFactors = TRUE)
head(dataframe)
m <- as.matrix(dataframe[, -1])
rownames(m) <- dataframe$Cytokine

#make heatmap
htmp <- Heatmap(
  m, 
  width = ncol(m)*unit(0.5,"cm"),
  height = ncol(m)*unit(7,"cm"),
  name = "log2FC",    
  col = colorRamp2(c(-2, 0, 2), c("#2166ac", "#FFFFFF", "#b2182b")),
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 10),
  
  show_column_names = TRUE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,

  heatmap_legend_param = list(
    legend_direction = "horizontal", 
    legend_width = unit(5, "cm")), 
)
draw(htmp, padding = unit(c(2, 2, 2, 20), "mm"), heatmap_legend_side="bottom",
     legend_grouping = "original")
