# Usage: Rscript figure_nonneuronal.R file_markers

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly=TRUE)
file_markers <- args[1]

markers_astro <- read.table(file_markers, header=TRUE) %>%
  filter(p_val_adj<0.05 & avg_log2FC>0 & pct.1>pct.2) %>%
  filter(cluster=="1") %>%
  .$gene

markers_oligo <- read.table(file_markers, header=TRUE) %>%
  filter(p_val_adj<0.05 & avg_log2FC>0 & pct.1>pct.2) %>%
  filter(cluster=="2") %>%
  .$gene

markers_nn <- read.table(file_markers, header=TRUE) %>%
  filter(p_val_adj<0.05 & avg_log2FC>0 & pct.1>pct.2) %>%
  filter(cluster=="6") %>%
  .$gene

sink("Markers of non-neuronal.txt")
print("Markers of non-neuronal in markers of Astrocytes")
table(markers_nn %in% markers_astro)

print("Markers of non-neuronal in markers of Oligodendrocytes")
table(markers_nn %in% markers_oligo)

print("Markers of non-neuronal in either astro or oligo")
table(markers_nn %in% c(markers_oligo, markers_astro))
sink(NULL)
