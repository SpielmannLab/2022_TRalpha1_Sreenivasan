# Makes plots for SI figure based on harmony-integrated OPC cluster
# Includes: FeaturePlots
# Note: separate rds file for the individual clusters were generated using separate_clusters.R script, and then they were individually integrated using integrate_clusters.R script (Harmony was preferred).

"Make - figure for OPC Characterisation

Usage: figure_opc.R --file_sc_obj_opc=<file>

Options:
  -h --help			Show this screen.
  --file_sc_obj_opc=<file>		The rds file of oligodendrocyte progenitor cells after integration
" -> doc

# Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(docopt))
suppressMessages(library(Matrix))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters
file_sc_obj_opc <- arguments$file_sc_obj_opc

message("file_sc_obj_opc: ", file_sc_obj_opc)

# read file
sc_obj <- readRDS(file_sc_obj_opc)

# Record which subcluster is being currently processed
coi <- sc_obj$main_cluster %>% unique

# Create Feature Plots for OPC for SI
# Features from LaManno extended data figure 7
features_lamanno <- c("Pdgfra","Lhfpl3","Bmp4","Neu4","Kcnj12","Tmem2","Rras2","Sema4f","Tmem141","Arap2","Erbb3","Opalin","Plekhh1","Mog","Gjb1","Mal")

p <- FeaturePlot(sc_obj, features=features_lamanno, cols=c("#bababa", "#ca0020"), keep.scale="all", pt.size=0.25) &
  theme_void() &
  theme(plot.title=element_text(hjust=0.5, face="italic")) &
  theme(legend.position="none")

filename <- paste0("FeaturePlot_cluster",coi,"_int_lamanno_markers.pdf")
ggsave(p + theme(legend.position="right"), file=filename,width=7,height=7)
filename <- paste0("FeaturePlot_cluster",coi,"_int_lamanno_markers.png")
ggsave(p, file=filename,width=7,height=7)
