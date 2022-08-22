"Make figure for the Neuron_1 cluster

Usage: figure_qc.R --file_sc_obj=<file> --file_functions=<file> --file_markers=<file>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		The rds file after clustering
  --file_markers=<file>		The tsv file containing marker genes for the cluster
  --file_functions=<file> The .R file containing all functions
" -> doc

# Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(docopt))
suppressMessages(library(RColorBrewer)) #to expand color palette

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters: Read in the output of the sc_multi_sample pipline
file_sc_obj <- arguments$file_sc_obj
file_functions <- arguments$file_functions
file_markers <- arguments$file_markers

message("file_sc_obj: ", file_sc_obj)
message("file_functions: ", file_functions)
message("file_markers: ", file_markers)

# load all functions defined in the r-script
source(file_functions)

# read file
sc_obj <- readRDS(file_sc_obj)

# First make heatmap for the DE genes identified by Seurat FindAllMarkers
markers <- read.table(file_markers, header=TRUE) %>%
  filter(avg_log2FC>0 & pct.1>pct.2 & p_val_adj<0.05) %>%
  group_by(cluster) %>%
  slice_head(n=10)
plot <- DoHeatmap(sc_obj, features=markers$gene)
filename=paste0("heatmap_markergenes_nolabel.pdf")
ggsave(plot, filename=filename, width=10, height=10)

# Custom settings for custom plots
markergenes <- c("Gad2", "Slc17a6", "Pmch", "Cartpt", "Pomc", "Avp",  "Prlr", "Pgr", "Npy", "Agrp", "Oxt", "Sst", "Hcrtr2")
# Rename idents
res <- 0.1 #This was the resolution at which the neurons were clustered
cluster_key <- sc_obj@meta.data %>%
  names %>%
  grep(pattern=paste0("snn_res.",res), value=TRUE)
plot <- DotPlot(sc_obj, features=markergenes, group.by=cluster_key) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
filename=paste0("dotplot_nolabel.pdf")
ggsave(plot, filename=filename, width=8, height=4)

# Name the clusters
Idents(sc_obj) <- sc_obj@meta.data[,cluster_key]
sc_obj <- RenameIdents(sc_obj, "0"= "GABA Neurons",
                              "1" = "Glutamatergic Neurons_1",
                              "2" = "Glutamatergic Neurons_2",
                              "3" = "Glutamatergic Neurons_3",
                              "4" = "MCH Neurons",
                              "5" = "CART Neurons_1",
                              "6" = "AVasopressin Neurons",
                              "7" = "POMC/CART Neurons",
                              "8" = "NPY/AgRP Neurons",
                              "9" = "Oxytocin+, AVP Neurons",
                              "10" = "Somatostatin+ Neurons",
                              "11" = "#11",
                              "12" = "Orexin+ Neurons",
                              "13" = "CART Neurons_2")
# reorder the plotting
Idents(sc_obj) <- factor(Idents(sc_obj), levels=levels(sc_obj)[c(1:6,14,7:13)])

plot <- DotPlot(sc_obj, features=markergenes) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Cluster") +
  xlab("Gene")
filename=paste0("dotplot_label.pdf")
ggsave(plot, filename=filename, width=10, height=4)

# Dotplot for campbell 2017 marker genes
features_campbell <-"Hdc\nSlc18a2\nMaob\nPenk\nGm8773\nChodl\nTh\nOnecut2\nNr4a2\nCalb2\nPnoc\nAp3b1\nHtr2c\nVcan\nSatb1\nOxt\nAngpt1\nFkbp11\nNxph1\nNr3c1\nBcl2\nCoch\nTrpc6\nOtof\nFam159b\nArhgap36\nIfi27\nGhrh\nGal\nSlc18a3\nTrh\nCox6a2\nH2-K1\nNr0b1\nAmigo2\nAgrp\nNpy\nSerpina3n\nAcvr1c\nFam159b\nTtr\nH2-Q2\nPomc\nCartpt\nNpy1r\nVip\nRgs16\nRora\nLmo1\nEnpp2\nKrt77\nCck\nNmu\nPtk2b\nNhlh2\nNr5a2\nRxfp1\nGfra1\nGlipr1\nPgr\nSytl4\nNpy5r\nLamp5\nCbln4\nSstr2\nSst\nPthlh\nCol25a1\nIcam5\nUnc13c\nSox14\nRxrg\nCrabp1\nVgll3\nHtr3b\nTbx19\nCd1d1\nUst\nQrfp\nC1ql3\nRprm\nSlc17a6\nAdcyap1\nFezf1\nSlc17a8\nNpnt\nPtgds\nSynpr\nAI593442\nPcsk2\nGpc3\nTmem163\nGabrq\n"
features_campbell <- strsplit(features_campbell, split="\n") %>% unlist() %>% unique()
plot <- DotPlot(sc_obj, features=features_campbell, cluster.idents=TRUE) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Cluster") +
  ylab("Gene")
filename=paste0("dotplot_label_clustered_campbelfeatures_neurons.pdf")
ggsave(plot, filename=filename, width=20, height=8)
plot <- DotPlot(sc_obj, features=features_campbell) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Cluster") +
  ylab("Gene")
filename=paste0("dotplot_label_campbelfeatures_neurons.pdf")
ggsave(plot, filename=filename, width=20, height=8)


# .......... UMAP
colorpalette <- colorRampPalette(brewer.pal(7, "Paired"))(length(levels(sc_obj)))
plot <- DimPlot(sc_obj, group.by="ident", reduction="umap", label=TRUE, repel=TRUE, pt.size=0.25) +
  scale_color_manual(values=rev(colorpalette)) +
  theme(legend.position="none")
filename=paste0("umap_clusters_res",res,".pdf")
ggsave(plot=plot, filename=filename, width=7, height=7)
plot <- DimPlot(sc_obj, group.by="ident", reduction="umap", pt.size=0.25) +
  scale_color_manual(values=rev(colorpalette)) +
  theme(legend.position="right")
filename=paste0("umap_clusters_res",res,".png")
ggsave(plot=plot, filename=filename, width=7, height=7)
plot <- DimPlot(sc_obj, group.by="ident", reduction="umap", pt.size=0.25) +
  scale_color_manual(values=rev(colorpalette)) +
  theme(legend.position="none")
filename=paste0("umap_clusters_res_nolegend",res,".png")
ggsave(plot=plot, filename=filename, width=7, height=7)

# ......... Calculate cellular composition and save plot
cluster_key <- "annotations"
sc_obj$annotations <- Idents(sc_obj)
metadata <- select(sc_obj@meta.data,c("genotype",all_of(cluster_key)))
# Use the function in cell_composition_calculation.R
output <- cell_composition(metadata=metadata, grouping_var="genotype", xvar=cluster_key, label_y_coord=30)

filename=paste0("cell_composition.pdf")
ggsave(plot=plot_grid(output$plot2 + theme_bw(), output$plot1 + theme_bw(), nrow=1), filename=filename, width=10, height=7)

# .......... UMAP by genotype
plot <- DimPlot(sc_obj, cols=c("#ca2027","#0272b0"), group.by="genotype", reduction="umap", label=TRUE, repel=TRUE, pt.size=0.25, shuffle=TRUE) +
  theme(legend.position="none")
filename=paste0("umap_by_genotype.pdf")
ggsave(plot=plot, filename=filename, width=7, height=7)
plot <- DimPlot(sc_obj, cols=c("#ca2027","#0272b0"), group.by="genotype", reduction="umap", pt.size=0.25, shuffle=TRUE) +
  theme(legend.position="none")
filename=paste0("umap_by_genotype.png")
ggsave(plot=plot, filename=filename, width=7, height=7)

# ***** UMAP by repeat
sc_obj$rep <- gsub(sc_obj$orig.ident, pattern="-\\w+", replacement="")

col_reps <- c("#4daf4a","#984ea3") # colorbrewer2.0, qualitative, 5th set.

plot <- DimPlot(sc_obj, cols=col_reps, group.by="rep", reduction="umap", label=TRUE, repel=TRUE, pt.size=0.25, shuffle=TRUE) +
  theme(legend.position="none")
filename=paste0("umap_by_repeat.pdf")
ggsave(plot=plot, filename=filename, width=7, height=7)
plot <- DimPlot(sc_obj, cols=col_reps, group.by="rep", reduction="umap", pt.size=0.25, shuffle=TRUE) +
  theme(legend.position="none")
filename=paste0("umap_by_repeat.png")
ggsave(plot=plot, filename=filename, width=7, height=7)
