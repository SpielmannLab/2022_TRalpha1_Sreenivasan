"Make figure for the Astrocytes_Tanycytes cluster

Usage: figure_astrotany.R --file_sc_obj=<file> --file_functions=<file> --file_markers=<file> --res=<value>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		The rds file after clustering
  --file_markers=<file>		The tsv file containing marker genes for the cluster
  --file_functions=<file> The .R file containing all functions
  --res=<value>           The resolution of clustering
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
res <- arguments$res

message("file_sc_obj: ", file_sc_obj)
message("file_functions: ", file_functions)
message("file_markers: ", file_markers)
message("res: ", res)

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
markergenes <- c("Agt", "Kcnd2", "Ntm", "Itih3", "Col23a1", "Gpr50", "Rax", "Frzb", "Efna5",  "Slit2", "Ch4dl1", "Igfbp2", "Islr", "Fa2h", "Prr5l", "St18")
# Rename idents

cluster_key <- sc_obj@meta.data %>%
  names %>%
  grep(pattern=paste0("snn_res.",res), value=TRUE)
plot <- DotPlot(sc_obj, features=markergenes, group.by=cluster_key) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
filename=paste0("dotplot_nolabel.pdf")
ggsave(plot, filename=filename, width=8, height=4)

# Make FeaturePlots
plot <- FeaturePlot(sc_obj, features=markergenes, split.by="genotype")
filename=paste0("featureplot.pdf")
ggsave(plot, filename=filename, width=8, height=32)

# Name the clusters
Idents(sc_obj) <- sc_obj@meta.data[,cluster_key]
sc_obj <- RenameIdents(sc_obj, "0"= "Astrocytes",
                              "1" = "Tanycytes",
                              "2" = "Unknown_1")

plot <- DotPlot(sc_obj, features=markergenes) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Cluster") +
  xlab("Gene")
filename=paste0("dotplot_label.pdf")
ggsave(plot, filename=filename, width=10, height=4)

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
