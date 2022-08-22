"Make figure for the tanycyte subcluster from the astrotany cluster

Usage: figure_tanycyte.R --file_sc_obj=<file> --file_functions=<file>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		The rds file after clustering
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

message("file_sc_obj: ", file_sc_obj)
message("file_functions: ", file_functions)

# load all functions defined in the r-script
source(file_functions)

# read file
sc_obj <- readRDS(file_sc_obj)

# Define the name of cluster of interest to make pdf filenames accordingly
coi <- "Tanycyte"

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

# Make featureplots for DE genes between WT and Mutant Tanycytes.
Idents(sc_obj) <- sc_obj$genotype
# Find marker genes
markers <- FindAllMarkers(sc_obj, logfc.threshold = 0.5, min.pct = 0.2, test.use = "wilcox", min.cells.group=0.5, only.pos=TRUE)

# save the markers
filename="tanycyte_markers_wtvstra1.tsv"
write.table(markers, file=filename, sep="\t")

# Make FeaturePlot
foi <- markers %>%
  filter(avg_log2FC>=0.5) %>%
  group_by(cluster) %>%
  slice_head(n=5) %>%
  .$gene

#Feature Plot in pdf
plot <- FeaturePlot(sc_obj, features=foi, split.by="genotype", keep.scale="all") &
  theme(legend.position="none",panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),   axis.ticks = element_blank(), panel.background = element_blank(), line=element_blank(), plot.title=element_text(family="sans", face="plain", hjust=0.5, size=15), aspect.ratio=1)
filename=paste0("FeaturePlot_cluster",coi,"_int_degenes.pdf")
ggsave(plot=plot, filename=filename, width=4*2, height=4*length(foi), limitsize=FALSE)

#Feature Plot in png
filename=paste0("FeaturePlot_cluster",coi,"_int_degenes.png")
ggsave(plot=plot, filename=filename, width=4*2, height=4*length(foi), limitsize=FALSE)  # size divided by 4 because it will be 1/4th of the size on the figure

# *** Make Violin Plots for the genes
foi <- rownames(markers)
plot <- VlnPlot(sc_obj, features=foi, group.by="genotype", ncol=length(foi), cols=c("#ca2027","#0272b0")) &
  theme(text=element_text(family="sans", size=12), plot.title=element_text(face="plain"), axis.title=element_blank(), axis.text=element_text(family="sans", size=12))
filename=paste0("cluster",coi,"_int_degenes.png")
ggsave(plot=plot, filename=filename, width=7*length(foi)/4, height=7/2) # size divided by 4 because it will be 1/4th of the size on the figure
