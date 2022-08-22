"Make - umaps for individual clusters

Usage: umap_separate_clusters_nf.R --file_sc_obj=<file>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		The rds file after clustering
" -> doc

# Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(docopt))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters: Read in the output of the sc_multi_sample pipline
file_sc_obj <- arguments$file_sc_obj

message("file_sc_obj: ", file_sc_obj)

# read file
sc_obj <- readRDS(file_sc_obj)

# Record which subcluster is being currently processed
coi <- sc_obj$main_cluster %>% unique

# ****** Make the plots
# .......... UMAP by genotype for each cluster
sc_obj$genotype <- gsub(sc_obj$orig.ident, pattern="\\w+-", replacement="")

plot <- DimPlot(sc_obj, cols=c("#ca2027","#0272b0"), group.by="genotype", reduction="umap", label=TRUE, repel=TRUE, pt.size=0.25, shuffle=TRUE) +
  theme(legend.position="none")
filename=paste0("cluster",coi,"_umap_by_genotype.pdf")

message("filename: ", filename)

ggsave(plot=plot, filename=filename, width=7, height=7)
plot <- DimPlot(sc_obj, cols=c("#ca2027","#0272b0"), group.by="genotype", reduction="umap", pt.size=0.25, shuffle=TRUE) + theme(legend.position="none")
filename=paste0("cluster",coi,"_umap_by_genotype.png")
ggsave(plot=plot, filename=filename, width=7, height=7)

# ***** UMAP by repeat for each cluster
sc_obj$rep <- gsub(sc_obj$orig.ident, pattern="-\\w+", replacement="")

col_reps <- c("#4daf4a","#984ea3") # colorbrewer2.0, qualitative, 5th set.

plot <- DimPlot(sc_obj, cols=col_reps, group.by="rep", reduction="umap", label=TRUE, repel=TRUE, pt.size=0.25, shuffle=TRUE) +
  theme(legend.position="none")
filename=paste0("cluster",coi,"_umap_by_repeat.pdf")
ggsave(plot=plot, filename=filename, width=7, height=7)
plot <- DimPlot(sc_obj, cols=col_reps, group.by="rep", reduction="umap", pt.size=0.25, shuffle=TRUE) +
  theme(legend.position="none")
filename=paste0("cluster",coi,"_umap_by_repeat.png")
ggsave(plot=plot, filename=filename, width=7, height=7)
