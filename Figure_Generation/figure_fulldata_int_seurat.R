"Make QC and Characterisation Figure for the paper

Usage: figure_fulldata_int.R --file_sc_obj=<file> --file_markers=<file> --file_functions=<file> --res=<value>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		The rds file after integration and clustering
  --file_markers=<file>		The tsv file containing marker genes for the cluster
  --file_functions=<file> The .R file containing all functions
  --res=<value>           The clustering resolution to be used.
" -> doc

# Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(docopt))
suppressMessages(library(ggridges))
suppressMessages(library(RColorBrewer))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters: Read in the output of the sc_multi_sample pipline
file_sc_obj <- arguments$file_sc_obj
file_markers <- arguments$file_markers
file_functions <- arguments$file_functions
res <- as.numeric(arguments$res)

message("file_sc_obj: ", file_sc_obj)
message("file_markers: ", file_markers)
message("file_functions: ", file_functions)
message("res: ", res)

# load all functions defined in the r-script
source(file_functions)

# read file
sc_obj <- readRDS(file_sc_obj)

# Get markers and rename idents
goi <- "Hcn1 Slc17a6 Gad1 Gad2 Aqp4 Agt Myo6 Mal Mag Mobp Pdgfra Cspg4 Hexb C1qa Cx3cr1 Snap25 Syt1 Igfbp2 Tbx18 Kcnj8" %>%
# "Agrp\nCart\nDdc\nGad1\nGad1\nGad2\nGad2\nHcn1\nHcrt\nNpy\nOxt\nPomc\nSlc17a6\nSlc17a6\nSlc32a1\nSnap25\nSnap25\nSyt1\nSyt1\nSyt2\nTrh\nAgt\nAgt\nAqp4\nGfap\nMyo6\nApoe\nApod\nCldn11\nEfnb3\nMag\nMal\nMbp\nMobp\nMog\nMog\nOpalin\nPlp1\nPtgds\nCspg4\nNeu4\nPdgfra\nSox10\nC1qa\nCsf1r\nCx3cr1\nHexb\nP2ry12\nTmem119\nAbcg2\nSlc47a1\nCcdc153\nCldn5\nAdgrf5\nCldn5\nKcnj8\nGfap\nGpr50\nRax\nIgfbp2\nTbx18\nVtn" %>%
  strsplit(,split=" ") %>%
  unlist() %>%
  unique() # These are the known genes

markers <- read.table(file_markers, header=TRUE) %>%
  filter(avg_log2FC>0 & pct.1>pct.2 & p_val_adj<0.05 & gene %in% goi) %>%
  filter(gene %in% rownames(sc_obj)) %>% # remove genes not in the "integrated"
  distinct(gene, .keep_all=TRUE)
markers$cluster <- factor(markers$cluster, levels=c(0,5,1:4,6))
markers <- markers %>% arrange(cluster)

# store annotations into a new metadata.
Idents(sc_obj) <- sc_obj$main_cluster
Idents(sc_obj) <- factor(Idents(sc_obj), levels=levels(sc_obj)[c(1,7,2,3,5,4,6)])

plot <- DotPlot(sc_obj, features=markers$gene) +
  scale_color_gradientn(colors=brewer.pal(n=9,"GnBu")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Cluster") +
  xlab("Gene")
filename=paste0("dotplot_label.pdf")
ggsave(plot, filename=filename, width=10, height=4)

# .......... UMAP
plot <- DimPlot(sc_obj, group.by="ident", reduction="umap", label=TRUE, repel=TRUE, pt.size=0.25) +
  scale_color_brewer(type="qual",aesthetics="col", palette=3, direction=-1) +
  theme(legend.position="none")
filename=paste0("umap_clusters_res",res,".pdf")
ggsave(plot=plot, filename=filename, width=7, height=7)
plot <- DimPlot(sc_obj, group.by="ident", reduction="umap", pt.size=0.25) +
  scale_color_brewer(type="qual",aesthetics="col", palette=3, direction=-1) +
  theme(legend.position="none")
filename=paste0("umap_clusters_res",res,".png")
ggsave(plot=plot, filename=filename, width=7, height=7)

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
col_reps <- c("#4daf4a","#984ea3") # colorbrewer2.0, qualitative, 5th set.

plot <- DimPlot(sc_obj, cols=col_reps, group.by="rep", reduction="umap", label=TRUE, repel=TRUE, pt.size=0.25, shuffle=TRUE) +
  theme(legend.position="none")
filename=paste0("umap_by_repeat.pdf")
ggsave(plot=plot, filename=filename, width=7, height=7)
plot <- DimPlot(sc_obj, cols=col_reps, group.by="rep", reduction="umap", pt.size=0.25, shuffle=TRUE) +
  theme(legend.position="none")
filename=paste0("umap_by_repeat.png")
ggsave(plot=plot, filename=filename, width=7, height=7)
