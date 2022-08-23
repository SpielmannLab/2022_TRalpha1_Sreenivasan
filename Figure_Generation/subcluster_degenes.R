"Cluster a merged seurat3 object using several resolution points. Adapted from the cluster.R and de_genes.R in scpipeline

Usage: figure_clusters_degenes.R --file_sc_obj=<file> --res=<value> --task=<value> --file_functions=<file>

Options:
  -h --help            Show this screen.
  --file_sc_obj=<file>    Optional. If !is.null, do not use jobname to guess scfilename.
  --res=<value>        String separated by comma with resolutions for clutering
  --task=<value>       integrate-harmony, integrate-seurat or merge
  --file_functions=<file> The .R file containing all functions
"-> doc

# Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(docopt))
suppressMessages(library(RColorBrewer)) #to expand color palette

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters
file_sc_obj <- arguments$file_sc_obj
ress <- as.numeric(unlist(strsplit(arguments$res, ",")))
task <- arguments$task
file_functions <- arguments$file_functions

message("file_sc_obj: ", file_sc_obj)
message("res: ", ress)
message("task: ", task)
message("file_functions: ", file_functions)

# load all functions defined in the r-script
source(file_functions)

if(is.null(task)){
	reduction <- 'pca'
} else if(task=='integrate-harmony') {
	reduction <- 'harmony'
} else {
	reduction <- 'pca'
}

# read file
sc_obj <- readRDS(file_sc_obj)
npcs <- ncol(sc_obj@reductions[[reduction]])

# Remove legacy clustering metadata
new_meta <- sc_obj@meta.data[!grepl("_snn_res", colnames(sc_obj@meta.data))]
rownames(new_meta) <- rownames(sc_obj@meta.data)
sc_obj@meta.data <- new_meta

# ***** Perform clustering at multiple resolutions and DE gene identification
sc_obj <- FindNeighbors(sc_obj, reduction=reduction, dims = seq(npcs))
markers <- list()
high_genes <- list()
colorpalette <- list()
for(res in ress){
  sc_obj <- FindClusters(sc_obj, resolution = res) #Find the clusters
  # Find the markers at each resolution
  markers[[as.character(res)]] <- FindAllMarkers(sc_obj, # find DE genes
    logfc.threshold = 0.5,
		min.pct = 0.5,
		test.use = "wilcox",
    only.pos=TRUE)
  high_genes[[as.character(res)]] <- FindHighExpressGenes(sc_obj, assay="RNA", slot="data", n=50) # find highly expressed genes
  colorpalette[[as.character(res)]] <- colorRampPalette(brewer.pal(7, "Paired"))(length(levels(sc_obj)))
}

cluster_metadata <- grep(colnames(sc_obj@meta.data), pattern="_snn_res", value=TRUE)
plots <- DimPlot(sc_obj, reduction = "umap", group.by = cluster_metadata, combine=FALSE, label=TRUE, repel=TRUE, pt.size=0.25)
for(i in 1:length(cluster_metadata)){
  plots[[i]] <- plots[[i]] +
    ggtitle(cluster_metadata[i]) +
    scale_color_manual(values=rev(colorpalette[[i]])) +
    theme(legend.position="none")
}

filename="umap_clusters.pdf"
ggsave(plot=plot_grid(plotlist=plots, ncol=2), filename=filename, width=7*2, height=7*ceiling(length(ress)/2))

file_sc_obj_clusters <- gsub(file_sc_obj, pattern=".rds", replacement="_clus.rds")
saveRDS(sc_obj, file = file_sc_obj_clusters)

for(i in names(markers)){
  write.table(markers[[i]], file=paste0("Markers_resolution_",i,".tsv"), sep="\t")
}

for(i in names(high_genes)){
  write.table(high_genes[[i]], file=paste0("Highly_Expressed_genes_resolution_",i,".tsv"), sep="\t")
}
