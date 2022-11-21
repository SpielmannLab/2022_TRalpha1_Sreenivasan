# Usage: Rscript hypoMap_integration.R file_sc_obj file_hm_obj file_functions

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(harmony))

# Get command arguments
args <- commandArgs(trailingOnly=TRUE)
file_hm_obj <- args[2]
file_sc_obj <- args[1]
file_functions <- args[3]

# source functions
source(file_functions)

# Read seurat objects
hm_obj <- readRDS(file_hm_obj)
# Load our data
sc_obj <- readRDS(file_sc_obj)

# Subset odc and opcs for integration
hm_o_cells <- hm_obj@meta.data %>% filter(Author_Class_Curated=="Oligodendrocytes" | Author_Class_Curated=="NG/OPC") %>% rownames()
hm_o_obj <- subset(hm_obj, cells=hm_o_cells)
hm_o_obj$orig.ident <- "HypoMap"
hm_o_obj$study <- "HypoMap"

# Subset ods and opcs from WT sample for integration
sc_o_cells <- sc_obj@meta.data %>% filter(main_cluster=="Oligodendrocytes" | main_cluster=="OPC") %>% filter(genotype=="wt") %>% rownames()
sc_o_obj <- subset(sc_obj, cells=sc_o_cells)
sc_o_obj$Dataset <- paste0("Sreenivasan10x_", sc_o_obj$rep)
sc_o_obj$Author_Class_Curated <- sc_o_obj$main_cluster
sc_o_obj$study <- "Sreenivasan10x"

# Choose number of npcs and nhvg
if(ncol(sc_o_obj) > 1000){
  npcs <- 20
  nhvg <- 4000
} else {
  npcs = 10
  nhvg <- 1000
}

# source perform_merge function from scpipeline
merge_obj <- perform_merge(list(sc_o_obj, hm_o_obj), nhvg=nhvg)

# Perform Harmony integration for different parameters - worked well with the sub-sampled dataset, with nclust=10 and only WT cells
nclust <- 2
merge_obj <- perform_harmony_integration(sc=merge_obj, npcs=npcs, nclust=nclust, group.by.vars="Dataset") %>%
    RunUMAP(dims=1:npcs, seed.use=111, reduction="harmony")

Idents(merge_obj) <- paste(merge_obj$study, merge_obj$Author_Class_Curated, sep="-")
d <- DimPlot(merge_obj, reduction="umap", shuffle=TRUE, raster=TRUE, pt.size=1, cols=c("#737373", "#bdbdbd", "#b2df8a", "#33a02c"), order=c("Sreenivasan10x-Oligodendrocytes", "Sreenivasan10x-OPC")) +
  theme(legend.position="right") +
  theme_void() +
  theme(aspect.ratio=1) +
  theme(plot.title=element_blank())
filename <- paste0("UmapPlot_hypmap_mittag_harmony_nclust",nclust,".pdf")
ggsave(d, file=filename, width=7,height=7)

features_lamanno <- c("Pdgfra","Lhfpl3","Bmp4","Neu4","Kcnj12","Tmem2","Rras2","Sema4f","Tmem141","Arap2","Erbb3","Opalin","Plekhh1","Mog","Gjb1","Mal")
p <- FeaturePlot(merge_obj, features=features_lamanno, cols=c("#bababa", "#ca0020"), pt.size=3, keep.scale="all", raster=TRUE) &
  theme(aspect.ratio=1) &
  theme_void() &
  theme(plot.title=element_text(hjust=0.5, face="italic")) &
  theme(legend.position="none")
filename <-paste0( "FeaturePlot_hypomap_mittag_harmony_nclust",nclust,"_lamanno_markers.pdf")
ggsave(p + theme(legend.position="right"), file=filename,width=10,height=10)
