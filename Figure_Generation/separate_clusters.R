"Generate Separate RDS file for individual clusters and run basic Seurat analysis

Usage: separate_clusters.R --file_sc_obj=<file> --res=<value>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		The rds file after clustering
  --res=<value>           The clustering resolution to be used.
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
res <- as.numeric(arguments$res)

message("file_sc_obj: ", file_sc_obj)
message("res: ", res)

# read file
sc_obj <- readRDS(file_sc_obj)

# Other parameters
assay <- "RNA" #The default assay

# Choose the metadata corresponding to the provided resolution
cluster_key <- sc_obj@meta.data %>%
  names %>%
  grep(pattern=paste0("snn_res.",res), value=TRUE)
Idents(sc_obj) <- sc_obj@meta.data[,cluster_key]
DefaultAssay(sc_obj) <- assay

# ***** subset clusters and make umap
# first define functions to run seurat analysis
seurat_analysis <- function(sc_obj=sc_obj_subset, assay=assay){
  # set npcs based on the number of cells
  if(ncol(sc_obj) > 1000){
          npcs <- 20
          nhvg <- 4000
          } else {
          npcs = 10
          nhvg <- 1000
          }

  # Analyse using Seurat
  set.seed(111)
  sc_obj <- sc_obj %>%
    NormalizeData(normalization.method="LogNormalize") %>%
    FindVariableFeatures(selection.method="vst", nfeatures=nhvg) %>%
    ScaleData() %>%
    RunPCA(npcs=npcs) %>%
    RunUMAP(dims=1:npcs, seed.use=111)
}

for(coi in levels(sc_obj)){
  sc_obj_subset <- subset(sc_obj, idents=coi) %>%
    seurat_analysis(assay=assay)

  file_sc_obj_subset <- gsub(file_sc_obj, pattern=".rds", replacement=paste0("_cluster",coi,".rds"))
  saveRDS(object=sc_obj_subset, file=file_sc_obj_subset)
}
