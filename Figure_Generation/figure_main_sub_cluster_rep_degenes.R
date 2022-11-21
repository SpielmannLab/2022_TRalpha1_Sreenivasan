"Collect sub cluster level annotations into the main dataset and generate pairwise genotype-based DE analysis for each main_sub_clusters

Usage: main_sub_cluster_collectAnnotations_n_separate.R --file_sc_obj=<file> --file_annotations_neurons=<file> --file_annotations_astrotany=<file>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		The rds file after clustering
  --file_annotations_neurons=<file>   tsv file containing main_cluster and sub_cluster annotations
  --file_annotations_astrotany=<file>  tsv file containing main_cluster and sub_cluster annotations
" -> doc

# Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))
suppressMessages(library(ggplot2))
suppressMessages(library(forcats))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters: Read in the output of the sc_multi_sample pipline
file_sc_obj <- arguments$file_sc_obj
file_annotations_neurons <- arguments$file_annotations_neurons
file_annotations_astrotany <- arguments$file_annotations_astrotany

message("file_sc_obj: ", file_sc_obj)
message("file_annotations_neurons: ", file_annotations_neurons)
message("file_annotations_astrotany: ", file_annotations_astrotany)

# read files
sc_obj <- readRDS(file_sc_obj)
annotations_neurons <- read.table(file=file_annotations_neurons, sep="\t")
annotations_astrotany <- read.table(file=file_annotations_astrotany, sep="\t")

# Choose the metadata corresponding to the provided resolution
sc_obj$main_sub_cluster <- sc_obj$main_cluster

annotations <- rbind(annotations_neurons, annotations_astrotany) %>%
  mutate(main_sub_cluster=paste0(main_cluster,"_",sub_cluster)) %>%
  select(main_sub_cluster)
sc_obj@meta.data[rownames(annotations), "main_sub_cluster"] <- as.character(annotations$main_sub_cluster)

# prepare for DE genes
Idents(sc_obj) <- paste0(sc_obj$rep, "_", sc_obj$main_sub_cluster)

# calculate DE genes
min_pct <- 0.1
main_sub_clusters <- unique(sc_obj$main_sub_cluster)
all_markers <- data.frame()
for(main_sub_cluster in main_sub_clusters){
  ident1 <- paste0("rep1_", main_sub_cluster)
  ident2 <- paste0("rep2_", main_sub_cluster)
  #compute markers for tra1
  markers_rep1 <- FindMarkers(sc_obj, ident.1=ident1, ident.2=ident2, logfc.threshold = 0.5, min.pct = min_pct, test.use = "wilcox", min.cells.group=0.5, only.pos=TRUE)
  if(NROW(markers_rep1)!=0){markers_rep1$rep <- "rep1"}
  #compute markers for wt
  markers_rep2 <- FindMarkers(sc_obj, ident.1=ident2, ident.2=ident1, logfc.threshold = 0.5, min.pct = min_pct, test.use = "wilcox", min.cells.group=0.5, only.pos=TRUE)
  if(NROW(markers_rep2)!=0){markers_rep2$rep <- "rep2"}

  markers <- rbind(markers_rep1, markers_rep2)
  if(NROW(markers)!=0){
    markers$main_sub_cluster <- main_sub_cluster
    markers$gene <- rownames(markers)
  }
  all_markers <- rbind(all_markers, markers)
}

#convert to factor for plotting sake
all_markers$main_sub_cluster <- factor(all_markers$main_sub_cluster, levels=main_sub_clusters)

# Save combined markerlist
filename <- paste0("markers_all_rep_min_pct_0.1_",
  DefaultAssay(sc_obj),
  ".tsv")
write.table(all_markers, file=filename, sep="\t")

# apply a bit stringent statistical filtering
all_markers <- all_markers %>%
  filter(avg_log2FC > 0.5 & p_val_adj<1E-20) %>%
  filter(main_sub_cluster!="Neurons_2" & main_sub_cluster!="Non-neuronal cells")

# plot the number of de genes per main_sub_cluster and save - in manuscript
filename <- paste0("no_of_rep_DEgenes_per_main_sub_cluster_",
  DefaultAssay(sc_obj),
  ".pdf")
p <- ggplot(data=all_markers, aes(x=fct_infreq(main_sub_cluster), fill=rep)) +
  geom_bar(position="stack") +
  scale_fill_manual(values=c("#4daf4a","#984ea3")) +
  scale_x_discrete(drop=FALSE) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    text=element_text(family="sans"),
    axis.text = element_text(size=12),
    legend.text = element_text(size=12)) +
  ylab("Number of positive DE genes") +
  xlab("Annotated clusters or sub-clusters") +
  ggtitle("log2FC > 0.5, p_val_adj < 1E-20, min.pct=10%")
ggsave(p, file=filename, width=10, height=7)

# Neurons_markers - not in manuscript
neurons_markers <- all_markers %>%
  filter(grepl(x=main_sub_cluster, pattern="Neurons_1"))

filename <- paste0("Output_Neuron_rep_markers_",
  DefaultAssay(sc_obj),
  ".txt")
sink(filename, append=TRUE)
print("Distinct neuronal de genes acros reps")
neurons_markers %>%
  .$gene %>%
  table
sink(NULL)

# Get a list of most frequent marker - not in manuscript
filename <- paste0("markers_frequent_rep_min_pct_0.1_",
  DefaultAssay(sc_obj),
  ".tsv")
frequent_degene <- all_markers %>%
  .$gene %>%
  table() %>%
  sort(decreasing=TRUE) %>%
  head(n=5) %>%
  names()
frequent_markers <- all_markers %>%
  filter(gene %in% frequent_degene) %>%
  arrange(gene)
write.table(frequent_markers, file=filename, sep="\t")
