# Makes plots for Figure 2 based on harmony-integrated ODC cluster (all other clusters are done in si_figure3_plots.R).
# Includes: UMAPs and Violin Plot
# Note: separate rds file for the individual clusters were generated using separate_clusters.R script, and then they were individually integrated using integrate_clusters.R script (Harmony was preferred).

"Make - figure for oligodendrocyte Characterisation

Usage: figure_oligo.R --file_sc_obj_oligo=<file> --file_sc_obj=<file>

Options:
  -h --help			Show this screen.
  --file_sc_obj_oligo=<file>		The rds file of oligodendrocytes after integration
  --file_sc_obj=<file>		The rds file after clustering of the entire dataset
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
file_sc_obj_oligo <- arguments$file_sc_obj_oligo
file_sc_obj <- arguments$file_sc_obj

message("file_sc_obj_oligo: ", file_sc_obj_oligo)
message("file_sc_obj: ", file_sc_obj)

# read file
sc_obj <- readRDS(file_sc_obj_oligo)

# Record which subcluster is being currently processed
coi <- sc_obj$main_cluster %>% unique

# Make featureplots for DE genes between WT and Mutant ODCs.
Idents(sc_obj) <- sc_obj$genotype
# Find marker genes
markers <- FindAllMarkers(sc_obj, logfc.threshold = 0.5, min.pct = 0.2, test.use = "wilcox", min.cells.group=0.5, only.pos=TRUE)

# save the markers
filename="oligo_markers_wtvstra1.tsv"
write.table(markers, file=filename, sep="\t")

# Make FeaturePlot
foi <- markers %>%
  filter(avg_log2FC>=1) %>%
  mutate(pct_diff=pct.1-pct.2) %>%
  group_by(cluster) %>%
  arrange(desc(pct_diff)) %>%
  slice_head(n=5) %>%
  .$gene

#Feature Plot in pdf
plot <- FeaturePlot(sc_obj, features=foi, keep.scale="all", ncol=5) &
  theme(legend.position="none",panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),   axis.ticks = element_blank(), panel.background = element_blank(), line=element_blank(), plot.title=element_text(family="sans", face="plain", hjust=0.5, size=15))
filename=paste0("FeaturePlot_cluster",coi,"_int_degenes.pdf")
ggsave(plot=plot, filename=filename, width=7*5/2, height=7*2/2) # size divided by 4 because it will be 1/4th of the size on the

#Feature Plot in png
filename=paste0("FeaturePlot_cluster",coi,"_int_degenes.png")
ggsave(plot=plot, filename=filename, width=7*5/2, height=7*2/2) # size divided by 4 because it will be 1/4th of the size on the figure

# Make FeaturePlot for hand-chosen genes
# override the foi - these are amongst the top genes.
foi <- c("Nckap5", "Col23a1", "Ccnd3", "Gsn", "Fbn2", "Tmem117", "Neat1", "Dpyd", "Hcn2", "Slco3a1")
#Feature Plot in pdf
plot <- FeaturePlot(sc_obj, features=foi, keep.scale="all", ncol=5) &
  theme(legend.position="none",panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),   axis.ticks = element_blank(), panel.background = element_blank(), line=element_blank(), plot.title=element_text(family="sans", face="plain", hjust=0.5, size=15))
filename=paste0("FeaturePlot_cluster",coi,"_int_degenes_chosen.pdf")
ggsave(plot=plot, filename=filename, width=7*5/2, height=7*2/2) # size divided by 4 because it will be 1/4th of the size on the
#Feature Plot in png
filename=paste0("FeaturePlot_cluster",coi,"_int_degenes_chosen.png")
ggsave(plot=plot, filename=filename, width=7*5/2, height=7*2/2) # size divided by 4 because it will be 1/4th of the size on the figure

# Create Feature Plots for ODC for SI
# Features from LaManno extended data figure 7
features_lamanno <- c("Pdgfra","Lhfpl3","Bmp4","Neu4","Kcnj12","Tmem2","Rras2","Sema4f","Tmem141","Arap2","Erbb3","Opalin","Plekhh1","Mog","Gjb1","Mal")

p <- FeaturePlot(sc_obj, features=features_lamanno, cols=c("lightgrey", "black"), keep.scale="all", pt.size=0.25) &
  theme_void() &
  theme(plot.title=element_text(hjust=0.5, face="italic")) &
  theme(legend.position="none")

filename <- paste0("FeaturePlot_cluster",coi,"_int_lamanno_markers.pdf")
ggsave(p + theme(legend.position="right"), file=filename,width=7,height=7)
filename <- paste0("FeaturePlot_cluster",coi,"_int_lamanno_markers.png")
ggsave(p, file=filename,width=7,height=7)

# *** Make Violin Plots for jens WB genes for the ODC cells
jens_markers <- c("Sox10", "Mbp", "Arsg", "Olig2", "Cnp", "Mog", "Plp1", "Cldn11", "Kcna1")
foi <- jens_markers
plot <- VlnPlot(sc_obj, features=foi, group.by="genotype", ncol=length(foi), cols=c("#ca2027","#0272b0")) &
  theme(text=element_text(family="sans", size=12), plot.title=element_text(face="plain"), axis.title=element_blank(), axis.text=element_text(family="sans", size=12))
filename=paste0("cluster",coi,"_int_jensViolin.png")
ggsave(plot=plot, filename=filename, width=7*length(foi)/4, height=7/2) # size divided by 4 because it will be 1/4th of the size on the figure

# Try to make boxplots
expr_data <- GetAssayData(sc_obj, slot="data", assay="RNA")[foi,] %>%
  t() %>%
  as.data.frame() %>%
  cbind(sc_obj$genotype) %>%
  rename('genotype'=`sc_obj$genotype`)
plot <- list()
for(gene in foi){
  plot[[gene]] <- ggplot(expr_data, aes_string("genotype", gene)) + geom_boxplot(fill=c("#ca2027","#0272b0")) + ggtitle(gene) + xlab("Genotype") + ylab("Expression") +
  theme_classic() +
  theme(text=element_text(family="sans", size=12), plot.title=element_text(face="plain", hjust=0.5))
}
filename=paste0("cluster",coi,"_int_jensBox.png")
ggsave(plot=plot_grid(plotlist=plot, ncol=length(plot)), filename=filename, width=7*length(foi)/4, height=7/2) # size divided by 4 because it will be 1/4th of the size on the figure

rm(sc_obj) # clear the object
# For the entire dataset, not just oligodendrocytes


# *** Make Violin Plots for jens WB genes for all cells
sc_obj <- readRDS(file_sc_obj) #Note the object after clustering somehow does not have metadata

sc_obj$genotype <- gsub(sc_obj$orig.ident, pattern="^\\w+-", replacement="")

# Make plot
jens_markers <- c("Sox10", "Mbp", "Arsg", "Olig2", "Cnp", "Mog", "Plp1", "Cldn11", "Kcna1")
foi <- jens_markers
plot <- VlnPlot(sc_obj, features=foi, group.by="genotype", ncol=length(foi), cols=c("#ca2027","#0272b0"), combine=FALSE)
plot <- lapply(plot, FUN=function(p){
  p+theme(text=element_text(family="sans", size=12), plot.title=element_text(face="plain"), axis.title=element_blank(), axis.text=element_text(family="sans", size=12), legend.position="none")
}) # change the theme of every plot
filename=paste0("allClusters_int_jensViolin.png")
ggsave(plot=plot_grid(plotlist=plot, ncol=length(plot)), filename=filename, width=7*length(foi)/4, height=7/2) # size divided by 4 because it will be 1/4th of the size on the figure
filename=paste0("allClusters_int_jensViolin.pdf")
ggsave(plot=plot_grid(plotlist=plot, ncol=length(plot)), filename=filename, width=7*length(foi)/4, height=7/2) # size divided by 4 because it will be 1/4th of the size on the figure

# Make boxplots
# Try to make boxplots
expr_data <- GetAssayData(sc_obj, slot="data", assay="RNA")[foi,] %>%
  t() %>%
  as.data.frame() %>%
  cbind(sc_obj$genotype) %>%
  rename('genotype'=`sc_obj$genotype`)
plot <- list()
for(gene in foi){
  plot[[gene]] <- ggplot(expr_data, aes_string("genotype", gene)) + geom_boxplot(fill=c("#ca2027","#0272b0")) + ggtitle(gene) + theme_classic()
}
filename=paste0("allClusters_int_jensBox.png")
ggsave(plot=plot_grid(plotlist=plot, ncol=length(plot)), filename=filename, width=7*length(foi)/4, height=7/2) # size divided by 4 because it will be 1/4th of the size on the figure
