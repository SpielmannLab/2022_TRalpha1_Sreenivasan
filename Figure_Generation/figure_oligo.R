# Makes plots for Figure 2 based on harmony-integrated ODC cluster (all other clusters are done in si_figure3_plots.R).
# Includes: UMAPs and Violin Plot
# Note: separate rds file for the individual clusters were generated using separate_clusters.R script, and then they were individually integrated using integrate_clusters.R script (Harmony was preferred).

"Make - figure for oligodendrocyte Characterisation

Usage: figure_oligo.R --file_sc_obj_oligo=<file>

Options:
  -h --help			Show this screen.
  --file_sc_obj_oligo=<file>		The rds file of oligodendrocytes after integration
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

message("file_sc_obj_oligo: ", file_sc_obj_oligo)

# read file
sc_obj <- readRDS(file_sc_obj_oligo)

# Record which subcluster is being currently processed
coi <- sc_obj$main_cluster %>% unique

# Make featureplots for DE genes between WT and Mutant ODCs.
Idents(sc_obj) <- sc_obj$genotype
# Find marker genes
min_pct <- 0.1
markers <- FindAllMarkers(sc_obj, logfc.threshold = 0.5, min.pct = min_pct, test.use = "wilcox", min.cells.group=0.5, only.pos=TRUE)

# save the markers
filename=paste0("oligo_markers_wtvstra1_min_pct_",min_pct,".tsv")
write.table(markers, file=filename, sep="\t")

# Make FeaturePlot
foi <- markers %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 1E-20) %>%
  mutate(pct_diff=pct.1-pct.2) %>%
  group_by(cluster) %>%
  arrange(desc(pct_diff)) %>%
  slice_head(n=5) %>%
  .$gene

#Feature Plot in pdf
plot <- FeaturePlot(sc_obj, features=foi, keep.scale="all", ncol=5, cols=c("#bababa", "#ca0020")) &
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
plot <- FeaturePlot(sc_obj, features=foi, keep.scale="all", ncol=5, cols=c("#bababa", "#ca0020")) &
  theme(legend.position="none", panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),   axis.ticks = element_blank(), panel.background = element_blank(), line=element_blank(), plot.title=element_text(family="sans", face="plain", hjust=0.5, size=15))
filename=paste0("FeaturePlot_cluster",coi,"_int_degenes_chosen.pdf")
ggsave(plot=plot + theme(legend.position="right"), filename=filename, width=7*5/2, height=7*2/2) # size divided by 4 because it will be 1/4th of the size on the
#Feature Plot in png
filename=paste0("FeaturePlot_cluster",coi,"_int_degenes_chosen.png")
ggsave(plot=plot, filename=filename, width=7*5/2, height=7*2/2) # size divided by 4 because it will be 1/4th of the size on the figure

# Violin plot of top 10 DE genes per genotype for Supplementary Figure
# Make FeaturePlot
foi_tra1 <- markers %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 1E-20) %>%
  mutate(pct_diff=pct.1-pct.2) %>%
  filter(cluster=="tra1") %>%
  arrange(desc(pct_diff)) %>%
  .$gene

foi_wt <- markers %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 1E-20) %>%
  mutate(pct_diff=pct.1-pct.2) %>%
  filter(cluster=="wt") %>%
  arrange(desc(pct_diff)) %>%
  .$gene

plot <- VlnPlot(sc_obj, features=foi_tra1, group.by="genotype", ncol=7, cols=c("#ca2027","#0272b0")) &
  theme(text=element_text(family="sans", size=12),
    plot.title=element_text(size=10, face="italic"),
    axis.title=element_blank(),
    axis.text=element_text(family="sans", size=12),
    legend.position="none",
    aspect.ratio=1) &
    scale_x_discrete(breaks=c("tra1","wt"), labels=c("TR\U03C3\U0031+m","WT"))
filename=paste0("ViolinPlot_cluster",coi,"_int_degenes_tra1.png")
ggsave(plot=plot, filename=filename, width=1.5*7, height=2*7) # size divided by 4 because it will be 1/4th of the size on the figure

plot <- VlnPlot(sc_obj, features=foi_wt[1:35], group.by="genotype", ncol=7, cols=c("#ca2027","#0272b0")) &
  theme(text=element_text(family="sans", size=12),
    plot.title=element_text(size=10, face="italic"),
    axis.title=element_blank(),
    axis.text=element_text(family="sans", size=12),
    legend.position="none",
    aspect.ratio=1)  &
    scale_x_discrete(breaks=c("tra1","wt"), labels=c("TR\U03C3\U0031+m","WT"))
filename=paste0("ViolinPlot_cluster",coi,"_int_degenes_wt_1.png")
ggsave(plot=plot, filename=filename, width=1.5*7, height=2*7)

plot <- VlnPlot(sc_obj, features=foi_wt[36:66], group.by="genotype", ncol=7, cols=c("#ca2027","#0272b0")) &
  theme(text=element_text(family="sans", size=12),
    plot.title=element_text(size=10, face="italic"),
    axis.title=element_blank(),
    axis.text=element_text(family="sans", size=12),
    legend.position="none",
    aspect.ratio=1) &
    scale_x_discrete(breaks=c("tra1","wt"), labels=c("TR\U03C3\U0031+m","WT"))
filename=paste0("ViolinPlot_cluster",coi,"_int_degenes_wt_2.png")
ggsave(plot=plot, filename=filename, width=1.5*7, height=2*7)

# Create Feature Plots for ODC for SI
# Features from LaManno extended data figure 7
features_lamanno <- c("Pdgfra","Lhfpl3","Bmp4","Neu4","Kcnj12","Tmem2","Rras2","Sema4f","Tmem141","Arap2","Erbb3","Opalin","Plekhh1","Mog","Gjb1","Mal")

p <- FeaturePlot(sc_obj, features=features_lamanno, cols=c("#bababa", "#ca0020"), keep.scale="all", pt.size=0.25) &
  theme_void() &
  theme(plot.title=element_text(hjust=0.5, face="italic")) &
  theme(legend.position="none")

filename <- paste0("FeaturePlot_cluster",coi,"_int_lamanno_markers.pdf")
ggsave(p + theme(legend.position="right"), file=filename,width=7,height=7)
filename <- paste0("FeaturePlot_cluster",coi,"_int_lamanno_markers.png")
ggsave(p, file=filename,width=7,height=7)

# *** Make Feature Plots for jens WB genes for the ODC cells for SI
jens_markers <- c("Sox10", "Mbp", "Arsg", "Olig2", "Cnp", "Mog", "Plp1","Cldn11")
plot <- FeaturePlot(sc_obj, features=jens_markers, keep.scale="all", ncol=4, cols=c("#bababa", "#ca0020"), raster=TRUE) &
  theme(legend.position="none", panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),   axis.ticks = element_blank(), panel.background = element_blank(), line=element_blank(), plot.title=element_text(family="sans", face="plain", hjust=0.5, size=15))
filename=paste0("FeaturePlot_cluster",coi,"_protein_markers.pdf")
ggsave(plot=plot + theme(legend.position="right"), filename=filename, width=7*4/2, height=7*2/2) # size divided by 4 because it will be 1/4th of the size on the
#Feature Plot in png
filename=paste0("FeaturePlot_cluster",coi,"_protein_markers.png")
ggsave(plot=plot, filename=filename, width=7*5/2, height=7*2/2) # size divided by 4 because it will be 1/4th of the size on the figure
