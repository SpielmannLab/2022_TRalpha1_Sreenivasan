"Make QC and Characterisation Figure for the paper

Usage: figure_qc.R --file_sc_obj=<file> --file_markers=<file> --file_functions=<file> --res=<value>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		The rds file after clustering
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

# Create new metadata called genotype and rep and save
sc_obj$genotype <- gsub(sc_obj$orig.ident, pattern="\\w+-", replacement="")
sc_obj$rep <- gsub(sc_obj$orig.ident, pattern="-\\w+", replacement="")
saveRDS(object=sc_obj, file=gsub(file_sc_obj, pattern=".rds", replacement="x.rds"))

# Doghnut plot of cell counts
cell_counts <- select(sc_obj@meta.data,c("orig.ident")) %>%
  group_by(orig.ident) %>%
  summarize(count=n()) %>%
  as.data.frame() %>%
  mutate(fraction=count/sum(count), ymax=cumsum(fraction), ymin=lag(ymax, default=0), pos=ymax-0.5*fraction) %>%
  mutate(genotype=gsub(orig.ident,pattern="^\\w+-",replacement="")) %>%
  mutate(rep=gsub(orig.ident,pattern="-\\w+$",replacement="")) %>%
  mutate(group=paste0(genotype,"-",rep))
cell_counts <- cell_counts[c(1,3,2,4),]

plot <- ggplot(cell_counts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=reorder(group,c(1,2,4,3)))) +
     geom_rect() +
     geom_text(aes(x=3.5,y=pos, label=paste0(group,"\nN: ",count)), colour="white", size=5, family="sans") +
     coord_polar(theta="y") +
     xlim(c(2, 4)) +
     scale_color_brewer(palette=5, aesthetics="fill", type="div") +
     theme_void() +
     theme(legend.position="none") +
     ggtitle(label=paste0(sum(cell_counts$count), " cells"))
filename=paste0("cell_counts_doughnutplot.pdf")
ggsave(plot=plot, filename=filename, width=7, height=7)

# Get markers and rename idents
goi <- "Hcn1 Slc17a6 Gad1 Gad2 Aqp4 Agt Myo6 Mal Mag Mobp Pdgfra Cspg4 Hexb C1qa Cx3cr1 Snap25 Syt1 Igfbp2 Tbx18 Kcnj8" %>%
# "Agrp\nCart\nDdc\nGad1\nGad1\nGad2\nGad2\nHcn1\nHcrt\nNpy\nOxt\nPomc\nSlc17a6\nSlc17a6\nSlc32a1\nSnap25\nSnap25\nSyt1\nSyt1\nSyt2\nTrh\nAgt\nAgt\nAqp4\nGfap\nMyo6\nApoe\nApod\nCldn11\nEfnb3\nMag\nMal\nMbp\nMobp\nMog\nMog\nOpalin\nPlp1\nPtgds\nCspg4\nNeu4\nPdgfra\nSox10\nC1qa\nCsf1r\nCx3cr1\nHexb\nP2ry12\nTmem119\nAbcg2\nSlc47a1\nCcdc153\nCldn5\nAdgrf5\nCldn5\nKcnj8\nGfap\nGpr50\nRax\nIgfbp2\nTbx18\nVtn" %>%
  strsplit(,split=" ") %>%
  unlist() %>%
  unique() # These are the known genes

markers <- read.table(file_markers, header=TRUE) %>%
  filter(avg_log2FC>0 & pct.1>pct.2 & p_val_adj<0.05 & gene %in% goi) %>%
  distinct(gene, .keep_all=TRUE)

# Reorder, so the Neurons_2 comes after Neurons_1
markers$cluster <- factor(markers$cluster, levels=c(0,5,1:4,6))
markers <- markers %>% arrange(cluster)

cluster_key <- sc_obj@meta.data %>%
  names %>%
  grep(pattern=paste0("snn_res.",res), value=TRUE)
plot <- DotPlot(sc_obj, features=markers$gene, group.by=cluster_key) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
filename=paste0("dotplot_nolabel.pdf")
ggsave(plot, filename=filename, width=15, height=8)

# Name the clusters
Idents(sc_obj) <- sc_obj@meta.data[,cluster_key]
sc_obj <- RenameIdents(sc_obj, "0"= "Neurons_1",
                              "1" = "Astrocytes",
                              "2" = "Oligodendrocytes",
                              "3" = "OPC",
                              "4" = "Microglia, Macrophages",
                              "5" = "Neurons_2",
                              "6" = "Non-neuronal cells")
# reorder the clusters plotting
Idents(sc_obj) <- factor(Idents(sc_obj), levels=levels(sc_obj)[c(1,6,2:5,7)])

plot <- DotPlot(sc_obj, features=markers$gene) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Cluster") +
  xlab("Gene")
filename=paste0("dotplot_label.pdf")
ggsave(plot, filename=filename, width=10, height=4)

plot <- DotPlot(sc_obj, features=markers$gene, cluster.idents=TRUE) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Cluster") +
  ylab("Gene")

filename=paste0("dotplot_label_clustered.pdf")
ggsave(plot, filename=filename, width=8, height=6)

# Plot features from Campbell et al
features_campbell <- c("Ptgds", "Car2", "Trf", "Plp1", "Mag","Cldn11","Hbb-bs","Hbb-bt", "Hba-a1", "Rgs5", "Acta2","Igfbp7", "Bcas1","Gpr17","Fyn","Pdgfra","Fabp7","Lhfpl3","Ctss","C1qb","C1qa","Dcn","Mgp","Apod","Ccdc153","Tmem212","Rarres2","Agt","Slc1a2","Slc1a3","Scn7a","Cldn10","Col25a1","6330403K07Rik","Nnat","Crym","Oxt","Pmch","Atp1b1","Vlp","Rgs16","Dlk1","Tac2","Prlr","Tmem35","Gal","Ghrh","Penk","Syt1","Slc18a2","Coch","Npy","Agrp","Crabp1","Ces1d","Epcam","Cyp2f2","Chga","Chgb","Scg2")
plot <- DotPlot(sc_obj, features=features_campbell, cluster.idents=TRUE) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Cluster") +
  ylab("Gene")
filename=paste0("dotplot_label_clustered_campbelfeatures.pdf")
ggsave(plot, filename=filename, width=20, height=6)
plot <- DotPlot(sc_obj, features=features_campbell) +
  scale_color_gradient(low="white", high="black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Cluster") +
  ylab("Gene")
filename=paste0("dotplot_label_campbelfeatures.pdf")
ggsave(plot, filename=filename, width=20, height=6)

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

# ......... Calculate cellular composition and save plot
cluster_key <- "annotations"
sc_obj$annotations <- Idents(sc_obj)
sc_obj$genotype <- gsub(sc_obj$orig.ident, pattern="^\\w+-", replacement="")
metadata <- select(sc_obj@meta.data,c("genotype",all_of(cluster_key)))
# Use the function in cell_composition_calculation.R
output <- cell_composition(metadata=metadata, grouping_var="genotype", xvar=cluster_key)

filename=paste0("cell_composition.pdf")
ggsave(plot=plot_grid(output$plot2 + theme_bw(), output$plot1 + theme_bw(), nrow=1), filename=filename, width=10, height=4)

# .... Make ridgeplots
plot <- ggplot(sc_obj@meta.data, aes(x=doublet_score, y=annotations, fill=annotations)) +
  geom_density_ridges() +
  ylab("") +
  xlab("Doublet score") +
  theme_ridges() +
  scale_color_brewer(type="qual",aesthetics="fill", palette=3, direction=-1)  +
  theme(legend.position = "none", axis.title.x=element_text(hjust=0.5), aspect.ratio=0.8)
filename=paste0("doublet_ridgeplot.pdf")
ggsave(plot=plot, filename=filename, width=7, height=7)

plot <- ggplot(sc_obj@meta.data, aes(x=pct_mt, y=annotations, fill=annotations)) +
  geom_density_ridges() +
  ylab("") +
  xlab("Percentage mitochondrial reads") +
  theme_ridges() +
  scale_color_brewer(type="qual",aesthetics="fill", palette=3, direction=-1)  +
  theme(legend.position = "none", axis.title.x=element_text(hjust=0.5), aspect.ratio=0.8)
filename=paste0("mitochondrial_ridgeplot.pdf")
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
