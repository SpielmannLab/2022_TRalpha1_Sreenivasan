# Usage: Rscript figure_pvalb_neurons.R file_sc_obj

suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(docopt))
suppressMessages(library(RColorBrewer))

# file_sc_obj <- "/data/humangen_mouse/hypthmsMittag/analysis/scpipeline_analysis/figures/figure_neurons/hypthmsMittag-merge_clusteringx_cluster0_int_clus.rds"
args <- commandArgs(trailingOnly=TRUE)
file_sc_obj <- args[1]
sc_obj <- readRDS(file_sc_obj)

pvalb_cells <- WhichCells(sc_obj, expression=Pvalb>0 & Syt2>0)
# pvalb_cells <- WhichCells(sc_obj, expression=Pvalb>0)

sc_obj$pvalb <- "Other"
sc_obj$pvalb[pvalb_cells] <- "Parvalbuminergic"
Idents(sc_obj) <- sc_obj$pvalb

pvalb_cellno <- sc_obj@meta.data[pvalb_cells,] %>%
  group_by(genotype, rep) %>%
  dplyr::summarize(pvalb_count=n()) %>%
  as.data.frame()

total_cellno <- sc_obj@meta.data %>%
  group_by(genotype, rep) %>%
  dplyr::summarize(total_count=n()) %>%
  as.data.frame()

cellno <- full_join(x=pvalb_cellno, y=total_cellno, by=c("genotype", "rep")) %>%
  mutate(pvalb_pt=pvalb_count/total_count*100)


sc_obj_pvalb <- subset(sc_obj, ident="Parvalbuminergic")
Idents(sc_obj_pvalb) <- sc_obj_pvalb$genotype
goi <- FindAllMarkers(sc_obj_pvalb, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.5) %>%
#  filter(p_val_adj < 1E-20) %>%
  group_by(cluster) %>%
  slice_head(n=3) %>%
  .$gene

# t.test
sink("statistics.txt")
tra1=cellno %>% filter(genotype=="tra1") %>% .$pvalb_pt
wt=cellno %>% filter(genotype=="wt") %>% .$pvalb_pt
t.test(x=tra1, y=wt)
sink(NULL)

#Make and save plots
filename="FeatureScatter_Pvalb_neurons.pdf"
plot <- FeatureScatter(sc_obj, feature1="Pvalb", feature2="Syt2", raster=TRUE, group.by="pvalb", cols=c("#bababa", "#ca0020"), pt.size=3) +
  theme_classic() +
  theme(aspect.ratio=1,
    text=element_text(family="sans"),
    axis.text = element_text(size=12),
    legend.text = element_text(size=12))
  # labs(title="Identification of Parvalbuminergic neurons", subtitle="Based on co-expression")
ggsave(plot, filename=filename, width=5, height=5)

filename <- "umap_Pvalb_neurons.pdf"
plot <- DimPlot(sc_obj, cols=c("#bababa", "#ca0020"), raster=TRUE, shuffle=TRUE, group.by="pvalb", order=c("Parvalbuminergic", "Other"), pt.size=3) +
    theme_classic() +
    theme(aspect.ratio=1,
      text=element_text(family="sans"),
      axis.text = element_text(size=12),
      legend.text = element_text(size=12))
ggsave(plot, filename=filename, width=5, height=5)

filename <- "cell_composition_Pvalb_neurons.pdf"
plot <- ggplot(data=cellno, aes(x=genotype, fill=rep, y=pvalb_pt, label=paste0(signif(pvalb_pt, digits=2), "%"))) +
  geom_col(position="stack", width=0.7) +
  theme_classic() +
  scale_fill_manual(values=c("#4daf4a","#984ea3")) +
  geom_text(position=position_stack(vjust=0.5)) +
  ylab("% of Pvalb+/Syt2+ Neurons") +
  xlab("Genotype") +
  theme(text=element_text(family="sans"),
    axis.text = element_text(size=12),
    legend.text = element_text(size=12))
ggsave(plot, filename=filename, width=4, height=5)


plot <- FeaturePlot(sc_obj, features=c("Pvalb", "Syt2"), keep.scale="all", ncol=2, raster=TRUE, cols=c("#bababa", "#ca0020"), pt.size=3) &
  theme(aspect.ratio=1,
      legend.position="none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      line=element_blank(),
      plot.title=element_text(family="sans", face="italic", hjust=0.5, size=15))

filename <- "featureplot_Pvalb_neurons.pdf"
ggsave(plot=plot + labs(colour="Expression") + theme(legend.position="right", legend.title=element_text(size=12)), filename=filename, width=12, height=5)

pdf("Neurons_Pvalb_other_plots.pdf")
ggplot(data=cellno, aes(x=genotype, fill=rep, y=pvalb_count)) +
  geom_col(position="stack") +
  theme_classic() +
  scale_fill_manual(values=c("#4daf4a","#984ea3")) +
  ylab("Number of Pvalb+/Syt2+ Neurons") +
  xlab("Genotype") +
  theme(text=element_text(family="sans"),
    axis.text = element_text(size=12),
    legend.text = element_text(size=12)) +
  ggtitle("Drop in number of Parvalbuminergic neurons in Tra1 mutant")
VlnPlot(sc_obj_pvalb, features=goi, group.by="genotype")
dev.off()


# Finding highly expressed genes
gene_expression <- AverageExpression(sc_obj,
  features=rownames(sc_obj),
  group.by=c("pvalb", "genotype"), assay="RNA")[["RNA"]]

gene_expression1 <- gene_expression %>%
  data.frame() %>%
  mutate(across(where(is.numeric), ~ signif(.x, digits=4))) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(desc(Parvalbuminergic_tra1))

write.table(gene_expression1, file="Highly_expressed_genes.tsv", sep="\t", row.names=FALSE)
