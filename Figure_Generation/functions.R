
#****** Calculate cellular composition and make plots - both my style and CX style *****
# Pass on the metadata object containing the columns grouping_var and xvar
# grouping_var is the name of the column containing the Sample information
# sub_trajectory is the name of the cluster or trajectory
# Value: Returns a list of dataframe and plot. Dataframe contains % cells / trajectory / sample normalized to total cells in that sample
cell_composition <- function(metadata, grouping_var="Mutant", xvar="sub_trajectory", label_y_coord=40){

# Calculate the total number of cells per sample type (or per grouping_var type)
count_total <- metadata %>%
                        group_by(across(all_of(grouping_var))) %>%
                        summarise(n=n()) %>%
                        data.frame() %>%
                        tibble::column_to_rownames(var=grouping_var)

# Generate a list of tibbles, each containing only one type of sample (defined by grouping_var)
count_split <- metadata %>%
                        group_by(across(all_of(grouping_var))) %>%
                        group_split()

cellno <- list() # Create an empty list to store the summarized cell count
for(i in seq(NROW(count_total))){
    cellno[[i]] <- count_split[[i]] %>%
                    group_by(across(all_of(xvar))) %>%
                    summarise(n=n())
    cellno[[i]][,grouping_var] <- rownames(count_total)[i]
    cellno[[i]]$n_percent <- cellno[[i]]$n/count_total[i,"n"]*100
    }

cellno1 <- cellno %>% bind_rows() %>% as.data.frame() #Combine list of tibbles into a single dataframe
cellno1$n_percent <- round(cellno1$n_percent, digits=2)

# Make my style plot
plot1 <- ggplot(cellno1,
  #aes(x=reorder(.data[[xvar]], desc(n_percent), sum),
  aes(x=.data[[xvar]], fill=.data[[grouping_var]], y=n_percent)) +
  geom_col(position="fill") +
  scale_fill_manual(values=c("#ca2027", "#0272b0")) +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  ylab("% cells") +
  coord_flip() +
  xlab(NULL)

cellno1 <- cellno1 %>%
  group_by(annotations) %>%
  summarize(n_pooled=sum(n)) %>%
  mutate(n_pooled_k=n_pooled/1000)

plot2 <- ggplot(cellno1,
  # aes(x=reorder(annotations, desc(n_pooled)), y=n_pooled_k)) +
  aes(x=annotations, y=n_pooled_k)) +
  geom_col(fill="black") +
  geom_text(aes(label=n_pooled, y=label_y_coord), angle=0, size=4, vjust=0, color="black") +
  xlab("Cluster") +
  ylab("Total cells per cluster x 1000") +
  coord_flip()

if(FALSE){
# Make a CX style plot to show the difference, but only if there are two grouping_vars. e.g., Mutant and WT
plot2 <- ggplot() + theme_void() # send a blank plot if more than two grouping_vars are present
if(length(unique(cellno1[,grouping_var]))==2){
  cellno2_1 <- cellno[[1]] %>% as.data.frame()
  cellno2_2 <- cellno[[2]] %>% as.data.frame()

  cellno2 <- merge(cellno2_1, cellno2_2, by=xvar,all.x=TRUE, all.y=TRUE) %>%
            mutate(n_percent_diff=n_percent.x/n_percent.y, n_total=n.x+n.y)

  plot2 <- ggplot(cellno2, aes(x=reorder(.data[[xvar]], n_percent_diff), y=n_percent_diff)) +
    geom_col(fill="cornsilk3") +
    scale_y_continuous(trans="log2") +
    geom_text(aes(label=n_total), angle=90, size=2, hjust=1, color="black") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
}
return(list(cellno=cellno, plot1=plot1, plot2=plot2)) #returns the ggplot object and the cell no dataframe
}


# Code copied from sc_pipeline
perform_standard_integration <- function(sc_list, reference, anchor.features, reduction="cca", dims=1:30, sample.tree, covars=NULL){
	anchors <- FindIntegrationAnchors(object.list = sc_list,
		reference = reference,
		anchor.features = anchor.features,
		normalization.method="LogNormalize",
		reduction=reduction,
		dims=dims,
		verbose=FALSE
		)
	sc <- IntegrateData(anchorset=anchors,
		normalization.method="LogNormalize",
		dims=dims,
		sample.tree=sample.tree,
		verbose=FALSE
		)
	sc <- ScaleData(sc,
		features = rownames(sc),
		assay = "integrated",
		vars.to.regress = covars,
		verbose=FALSE
		)
	return(sc)
}

# code copied from sc_pipeline
perform_harmony_integration <- function(sc, npcs, nclust, group.by.vars){
	sc <- RunPCA(sc,
		features = VariableFeatures(object = sc),
		assay = "RNA",
		npcs = npcs,
		approx=FALSE,
		verbose=FALSE
		)
	sc <- RunHarmony(object = sc,
		group.by.vars=group.by.vars,
		assay.use="RNA",reduction.save="harmony",
		nclust=nclust
		)
	return(sc)
}

seurat_to_anndata <- function(sc_obj=sc_obj, assay="RNA", filename="anndata.h5ad"){
	# Saves an AnnData object from the Seurat object containing the following:
	# Normalized "Data" slot in the assay of choice
	# Assay names "spliced" and "unspliced" containing exonic and intronic counts

	# Add libraries
	suppressMessages(library(reticulate))
	suppressMessages(library(Seurat))
	suppressMessages(library(anndata))
	suppressMessages(library(Matrix))

	# Import python libraries using the reticulate package
	# use_condaenv("scVelo", required = TRUE)
	anndata <- import("anndata")

  counts_matrix <- GetAssayData(sc_obj, assay=assay, slot="data") # use data slot, because scVelo does not normalize it
	spliced_matrix <- GetAssayData(sc_obj, assay='spliced', slot='counts')
	unspliced_matrix <- GetAssayData(sc_obj, assay='unspliced', slot='counts')

  # get genes that are common to all three count matrices and conver to dataframe
  genes <- rownames(counts_matrix) %>%
            data.frame()
  rownames(genes) <- genes[,1]
  colnames(genes) <- "id"
  cells <- colnames(counts_matrix)

	# the package anndata is need for the following creation of AnnData object in R
	adata <- AnnData(X=t(counts_matrix[genes$id,]),
			obs=sc_obj@meta.data,
			var=genes,
			layers=list(spliced=t(spliced_matrix[genes$id,]),
				    unspliced=t(unspliced_matrix[genes$id,])),
			obsm=list(X_pca=sc_obj@reductions$pca@cell.embeddings,
				    X_umap=sc_obj@reductions$umap@cell.embeddings))

	# Write to file
	adata$write_h5ad(filename)
	message("Saved AnnData object as: ", filename)
}


FindHighExpressGenes <- function(sc_obj, assay="RNA", slot="data", n=50){
# A custom function to find the top expressed genes per Ident. Not DE genes
# n is the number of highly expressed genes to be found

  # extract data from sc_obj
  grouping <- Idents(sc_obj)
  groups <- levels(sc_obj)
  matrix <- GetAssayData(sc_obj, assay=assay, slot=slot)

  # Define a subfunction to find highly expressed genes per cluster
  find_high_genes_per_group <- function(x){
    cells_in_group <- grouping %in% x
    matrix_sub <- matrix[, cells_in_group]
    top_genes <- Matrix::rowSums(matrix_sub) %>%
      sort(decreasing=TRUE) %>%
      head(n=n)
    return=data.frame(expression=top_genes, cluster=as.numeric(x), gene=names(top_genes))
  }

  high_genes_list <- lapply(X=groups, FUN=find_high_genes_per_group)
  high_genes <- do.call(rbind, high_genes_list)
  return(high_genes)
}
