" Merge or integrate two or more Seurat objects. Merge only using RNA assay. Integration using both RNA or SCT assays. Integration algorithm can be Seurat or Harmony

Usage: merge_or_integrate.R --task=<value> --method=<value> --project=<value> --nhvg=<value> --ncores=<value> [--npcs=<value>] [--reference=<file> <file>...] [--integrateBy=<value>] [--sample.tree=<value>] [--nclust=<value>] [--covars=<value>] [<file> <file>...]

Options:
-h --help               	Show this screen.
--task=<value>			Choose: merge, integrate-seurat, integrate-harmony
--method=<value>		Choose: standard or SCT
--project=<value>		Name of the combined dataset, and the output *.rds file
--nhvg=<value>			Number of relevant features. Used for all three tasks
--npcs=<value>			Required only for integrate-seurat or integrate-harmony. Number of PCs.
--reference=<value>		[optional] The # of the sample(s) in the list which needs to be used as the reference during integration. Sample separated
--integrateBy=<value>	Required only for integrate-harmony task. Can be <orig.ident>, <group>, or <replicate>.
--sample.tree=<value>		[optional, Advanced syntax] Specify the order of integration based on the # in the list.
--nclust=<value>		Required only for integrate-harmony
--ncores=<value>	        Number of cores to use for parallelisation using future package
--covars=<value>		[optional] Which metadata should be regressed out? Comma separated. Used for all three tasks

"-> doc

# --- Define functions
perform_merge <- function(sc_list, project="MyProject", nhvg, covars=NULL){
	sc <- merge(x=sc_list[[1]],
		y=sc_list[-1],
		project=project,
		merge.data=TRUE
		)
	sc <- FindVariableFeatures(sc,
		assay = "RNA",
		selection.method = "vst",
		nfeatures = nhvg,
		verbose=FALSE
		)
	sc <- ScaleData(sc,
		features = rownames(sc),
		assay = "RNA",
		vars.to.regress = covars,
		verbose=FALSE
		)
	return(sc)
}

perform_standard_integration <- function(sc_list, reference, anchor.features, reduction="cca", dims=1:30, sample.tree, covars){
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

perform_SCT_integration <- function(sc_list, reference, anchor.features, reduction="cca", dims=1:30, sample.tree){
	features <- SelectIntegrationFeatures(object.list = sc_list,
		nfeatures=anchor.features
		)
	sc_list <- PrepSCTIntegration(object.list = sc_list,
		anchor.features=features,
		verbose=FALSE
		)
	anchors <- FindIntegrationAnchors(object.list = sc_list,
		reference = reference,
		anchor.features = features,
		normalization.method="SCT",
		reduction=reduction,
		dims=dims,
		verbose=FALSE
		)
	sc <- IntegrateData(anchorset = anchors,
		normalization.method = "SCT",
		dims = dims,
		sample.tree = sample.tree,
		verbose=FALSE
		)
	return(sc)
}

perform_harmony_integration <- function(sc, npcs, nclust, integrateBy){
	sc <- RunPCA(sc,
		features = VariableFeatures(object = sc),
		assay = "RNA",
		npcs = npcs,
		approx=FALSE,
		verbose=FALSE
		)
	sc <- RunHarmony(object = sc,
		group.by.vars=integrateBy,
		assay.use="RNA",reduction.save="harmony",
		nclust=nclust
		)
	return(sc)
}

# --- Load all the libraries
suppressMessages(library(docopt))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(future))
suppressMessages(library(harmony))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Parameters
sc_file_list <- arguments$file
task <- arguments$task
method <- arguments$method
project <- arguments$project
ncores <- as.numeric(arguments$ncores)
nhvg <- as.numeric(arguments$nhvg)

if(is.null(arguments$npcs) || arguments$npcs=='NULL' || arguments$npcs =='null') {
	npcs <- NULL
} else {
	npcs <- as.numeric(arguments$npcs)
}

if(is.null(arguments$reference) || arguments$reference=='NULL' || arguments$reference =='null') {
	reference <- NULL
} else {
	reference <- arguments$reference %>% strsplit(split=",") %>% unlist() %>% as.numeric()
}

if(is.null(arguments$integrateBy) || arguments$integrateBy=='NULL' || arguments$integrateBy =='null') {
	integrateBy <- NULL
} else {
	integrateBy <- arguments$integrateBy
}

if(is.null(arguments$sample.tree) || arguments$sample.tree=='NULL' || arguments$sample.tree =='null') {
	sample.tree <- NULL
} else {
	sample.tree <- as.numeric(arguments$sample.tree)
}

if(is.null(arguments$nclust) || arguments$nclust=='NULL' || arguments$nclust =='null') {
	nclust <- NULL
} else {
	nclust <- as.numeric(arguments$nclust)
}

if(is.null(arguments$covars) || arguments$covars=='NULL' || arguments$covars =='null') {
	covars <- NULL
} else {
	covars <- arguments$covars %>% strsplit(split=",") %>% unlist()
}

message("sc_file_list: ", sc_file_list)
message("task: ", task)
message("method: ", method)
message("project: ", project)
message("ncores: ", ncores)
message("nhvg: ", nhvg)
message("npcs: ", npcs)
message("reference: ", reference)
message("integrateBy: ", integrateBy)
message("sample.tree: ", sample.tree)
message("nclust: ", nclust)
message("covars: ", covars)

if(task=="integrate-harmony"){if(any(is.null(nclust),is.null(npcs),is.null(integrateBy))){stop("nclust, npcs, and integrateBy are required for task integrate-harmony")}}
if(task=="integrate-seurat"){if(is.null(npcs)){stop("npcs is required for task integrate-seurat")}}

# --- Setting multicore env for some Seurat applications
plan("multiprocess", workers = ncores)
plan()

# setting parameters of the future package for memory?!
l <-25600*1024^2
options(future.globals.maxSize=l)

# --- Run

#get the sample files and read into a list
message("Reading the input files")
sc_list <- lapply(sc_file_list, FUN=readRDS)
message("Finished reading the files")



# Get the assay and normalization methods
if(method=="standard"){
	assay <- "RNA"
	normalization.method="LogNormalize"
	} else if(method=="SCT") {
		assay <- "SCT"
		normalization.method="SCT"
		} else {
			stop("No appropriate method specified. Choose method=\"RNA\" or \"SCT\"")
			}

# Enforce the default assay to all the seurat objects in the list
sc_list <- lapply(sc_list, FUN=function(sc){
		DefaultAssay(sc) <- assay
		return(sc)
		}
	)
print(sc_list)
#### If the task is to only merge:
if(task=="merge" | task=="integrate-harmony"){

	# Carry out merge only for standard workflow, i.e., assay is "RNA"
	if(!assay=="RNA"){
		stop("Merging or Harmony is only possible with the standard workflow. Please choose method=\"standard\"")
	}
	message("Performing merging on RNA assay")
	sc <- perform_merge(sc_list=sc_list,
		project=project,
		nhvg=nhvg,
		covars=covars
		)
	message("The merged data can be found in the assay \"RNA\" in slots \"data\" and \"scale.data\"")
}

#### If the task is to integrate
if(task=="integrate-seurat"){

	# This can be either Standard Workflow or SCT
	if(assay=="RNA"){
		message("Performing standard integration on RNA assay")

		sc <- perform_standard_integration(sc_list,
			reference=reference,
			anchor.features=nhvg,
			reduction="cca",
			dims=1:npcs,
			sample.tree=sample.tree,
			covars=covars
			)
		message("The integrated data can be found in the assay \"integrated\" in slots \"data\" and \"scale.data\"")
	}

	if(assay=="SCT"){
		message("Performing SCT integration on SCT assay")
		sc <- perform_SCT_integration(sc_list,
			reference=reference,
			anchor.features=nhvg,
			reduction="cca",
			dims=1:npcs,
			sample.tree=sample.tree
			)
		message("The integrated data can be found in the assay \"integrated\" in slot \"scale.data\"")
	}
}

# If the task is to integrate with harmony
if(task=="integrate-harmony"){
	message("Performing Harmony integration on RNA assay")
	sc <- perform_harmony_integration(sc,
		npcs=npcs,
		nclust=nclust,
		integrateBy=integrateBy
		)
	message("The integrated data in found in the reduction \"harmony\"")
}

message("Saving the output file")
saveRDS(object=sc, file=paste0(project,"-", task, "_normalized.rds"))
message("Finished saving the output file")
