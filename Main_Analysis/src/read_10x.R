"CellRanger count matrix to Seurat Object

Usage: read_10x.R --infolder=<folder> --samplename=<value>
Options:
  -h --help			Show this screen.
  --infolder=<file>		Path to containg cell ranger barcodes, cells, and matrix
  --samplename=<value>		Name of the sample

" -> doc

# --- Load all the libraries
suppressMessages(library(Seurat))
suppressMessages(library(docopt))

# --- Read all the arguments passed
arguments<-docopt(doc, quoted_args=TRUE)

# --- Parameters
infolder <- arguments$infolder
samplename <- arguments$samplename

message("infolder: ", infolder)
message("samplename: ",samplename)

# --- Run
sc_data <- Read10X(data.dir=infolder) # Read the data
sc <- CreateSeuratObject(counts=sc_data,project=samplename) # Create Seurat Object

print(sc) #preview the created Seurat Object

# --- Save the output rds file
filename <- paste0(samplename,"_10x.rds")

saveRDS(sc, file = filename)
message("The cell ranger count matrix read and seurat object saved as:", filename)
