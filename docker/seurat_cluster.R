suppressMessages(suppressWarnings(library(c("singleCellTK", "Matrix"))))

# args from command line:
args <- commandArgs(TRUE)
RAW_COUNT_MATRIX <- args[1]
OUTPUT_UMAP_BASE <- 'umap_matrix.tsv'

# change the working directory to co-locate with the counts file:
working_dir <- dirname(RAW_COUNT_MATRIX)
setwd(working_dir)

# Import counts as a data.frame
# expects the data frame write output as:
# write.table(
#     counts(sce),
#     file = filename,
#     col.names = TRUE,
#     row.names = TRUE,
#     sep = "\t",
#     quote = FALSE
# )
counts <- read.table(
    file = RAW_COUNT_MATRIX,
    sep = "\t",
    row.names = 1
)
counts <- as(as.matrix(counts), "sparseMatrix")

# Create an SCE object from the counts
sce <- SingleCellExperiment(
    assays=list(counts=counts)
)

# Pre-process the SCE object
sce <- scater_logNormCounts(sce, "logcounts")
# Go through the Seurat curated workflow to get basic 
# scaled assay and dimension reduction
# Because Seurat is a 3rd party package, it requires manually calling
# much of the necessary pre-processing.
sce <- seuratNormalizeData(inSCE = sce, useAssay = "counts")
sce <- seuratFindHVG(inSCE = sce, useAssay = "seuratNormData")
sce <- seuratScaleData(inSCE = sce, useAssay = "seuratNormData")
sce <- seuratPCA(inSCE = sce, useAssay = "seuratScaledData")
sce <- seuratRunUMAP(sce)

# Find the clusters with Seurat
sce <- seuratFindClusters(
    inSCE = sce, 
    useAssay = "seuratScaledData",
    useReduction = "pca",
    dims = 10, # How many components to use for clusters
    algorithm = "louvain", # probably best to leave alone
    resolution = 0.8
)

# Create a data frame from the Seurat factors in the SCE object
# cell_barcode directly referencing the counts colnames is a bit of a hack
# No guarantee order was preserved in the creation of the 
# SingleCellExperiment() object.
df.seurat <- data.frame(
    cell_barcode = as.vector(colnames(counts)),
    seurat_cluster = as.vector(sce$Seurat_louvain_Resolution0.8)
)

# Write to file
output_filename <- paste(
    working_dir, 
    OUTPUT_UMAP_BASE, 
    sep='/'
)
write.table(
    df.seurat, 
    output_filename, 
    sep='\t', 
    quote=F, 
    row.names = FALSE
)