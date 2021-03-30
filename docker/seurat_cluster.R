suppressMessages(suppressWarnings(library("singleCellTK")))
suppressMessages(suppressWarnings(library("Matrix")))

# args from command line:
args <- commandArgs(TRUE)
RAW_COUNT_MATRIX <- args[1]
DIMS <- as.integer(args[2])
OUTPUT_CLUSTER_MAPPING <- 'seurat_clusters.tsv'

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
cnts <- read.table(
    file = RAW_COUNT_MATRIX,
    sep = "\t",
    row.names = 1
)
cnts <- as(as.matrix(cnts), "sparseMatrix")

# Create an SCE object from the counts
sce <- SingleCellExperiment(
    assays=list(counts=cnts)
)

# Pre-process the SCE object
# Go through the Seurat curated workflow to get basic 
# scaled assay and dimension reduction
# Because Seurat is a 3rd party package, it requires manually calling
# much of the necessary pre-processing.
sce <- seuratNormalizeData(inSCE = sce, useAssay = "counts")
sce <- seuratFindHVG(inSCE = sce, useAssay = "seuratNormData")
sce <- seuratSCTransform(inSCE = sce, normAssayName = "SCTCounts", useAssay = "counts")
sce <- seuratPCA(inSCE = sce, useAssay = "SCTCounts", nPCs=DIMS)

# Find the clusters with Seurat
sce <- seuratFindClusters(
    inSCE = sce, 
    useAssay = "SCTCounts",
    useReduction = "pca",
    dims = DIMS, # How many components to use for clusters
    algorithm = "louvain", # probably best to leave alone
    resolution = 0.8
)

# Create a data frame from the Seurat factors in the SCE object
# cell_barcode directly referencing the counts colnames is a bit of a hack
# No guarantee order was preserved in the creation of the 
# SingleCellExperiment() object.
df.seurat <- data.frame(
    cell_barcode = as.vector(colnames(cnts)),
    seurat_cluster = as.vector(sce$Seurat_louvain_Resolution0.8)
)

output_filename <- paste(working_dir, OUTPUT_CLUSTER_MAPPING, sep='/')
write.table(
    df.seurat, 
    output_filename, 
    sep='\t', 
    quote=F, 
    row.names = FALSE
)

# to work with MEV, need to create an outputs file
json_str = paste0('{"seurat_clusters":"', output_filename, '"}')
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)
