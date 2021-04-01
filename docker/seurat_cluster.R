suppressMessages(suppressWarnings(library("singleCellTK")))
suppressMessages(suppressWarnings(library("Matrix")))
suppressMessages(suppressWarnings(library("dplyr")))
suppressMessages(suppressWarnings(library("rjson", character.only=T, warn.conflicts = F, quietly = T)))

# args from command line:
args <- commandArgs(TRUE)
RAW_COUNT_MATRIX <- args[1]
DIMS <- as.integer(args[2])
OUTPUT_CLUSTER_MAPPING <- args[3]
CLUSTER_COUNTS <- args[4]

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
    row.names = 1,
    header=T
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

write.table(
    df.seurat, 
    OUTPUT_CLUSTER_MAPPING, 
    sep='\t', 
    quote=F, 
    row.names = FALSE
)

# count the number in each cluster for summary
count.df = df.seurat %>% count(seurat_cluster)
j = toJSON(setNames(count.df$n, count.df$seurat_cluster))
write(j, CLUSTER_COUNTS)

