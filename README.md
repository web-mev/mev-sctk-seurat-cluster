# mev-sctk-seurat-cluster

This repository contains a WDL-format Cromwell-compatible workflow for executing clustering single-cell RNA-seq data using the Seurat package (https://satijalab.org/seurat/) exposed via the Single-Cell Toolkit (https://github.com/compbiomed/singleCellTK).


The outputs include:
- a tab-delimited file giving the cell barcodes assigned to the clusters
- a JSON-format count summary of how many cells were assigned to each cluster.

---

To use, simply fill in the the `inputs.json` with the path to the single-cell counts file and the number of PCA-dimensions (for dimensional reduction-- 50 is a reasonable default). Then submit to a Cromwell runner. 

Alternatively, you can pull the docker image (https://github.com/web-mev/mev-sctk-seurat-cluster/pkgs/container/mev-sctk-seurat-cluster), start the container, and run: 

```
Rscript /opt/software/seurat_cluster.R \
    <path to raw counts tab-delimited file> \
    <number of PCA dimensions, as integer> \
    <output file name for cluster assignments \
    <output file name for cluster count summary>
```