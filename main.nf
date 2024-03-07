process run_seurat {

    tag "Run SCTK Seurat clustering"
    publishDir "${params.output_dir}/SctkSeuratCluster.seurat_output", mode:"copy", pattern:"${output_name}"
    publishDir "${params.output_dir}/SctkSeuratCluster.cluster_counts", mode:"copy", pattern:"${edge_cluster_counts_name}"
    container "ghcr.io/web-mev/mev-sctk-seurat-cluster"
    cpus 2
    memory '12 GB'

    input:
        path raw_counts

    output:
        path "${output_name}"
        path "${edge_cluster_counts_name}"

    script:
        output_name = 'seurat_clusters.tsv'
        cluster_counts_name = 'cluster_counts.json'
        """
        Rscript /opt/software/seurat_cluster.R \
            ${raw_counts} \
            ${params.pca_dims} \
            ${output_name} \
            ${cluster_counts_name}

        """
}

workflow {
    run_seurat(params.raw_counts)
}