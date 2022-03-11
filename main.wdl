workflow SctkSeuratCluster {
    
    # An integer matrix of counts
    File raw_counts

    # The number of PCA dimensions to use PRIOR to clustering.
    Int pca_dims

    call runSeurat {
        input:
            raw_counts = raw_counts,
            pca_dims = pca_dims
    }

    output {
        File seurat_output = runSeurat.fout
        File cluster_counts = runSeurat.cluster_counts
    }
}

task runSeurat {
    File raw_counts
    Int pca_dims

    String output_name = 'seurat_cluster.tsv'
    String cluster_counts_name = 'cluster_counts.json'

    Int disk_size = 20

    command {
        Rscript /opt/software/seurat_cluster.R ${raw_counts} ${pca_dims} ${output_name} ${cluster_counts_name}
    }

    output {
        File fout = "${output_name}"
        File cluster_counts = "${cluster_counts_name}"
    }

    runtime {
        docker: "ghcr.io/web-mev/mev-sctk-seurat-cluster"
        cpu: 2
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
