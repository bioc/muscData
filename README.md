# muscData
Multi-sample multi-group scRNA-seq data

The `muscData` package is aimed at providing a set of publicly available single-cell RNA sequencing (scRNA-seq) datasets with complex experimental designs, i.e., datasets that contain multiple samples (e.g., individuals) measured across multiple experimental conditions (e.g., treatments), formatted into `SingleCellExperiment` (SCE) Bioconductor objects. All provided SCEs contain *unfiltered* raw counts, and any gene and cell metadata available from the original data source.