library(DropletUtils)
library(dplyr)
library(matrixStats)
library(muscat)
library(readxl)
library(R.utils)
library(scater)
library(scds)
library(Seurat)
library(SingleCellExperiment)

# download & unpack data from figshare
url <- "https://ndownloader.figshare.com/articles/8976473/versions/1"
dir <- "Crowell19"; tar <- "foo.tar"
download.file(url, destfile = tar)
untar(tar, exdir = dir)
fns <- list.files(dir, ".zip", full.names = TRUE)
for (fn in fns) untar(fn, exdir = dir)
dirs <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
for (fn in list.files(dirs, ".gz", full.names = TRUE)) gunzip(fn)

# prep. data -------------------------------------------------------------------

# load raw counts
fastq_dirs <- list.dirs(dir, recursive = FALSE, full.names = TRUE)
names(fastq_dirs) <- basename(fastq_dirs)
sce <- read10xCounts(fastq_dirs)

# rename row/colData colnames & SCE dimnames
names(rowData(sce)) <- c("ENSEMBL", "SYMBOL")
names(colData(sce)) <- c("sample_id", "barcode")
sce$sample_id <- factor(basename(sce$sample_id))
dimnames(sce) <- list(
    with(rowData(sce), paste(ENSEMBL, SYMBOL, sep = ".")),
    with(colData(sce), paste(barcode, sample_id, sep = ".")))

# load metadata
md <- read_excel(file.path(dir, "metadata.xlsx"))
m <- match(sce$sample_id, md$`Sample ID`)
sce$group_id <- md$Characteristics[m]

# preprocessing ----------------------------------------------------------------

# remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]

# split SCE by sample
cs_by_s <- split(colnames(sce), sce$sample_id)
sce_by_s <- lapply(cs_by_s, function(cs) sce[, cs])

# run 'scds'
sce_by_s <- lapply(sce_by_s, function(u) 
    cxds_bcds_hybrid(bcds(cxds(u))))

# remove doublets
sce_by_s <- lapply(sce_by_s, function(u) {
    # compute expected nb. of doublets (10x)
    n_dbl <- ceiling(0.01 * ncol(u)^2 / 1e3)
    # remove 'n_dbl' cells w/ highest doublet score
    o <- order(u$hybrid_score, decreasing = TRUE)
    u[, -o[seq_len(n_dbl)]]
})

# merge back into single SCE
sce <- do.call("cbind", sce_by_s)

# calculate QC Metrics
mito <- grep("mt-", rownames(sce), value = TRUE)
sce <- addPerCellQC(sce, subsets = list(Mt = mito))

# filtering; get sample-specific outliers
cols <- c("sum", "detected", "subsets_Mt_percent")
type <- c("both", "both", "higher")
log <- c(TRUE, TRUE, FALSE)

drop_cols <- paste0(cols, "_drop")
for (i in seq_along(cols))
    colData(sce)[[drop_cols[i]]] <- isOutlier(sce[[cols[i]]], 
        nmads = 2.5, type = type[i], log = log[i], batch = sce$sample_id)

# drop outlier cells
ol <- rowAnys(as.matrix(colData(sce)[drop_cols]))
sce <- sce[, !ol]

# require count > 1 in at least 20 cells
sce <- sce[rowSums(counts(sce) > 1) >= 20, ]

# intergation, dimensionality reduction & clustering ---------------------------

# increase future's maximum allowed size of objects
# to be exported from default of 500 MB to 2 GB
options(future.globals.maxSize = 2048 * 1024 ^ 2)

# create 'SeuratObject'
so <- CreateSeuratObject(
    counts = counts(sce),
    meta.data = data.frame(colData(sce)),
    project = "Crowell19_LPS")

# split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so <- lapply(cells_by_sample, function(i) subset(so, cells = i))

# normalize, find variable genes, and scale
so <- lapply(so, NormalizeData, verbose = FALSE)
so <- lapply(so, FindVariableFeatures, nfeatures = 2e3, 
    selection.method = "vst", verbose = FALSE)
so <- lapply(so, ScaleData, verbose = FALSE)

# find anchors & integrate
as <- FindIntegrationAnchors(so, verbose = FALSE)
so <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE)

# scale integrated data
DefaultAssay(so) <- "integrated"
so <- ScaleData(so, verbose = FALSE)

# dimensionality reduction
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- RunTSNE(so, reduction = "pca", dims = seq_len(20),
    seed.use = 1, do.fast = TRUE, verbose = FALSE)
so <- RunUMAP(so, reduction = "pca", dims = seq_len(20),
    seed.use = 1, verbose = FALSE)

# clustering
so <- FindNeighbors(so, reduction = "pca", dims = seq_len(20), verbose = FALSE)
so <- FindClusters(so, resolution = 0.2, random.seed = 1, verbose = FALSE)

# convert 'SeuratObject' to SCE
sce <- as.SingleCellExperiment(so, assay = "RNA")

# cell type annoation ----------------------------------------------------------

# set cluster IDs to resolution 0.2 clustering
sce$cluster_id <- FetchData(so, "integrated_snn_res.0.2")[, 1]

# annotation
anno <- list(
    "unassigned" = 18,
    "Astrocytes" = 3,
    "Endothelial" = 13,
    "Microglia" = c(19, 21),
    "Oligodendrocytes" = 4,
    "OPC" = 14, 
    "CPE cells" = 20,
    "Excit. Neuron" = c(0, 1, 2, 5, 6, 7, 8, 11, 15, 17),
    "Inhib. Neuron" = c(9, 10, 12, 16))

m <- match(sce$cluster_id, unlist(anno))
ns <- vapply(anno, length, numeric(1))
lab <- rep.int(names(anno), ns)[m]
sce$cluster_id <- factor(lab, levels = names(anno)[-1])

# remove unassigned
sce <- sce[, !is.na(sce$cluster_id)]

# compute log-library-size normalized counts
sce <- logNormCounts(sce)

# make WT reference group & rename
sce$group_id <- factor(sce$group_id, 
    levels = c("WT", "LPS"), 
    labels = c("Vehicle", "LPS"))

# reorder sample levels
m <- match(levels(sce$sample_id), sce$sample_id)
o <- order(sce$group_id[m])
sce$sample_id <- factor(sce$sample_id, 
    levels = levels(sce$sample_id)[o])

# separate ensembl IDs & gene symbols
ss <- strsplit(rownames(sce), ".", fixed=TRUE)
rowData(sce)$ENSEMBL <- sapply(ss, .subset, 1)
rowData(sce)$SYMBOL <- sapply(ss, .subset, 2)

# prep. SCE for 'muscat'
sce <- prepSCE(sce)
save(sce, file = "Crowell19_4vs4.Rda")
