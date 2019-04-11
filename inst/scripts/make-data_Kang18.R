# description: 
#   10x droplet-based scRNA-seq PBMC data from 8 Lupus patients
#   before and after 6h-treatment with INF-beta (16 samples in total).
#
# availability: 
#   The original data is deposited in the Gene Expression Ombnibus (GEO) 
#   under accession number GSE96583.
#
# reference: 
#   Kang et al., 2019: "Multiplexed droplet single-cell RNA-sequencing 
#   using natural genetic variation". Nature Biotechnology (36), 89–94
#
# link to reference: https://www.ncbi.nlm.nih.gov/pubmed/29227470
# link to raw data:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583
#
# Helena L. Crowell; last modified: April 7th, 2019
# ------------------------------------------------------------------------------

# download data
geo <- "GSE96583"
raw <- getGEOSuppFiles(geo, filter_regex = "(RAW)|(batch2)")

# unpack raw data
raw_fn <- "GSE96583_RAW.tar"
mtx_fns <- c("GSM2560248_2.1.mtx.gz", "GSM2560249_2.2.mtx.gz")
untar(file.path(geo, raw_fn), exdir = geo, files = mtx_fns)

# load cell & gene metadata
cell_md <- read.delim(file.path(geo, "GSE96583_batch2.total.tsne.df.tsv.gz"))
gene_md <- read.delim(file.path(geo, "GSE96583_batch2.genes.tsv.gz"),
    header = FALSE, col.names = c("gene", "feature"))

# load counts
counts <- lapply(file.path(geo, mtx_fns), readMM)
counts <- do.call("cbind", counts)
counts <- as(counts, "dgCMatrix")
dimnames(counts) <- lapply(list(gene_md, cell_md), rownames)

# pull reduced dimensions
cols <- c("tsne1", "tsne2")
tsne <- as.matrix(select(cell_md, cols))
cell_md <- select(cell_md, -cols)

# construct SCE
sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = DataFrame(cell_md),
    rowData = DataFrame(gene_md),
    reducedDims = SimpleList(TSNE = tsne))

# save object & remove raw data
save(sce, file = "Kang18_8vs8.Rda")
unlink(geo, recursive = TRUE)