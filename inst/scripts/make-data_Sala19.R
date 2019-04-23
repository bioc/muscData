# description: 
#   10x scRNA-seq data from E8.5 mouse chimeric embryos in which Tal1-knockout
#   (KO) embryonic stem cells were injected into a wild-type (WT) blastocyst.
#
# availability: 
#   The original data is made available by the 
#   Cancer Research UK Cambridge Institute (see link below).
#
# reference: 
#   Pijuan-Sala et al., 2019: "A single-cell molecular map of mouse 
#   gastrulation and early organogenesis". Nature (566), 490–495
#
# link to reference: https://www.ncbi.nlm.nih.gov/pubmed/30787436
# link to raw data:  https://content.cruk.cam.ac.uk/jmlab/chimera_tal1_data
#
# Helena L. Crowell; last modified: April 7th, 2019
# ------------------------------------------------------------------------------

# create temporary directory
tmp_dir <- "tmp_dir"
dir.create(file.path(tmp_dir))

# download data
url <- "https://content.cruk.cam.ac.uk"
dir <- file.path(url, "jmlab", "chimera_tal1_data")
download.file(
    url = file.path(dir, "raw_counts.mtx.gz"), 
    destfile = file.path(tmp_dir, "raw_counts"))
download.file(
    url = file.path(dir, "meta.tab.gz"), 
    destfile = file.path(tmp_dir, "cell_md"))
download.file(
    url = file.path(dir, "genes.tsv.gz"), 
    destfile = file.path(tmp_dir, "gene_md"))

# load cell & gene metadata
cell_md <- read.delim(file.path(tmp_dir, "cell_md"), 
    stringsAsFactors = FALSE)
gene_md <- read.delim(file.path(tmp_dir, "gene_md"), 
    header = FALSE, stringsAsFactors = FALSE,
    col.names = c("ENSEMBL", "SYMBOL"))

# load raw counts
counts <- readMM(file.path(tmp_dir, "raw_counts"))
counts <- as(counts, "dgCMatrix") 
rownames(counts) <- with(gene_md, uniquifyFeatureNames(ENSEMBL, SYMBOL))

# construct SCE
sce <- SingleCellExperiment(
    assays = list(counts = counts),
    rowData = gene_md, colData = cell_md)

# save object & remove raw data
save(sce, file = "Sala19_2vs2.Rda")
unlink("tmp_dir", recursive = TRUE)
