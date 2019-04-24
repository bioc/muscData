ids <- c("Kang18_8vs8", "Sala19_2vs2")

Kang18_8vs8 <- data.frame(
    stringsAsFactors = FALSE,
    Title = "Kang18_8vs8",
    Description = paste(
    "Droplet-based scRNA-seq PBMC data from 8 Lupus patients",
    "before and after 6h-treatment with INF-beta.",
    "Represented as a SingleCellExperiment;",
    "derived from GEO accession number GSE96583."),
    BiocVersion = "3.9",
    Genome = NA,
    SourceType = "tar.gz",
    SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583",
    SourceVersion = "Mar 27 2019",
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based  = NA,
    DataProvider = "GEO",
    Maintainer = "Helena L. Crowell <helena.crowell@uzh.ch>",
    RDataClass = "SingleCellExperiment",
    DispatchClass = "Rda",
    RDataPath = "muscData/Kang18_8vs8.Rda")

Sala19_2vs2 <- data.frame(
    stringsAsFactors = FALSE,
    Title = "Sala19_2vs2",
    Description = paste(
        "10x scRNA-seq data from chimeric mice embryos", 
        "in which Tal1-knockout (KO) embryonic stem cells", 
        "were injected into a wild-type (WT) blastocyst.", 
        "Represented as a SingleCellExperiment; derived from data",
        "made available by the Cancer Research UK Cambridge Institute."),
    BiocVersion = "3.9",
    Genome = NA,
    SourceType = "tar.gz",
    SourceUrl = "https://content.cruk.cam.ac.uk/jmlab/chimera_tal1_data",
    SourceVersion = "Feb 02 2018",
    Species = "Mus musculus",
    TaxonomyId = "10090",
    Coordinate_1_based  = NA,
    DataProvider = "CRUK",
    Maintainer = "Helena L. Crowell <helena.crowell@uzh.ch>",
    RDataClass = "SingleCellExperiment",
    DispatchClass = "Rda",
    RDataPath = "muscData/Sala19_2vs2.Rda")

# write to .csv
df <- do.call("rbind", lapply(ids, get))
write.csv(df, file = "../extdata/metadata.csv", row.names = FALSE)
