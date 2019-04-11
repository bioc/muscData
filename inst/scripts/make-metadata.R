ids <- c("Kang18_8vs8")

Kang18_8vs8 <- data.frame(
    stringsAsFactors = FALSE,
    Title = "Kang18_8vs8",
    Description = paste(
    "Droplet-based scRNA-seq PBMC data from 8 Lupus patients",
    "before and after 6h-treatment with INF-beta.",
    "Represented as a SingleCellExperiment;",
    "derived from GEO accession number GSE96583."),
    BiocVersion = "3.8",
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

# write to .csv
df <- bind_rows(lapply(ids, get))
write.csv(df, file = "../extdata/metadata.csv", row.names = FALSE)
