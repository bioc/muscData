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

Crowell19_4vs4 <- data.frame(
    stringsAsFactors = FALSE,
    Title = "Crowell19_4vs4",
    Description = paste(
        "Single-nuclei RNA-seq data of 8 CD-1 male mice (age 11 weeks)",
        "split into 2 groups with 4 animals each: vehicle and",
        "peripherally lipopolysaccharaide (LPS) treated mice;",
        "derive from Figshare DOI:10.6084/m9.figshare.8976473.v1."),
    BiocVersion = "3.11",
    Genome = NA,
    SourceType = "Zip",
    SourceUrl = "https://doi.org/10.6084/m9.figshare.8976473.v1",
    SourceVersion = "July 22 2019",
    Species = "Mus musculus",
    TaxonomyId = "10090",
    Coordinate_1_based  = NA,
    DataProvider = "F. Hoffmann-La Roche Ltd.",
    Maintainer = "Helena L. Crowell <helena.crowell@uzh.ch>",
    RDataClass = "SingleCellExperiment",
    DispatchClass = "Rda",
    RDataPath = "muscData/Crowell19_4vs4.Rda")

# write to .csv
ids <- c("Kang18_8vs8", "Crowell19_4vs4")
df <- do.call("rbind", lapply(ids, get))
write.csv(df, file = "../extdata/metadata.csv", row.names = FALSE)
