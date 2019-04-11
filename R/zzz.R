.onLoad <- function(libname, pkgname) {
    fn <- system.file("extdata", "metadata.csv", package = "muscData")
    ts <- read.csv(fn, stringsAsFactors = FALSE)$Title
    createHubAccessors(pkgname, ts)
}