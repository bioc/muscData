\name{Kang18_8vs8}
\title{\code{Kang18_8vs8} dataset}
\alias{Kang18}
\alias{Kang18_8vs8}

\description{
\itemize{
\item{10x droplet-based scRNA-seq PBMC data from 8 Lupus 
patients before and after 6h-treatment with INF-beta.}
\item{The dataset contains 35635 features across 29065 cells, 
and no filtering or preprocessing has been applied.}
\item{The original data is deposited in the Gene Expression Ombnibus (GEO) 
under accession number GSE96583.}}}

\usage{Kang18_8vs8(metadata = FALSE)}

\arguments{
\item{metadata}{
    \code{logical} value indicating whether ExperimentHub metadata 
    (describing the overall dataset) should be returned only, or 
    if the whole dataset should be loaded. Defaults to \code{FALSE}.}}

\details{
\describe{
    \item{link to reference:}{
        https://www.ncbi.nlm.nih.gov/pubmed/29227470}
    \item{link to raw data:}{
        https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583}}}

\value{an object of class \code{SingleCellExperiment}.}

\examples{
# load metadata only
Kang18_8vs8(metadata = TRUE)

# load SingleCellExperiment
Kang18_8vs8(metadata = FALSE)
}

\author{Helena L. Crowell \email{helena.crowell@uzh.ch}} 

\references{
Kang et al. (2019).
"Multiplexed droplet single-cell RNA-sequencing using natural genetic variation",
\emph{Nature Biotechnology} \bold{36}, 89-94.}