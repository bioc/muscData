\name{Sala19_2vs2}
\title{\code{Sala19_2vs2} dataset}
\alias{Sala19}
\alias{Sala19_2vs2}

\description{
\itemize{
\item{10x scRNA-seq data from E8.5 mouse chimeric embryos 
in which Tal1-knockout (KO) embryonic stem cells were injected 
into a wild-type (WT) blastocyst.}
\item{The dataset contains 29453  features across 56122 cells, 
and no filtering or preprocessing has been applied.}
\item{The original data is made available by the 
Cancer Research UK Cambridge Institute (see link below).}}}

\usage{Sala19_2vs2(metadata = FALSE)}

\arguments{
\item{metadata}{
    \code{logical} value indicating whether ExperimentHub metadata 
    (describing the overall dataset) should be returned only, or 
    if the whole dataset should be loaded. Defaults to \code{FALSE}.}}

\details{
\describe{
    \item{link to reference:}{
        https://www.ncbi.nlm.nih.gov/pubmed/30787436}
    \item{link to raw data:}{
        https://content.cruk.cam.ac.uk/jmlab/chimera_tal1_data}}}

\value{an object of class \code{SingleCellExperiment}.}

\examples{
# load metadata only
Sala19_2vs2(metadata = TRUE)

# load SingleCellExperiment
Sala19_2vs2(metadata = FALSE)}

\author{Helena L. Crowell \email{helena.crowell@uzh.ch}} 

\references{
Pijuan-Sala et al. (2019).
"A single-cell molecular map of mouse gastrulation and early organogenesis",
\emph{Nature} \bold{566}, 490-495.}