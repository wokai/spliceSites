\name{addGeneAligns}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{addGeneAligns}
\alias{addGeneAligns-method}
\alias{addGeneAligns,gapSites-method}
%
\title{Reads a bamRange object for a given \code{bamReader},
        \code{refGenome} and gene name.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{Locates gene in genome via \code{refGenome} and reads a
    \code{bamRange} from the determined region.}
\usage{
    addGeneAligns(x)
}
\arguments{
    \item{x}{gapSites. The result contains a copy of the passed object. }
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{ The function adds a \code{gene_aligns} column 
    to the contained \code{data.frame}.
}
%
\value{gapSites}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# A) Read gapSites
bam <- system.file("extdata", "rna_fem.bam", package="spliceSites")
reader <- bamReader(bam, idx=TRUE)
ga <- alignGapList(reader)
bamClose(reader)

# B) Annotate
ucf <- system.file("extdata", "uc_small_junc.RData", package="spliceSites")
ucj <- loadGenome(ucf)
annotation(ga) <- annotate(ga, ucj)

# C) align part
gal <- addGeneAligns(ga)
gal
}
\keyword{addGeneAligns}
\keyword{gapAligns}
