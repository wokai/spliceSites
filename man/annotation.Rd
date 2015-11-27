\name{annotation}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{annotate}
\alias{annotate-methods}
\alias{annotate,gapSites,refJunctions-method}
%
\alias{annotation}
\alias{annotation-methods}
\alias{annotation,gapSites-method}
%
\alias{annotation<-}
\alias{annotation<--methods}
\alias{annotation<-,gapSites,ANY-method}
\alias{annotation<-,gapSites,refJunctions-method}
%
\title{Annotation functions for \code{gapSites} objects}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Desription
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The \code{annotate} function takes a \code{gapSites}
    and a \code{refGenome} object and returns a list
    which additionally contains a 'class' attribute 'annotationResult'.
    The object is intended as input for the \code{annotation} member
    function of class \code{gapSites}.
    The \code{annotation} member functions act as writing and reading accessor
    for annotation data inside \code{gapSites} objects.}
%
\usage{annotate(object, junc)}
%
\arguments{
    \item{object}{[gapSites]. Align-gap data for which annotations are 
        provided via overlap.}
    %
    \item{junc}{[refJunctions]. Object which provides annotated splice site
        positions.}
}
\value{annAlignGaps}
\details{The \code{annotation} reading accessor takes a \code{gapSites}
    object and returns a \code{annAlignGaps} object.
    The \code{annotation} writing accessor takes a \code{gapSites}
    and a \code{annotationResult} object and copies the contained table
    into the \code{annotation} slot of the \code{gapSites} object.
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# A) Create gapSites object
bam <- system.file("extdata", "rna_fem.bam", package="spliceSites")
reader <- bamReader(bam[1], idx=TRUE)
ga <- alignGapList(reader)
bamClose(reader)

# B) Read refGenome object
ucf <- system.file("extdata", "uc_small_junc.RData", package="spliceSites")
ucj <- loadGenome(ucf)

# C) Add annotation data
annotation(ga) <- annotate(ga, ucj)
}
\keyword{gapSites}
\keyword{annotation}
\keyword{annotate}
\keyword{refGenome}
