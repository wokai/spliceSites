\name{annotation}
\alias{annotate}
\alias{annotate-methods}
\alias{annotate,gapSites-method}
\alias{annotation}
\alias{annotation-methods}
\alias{annotation,gapSites-method}
\alias{annotation<-}
\alias{annotation<--methods}
\alias{annotation<-,gapSites,ANY-method}
\title{Annotation functions for \code{gapSites} objects}
\description{The \code{annotate} function takes a \code{gapSites} and a \code{refGenome} object and returns a list which additionally contains a 'class' attribute 'annotationResult'. The object is intended to be solely used as input for the \code{annotation} member function of class \code{gapSites}. The \code{annotation} member functions act as writing and reading accessor for annotation data inside \code{gapSites} objects. }
\usage{
annotate(object,genome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{gapSites. Align-gap data for which annotations are provided via overlap.}
  \item{genome}{refGenome. Genome object which provides genomic positions for annotated features.}
}
\value{gapSites}
\details{
The \code{annotation} reading accessor takes a \code{gapSites} object and returns a \code{annAlignGaps} object. The \code{annotation} writing accessor takes a \code{gapSites} and a \code{annotationResult} object and copies two contained tables (\code{left} and \code{right}) into the annotation slots (\code{lann} and \code{rann}) of the \code{gapSites} object.
}
\author{
Wolfgang Kaisers
}
\examples{
# A) Create gapSites object
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam[1],idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)

# B) Read refGenome object
ucf<-system.file("extdata","uc_small.RData",package="spliceSites")
uc<-loadGenome(ucf)

# C) Add annotation data
annotation(ga)<-annotate(ga,uc)
}
\keyword{gapSites}
\keyword{annotation}
\keyword{annotate}
\keyword{refGenome}
