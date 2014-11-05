\name{xJuncStrand}
\alias{lJuncStrand}
\alias{lJuncStrand-methods}
\alias{lJuncStrand,gapSites-method}
\alias{lrJuncStrand}
\alias{lrJuncStrand-methods}
\alias{lrJuncStrand,gapSites-method}
\alias{rJuncStrand}
\alias{rJuncStrand-methods}
\alias{rJuncStrand,gapSites-method}
\title{xJuncStrand methods: lJuncStrand, rJuncStrand, lrJuncStrand}
\description{
The term 'xJuncStrand' envelopes three functions: lJuncStrand, rJuncStrand and lrJuncStrand. All three functions take a \code{gapSites} object and return ranges which are restricted around align-gap (exon-intron) boundaries. The functions \code{lJuncStrand} and \code{rJuncStrand} return \code{cRanges} objects, the \code{lrJuncStrand} function returns a \code{gapSites} object. The resulting objects contain strand information which is copied from the input objects.
}
\usage{lJuncStrand(x,featlen,gaplen,...)}
\arguments{
  \item{x}{gapSites. Object from which the lJuncStrand values are calculated.}
  \item{featlen}{Numeric. Number of nucleotides on feature (exon) side of boundary. }
  \item{gaplen}{Numeric. Number of nucleotides on gap (intron) side of boundary. }
  \item{...}{Optional arguments passed additionally to the function (currently unused).}
}
\value{cRanges}
\details{The functions are intended to provide position information which crosses exon-intron boundaries. Added DNA sequences can be used to produce seqlogos. The functions are intended to be used in advance of \code{xCodons} functions. Later on added AA sequences can be used to search for proteins where intronic sequences are retained.}
\author{
Wolfgang Kaisers
}
\examples{
# A) Create gapSites object
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam[1],idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
ga

# B) Extract JuncStrandtion data
lj<-lJuncStrand(ga,featlen=6,gaplen=6)
ljm<-lJuncStrand(ga,featlen=6,gaplen=6)
rj<-rJuncStrand(ga,featlen=6,gaplen=6)
rjm<-rJuncStrand(ga,featlen=6,gaplen=6)
lrj<-lrJuncStrand(ga,lfeatlen=6,rfeatlen=6)
lrjm<-lrJuncStrand(ga,lfeatlen=6,rfeatlen=6)
}
\keyword{lCodons}
\keyword{rCodons}
\keyword{stop-codon}
