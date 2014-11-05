\name{trim}
\alias{trim_left}
\alias{trim_left-methods}
\alias{trim_left,gapSites-method}
\alias{trim_right}
\alias{trim_right-methods}
\alias{trim_right,gapSites-method}
\alias{resize_left}
\alias{resize_left-methods}
\alias{resize_left,gapSites-method}
\alias{resize_right}
\alias{resize_right-methods}
\alias{resize_right,gapSites-method}
\title{trim and resize methods: trim_left, trim_right, resize_left, resize_right}
\description{
The trim and resize functions change number of nucleotides contained in align-gap features (exonic). Trim functions cut feature sizes down to \code{maxlen}. Resize functions reset all sizes to a fixed value. The functions operate directly on the passed objects. There is no return value.
}
\usage{trim_left(x,maxlen)}
\arguments{
  \item{x}{gapSites. Object from which the lJunc values are calculated.}
  \item{maxlen}{Numeric. Maximum number of nucleotides on feature (exon) side of boundary.}
}
\value{None.}
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
# B) Trim 
trim_left(ga,3)
trim_right(ga,2)
ga
# C) Resize
resize_left(ga,5)
resize_right(ga,6)
ga
}
\keyword{lCodons}
\keyword{rCodons}
\keyword{stop-codon}
