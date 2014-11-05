\name{addMaxEnt}
\alias{addMaxEnt}
\alias{addMaxEnt-methods}
\alias{addMaxEnt,gapSites-method}
\alias{setMeStrand}
\alias{setMeStrand-methods}
\alias{setMeStrand,gapSites-method}
\alias{getMeStrand}
\alias{getMeStrand-methods}
\alias{getMeStrand,gapSites-method}
\title{addMaxEnt: Extract subset of data contained in given range given object. }
\description{addMaxEnt adds new columns to object data which contain MaxEnt-Score derived values. \code{mxe_ps5} contains \code{score5} values for left align-gap (exon-intron) boundary (i.e. assumed to reside on '+'-strand. \code{mxe_ps3} contains \code{score3} (maxent) values for right align-gap (intron-exon) boundary (i.e. assumed to reside on '+'-strand). \code{mxe_ms5} contains \code{score5} values for right align-gap (exon-intron) boundary on \code{reverseComplement} transformed sequence (i.e. assumed to reside on '-'-strand). \code{mxe_ms3} contains \code{score3} values for left align-gap (intron-exon) boundary on \code{reverseComplement} transformed sequence (i.e. assumed to reside on '-'-strand). From these values, \code{s3strand},  \code{s5strand} and \code{meStrand} are derived: \code{s3strand} is '+' when \code{mxe_ps5} >= \code{mxe_ms5} and '-' otherwise; \code{s3strand} is '+' when \code{mxe_ps3} >= \code{mxe_ms3} and '-' otherwise. \code{meStrand} equals \code{s5strand} when \code{s5strand}=\code{s3strand} and '*' otherwise. The function \code{setMeStrand} copies existing \code{meStrand} values into \code{strand} column (and throws an error when \code{meStrand} does not exist).}
\usage{
addMaxEnt(x,dna,maxent,digits=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{gapSites.}
  \item{dna}{DNAStringSet. Reference sequence identifier.}
  \item{maxent}{maxEnt. Contains score table which are internally used by \code{score3} and \code{score5} methods.}
  \item{digits}{Numeric. Default value: 1. Internally calculated maxent scores are rounded to given number of decimal places.}
}
\value{gapSites}
\author{Wolfgang Kaisers}
\examples{
# A) Read gapSites
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
# B) Load DNA
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)
# C) maxEnt
mes<-load.maxEnt()
gae<-addMaxEnt(ga,dna_small,mes)
getMeStrand(gae)
sae<-setMeStrand(gae)
}
\keyword{addMaxEnt}
\keyword{maxEnt}
