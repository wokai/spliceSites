\name{dnaRanges}
\alias{dnaRanges}
\alias{dnaRanges-methods}
\alias{dnaRanges,cRanges-method}
\title{Reads a bamRange object for a given \code{bamReader}, \code{refGenome} and gene name.}
\description{Locates gene in genome via \code{refGenome} and reads a \code{bamRange} from the determined region.}
\usage{
dnaRanges(x,dnaset,useStrand=TRUE,removeUnknownStrand=TRUE,verbose=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{cRanges. Range-data will be copied from this object.}
  \item{dnaset}{DNAStringSet. Contains the reference sequence from which the DNA-sequence is extracted. }
  \item{useStrand}{logical. When TRUE, sequences for which strand='-' are reverse-complemented.}
  \item{removeUnknownStrand}{logical. When TRUE, sequences for which strand='-' are removed.}
  \item{verbose}{logical. Determines amount of console output during routine runtime.}
  \item{...}{Optional additional arguments (currently unused).}
}
\value{cdRanges}
\author{
Wolfgang Kaisers
}
\examples{
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
lj<-lJunc(ga,featlen=6,gaplen=6,strand='+')
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)
ljd<-dnaRanges(lj,dna_small)
seqlogo(ljd)
}
\keyword{dnaRanges}