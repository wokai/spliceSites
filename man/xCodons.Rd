\name{xCodons}
\alias{lCodons}
\alias{lCodons-methods}
\alias{lCodons,cRanges-method}
\alias{rCodons}
\alias{rCodons-methods}
\alias{rCodons,cRanges-method}
\title{
xCodon methods
}
\description{
The xCodon functions work on \code{cRanges} objects. On the contained ranges the function can have two effects: an upstream frame-shift of 0 to 2 positions and a downstream trim to full codons (i.e. (end-start+1)\%\%3==0). The \code{lCodon} function acts in '+' strand mode (left frame-shift, right truncation) and the \code{rCodon} function acts in '-' strand mode (right frame-shift, left truncation).
}
\usage{lCodons(x,frame=1,keepStrand=TRUE)}
\arguments{
  \item{x}{cRanges. Object in which codon positions are calculated}
  \item{frame}{Numeric. Default is 1. Accepted values are 1,2 or 3. The value causes a frame-shift of size (frame-1).}
  \item{keepStrand}{Logical. Default is TRUE. When FALSE, \code{lCodons} overwrites strand entries by '+' and \code{rCodons} overwrites strand entries by '-'. }
}
\details{The function causes an upstream frameshift and a downstream truncation. \code{lCodon} works with '+'-strand view (left-to-right) and \code{rCodon} works with '-'-strand view (right to left). The underlying rationale is: The \code{cRanges} object contains ranges around exon-intron boundaries. The boundary itself is marked by the position value. The functions calculate genomic ranges which can be supplemented by the reference DNA-sequence which then can readily be translated into amino-acid sequences. The different values for \code{frame} and \code{keepStrand} are used to produce all six putative amino-acid sequences for this exon-intron boundary.  }
\author{
Wolfgang Kaisers
}
\examples{
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)
ucf<-system.file("extdata","uc_small.RData",package="spliceSites")
uc<-loadGenome(ucf)

lj<-lJunc(ga,featlen=21,gaplen=21,strand='+')
ljc<-lCodons(lj,frame=1,keepStrand=TRUE)
rj<-rJunc(ga,featlen=21,gaplen=21,strand='-')
rjc<-rCodons(rj,frame=1,keepStrand=TRUE)
}
\keyword{lCodons}
\keyword{rCodons}
\keyword{stop-codon}