\name{seqlogo}
\alias{seqlogo}
\alias{seqlogo-methods}
\alias{seqlogo,cdRanges-method}
\title{seqlogo: Plotting sequence logo for 'cdRanges' and 'dnaGapSites' objects.}
\description{The function produces a sequence logo plot based on the contained sequences.}
\usage{
seqlogo(x,strand="+",useStrand=TRUE,...)
}
\arguments{
  \item{x}{\code{cdRanges} or \code{dnaGapSites} Object.}
  \item{strand}{Character. Determines the subset for which the seqlogo is plotted. This option is only used when useStrand is given as 'TRUE'.}
  \item{useStrand}{Logical. Determines whether the given strand information is used. For \code{useStrand=FALSE} the plot is made up from all contained sequences.}
  \item{...}{(Currently unused)}
}
\details{The function fails with an error message when the dataset does not contain any records with the given strand (except useStrand=FALSE).}
\value{None}
\author{Wolfgang Kaisers}
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
\keyword{seqlogo}