\name{trypsinCleave}
\alias{trypsinCleave}
\alias{trypsinCleave-methods}
\alias{trypsinCleave,caRanges-method}
\alias{trypsinCleave,aaGapSites-method}
\title{
trypsinCleave method
}
\description{
trypsinCleaves amino acid sequences and returns the fragment which contains the position described by position entry in data.frame.
}
\usage{
trypsinCleave(x, minLen=5, ...)
}
\arguments{
  \item{x}{caRanges (aaGapSites). Object in which amino-acid sequences are to be truncated.}
  \item{minLen}{Numeric. Default is 5. Data sets where the remaining sequence fragment is shorter than \code{minLen} are excluded. }
  \item{...}{Additional arguments which may be passed to the routine (currently unused).}
}
\details{The routine implements the "Keil"-rule, where sites are described by the regex  "[RK](?!P)". The cut position is between [RK] and the following character.  The sequence fragment which contains the exon-intron boundary (depicted by \code{position}) is returned. Dependent numeric values are recalculated. The returned sequence ends on "[RK]" unless the returned fragment is a sequence suffix.}
\value{Same class as given object.}
\author{
Wolfgang Kaisers
}
\examples{
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
ga
lj<-lJunc(ga,featlen=21,gaplen=21,strand='+')
ljc<-lCodons(lj,frame=1,keepStrand=TRUE)
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)
cdr<-dnaRanges(ljc,dna_small)
ar<-translate(cdr)
tra<-truncateSeq(ar)
tyc<-trypsinCleave(tra)
}
\keyword{trypsinCleave}
