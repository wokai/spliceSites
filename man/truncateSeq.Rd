\name{truncateSeq}
\alias{truncateSeq}
\alias{truncateSeq-methods}
\alias{truncateSeq,caRanges-method}
\alias{truncateSeq,aaGapSites-method}
\title{
truncateSeq method
}
\description{
truncateSeqs amino acid sequences at positions depicted by '*' (stop-codon).
}
\usage{

truncateSeq(x,rme=TRUE,trunc=42L)
}
\arguments{
  \item{x}{caRanges. Object in which amino-acid sequences are to be truncated.}
  \item{rme}{Logical. Default is TRUE. When TRUE, sites with resulting empty sequence (i.e. stop-codon upstream of the splice position) are removed from dataset.}
  \item{trunc}{Integer. ASCII code for character at which truncation should occur. Default value is 42='*' (stop-codon). }
}
\details{The function truncateSeqs the contained amino acid sequences. When the stop-codon is found on the left side of position, the function returns an empty sequence for that site. The position values for these records are also set to 0.}
\value{Object of same class as input.}
\author{
Wolfgang Kaisers
}
\examples{
# A) Read gap-sites from BAM-file
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
# B) Load DNA sequence
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)
ucf<-system.file("extdata","uc_small.RData",package="spliceSites")
uc<-loadGenome(ucf)
# C) Calculate codon frame data and add DNA
lj<-lJunc(ga,featlen=21,gaplen=21,strand='+')
ljc<-lCodons(lj,frame=1,keepStrand=TRUE)
cdr<-dnaRanges(ljc,dna_small)
# D) Translate DNA to amino acid and truncate
ar<-translate(cdr)
tra<-truncateSeq(ar)
}
\keyword{truncateSeq}
\keyword{stop-codon}
