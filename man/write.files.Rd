\name{write.files}
\alias{write.files}
\alias{write.files-methods}
\title{
write.files
}
\description{
Writes table data and sequence in separate files.
}
\usage{
write.files(x, path, filename,...)
}
\arguments{
  \item{x}{caRanges or aaGapSites object for which data is written.}
  \item{path}{Path for writing files.}
  \item{filename}{ Basic filename to which suffixes are added.}
  \item{...}{Other arguments passed to "write.table".}
}
\details{There are two files written: A text file with tabulated values from data.frame (separated by ";") and a fasta file which contains the stored dna sequence.}
\value{None.}
\author{
Wolfgang Kaisers
}
\note{The function tries to create directory 'path' when it does not exist.}
\examples{
# A) Read gap-sites from BAM-files
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
# B) Load DNA sequence
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)
# C) Add DNA sequence
lj<-lJunc(ga,featlen=21,gaplen=21,strand='+')
ljc<-lCodons(lj,frame=1,keepStrand=TRUE)
cdr<-dnaRanges(ljc,dna_small)
# D) Translate DNA to amino-acid
ar<-translate(cdr)
# E) Write "ar.csv" and "ar.fa"
# write.files(ar,".","ar")
}
\keyword{write.files}
