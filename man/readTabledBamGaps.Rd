\name{readTabledBamGaps}
\alias{readTabledBamGaps}
\title{
readTabledBamGaps function
}
\description{
readTabledBamGaps 
}
\usage{
readTabledBamGaps(infiles,idxInfiles=paste(infiles,".bai",sep=""),prof,rpmg=TRUE)
}
\arguments{
  \item{infiles}{character. Names of BAM-files.}
  \item{idxInfiles}{character. Names of BAM-index files. When given index file is not found, the function attempts to create a BAM-index file whith the depicted name.}
  \item{prof}{data.frame. Contains group affiliations for each BAM-file. Each column describes an entity by which values are grouped. The row-number in \code{prof} must be equal to the number of given BAM-files. The order of BAM infiles and prof defines the group classification for each BAM file. All \code{prof} columns must be factors.}
  \item{rpmg}{logical. When TRUE, there will be group specific rpmg align-rates be added to the result table}
}
\value{gapSites}
\details{The function reads gap-align data from all given BAM-files. For each factor level, the number of probes and aligns are counted. When \code{gptm=TRUE} also the gptm values are written for each group. The result table contains for each prof factor level 2 (or 3) extra columns.}
\author{
Wolfgang Kaisers
}
\examples{
bam<-character(2)
bam[1]<-system.file("extdata","rna_fem.bam",package="spliceSites")
bam[2]<-system.file("extdata","rna_mal.bam",package="spliceSites")
prof<-data.frame(gender=c("f","m"))
rtbg<-readTabledBamGaps(bam,prof=prof,rpmg=TRUE)
rtbg
getProfile(rtbg)
}
\keyword{readTabledBamGaps}
\keyword{keyProfiler}
