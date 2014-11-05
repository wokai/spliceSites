\name{readMergedBamGaps}
\alias{readMergedBamGaps}
\title{Reads an object of type \code{gapSites} using a vector of BAM file names.}
\description{The function takes a vector of BAM-file names and corresponding BAM-index file names. For each given filename, the BAM-file will be opened. The functions uses the \code{bamGapList} function (rbamtools) to obtain a \code{data.frame} from an \code{bamReader}. Values for 'gptm' and 'rpmg' are added. Both are rounded to the number of given digits. The function tests for open connection to BAM-file and for initialized index.}
\usage{
readMergedBamGaps(infiles,idxInfiles=paste(infiles,".bai",sep=""),digits=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infiles}{character. Name of BAM-files to be opened.}
  \item{idxInfiles}{character. Name of corresponding BAM-index files. Default: paste(infiles,".bai",sep="")}
  \item{digits}{numeric. \code{gptm} and \code{rpmg} values will be rounded to the number of decimal places given.}
}
\value{gapSites}
\author{
Wolfgang Kaisers
}
\examples{
bam<-character(2)
bam[1]<-system.file("extdata","rna_fem.bam",package="spliceSites")
bam[2]<-system.file("extdata","rna_mal.bam",package="spliceSites")
mbg<-readMergedBamGaps(bam)
}
\keyword{gapSites}
\keyword{readMergedBamGaps}
