\name{countByGeneName}
\alias{countByGeneName}
\title{Reads align number for selected gene from multiple BAM-files.}
\description{Opens multiple BAM-files and reads aligns for selected gene for each file. The function counts the tag-selected value which either is a BAM-cigar operation (like "N" or "M") or the total number of aligns.}
\usage{
countByGeneName(object,infiles,idxInfiles=paste(infiles,".bai",sep=""),gene,tag="N")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class "refGenome"}
  \item{infiles}{Vector of BAM-files}
  \item{idxInfiles}{(Optional) Vector of BAM-index files.}
  \item{gene}{Gene name}
  \item{tag}{Character. Passed to (rbamtools) 'bamCountAll' function. Default value is "N". Other accepted values include "nAligns","M","I","D".}
}
\value{Numeric vector. Length equals number of BAM-input files.}
\details{
\code{countByGeneName} first uses the \code{extractByGeneName} and \code{getGenePositions} from 'refGenome' in order to calculate coordinates from the given gene name. Then for each given BAM-file name, the functions calls the \code{bamCount} function and returns a vector with a count value for each given file. Internally \code{countByGeneName} also checks for existing BAM-index file and tries to create index files which do not exist. 
}
\author{
Wolfgang Kaisers
}
\examples{
# A) Read filenames
ucf<-system.file("extdata","uc_small.RData",package="spliceSites")
uc<-loadGenome(ucf)
bam<-character(2)
bam[1]<-system.file("extdata","rna_fem.bam",package="spliceSites")
bam[2]<-system.file("extdata","rna_mal.bam",package="spliceSites")
# B) count
countByGeneName(uc,bam,gene="WASH7P",tag="N")
countByGeneName(uc,bam,gene="WASH7P",tag="nAligns")
}
\keyword{gapSites}
\keyword{bamCount}
