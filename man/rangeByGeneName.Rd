\name{rangeByGeneName}
\alias{rangeByGeneName}
\title{Reads a bamRange object for a given \code{bamReader}, \code{refGenome} and gene name.}
\description{Locates gene in genome via \code{refGenome} and reads a \code{bamRange} from the determined region.}
\usage{
rangeByGeneName(reader,genome,gene,complex=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reader}{Object of class (rbamtools) \code{bamReader}. The reader must pass \code{isOpen} and \code{index.initialized} test.}
  \item{genome}{Object of class (refgenome) \code{refGenome}.}
  \item{gene}{Single gene name (character)}
  \item{complex}{Logical. Passed to 'bamRange' function. When TRUE, only aligns with nCigar>1 are counted.}
}
\value{bamRange}
\author{Wolfgang Kaisers}
\examples{
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ucf<-system.file("extdata","uc_small.RData",package="spliceSites")
uc<-loadGenome(ucf)
range<-rangeByGeneName(reader,uc,"WASH7P")
size(range)
}
\keyword{rangeByGeneName}
\keyword{refGenome}
\keyword{bamRange}
