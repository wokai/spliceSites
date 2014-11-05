\name{getGapSites}
\alias{getGapSites}
\alias{alignGapList}
\title{Read \code{gapSites}}
\description{
\code{getGapSites} and \code{alignGapList} read gap-site data from single BAM-files (given as \code{bamReader}) and return a \code{gapSites} object. \code{getGapSites} reads data for one seqid (given as 1-based numeric value). \code{alignGapList} reads the whole BAM-file. The functions test for opened reader and initialized index.
}
\usage{
getGapSites(reader,seqid,startid=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reader}{bamReader (rbamtools). An opened instance of \code{bamReader} with initialized index.}
  \item{seqid}{Numeric. 1-based index of reference sequence for which gap-sites are to be read.}
  \item{startid}{Numeric. Default: 1. Determines start value for \code{id} column from which the values are ascending ennumerated. \code{startid} greater than 1 allow to produce unique values over multiple BAM-files.}
}
\value{gapSites}
\details{\code{getGapSites}  internally calls \code{rbamtools::gapList}. \code{alignGapList} internally calls \code{rbamtools::bamGapList}. 'nProbes' values are set to 1.
}
\author{Wolfgang Kaisers}
\examples{
bam<-system.file("extdata", "rna_fem.bam", package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
gal<-getGapSites(reader,1,startid=10)
gal
gal<-alignGapList(reader)
gal
}
\keyword{gapSites}
\keyword{alignGapList}
\keyword{bamReader}
