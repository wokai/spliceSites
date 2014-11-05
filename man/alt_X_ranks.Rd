\name{alt_X_ranks}
\alias{alt_left_ranks}
\alias{alt_left_ranks-methods}
\alias{alt_left_ranks,gapSites-method}
\alias{alt_right_ranks}
\alias{alt_right_ranks-method}
\alias{alt_right_ranks,gapSites-method}
\title{alt_X_ranks functions}
\alias{alt_ranks}
\alias{alt_ranks-methods}
\alias{alt_ranks,gapSites-method}
\alias{plot_diff_ranks}
\alias{plot_diff_ranks-methods}
\alias{plot_diff_ranks,gapSites-method}
\description{alt_X_ranks covers the functions \code{alt_left_ranks} and \code{alt_right_ranks}. Both functions identify alternative splice-sites. \code{alt_left_ranks} finds sites which share the same \code{rstart} value (on the same seqid). \code{alt_right_ranks} finds sites which share the same \code{lend} value (on the same seqid). \code{alt_ranks} combines the results of both functions together wich \code{seqid}, \code{lend} and \code{rstart} values in one table.}
\usage{
alt_left_ranks(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{gapSites. Object for which alternative ranks are calculated}
}
\details{The function \code{alt_left_ranks} groups align-gaps (splice-sites) which share identical \code{rstart} position and have different \code{lend} position. Each Group is assigned a unique \code{alt_id} (integer value beginning from 1). The first column in the returned data.frame is an id-column which facilitates table merging with the source table. The result has the same number of rows as the source and the id-column.
}
\value{data.frame. The table contains the columns nr_alt, alt_id, id, diff_ranks and gap_diff.}
\author{
Wolfgang Kaisers
}
\examples{
# A) Read gapSites
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
# B) alt_ranks
alr<-alt_left_ranks(ga)
ar<-alt_ranks(ga)
}
\keyword{addGeneAlignPart}
\keyword{gapSites}
