\name{extractRange}
\alias{extractRange}
\alias{extractRange-methods}
\alias{extractRange,gapSites-method}
\alias{extractRange,cRanges-method}
\title{extractRange: Extract subset from object where records lie in given range.}
\description{Searches in object for data which lie inside the given range and returns an object of same type containing extracted data. }
\usage{
extractRange(object,seqid,start,end)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{gapSites or cRanges (or derived). Object inside which the data is searched for.}
  \item{seqid}{character. Reference sequence identifier.}
  \item{start}{numeric. Start position of given range. }
  \item{end}{numeric. End position of given range.}
}
\value{Same type as object}
\author{Wolfgang Kaisers}
\examples{
# A) Read gapSites
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
# B) Load refGenome
ucf<-system.file("extdata","uc_small.RData",package="spliceSites")
uc<-loadGenome(ucf)
# C) For gapSites
extractRange(ga,seqid="chr1",start=14000,end=30000)
# D) For cRanges
lj<-lJunc(ga,featlen=3,gaplen=6,strand='+')
extractRange(lj,seqid="chr1",start=14000,end=30000)
}
\keyword{extractRange}
\keyword{gapSites}
\keyword{cRanges}
\keyword{cdRanges}
\keyword{caRanges}
