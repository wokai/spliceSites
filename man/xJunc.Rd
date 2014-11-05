\name{xJunc}
\alias{lJunc}
\alias{lJunc-methods}
\alias{lJunc,gapSites-method}
\alias{lrJunc}
\alias{lrJunc-methods}
\alias{lrJunc,gapSites-method}
\alias{rJunc}
\alias{rJunc-methods}
\alias{rJunc,gapSites-method}
\title{xJunc methods: lJunc, rJunc, lrJunc}
\description{
The term 'xJunc' envelopes three functions: lJunc, rJunc and lrJunc. All three functions take a \code{gapSites} object and return ranges which are restricted around align-gap (exon-intron) boundaries. The functions \code{lJunc} and \code{rJunc} return \code{cRanges} objects, the \code{lrJunc} function returns a \code{gapSites} object.
}
\usage{lJunc(x,featlen,gaplen,unique=FALSE,strand,...)}
\arguments{
  \item{x}{gapSites. Object from which the lJunc values are calculated.}
  \item{featlen}{Numeric. Number of nucleotides on feature (exon) side of boundary. }
  \item{gaplen}{Numeric. Number of nucleotides on gap (intron) side of boundary. }
  \item{unique}{Logical. Default is 'FALSE'. When 'TRUE', the function removes duplicate entries which can be due to alternative splice events.}  
  \item{strand}{Character. Mandatory. All strand entries are set to the given value. }
  \item{...}{Optional arguments passed additionally to the function (currently unused).}
}
\value{cRanges}
\details{The functions are intended to provide position information which crosses exon-intron boundaries. Added DNA sequences can be used to produce seqlogos. The functions are intended to be used in advance of \code{xCodons} functions. Later on added AA sequences can be used to search for proteins where intronic sequences are retained.}
\author{
Wolfgang Kaisers
}
\examples{
# A) Create gapSites object
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam[1],idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
ga

# B) Extract junction data
lj<-lJunc(ga,featlen=6,gaplen=6,strand='+')
ljm<-lJunc(ga,featlen=6,gaplen=6,strand='-')
rj<-rJunc(ga,featlen=6,gaplen=6,strand='+')
rjm<-rJunc(ga,featlen=6,gaplen=6,strand='-')
lrj<-lrJunc(ga,lfeatlen=6,rfeatlen=6,strand='+')
lrjm<-lrJunc(ga,lfeatlen=6,rfeatlen=6,strand='-')
}
\keyword{lCodons}
\keyword{rCodons}
\keyword{stop-codon}
