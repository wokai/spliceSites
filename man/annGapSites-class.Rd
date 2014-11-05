\name{annGapSites-class}
\Rdversion{1.1}
\docType{class}
\alias{annGapSites-class}
\alias{plot_diff}
\alias{plot_diff-methods}
\alias{plot_diff,annGapSites-method}
\title{Class \code{"annGapSites"}}
\description{
Contains data from align gaps together with annotation data (and optional data about alternative splice positions). Objects of this clasa are returned from the \code{annotation} member function for class \code{gapSites}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{annotation} on \code{gapSites} objects.
}
\section{Slots}{
  \describe{
    \item{\code{nAligns}:}{Object of class \code{"numeric"} Total number of aligns. }
    \item{\code{nAlignGaps}:}{Object of class \code{"numeric"} Total number of gapped aligns. }
    \item{\code{dt}:}{\code{"data.frame"}. Contains gap-positions, annotation data and optional alternative position data. }
    \item{\code{lann}:}{\code{"data.frame"}. Contains annotation data for left side. }
    \item{\code{rann}:}{\code{"data.frame"}. Contains annotation data for right side. }
    \item{\code{profile}:}{\code{"data.frame"}. Contains descriptive data for source probes (BAM-files).}
  }
}
\section{Extends}{
Class \code{"\linkS4class{gapSites}"}, directly.
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "annGapSites")}: Returns the contained data. }
	 }
}
\details{
\code{plot_diff} plots tabled distance between inner gap-site border and annotated exon-intron boundaries.
}
\author{
Wolfgang Kaisers
}
\examples{
# A) Read gapSites from BAM
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
# B) Load annotation data
ucf<-system.file("extdata","uc_small.RData",package="spliceSites")
uc<-loadGenome(ucf)
# C) Add Annotation
annotation(ga)<-annotate(ga,uc)
# D) Retrieve annotation
aga<-annotation(ga)
aga
# D) plot_diff
aga<-annotation(ga)
plot_diff(aga)
}
\keyword{annGapSites}
