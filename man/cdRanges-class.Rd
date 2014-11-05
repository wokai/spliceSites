\name{cdRanges-class}
\Rdversion{1.1}
\docType{class}
\alias{cdRanges-class}
\alias{c,cdRanges-method}
\alias{head,cdRanges-method}
\alias{initialize,cdRanges-method}
\alias{getSequence}
\alias{getSequence-methods}
\alias{getSequence,cdRanges-method}
\alias{seqlogo}
\alias{seqlogo-methods}
\alias{seqlogo,cdRanges-method}
\alias{translate}
\alias{translate-methods}
\alias{translate,cdRanges-method}

\title{Class \code{"cdRanges"}}
\description{
\code{"cdRanges"} Objects that contain centered Ranges (exon-intron junctions) and dna-sequences.}
\section{Objects from the Class}{
Objects are usually created from \code{"cRanges"} with the function \code{"dnaRanges"}.
}
\section{Slots}{
  \describe{
    \item{\code{dt}:}{Object of class \code{"data.frame"}. Contains the columns \code{"seqid"},\code{"start"},\code{"end"},\code{"strand"},\code{"position"},\code{"id"},\code{"frame"}. }  
    \item{\code{seq}:}{Object of class \code{"DNAStringSet"}. Contains the dna-sequence of ranges described in dt.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{cRanges}"}, directly.
}
\section{Methods}{
  \describe{
    \item{c}{\code{signature(x = "cdRanges")}: Generic combining for \code{cdRanges} objects.}
    \item{getSequence}{\code{signature(x="cdRanges")}: Returns contained sequence (DNAStringSet).}
    \item{head}{\code{signature(x = "cdRanges")}: Prints first items from object. }
    \item{initialize}{\code{signature(.Object = "cdRanges")}: Create an instance of class using \code{new}.}
    \item{seqlogo}{\code{signature(x = "cdRanges")}: Show a seqlogo of contained sequences }
    \item{translate}{\code{signature(x = "cdRanges")}: Translates dna-sequence into amino-acid-sequence. Returns an object of class \code{"caRanges"}.}
	 }
}
\author{Wolfgang Kaisers}
\seealso{cRanges}
\examples{
# A) Read gapSites object
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
ga
# B) Create cRanges object
lj<-lJunc(ga,featlen=21,gaplen=21,strand='+')
ljc<-lCodons(lj,frame=1,keepStrand=TRUE)
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)
# C) Add DNA sequence
cdr<-dnaRanges(ljc,dna_small)
# D) seqLogo ...
seqlogo(cdr)
}
\keyword{classes}
\keyword{cdRanges}