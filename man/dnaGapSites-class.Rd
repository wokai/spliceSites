\name{dnaGapSites-class}
\Rdversion{1.1}
\docType{class}
\alias{dnaGapSites-class}
\alias{head,dnaGapSites-method}
\alias{seqlogo,dnaGapSites-method}
\alias{show,dnaGapSites-method}
\alias{translate,dnaGapSites-method}

\title{Class \code{"dnaGapSites"}}
\description{dnaGapSites contains all data which is stored in objects of class \code{"gapSites"} plus additional DNA sequences in the \code{"seq"} slot.}
\section{Objects from the Class}{
Objects are usually created from \code{gapSites} via \code{dnaGapSites}.
}
\section{Slots}{
  \describe{
    \item{\code{seq}:}{\code{"DNAStringSet"}. Contains DNA sequence. }
    \item{\code{nAligns}:}{code{"numeric"}. Contains total number of aligns. }
    \item{\code{nAlignGaps}:}{\code{"numeric"}. Contains total number of align gaps. }
    \item{\code{dt}:}{code{"data.frame"}. Contains data on gap-sites. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{gapSites}"}, directly.
}
\section{Methods}{
  \describe{
    \item{head}{\code{signature(x = "dnaGapSites")}: Returns head of dt.}
    \item{seqlogo}{\code{signature(x = "dnaGapSites")}: Prints seq-logo of stored dna-sequence.}
    \item{show}{\code{signature(object = "dnaGapSites")}: Prints head of dt.}
    \item{translate}{\code{signature(x = "dnaGapSites")}: Returns an object of class aaalignGaps by translating seq into amino acids. }
	 }
}
\author{
Wolfgang Kaisers
}

\seealso{gapSites}
\examples{
# A) Read gapSites
bam<-system.file("extdata", "rna_fem.bam", package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
# B) Load DNA sequence
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)
# C 1) Add DNA
dga<-dnaGapSites(ga,dna_small)
dga
# C 2) Calculate codon positions
lrj<-lrJunc(ga,lfeatlen=6,rfeatlen=6,strand='+')
lrc<-lrCodons(lrj,frame=1,strand='+')
# D) Add DNA sequence and translate
lrd<-dnaGapSites(lrc,dna_small)
lra<-translate(lrd)
lra
}
\keyword{classes}
