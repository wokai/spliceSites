\name{caRanges-class}
\Rdversion{1.1}
\docType{class}
\alias{caRanges-class}
\alias{c,caRanges-method}
\alias{getSequence,caRanges-method}
\alias{head,caRanges-method}
\alias{write.files,caRanges-method}

\title{Class \code{"caRanges"}}
\description{\code{"caRanges"} Objects that contain a centered genomic range and amino acid sequences.}
\section{Objects from the Class}{
Objects are usually created from objects of class \code{"cdRanges"} by the \code{"translate"} function.
}
\section{Slots}{
  \describe{
    \item{\code{dt}:}{Object of class \code{"data.frame"}. Contains the columns \code{"seqid"},\code{"start"},\code{"end"},\code{"strand"},\code{"position"},\code{"id"},\code{"frame"}}  
    \item{\code{seq}:}{Object of class \code{"AAStringSet"}. Contains amino-acid-sequence of ranges described in dt.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{cRanges}"}, directly.
}
\section{Methods}{
  \describe{
    \item{c}{\code{signature(x = "caRanges")}: Generic combining for \code{caRanges} objects.}
    \item{getSequence}{\code{signature(x="caRanges")}: Returns contained sequence (DNAStringSet).}
    \item{head}{\code{signature(x = "aaGapAligns")}: Returns the first lines of object. }
    \item{show}{\code{signature(object = "aaGapAligns")}: Returns the last lines of object. }
    \item{truncateSeq}{\code{signature(x="caRanges",rme=TRUE,trunc=42L)}: Truncates contained sequence when character (given by ASCII code in \code{trunc}). The default (42L) encodes for character '*' which indicates stop-codon. } 
    \item{trypsinCleave}{\code{signature(x = "caRanges",minLen = 5)}: Performs in silico trypsinization of contained sequence. The sequence fragment which contains the (position depicted) exon-intron boundary is returned. Datasets for which the truncated sequence is shorter than minLen are excluded. }
    \item{write.files}{\code{signature(x = "caRanges")}: Exports contained data into "csv" file. }
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
# D) Translate into AA sequence
ar<-translate(cdr)
# E) Truncate and cleave...
tra<-truncateSeq(ar)
tyc<-trypsinCleave(tra)
}
\keyword{classes}
\keyword{caRanges}