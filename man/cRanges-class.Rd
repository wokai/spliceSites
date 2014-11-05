\name{cRanges-class}
\Rdversion{1.1}
\docType{class}
\alias{cRanges-class}
\alias{as.data.frame.cRanges}
%\alias{[,cRanges-method}
\alias{c,cRanges-method}
\alias{count}
\alias{count-methods}
\alias{count,cRanges-method}
\alias{dim,cRanges-method}
\alias{end,cRanges-method}
\alias{head,cRanges-method}
\alias{id}
\alias{id-methods}
\alias{id,cRanges-method}
\alias{initialize,cRanges-method}
\alias{seqid,cRanges-method}
\alias{show,cRanges-method}
\alias{sortTable}
\alias{sortTable,cRanges-method}
\alias{start,cRanges-method}
\alias{strand}
\alias{strand-methods}
\alias{strand,cRanges-method}
\alias{width,cRanges-method}

\title{Class \code{"cRanges"}: Centered ranges.}
\description{
\code{"cRanges"} Objects that contain centered genomic ranges. The center position marks a prominent position inside the range, generally an exon-intron junction. Position values represent the 0-based position of last exon nucleotide.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cRanges", seqnames, start, end, width, strand, position, id)}.
}
\section{Slots}{
  \describe{
    \item{\code{dt}:}{Object of class \code{"data.frame"}. The data.frame contains the columns id, seqnames, start, end, width, strand and position. Each row contains data for one centered range. }
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "cRanges")}: Returns a copy of the contained data inside a \code{data.frame} object. }
%    \item{[}{\code{signature(x = "cRanges")}: Indexed access to data}
    \item{c}{\code{signature(x = "cRanges")}: Generic combining for \code{cRanges} objects.}
    \item{count}{\code{signature(x = "cRanges")}: Returns the number of contained ranges (number of rows). }
    \item{dim}{\code{signature(x = "cRanges")}: Returns the dim of the contained data.frame. }
    \item{dnaRanges}{\code{signature(x = "cRanges", dnaset="DNAStringSet",useStrand="logical",removeUnknownStrand=logical)}: Takes a \code{cRanges} object and a \code{DNAStringSet} (a reference sequence) and adds the appropriate DNA sequence to the genomic ranges. Returns a \code{cdRanges} object. }
    \item{end}{\code{signature(x = "cRanges")}: Returns \code{end} column of data.}
    \item{head}{\code{signature(x = "cRanges",n="numeric",digits="numeric")}: Returns first n (default: n=6) lines of contained data.frame.}
    \item{id}{\code{signature(x = "cRanges")}: Returns id column from contained data.frame. }
    \item{initialize}{\code{signature(.Object = "cRanges")}: Generic class initialisation method. }
    \item{lCodons}{\code{signature(x = "cRanges",frame="numeric", keepStrand="logical")}: Returns cRanges object which represents ranges truncated to codon size. When 'keepStrand' is set to FALSE, strand is set to '+'. The intention is that appended DNA sequences which then can be translated into amino acids.}
    \item{rCodons}{\code{signature(x = "cRanges",frame="numeric", keepStrand="logical")}: Returns cRanges object which represents ranges truncated to codon size. When 'keepStrand' is set to FALSE, strand is set to '+'. The intention is that appended DNA sequences which then can be translated into amino acids.}
    \item{seqid}{\code{signature(x = "cRanges")}: Returns vector with seqid's. }
    \item{show}{\code{signature(object = "cRanges")}: Generic print function. }
    \item{sortTable}{\code{signature(x="cRanges")}: Sort contained tables by \code{seqid},\code{lend} and \code{rstart}.}
    \item{start}{\code{signature(x = "cRanges")}: Returns start column from contained data.frame. }
    \item{strand}{\code{signature(x = "cRanges")}: Returns strand column from contained data.frame.}
    \item{width}{\code{signature(x = "cRanges")}: Returns width of contained ranges (=end-start+1).}
	 }
}
\author{Wolfgang Kaisers}
\seealso{gapRanges}
\examples{
# A) Create cRanges object from scratch
sq<-factor(c(1,1,2,2,3,3),labels=c("chr1","chr2","chr3"))
st<-c(100,200,100,300,100,400)
en<-c(120,210,110,310,110,410)
pos<-c(2,3,4,5,6,7)
cr<-new("cRanges",seqid=sq,start=st,end=en,position=pos)
cr
seqid(cr)
start(cr)
end(cr)
width(cr)
strand(cr)
id(cr)
lCodons(cr,frame=1,keepStrand=TRUE)
lCodons(cr,frame=1,keepStrand=FALSE)
lCodons(cr,frame=2,keepStrand=TRUE)
rCodons(cr,frame=1,keepStrand=FALSE)
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# B) Intended way to create a cRanges object from BAM data
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
ga
lj<-lJunc(ga,featlen=3,gaplen=6,strand='+')
lj
# C) ...
table(strand(lj))
}
\keyword{classes}
\keyword{cRanges}