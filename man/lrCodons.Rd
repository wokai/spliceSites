\name{lrCodons}
\alias{lrCodons}
\alias{lrCodons-methods}
\alias{lrCodons,gapSites-method}
\title{
lrCodon methods
}
\description{
The \code{lrCodon} function works on \code{gapSites} objects. \code{gapSites} manage data on align-gaps which represent data on RNA splice sites. On the contained ranges the function can have two effects: an upstream frame-shift of 0 to 2 positions and a downstream trim to full codons (i.e. (end-start+1)\%\%3==0). The strand argument controls direction of effects: '+' strand mode means left frame-shift and right truncation. '-' strand mode means right frame-shift and left truncation. 
}
\usage{lrCodons(x,frame=1L,strand="+")}
\arguments{
  \item{x}{gapSites. Object in which codon positions are calculated}
  \item{frame}{Numeric. Default is 1. Accepted values are 1,2 or 3. The value causes a frame-shift of size (frame-1).}
  \item{strand}{Character or numeric. Default is '+' which is equivalent to 0. Any other value will be interpreted as '-' which is equivalent to 1. }
}
\details{The function causes an upstream frameshift and a downstream truncation. \code{gapSites} objects contain data on gap aligns which represent a related pair of exon-intron boundaries. The returned object is of the same class as the input. Supplemented DNA sequence \code{gapSites} objects will omit introns and will represent the 'spliced' DNA around the splice site. \code{lrCodon} function is intended to shift coordinates, so that the resulting DNA-sequence can readily be translated in a putative amino-acid sequence which contains the splice-site.}
\author{
Wolfgang Kaisers
}
\examples{
# A) Create gapSites object
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)

# B) lr-Junctions for '+' -strand
lrj<-lrJunc(ga,lfeatlen=6,rfeatlen=6,strand='+')
lr1<-lrCodons(lrj,frame=1)
lr2<-lrCodons(lrj,frame=2)
lr3<-lrCodons(lrj,frame=3)
lr<-c(lr1,lr2,lr3)

# C) lr-Junctions for '-' -strand
lrj<-lrJunc(ga,lfeatlen=6,rfeatlen=6,strand='-')
lr1<-lrCodons(lrj,frame=1)
lr2<-lrCodons(lrj,frame=2)
lr3<-lrCodons(lrj,frame=3)
lr<-c(lr1,lr2,lr3)

}
\keyword{lrCodons}
