\name{gapSites}
\alias{gapSites}
\alias{dnaGapSites}
\alias{dnaGapSites-methods}
\alias{dnaGapSites,gapSites-method}
\title{Creating 'gapSites' and 'dnaGapSites' objects.}
\description{\code{gapSites} creates objects of class \code{gapSites} from scratch. \code{dnaGapSites} creates objects of class \code{dnaGapSites} from \code{gapSites} objects.}
\usage{
gapSites(seqid=factor(),lstart=integer(),lend=integer(),
                    rstart=integer(),rend=integer(),gaplen,strand,
                    nr_aligns=1,nAligns=sum(nr_aligns),
                    nAlignGaps=sum(nr_aligns),nProbes=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqid}{Character.Identifies reference sequence.}
  \item{lstart}{Coordinates for start of left range.}
  \item{lend}{Coordinates for end of left range. Usually exon-intron boundary.}
  \item{rstart}{Coordinates for start of right range. Usually exon-intron boundary.}
  \item{rend}{Coordinates for end of right range.}
  \item{gaplen}{Length of enclosed gap. Should equal rstart-lend-1.}
  \item{strand}{+ or - or * (for unknown). Default: '*'.}
  \item{nr_aligns}{Number of gapped aligns which have the same exon-intron boundaries (lend and rstart)}
  \item{nAligns}{Total number of aligns for probeset.}
  \item{nAlignGaps}{Total lnumber of gapped aligns for probeset.}
  \item{nProbes}{ Numeric. Number of probes in which this gapped position is present.}
}
\value{An object of class 'gapSites'.}
\details{
The intended way to create a \code{gapSites} object is to use the \code{alignGapList} function which in turn calls the (rbamtools) \code{bamGapList} function. When a BAM file almoust exclusively contains gapped aligns which sometimes are multiply gapped, possibly the 'nAlignGaps' value is greater than the 'nAligns'. When reading BAM files which contain the complete date of an alignment, usually the 'nAlignGaps' value is about $1/3$ of the 'nAligns' value.
}
\author{
Wolfgang Kaisers
}
\examples{
# A) Construct source data from scratch
seqid<-c("chr1","chr1","chr2","chr2","chr2")
lstart<-c(900, 1900,900 ,900, 1900)
lend  <-c(1000,2000,1000,1000,2000)
rstart<-c(1100,2100,1100,1200,2100)
rend  <-c(1200,2200,1200,1300,2200)
nr_aligns<-c(10,20,30,40,10)

# B) Construct gapSites object
ga<-gapSites(seqid,lstart,lend,rstart,rend,nr_aligns=nr_aligns)
ga

# C) Use gapSites accessors
seqid(ga)
lend(ga)
rstart(ga)
strand(ga)
gptm(ga)
rpmg(ga)
nAligns(ga)
nAlignGaps(ga)

# D) Create 
bam<-system.file("extdata", "rna_fem.bam", package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
ga
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)
dga<-dnaGapSites(ga,dna_small)

}
\keyword{gapSites}
