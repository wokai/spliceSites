\name{addGeneAlignPart}
\alias{addGeneAlignPart}
\alias{addGeneAlignPart-method}
\alias{addGeneAlignPart,gapSites-method}
\title{Reads a bamRange object for a given \code{bamReader}, \code{refGenome} and gene name.}
\description{Locates gene in genome via \code{refGenome} and reads a \code{bamRange} from the determined region.}
\usage{
addGeneAlignPart(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{gapSites. The result contains a copy of the passed object. }
}
\details{
The function adds a \code{align_part} column to contained \code{data.frame}. From the contained align-gaps the ones with existing (left or right) gene_id and equal left_gene_id and right_gene_id are extracted. For each gene_id, the gene-specific fraction of aligns is added. The value is intended to serve as limit for exclusion of relative underrepresented align-gaps. 
}
\value{gapSites}
\author{
Wolfgang Kaisers
}
\examples{
# A) Read gapSites
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
# B) Annotate
ucf<-system.file("extdata","uc_small.RData",package="spliceSites")
uc<-loadGenome(ucf)
annotation(ga)<-annotate(ga,uc)
# C) align part
gal<-addGeneAlignPart(ga)
gal
}
\keyword{addGeneAlignPart}
\keyword{gapAligns}
