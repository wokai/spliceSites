\name{extractByGeneName}
\alias{extractByGeneName}
\alias{extractByGeneName-methods}
\alias{extractByGeneName,cRanges-method}
\alias{extractByGeneName,gapSites-method}
\title{extractByGeneName: Extract subset for sites which lie in range(s) defined by given gene list.}
\description{The function takes a 'cRanges' object (or derived) and searches inside of given 'refGenome' object for gene names. From identified gene-name matches genomic target regions can be defined for wich in turn the contained sites are extracted.}
\usage{
extractByGeneName(object,geneNames,src,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{gapSites or cRanges (or derived). Object inside which the data is searched for.}
  \item{geneNames}{Character. Vector of gene names for which data is to be extracted.}
  \item{src}{refGenome. Contains gene annotation (for conversion of gene-name to genomic coordinates). }
  \item{...}{(currently unused)}
}
\value{Same type as object}
\details{The function internally calls 'extractByGeneName' on 'refGenome'. This function also prints out non matching gene names. On the result, the function calls 'getGenePositions' from which the genomic regions can be extracted. For each gene, data is extracted via 'extractRange' and the resulting objects are then concatenated.}
\author{Wolfgang Kaisers}
\examples{
# A) Read gapSites from BAM
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)
# B) Load DNAStringSet
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)
# C) load refGenome
ucf<-system.file("extdata","uc_small.RData",package="spliceSites")
uc<-loadGenome(ucf)
# D) For cRanges
lj<-lJunc(ga,featlen=6,gaplen=3,strand='+')
ljw<-extractByGeneName(lj,geneNames="WASH7P",src=uc)
# E) For cdRanges
ljc<-lCodons(lj,frame=2)
ljcd<-dnaRanges(ljc,dna_small)
ljcdw<-extractByGeneName(ljcd,geneNames="WASH7P",src=uc)
# F) For caRanges
ljca<-translate(ljcd)
ljcaw<-extractByGeneName(ljca,geneNames="WASH7P",src=uc)
# G) For gapSites
lrj<-lrJunc(ga,lfeatlen=6,rfeatlen=6,strand='+')
lrjw<-extractByGeneName(lrj,geneNames="WASH7P",src=uc)
}
\keyword{extractByGeneName}
\keyword{gapSites}
\keyword{cRanges}
\keyword{cdRanges}
\keyword{caRanges}
