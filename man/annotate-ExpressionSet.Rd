\name{annotate-ExpressionSet}
\alias{annotate,ExpressionSet-method}
\title{Adds annotation data to existing ExpressionSet (created by \code{readExpSet})}
\description{Reads \code{featureData} from incoming Expression set which should contain range data on embedding exons for gap-sites. The annotate function then overlaps the ranges with given annotation data. The result of overlapping is written into a AnnotatedDataFrame.}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ExpressionSet}
  \item{genome}{refGenome}
}
\value{AnnotatedDataFrame}
\author{
Wolfgang Kaisers
}
\examples{
# A) Names of BAM-files
bam<-character(2)
bam[1]<-system.file("extdata","rna_fem.bam",package="spliceSites")
bam[2]<-system.file("extdata","rna_mal.bam",package="spliceSites")
# B) Experiment Profile
prof<-data.frame(gender=c("f","m"))
meta<-data.frame(labelDescription=names(prof),row.names=names(prof))
pd<-new("AnnotatedDataFrame",data=prof,varMetadata=meta)
# C) Read ExpressionSet
es<-readExpSet(bam,phenoData=pd)
# D) Annotate ExpressionSet
ucf<-system.file("extdata","uc_small.RData",package="spliceSites")
uc<-loadGenome(ucf)
ann<-annotate(es,uc)
}
\keyword{annotate}
\keyword{ExpressionSet}