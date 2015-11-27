\name{uniqueJuncAnn}
\alias{uniqueJuncAnn}
\alias{uniqueJuncAnn-methods}
\alias{uniqueJuncAnn,ExpressionSet,refJunctions-method}
\title{
uniqueJuncAnn method for ExpressionSet
}
\description{
uniqueJuncAnns adds annotation data to ExpressionSet and removes not-matching sites.
}
\usage{
uniqueJuncAnn(object,junc,ann=TRUE,...)
}
\arguments{
  \item{object}{\code{ExpressionSet}. Object containing gap-site expression data.}
  \item{junc}{\code{refJunctions}. Object containing splice-junction sites.}
  \item{ann}{\code{logical}. Default: \code{TRUE}. When \code{TRUE} the unannotated sites are removed, otherwise the annotated sites are removed.}
  \item{...}{Unused.}
}
\value{ExpresssionSet}
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
ucj<-getSpliceTable(uc)

# E) Extract unique annotated junction sites.
uja<-uniqueJuncAnn(es,ucj)

}
\keyword{uniqueJuncAnn}
\keyword{stop-codon}
