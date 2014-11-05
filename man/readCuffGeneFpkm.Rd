\name{readCuffGeneFpkm}
\alias{readCuffGeneFpkm}
\title{Reads FPKM values into ExpressionSet.}
\description{Opens fpkm_tracking files and collects FPKM values into ExpressionSet. The function is intended to work whith \code{genes.fpkm_tracking} files. In order to get unique gene identifier, the contained values are grouped and for each gene the maximum FPKM values is selected. There should only be a few hundred multiple occurring genes and the maximum value should give a (slight) underestimation of the real value.}
\usage{
readCuffGeneFpkm(cuff,phenoData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cuff}{Vector of cufflinks files}
  \item{phenoData}{AnnotatedDataFrame}
}
\value{ExpressionSet}
\author{
Wolfgang Kaisers
}
\examples{
n<-10
cuff<-system.file("extdata","cuff_files",
  paste(1:n,"genes","fpkm_tracking",sep="."),package="spliceSites")
  
gr<-system.file("extdata","cuff_files","groups.csv",package="spliceSites")
groups<-read.table(gr,sep="\t",header=TRUE)
meta<-data.frame(labelDescription=c("gender","age-group","location"),
    row.names=c("gen","agg","loc"))
phenoData<-new("AnnotatedDataFrame",data=groups,varMetadata=meta)
exset<-readCuffGeneFpkm(cuff,phenoData)
}
\keyword{readCuffGeneFpkm}
\keyword{ExpressionSet}