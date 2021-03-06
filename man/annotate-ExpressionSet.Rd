\name{annotate-ExpressionSet}
%
\alias{annotate,ExpressionSet,refJunctions-method}
%
\title{Adds annotation data to existing ExpressionSet
(created by \code{readExpSet})}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{Reads \code{featureData} from incoming Expression set which 
    should contain range data on embedding exons for gap-sites. 
    The annotate function then overlaps the ranges with given annotation data.
    The result of overlapping is written into a AnnotatedDataFrame.}
\arguments{
    \item{object}{ExpressionSet}
    \item{genome}{refGenome}
}
\value{AnnotatedDataFrame}
\author{
Wolfgang Kaisers
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# A) Names of BAM-files
bam <- character(2)
bam[1] <- system.file("extdata", "rna_fem.bam", package="spliceSites")
bam[2] <- system.file("extdata", "rna_mal.bam", package="spliceSites")

# B) Experiment Profile
prof <- data.frame(gender=c("f", "m"))
meta <- data.frame(labelDescription=names(prof), row.names=names(prof))
pd<-new("AnnotatedDataFrame", data=prof, varMetadata=meta)

# C) Read ExpressionSet
es <- readExpSet(bam, phenoData=pd)

# D) Annotate ExpressionSet
ucf <- system.file("extdata", "uc_small.RData", package="spliceSites")
uc <- loadGenome(ucf)
juc <- getSpliceTable(uc)
ann <- annotate(es, juc)
}
\keyword{annotate}
\keyword{ExpressionSet}
