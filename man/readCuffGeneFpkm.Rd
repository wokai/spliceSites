\name{readCuffGeneFpkm}
\alias{readCuffGeneFpkm}
\title{Reads FPKM values into ExpressionSet.}
%
\description{Opens fpkm_tracking files and collects FPKM values into
ExpressionSet.
The function is intended to work whith \code{genes.fpkm_tracking} files.
In order to get unique gene identifier, the contained values are grouped
and for each gene the maximum FPKM values is selected.
There should only be a few hundred multiple occurring genes and the maximum
value should give a (slight) underestimation of the real value.}
%
\usage{
readCuffGeneFpkm(cuff, phenoData, summ="max")
}
%
\arguments{
    \item{cuff}{\code{character}: Vector of cufflinks files}
    %
    \item{phenoData}{\code{AnnotatedDataFrame}: 
        Requirement for construction of an \code{ExpressionSet}}
    %
    \item{summ}{\code{character}: Must be either 'max' or 'sum'.
    A handful of tracking id's occur multiple times due to multiple transcripts
    which partially are non-overlapping.
    The summ (summarize) Argument determines the way the multiplets are
    handled.}
}
\value{ExpressionSet}
\author{
Wolfgang Kaisers
}
\examples{
n<-10
cuff <- system.file("extdata", "cuff_files",
    paste(1:n, "genes", "fpkm_tracking", sep="."), package="spliceSites")

## Create Pheno - data
gr <- system.file("extdata", "cuff_files", "groups.csv", package="spliceSites")
groups <- read.table(gr, sep="\t", header=TRUE)
meta <- data.frame(labelDescription=c("gender", "age-group", "location"),
            row.names=c("gen", "agg", "loc"))
phenoData <- new("AnnotatedDataFrame", data=groups, varMetadata=meta)

## Read ExpressionSet
exset <- readCuffGeneFpkm(cuff, phenoData)
}
\keyword{readCuffGeneFpkm}
\keyword{ExpressionSet}
