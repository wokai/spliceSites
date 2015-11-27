\name{readExpSet}
%
\alias{readExpSet}
%
\title{Reads align number or gptm or rpmg value from all given BAM-files
and all identified align gaps into ExpressionSet.}
%
\description{Opens multiple BAM-files and reads aligns for selected gene
    for each file.
    Number of alignes is counted.}
\usage{readExpSet(bam, idx, val="nAligns", phenoData, expData)}
%
\arguments{
    \item{bam}{Vector of BAM-files}
    %
    \item{idx}{Vector of index files (optional)}
    %
    \item{val}{"gptm", "rpmg" or "nAligns". 
        Value type which is written to ExpressionSet matrix
        (nAligns = read count).}
    %
    \item{phenoData}{\code{AnnotatedDataFrame}. Each BAM-file must correspond
        to one identifier.}
    %
    \item{expData}{\code{MIAME}. Optional. 
        Experiment data which can be added to ExpressionSet}
}
%
\value{ExpressionSet}
%
\author{Wolfgang Kaisers}
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
pd <- new("AnnotatedDataFrame", data=prof, varMetadata=meta)

# C) Read ExpressionSet
es <- readExpSet(bam, phenoData=pd)

# D) Annotate ExpressionSet
ucf <- system.file("extdata", "uc_small_junc.RData", package="spliceSites")
juc <- loadGenome(ucf)
ann <- annotate(es, juc)
phenoData(es) <- ann
}
\keyword{readExpSet}
\keyword{ExpressionSet}
