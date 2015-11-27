\name{plotGeneAlignDepth}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{plotGeneAlignDepth}
\alias{plotGeneAlignDepth-methods}
\alias{plotGeneAlignDepth,bamReader-method}
%su
\title{
plotGeneAlignDepth: Plots of read alignment depth for genetic regions
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
The function takes a bamReader and a refGenome object together with 
a gene name and an optional transcript name and plots the read alignment
depth for the regtion of the gene in the opened BAM file. When transcript
data is present, the exonic ranges are added as rectangles on a 
chromosomal line.
}
%
\usage{
plotGeneAlignDepth(reader, genome, 
            gene=NULL, transcript=NULL,
            log="y", cex.main=2,
            col="grey50", fill="grey90", grid=TRUE, 
            box.col="grey20", box.border="grey80")
}
%
\arguments{
    \item{reader}{bamReader (rbamtools). Must be opened and 
        have initialized index.}
    \item{genome}{refGenome. Object which contains genomic annotatin data.}
    \item{gene}{character. Name of one single gene. }
    \item{transcript}{character (optional). Name of one single transcript.}
    \item{log}{character. Name of one single gene. }
    \item{cex.main}{numeric. Determines size of main title. }
    \item{col}{color. A color for align depth line. }
    \item{fill}{color. A color for the interior of align depth area. }
    \item{grid}{logical. When TRUE, a grid is drawn. }
    \item{box.col}{color. A color for the interior of exon rectangles. }
    \item{box.border}{color. A color for the border of exon rectangles. }
}
%
\details{The function checks for opened bamReader and initialized index.
When transcript name is given, the function will plot the positions
of the transcript beneath the alignment depth.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Example
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# Open bamReader
bam <- system.file("extdata","rna_fem.bam",package="spliceSites")
reader <- bamReader(bam, idx=TRUE)
# Load annotation data
ucf <- system.file("extdata", "uc_small.RData", package="spliceSites")
uc <- loadGenome(ucf)
plotGeneAlignDepth(reader, uc, gene="WASH7P", transcript="uc001aac.4")
}
\keyword{plotGeneAlignDepth}
