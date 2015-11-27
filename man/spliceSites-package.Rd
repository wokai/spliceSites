\name{spliceSites-package}
%
\alias{spliceSites-package}
\alias{spliceSites}
%
\docType{package}
%
\title{Calculate information on splice-sites from gapped alignments
    in RNA-seq data.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The package defines 'cRanges' the (centered ranges) class
    which represents a genomic range that contains
    a highlighted position (center): 
    This will usually be the boundary between an exon and an intron.
    The second defined type is the class 'gapSites' which represents
    two exonic regions divided by a gap (usually an intron).
    There are subclasses which additionally contain DNA or AA sequences.
}
%
\details{
\tabular{ll}{
Package: \tab spliceSites\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-10-28\cr
License: \tab GPL-2 \cr
Depends: \tab methods,rbamtools,refGenome,Biobase,BiocGenerics,Biostrings,seqLogo\cr
}
%
}
\author{Wolfgang Kaisers
Maintainer: Wolfgang Kaisers <kaisers@med.uni-duesseldorf.de>
}
%
\references{
Yeo G, Burge CB
Maximum entropy modeling of short sequence motifs
with applications to RNA splicing signals.
J Comput Biol 2004; 11(2-3):377-94
\url{http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html}
}
%
\keyword{ package }
\seealso{
 \code{\link[rbamtools:rbamtools-package]{rbamtools}}
 \code{\link[refGenome:refGenome-package]{refGenome}}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
bam <- system.file("extdata", "rna_fem.bam", package="spliceSites")
reader <- bamReader(bam, idx=TRUE)
ga <- alignGapList(reader)
bamClose(reader)
dnafile <- system.file("extdata", "dna_small.RData", package="spliceSites")
load(dnafile)
ucf <- system.file("extdata","uc_small_junc.RData", package="spliceSites")
ucj <- loadGenome(ucf)
annotation(ga) <- annotate(ga, ucj)
ga
}
