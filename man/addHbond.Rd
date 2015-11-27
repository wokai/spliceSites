\name{addHbond}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{addHbond}
\alias{addHbond-methods}
\alias{addHbond,gapSites-method}
\alias{addHbond,cdRanges-method}
%
\title{Class \code{"hbond"}: Provides data and functions for calculation of
    HBond scores for 5' splice-sites.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{The \code{addHbond} methods add HBond scores to \code{gapSites} 
and \code{cdRanges} objects.
HBond scores provide a measure for the capability of a 5' splice-site to
form H-bonds with the U1 snRNA. 
The function requires at least 3 exon nucleotides and 8 intron nucleotides.
The first two intron nucleotides are expected to be 'GT' 
(for other values the returned score will be 0).
The routine equally accepts upper and lower case characters.}
\usage{
addHbond(x,dna)
}
\arguments{
    \item{x}{\code{gapSites}. The object to which HBond scores are added.}
    %
    \item{dna}{\code{DNAStringSet}. Reference sequence identifier.}
}
\references{
http://www.uni-duesseldorf.de/rna/html/hbond_score.php
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{In \code{cdRanges} objects, the function adds a \code{hbond} column.
In \code{gapSites} objects, the function adds a \code{lhbond} (left side) and
a \code{rhbond} (right side) column.
The \code{lhbond} values always assume '+'-strand (because HBond works on
the 5' side).
The \code{rhbond} values always assume '-'-strand.
Therefore, there will be discrepancies in the output of
\code{write.annDNA.tables} because the \code{leftseq} and \code{rightseq}
sequences are reverse-complemented according to the \code{strand} column: 
The xhbond may be > 0 without GT at position 4 (but with AC at position 7).}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# A) Read gapSites
bam<-system.file("extdata","rna_fem.bam",package="spliceSites")
reader<-bamReader(bam,idx=TRUE)
ga<-alignGapList(reader)
bamClose(reader)

# B) Load DNA
dnafile<-system.file("extdata","dna_small.RData",package="spliceSites")
load(dnafile)

# C) HBond
gab<-addHbond(ga,dna_small)

# D) cdRanges
lj<-lJunc(ga,featlen=3,gaplen=8,strand='+')
ljd<-dnaRanges(lj,dna_small)
ljdh<-addHbond(ljd)
}
\keyword{hbond}
