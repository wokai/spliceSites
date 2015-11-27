\name{aaGapSites-class}
\Rdversion{1.1}
\docType{class}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{aaGapSites-class}
\alias{head,aaGapSites-method}
\alias{show,aaGapSites-method}
\alias{write.files,aaGapSites-method}
%
\title{Class \code{"aaGapSites"}}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
Contains gapAligns data and a AAStringSet.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("aaGapSites", ...)}.
}
\section{Slots}{
    \describe{
        \item{\code{seq}:}{\code{"AAStringSet"}:
            Contains amino acid sequences.}
        %
        \item{\code{nAligns}:}{\code{"numeric"}:
            Contains total number of aligns.}
        %
        \item{\code{nAlignGaps}:}{\code{"numeric"}:
            Contains total number of align gaps.}
        %
        \item{\code{dt}:}{\code{"data.frame"}: 
            Contains data for all gap sites.}
    }
}
\section{Extends}{
Class \code{"\linkS4class{gapSites}"}, directly.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
    \describe{
        \item{head}{\code{signature(x = "aaGapSites")}:
            Returns the first lines of object.}
        %
        \item{show}{\code{signature(object = "aaGapSites")}:
            Returns the last lines of object.}
        %
        \item{truncateSeq}{\code{signature(x="caRanges",rme=TRUE,trunc=42L)}: 
            Truncates contained sequence when character (given by ASCII code
            in \code{trunc}). The default (42L) encodes for character '*'
            which indicates stop-codon.}
        %
        \item{trypsinCleave}{\code{signature(x = "caRanges",minLen = 5)}:
            Performs in silico trypsinization of contained sequence.
            The sequence fragment which contains the (position depicted)
            exon-intron boundary is returned. Datasets for which the truncated 
            sequence is shorter than minLen are excluded.}
        %
        \item{write.files}{\code{signature(x = "caRanges")}: 
            Exports contained data into "csv" file.}
}}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
# A) Read gap-sites from BAM-file
bam <- system.file("extdata", "rna_fem.bam", package="spliceSites")
reader <- bamReader(bam, idx=TRUE)
ga <- alignGapList(reader)
bamClose(reader)

# B) Load reference dna
dnafile <- system.file("extdata", "dna_small.RData", package="spliceSites")
load(dnafile)

# C) Calculate cross junctional ranges
lrj <- lrJunc(ga, lfeatlen=6, rfeatlen=6, strand='+')
lr1 <- lrCodons(lrj, frame=1, strand='+')
lr2 <- lrCodons(lrj, frame=2, strand='+')
lr3 <- lrCodons(lrj, frame=3, strand='+')
lr <- c(lr1, lr2, lr3)

# D) Add DNA-sequence
lrd <- dnaGapSites(lr, dna_small)

# E) Translate DNA to amino acid
lra <- translate(lrd)
}
\keyword{classes}
