\name{gapSites-class}
\Rdversion{1.1}
\docType{class}
\alias{gapSites-class}
\alias{as.data.frame.gapSites}
%\alias{[,gapSites-method}
\alias{c,gapSites-method}
\alias{dim,gapSites-method}
\alias{do_group_align_data}
\alias{getAnnStrand}
\alias{getAnnStrand-methods}
\alias{getAnnStrand,gapSites-method}
\alias{getProfile}
\alias{getProfile-methods}
\alias{getProfile,gapSites-method}
\alias{gptm}
\alias{gptm-methods}
\alias{gptm,gapSites-method}
\alias{head,gapSites-method}
\alias{lend}
\alias{lend-methods}
\alias{lend,gapSites-method}
\alias{lstart}
\alias{lstart-methods}
\alias{lstart,gapSites-method}
\alias{merge.gapSites}
\alias{nAligns}
\alias{nAligns-methods}
\alias{nAligns,gapSites-method}
\alias{nAlignGaps}
\alias{nAlignGaps-methods}
\alias{nAlignGaps,gapSites-method}
\alias{overlap_genome}
\alias{rend}
\alias{rend-methods}
\alias{rend,gapSites-method}
\alias{rstart}
\alias{rstart-methods}
\alias{rstart,gapSites-method}
\alias{seqid}
\alias{seqid-methods}
\alias{seqid,gapSites-method}
\alias{rpmg}
\alias{rpmg-methods}
\alias{rpmg,gapSites-method}
\alias{show,gapSites-method}
\alias{sortTable,gapSites-method}
\alias{seqid,gapSites-method}
\alias{strand,gapSites-method}
\alias{strand<-,gapSites-method}
\alias{write.annDNA.tables}
\alias{write.annDNA.tables-methods}
\alias{write.annDNA.tables,gapSites-method}

\title{Class \code{"gapSites"}}
\description{
Contains data from align gaps. \code{"getalignGaps(reader,seqid)"} reads gapped alignments for the specified seqid from a BAM file (via rbamtools) into an object of class \code{"gapSites"}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{alignGapList(reader)}.
}
\section{Slots}{
  \describe{
    \item{\code{nAligns}:}{Object of class \code{"numeric"} Total number of aligns in alignment. }
    \item{\code{nAlignGaps}:}{Object of class \code{"numeric"} Total number of gapped aligns in alignment. }
    \item{\code{dt}:}{Object of class \code{"data.frame"} Table containing basic data for object.}
    \item{\code{lann}:}{Object of class \code{"dataFrameOrNULL"} Optional data.frame containing annotation overlap data for left side.}
    \item{\code{rann}:}{Object of class \code{"dataFrameOrNULL"} Optional data.frame containing annotation overlap data for right side.}
    \item{\code{profile}:}{\code{dataFrameOrNULL} Optional. Contains probe information (Name of BAM-file, group affiliation, number of sites).}
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "gapSites")}: Returns copy of contained data.frame. }
    \item{c}{\code{signature(x = "gapSites")}: Specialisation of generic combine function. }
    \item{dim}{\code{signature(x = "gapSites")}: Specialisation of generic dim function. }
    \item{dnaGapSites}{\code{signature(x = "gapSites")}: Create \code{dnaGapSites} object by adding DNA sequences. }
    \item{getAnnStrand}{\code{signature(x)}: Return strand vector based on annotation content.}
    \item{getProfile}{\code{signature(x)}: Return profile table (data.frame) which contains BAM-file names, group affiliation and number of Sites.}
    \item{gptm}{\code{signature(x = "gapSites")}: Reading accessor for gptm values. }
    \item{head}{\code{signature(x = "gapSites")}: Specialisation of generic head function. }
    \item{lrCodons}{\code{signature(x = "gapSites")}: }
    \item{lJunc}{\code{signature(x = "gapSites",featlen="numeric",gaplen="numeric",keepStrand="logical" (FALSE),unique="logical" (FALSE))}: featlen: Number nucleotides of feature (=exon). gaplen: Number of nucleotides of gap (=intron). keepStrand: Values for strand are copied from argument, otherwise all positions are marked as "+". unique: Multiple identical positions (arising from alternative splice sites on the right side) are collapsed to one line (number of sites is counted in "mult").  Position: 0-based position of last exon nucleotide in DNA sequence.}
    \item{lrJunc}{\code{signature(x = "gapSites")}: ... }
    \item{addGeneAlignPart}{\code{signature(object="gapSites")}: Addes genewise partition of nAligns numbers. Annotation tables must be present. Otherwise an error occurs. }
    \item{merge}{\code{signature(x = "gapSites", y = "ANY")}: Specialisation of generic merge (data.frame) function. }
    \item{nAligns}{\code{signature(object = "gapSites")}: Reading accessor for nAligns value.}
    \item{nAlignGaps}{\code{signature(object = "gapSites")}: Reading accessor for nAlignGaps value. }
    \item{rpmg}{\code{signature(x = "gapSites")}: Reading accessor for rpmg values. }
    \item{show}{\code{signature(object = "gapSites")}: Specialisation of generic show function. }
    \item{sortTable}{\code{signature(x="gapSites")}: Sorts all contained tables by seqid, lend and rstart.}
    \item{write.annDNA.tables}{\code{signature(x="gapSites")}: dnaset,filename,featlen=3,gaplen=8,sep=";",dec=".",row.names=FALSE. Writes csv file with gap-positions, annotations and dna-sequence. }
	 }
}
\author{Wolfgang Kaisers}
\seealso{dnaGapSites}
\examples{
bam<-character(2)
bam[1]<-system.file("extdata","rna_fem.bam",package="spliceSites")
bam[2]<-system.file("extdata","rna_mal.bam",package="spliceSites")
reader<-bamReader(bam[1],idx=TRUE)
agl<-alignGapList(reader)
agl
bamClose(reader)

mbs<-readMergedBamGaps(bam)
mbs
getProfile(mbs)
}
\keyword{classes}
\keyword{gapSites}
