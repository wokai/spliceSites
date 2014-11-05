\name{hbond-class}
\Rdversion{1.1}
\docType{class}
\alias{hbond-class}
\alias{load.hbond}
\alias{basedir,hbond-method}
\alias{basedir<-,hbond-method}
\alias{hbond}
\alias{hbond-methods}
\alias{hbond,hbond-method}
\title{Class \code{"hbond"}}
\description{Provides methods and data for calculation of HBond 5'splice-site scores. HBond scores provide a measure for the capability of a 5' splice-site to form an RNA duplex with the 5' end of the U1 snRNA. The function requires at least 3 exon nucleotides and 8 intron nucleotides. The \code{hbond} function takes a vector DNA sequences and a vector of position (\code{pos}) values. The position values represent the 1-based position of the last exon nucleotide. Therefore all position values must be >=3 and the sequence length must be >= pos+8.}
\section{Creation of \code{hbond} objects}{
Objects can be created by \code{load.hbond()}.
}
\section{Slots}{
  \describe{
    \item{\code{ev}:}{Object of class \code{"environment"} Contains external score data. }
    \item{\code{basedir}:}{Object of class \code{"character"} Directory from which external data is restored. }
    }
}
\section{Methods}{
  \describe{
    \item{basedir}{\code{signature(x = "hbond")}: Returns basedir value.}
    \item{basedir<-}{\code{signature(x = "hbond",value="character")}: Sets basedir value.}
    \item{hbond}{\code{signature(x = "hbond",seq="character",pos="integer")}: Calculates score5 value for seq at given position. }
}}
\details{
The first two intron nucleotides must be 'GT'. Otherwise the returned value is 0. All other sequence characters must be in "ATCG" (capitalization does not matter). When any other character (such as N) is found, the function also returns 0. 
}

\references{
http://www.uni-duesseldorf.de/rna/html/hbond_score.php
}
\author{Wolfgang Kaisers}
\examples{
hb<-load.hbond()
seq<-c("CAGGTGAGTTC","ATGCTGGAGAA","AGGGTGCGGGC","AAGGTAACGTC","AAGGTGAGTTC")
hbond(hb,seq,3)
}
\keyword{classes}
