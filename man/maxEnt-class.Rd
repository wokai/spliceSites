\name{maxEnt-class}
\Rdversion{1.1}
\docType{class}
\alias{maxEnt-class}
\alias{load.maxEnt}
\alias{saveMaxEnt}
\alias{saveMaxEnt-methods}
\alias{saveMaxEnt,maxEnt-method}
\alias{basedir,maxEnt-method}
\alias{basedir<-,maxEnt-method}
\alias{score5}
\alias{score5-methods}
\alias{score5,maxEnt-method}
\alias{score3}
\alias{score3-methods}
\alias{score3,maxEnt-method}
\alias{scoreSeq5}
\alias{scoreSeq5-methods}
\alias{scoreSeq5,maxEnt-method}
\alias{scoreSeq3}
\alias{scoreSeq3-methods}
\alias{scoreSeq3,maxEnt-method}
\title{Class \code{"maxEnt"}}
\description{Provides methods for calculation of Splice-site scores. Both functions (\code{score5} and \code{score3}) are intended to work on the '+' strand. \code{score5} scores the 5' side (i.e. the splice donor, left) and the \code{score3} scores the 3' side (i.e. the splice acceptor, right).}
\section{Creation of \code{maxEnt} objects}{
Objects can be created by \code{load.maxEnt()}.
}
\section{Slots}{
  \describe{
    \item{\code{ev}:}{Object of class \code{"environment"} Contains external score data. }
    \item{\code{basedir}:}{Object of class \code{"character"} Directory from which external data is restored. }
    }
}
\section{Methods}{
  \describe{
    \item{basedir}{\code{signature(x = "maxEnt")}: Returns basedir value.}
    \item{basedir<-}{\code{signature(x = "maxEnt",value="character")}: Sets basedir value.}
    \item{score5}{\code{signature(x = "maxEnt",seq="character",pos="integer")}: Calculates score5 value for seq at given position. }
    \item{scoreSeq5}{\code{signature(x="maxEnt",seq="character",frame="integer")}: Calculates score5 values for a single sequence and a series of positions (frame).}
    \item{score3}{\code{signature(x = "maxEnt",seq="character",pos="integer",which="character")}: Calculates score3 value for seq at given position. Accepted values for which are: "ent","wmm" and "emm".}
    \item{scoreSeq3}{\code{signature(x = "maxEnt",seq="character",frame="integer",which="character")}: Calculates score3 values for a single sequence and a series of positions (frame). Accepted values for which are: "ent","wmm" and "emm".}
}}
\author{Wolfgang Kaisers}
\examples{
mes<-load.maxEnt()
score5(mes,"CCGGGTAAGAA",4)                  # 9.844127
score3(mes,"CTCTACTACTATCTATCTAGATC",pos=20) # 6.706947

# scoreSeq functions
sq5<-scoreSeq5(mes,seq="ACGGTAAGTCAGGTAAGT")
sq3<-scoreSeq3(mes,seq="TTTATTTTTCTCACTTTTAGAGACTTCATTCTTTCTCAAATAGGTT")

}
\keyword{classes}
