\name{keyProfiler-class}
\Rdversion{1.1}
\docType{class}
\alias{keyProfiler-class}
\alias{addKeyTable}
\alias{addKeyTable-methods}
\alias{addKeyTable,keyProfiler-method}
\alias{appendKeyTable}
\alias{appendKeyTable-methods}
\alias{appendKeyTable,keyProfiler-method}
\alias{getKeyTable}
\alias{getKeyTable-methods}
\alias{getKeyTable,keyProfiler-method}
\alias{initialize,keyProfiler-method}
\title{Class \code{"keyProfiler"}}
\description{
Internal class that counts occurrence of profile factors (e.g. gender male and female) successively for added key-tables. The columns of the key-tables define the groups (e.g. genomic positions: seqid, start, end) for each all profile factors are counted.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("annAligns", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{ev}:}{Environment: contains the main data of each object. The environment contains the data.frames 'dtb' (key-tabled profiles) and 'prof' (profiles: a table that contains the profile definition for each added key-table) as well as 'groupExpr', an unevaluated Expression which does the data.frame-grouping after addition of a new key-table.}
    \item{\code{unique}:}{Logical: When true, there can be maximal one table added for each indexed profile }
    \item{\code{counted}:}{Logical: Stores the information which profile already has been counted. Is only used when 'unique' is 'TRUE'.}
        \item{\code{useValues}:}{Logical: When TRUE, the object tables the values given together with each key-table, otherwise the profiles are simply counted.}
  }
}
\section{Methods}{
  \describe{
    \item{addKeyTable}{\code{signature(x = "keyProfiler",keyTable="data.frame",index="numeric",values="numeric")}: Adds keyed data to key-table and counts values according to profile (which is defined by index via profile table). }
    \item{getKeyTable}{\code{signature(x = "keyProfiler")}: Returns key-table. }
    \item{appendKeyTable}{\code{signature(x = "keyProfiler",keytable="data.frame",prefix="character",valFactor="numeric",rateFactor="numeric",digits="numeric")}: cbinds internal key-table to keytable-argument. A prefix can be added to column-names. A given valFactor is multiplied with the counted values. A given rateFactor causes counted values to be converted into rates (i.e. divided by column-sums and multiplied with rateFactor value. Values are rounded when a digits argument is provided.) }
  }
}
\author{
Wolfgang Kaisers
}
\examples{
# Loads profile, position data (key) and aggregated values (ku) data.frames
load(system.file("extdata", "key.RData", package="spliceSites"))
# Group positions
kpc<-new("keyProfiler",keyTable=key1[,c("seqid","lend","rstart")],prof=prof)
addKeyTable(kpc,keyTable=key2[,c("seqid","lend","rstart")],
                index=2,values=key2$nAligns)
addKeyTable(kpc,keyTable=key3[,c("seqid","lend","rstart")],
                index=4,values=key3$nAligns)
cp<-appendKeyTable(kpc,ku,prefix="c.")
}
\keyword{keyProfiler}
