\name{silic_tryp}
\alias{silic_tryp}
\title{
silic_tryp function
}
\description{
silic_tryp performs silicon trypsination and returns the fragments to which the position coordinate points. The position value is corrected so that it afterwards points to the same amino-acid as before.
}
\usage{
silic_tryp(seq,pos,id)
}
\arguments{
  \item{seq}{Character. Amino-acid sequences which are to be truncated. }
  \item{pos}{Numeric. Points to an amino-acid inside the sequence. }
  \item{id}{Numeric. An identifier which is copied to the result table. }
}
\details{The routine implements the "Keil"-rule, where sites are described by the regex  "[RK](?!P)". The cut position is between [RK] and the following character. The sequence fragment which contains the exon-intron boundary (depicted by \code{position}) is returned. Dependent numeric values are recalculated.}
\value{data.frame}
\author{
Wolfgang Kaisers
}
\examples{
silic_tryp(seq="AXKUEMRFG",pos=4)
}
\keyword{silic_tryp}