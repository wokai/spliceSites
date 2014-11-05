\name{truncate_seq}
\alias{truncate_seq}
\title{
truncate_seq function
}
\description{
truncateSeqs amino acid sequences at positions depicted by '*' (stop-codon).
}
\usage{
truncate_seq(seq,pos,id,rme=TRUE,trunc=42L)
}
\arguments{
  \item{seq}{Character. Amino-acid sequences which are to be truncated. }
  \item{pos}{Numeric. Points to an amino-acid inside the sequence. }
  \item{id}{Numeric. An identifier which is copied to the result table. }
  \item{rme}{Logical. Empty sequences are removed when set to TRUE.}
  \item{trunc}{Integer. ASCII code for character at which truncation should occur. Default value is 42='*' (stop-codon).}
}
\details{The function truncateSeqs the contained amino acid sequences. When the stop-codon is found on the left side of position, the function returns an empty sequence for that site. The position values for these records are also set to 0.}
\value{data.frame}
\author{
Wolfgang Kaisers
}
\examples{
truncate_seq(seq="ARPX*QR",pos=3)
}
\keyword{truncate_seq}