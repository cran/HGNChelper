\name{hgnc.table}
\alias{hgnc.table}
\docType{data}
\title{
All current and withdrawn HGNC symbols.
}
\description{
All current and withdrawn symbols from genenames.org.
}
\usage{
data(hgnc.table)
}
\format{
  A dataframe with the first column providing a gene symbol or known alias
  (including withdrawn symbols), second column
  providing the approved HGNC symbol.
}
\source{
Extracted from genenames.org on July 21, 2012.
}
\examples{
data(hgnc.table)
head(hgnc.table)
}
\keyword{datasets}
