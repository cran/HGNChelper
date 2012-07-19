\name{hgnc.symbols}
\alias{hgnc.symbols}
\docType{data}
\title{
All current and withdrawn HGNC symbols.
}
\description{
All current and withdrawn symbols from genenames.org.
}
\usage{
data(hgnc.symbols)
}
\format{
  A dataframe with the first column providing HGNC ID, second column
  providing the HGNC symbol (including withdrawn symbols).
}
\source{
Extracted from genenames.org on July 19, 2012.
}
\examples{
data(hgnc.symbols)
head(hgnc.symbols)
}
\keyword{datasets}
