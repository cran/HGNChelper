\name{rToSymbol}
\alias{rToSymbol}
\title{function to reverse the conversion made by symbolToR}
\description{This function reverses the actions of the symbolToR function.}
\usage{rToSymbol(x)}
\arguments{
  \item{x}{the character vector returned by the symbolToR function.}
}

\value{a character vector of HGNC gene symbols, which are not in general valid R names.}

\examples{
###See that transformation is reversible:
data(hgnc.symbols)
hgnc.symbols <- hgnc.symbols[ ,2]
hgnc.symbols <- hgnc.symbols[!grepl("withdrawn", hgnc.symbols)]
if( !identical(all.equal(hgnc.symbols, rToSymbol(symbolToR(hgnc.symbols))), TRUE))
    stop("HGNC mapping was not reversible.")
}
\author{Levi Waldron}