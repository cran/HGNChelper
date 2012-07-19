\name{symbolToR}
\alias{symbolToR}
\title{function to *reversibly* convert HGNC gene symbols to valid R names.}
\description{This function reversibly converts HGNC gene symbols to valid R
names by prepending "symbol.", and making the following
substitutions: "-" to "hyphen", "@" to "ampersand", and "/" to
"forwardslash".}
\usage{symbolToR(x)}
\arguments{
  \item{x}{vector of HGNC symbols}
}

\value{a vector of valid R names, of the same length as x, which can be
converted to the same HGNC symbols using the rToSymbol function.}

\examples{
###See that transformation is reversible:
data(hgnc.symbols)
hgnc.symbols <- hgnc.symbols[ ,2]
hgnc.symbols <- hgnc.symbols[!grepl("withdrawn", hgnc.symbols)]
if( !identical(all.equal(hgnc.symbols, rToSymbol(symbolToR(hgnc.symbols))), TRUE))
    stop("HGNC mapping was not reversible.")
}

\author{Levi Waldron}
