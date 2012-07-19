\name{affyToR}
\alias{affyToR}
\title{function to convert Affymetrix probeset identifiers to valid R names}
\description{ Simply prepends "affy." to the input vector. }
\usage{affyToR(x)}
\arguments{
  \item{x}{vector of Affymetrix probeset identifiers, or any identifier which may with a digit.}
}

\value{a character vector that is simply x with "affy." prepended to each value.}

\author{Levi Waldron}





