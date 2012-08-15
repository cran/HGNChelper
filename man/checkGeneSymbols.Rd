\name{checkGeneSymbols}
\alias{checkGeneSymbols}
\title{function to identify outdated or Excel-mogrified gene symbols}
\description{This function identifies gene symbols which are outdated or may have been
mogrified by Excel or other spreadsheet programs.  If output is
assigned to a variable, it returns a data.frame of the same number of
rows as the input, with a second column indicating whether the symbols
are valid and a third column with a corrected gene list.}
\usage{checkGeneSymbols(x)}
\arguments{
  \item{x}{Vector of gene symbols to check for mogrified or outdated values}
}

\value{The function will return a data.frame of the same number of rows as the input,
with corrections possible from hgnc.table. }

\author{Levi Waldron and Markus Riester}




\examples{library(HGNChelper)

x = c("FN1", "TP53", "UNKNOWNGENE","7-Sep", "1-Mar") 

res <- checkGeneSymbols(x)

stopifnot(sum(res[,1] != x) == 0)
stopifnot(sum(res[,2] != c(TRUE,TRUE,FALSE,FALSE,FALSE)) == 0)
stopifnot( all.equal(as.character(res[,3]), c("FN1","TP53",NA,"SEPT7","MARCH1 /// MARC1")) 
)

}
