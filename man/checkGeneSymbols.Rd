\name{checkGeneSymbols}
\alias{checkGeneSymbols}
\seealso{hgnc.table}
\title{function to identify outdated or Excel-mogrified gene symbols}
\description{This function identifies gene symbols which are outdated or may have been
mogrified by Excel or other spreadsheet programs.  If output is
assigned to a variable, it returns a data.frame of the same number of
rows as the input, with a second column indicating whether the symbols
are valid and a third column with a corrected gene list.}
\usage{checkGeneSymbols(x, unmapped.as.na=TRUE, hgnc.table=NULL)}
\arguments{
  \item{x}{Vector of gene symbols to check for mogrified or outdated
    values}
  \item{unmapped.as.na}{If TRUE, unmapped symbols will appear as NA in the
    Suggested.Symbol column.  If FALSE, the original unmapped symbol
    will be kept when no correction can be found.}
  \item{hgnc.table}{If hgnc.table is a data.frame with colnames(hgnc.table) identical
to c("Symbol", "Approved.Symbol"), it is used to correct gene symbols in
x.  Otherwise, the default table data("hgnc.table",
package="HGNChelper") is used.  The function used for creating this
the hgnc.table dataframe can be found in the inst/hgncLookup.R file, in
the source of this package.  By default the hgnc.table dataframe shipped
with this package is used (see ?hgnc.table)}
}

\value{The function will return a data.frame of the same number of rows as the input,
with corrections possible from hgnc.table. }

\author{Levi Waldron and Markus Riester}

\examples{library(HGNChelper)

x = c("FN1", "TP53", "UNKNOWNGENE","7-Sep", "1-Mar") 

res <- checkGeneSymbols(x)
res

stopifnot(sum(res[,1] != x) == 0)
stopifnot(sum(res[,2] != c(TRUE,TRUE,FALSE,FALSE,FALSE)) == 0)
stopifnot( all.equal(as.character(res[,3]), c("FN1","TP53",NA,"SEPT7","MARCH1 /// MARC1")) 
)

}
