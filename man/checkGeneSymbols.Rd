% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkGeneSymbols.R
\name{checkGeneSymbols}
\alias{checkGeneSymbols}
\title{Title Identify outdated or Excel-mogrified gene symbols}
\usage{
checkGeneSymbols(x, unmapped.as.na = TRUE, map = NULL,
  species = "human")
}
\arguments{
\item{x}{Vector of gene symbols to check for mogrified or outdated values}

\item{unmapped.as.na}{If TRUE (default), unmapped symbols will appear as NA in the
Suggested.Symbol column.  If FALSE, the original unmapped symbol
will be kept.}

\item{map}{Specify if you do not want to use the default maps provided 
by setting species equal to "mouse" or "human". 
map can be any other data.frame with colnames(map) identical
to c("Symbol", "Approved.Symbol"). The default maps can be updated 
by running the interactive example below.}

\item{species}{A character vector of length 1, either "human" (default) or "mouse". 
If NULL, or anything other than "human" or "mouse", then the map
argument must be provided.}
}
\value{
The function will return a data.frame of the same number of rows as the input,
with corrections possible from map.
}
\description{
This function identifies gene symbols which are outdated or may have been
 mogrified by Excel or other spreadsheet programs.  If output is
 assigned to a variable, it returns a data.frame of the same number of
 rows as the input, with a second column indicating whether the symbols
 are valid and a third column with a corrected gene list.
}
\examples{
library(HGNChelper)
human = c("FN1", "TP53", "UNKNOWNGENE","7-Sep", "9/7", "1-Mar", "Oct4", "4-Oct",
      "OCT4-PG4", "C19ORF71", "C19orf71")
checkGeneSymbols(human)
## mouse
mouse <- c("1-Feb", "Pzp", "A2m")
checkGeneSymbols(mouse, species="mouse")
if (interactive()){
  ##Run checkGeneSymbols with a brand-new map downloaded from HGNC:
  source(system.file("hgncLookup.R", package = "HGNChelper"))
  ## You should save this if you are going to use it multiple times,
  ## then load it from file rather than burdening HGNC's servers.
  save(hgnc.table, file="hgnc.table.rda", compress="bzip2")
  load("hgnc.table.rda")
  checkGeneSymbols(human, species=NULL, map=hgnc.table)
  checkGeneSymbols(human, species=NULL, map=mouse.table)
}
}
\seealso{
\code{\link{mouse.table}} for the mouse lookup table, 
\code{\link{hgnc.table}} for the human lookup table
}
