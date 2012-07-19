\name{findExcelGeneSymbols}
\alias{findExcelGeneSymbols}
\title{function to identify Excel-mogrified gene symbols}
\description{This function identifies gene symbols which may have been
mogrified by Excel or other spreadsheet programs.  If output is
assigned to a variable, it returns a vector of the same length
where symbols which could be mapped have been mapped.}
\usage{findExcelGeneSymbols(x, mog.map, regex)}
\arguments{
  \item{x}{Vector of gene symbols to check for mogrified values}
  \item{mog.map}{(optional) map of known mogrifications.  This should be
a dataframe with two columns: original and mogrified, containing
the correct and incorrect symbols, respectively.  The default map is
provided in extdata/mog_map.csv (see examples)
}
  \item{regex}{(optional) regular expression, recognized by the base::grep function, which
is called with ignore.case=TRUE.  It is not necessary for all
matches to have a corresponding entry in mog.map$mogrified;
values in x which are matched by this regex but are not found in
mog.map$mogrified simply will not be corrected.  The default regex is:

"[0-9]\\-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)|[0-9]\\.[0-9][0-9]E\\+[[0-9][0-9]"

This regex is based on that provided by Zeeberg et al., Mistaken
Identifiers: Gene name errors can be introduced inadvertently when using
Excel in bioinformatics.  BMC Bioinformatics 2004, 5:80.}
}

\value{if the return value of the function is assigned to a variable, the
function will return a vector of the same length as the input,
with corrections possible from mog.map made.}

\author{Levi Waldron}

\examples{
### See the default mogrification map:
read.csv(system.file("extdata/mog_map.csv", package = "HGNChelper"), as.is=TRUE)
}




