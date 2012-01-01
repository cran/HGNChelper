checkGeneSymbols <-
function  #function to identify outdated or Excel-mogrified gene symbols
### This function identifies gene symbols which are outdated or may have been
### mogrified by Excel or other spreadsheet programs.  If output is
### assigned to a variable, it returns a data.frame of the same number of
### rows as the input, with a second column indicating whether the symbols
### are valid and a third column with a corrected gene list.
(x,
 ### Vector of gene symbols to check for mogrified or outdated values
 unmapped.as.na=TRUE
 ### If TRUE, unmapped symbols will appear as NA in the
 ### Suggested.Symbol column.  If FALSE, the original unmapped symbol
 ### will be kept.
 ){
    if(class(x) != "character"){
        x <- as.character(x)
        warning("coercing x to character.")
    }
    data(hgnc.table)
    approved <- x %in% hgnc.table$Approved.Symbol
    ##change to upper case, then change orfs back to lower case:
    x.casecorrected <- toupper(x)
    x.casecorrected <- sub("(.*C[0-9XY]+)ORF(.+)", "\\1orf\\2", x.casecorrected)
    if (!all(x == x.casecorrected))
        warning("Some lower-case letters were found and converted to upper-case.
                 HGNChelper is intended for human symbols only, which should be all
                 upper-case except for open reading frames (orf).")
    alias <- x.casecorrected %in% hgnc.table$Symbol
    df <- data.frame(x=x,
                     Approved=approved,
                     Suggested.Symbol=sapply(1:length(x), function(i)
                     ifelse(approved[i], x[i], ifelse(alias[i],paste(hgnc.table$Approved.Symbol[x.casecorrected[i] == hgnc.table$Symbol], 
                     collapse=" /// "), NA))),
                     stringsAsFactors=FALSE)
    if(!unmapped.as.na){
        df[is.na(df$Suggested.Symbol), "Suggested.Symbol"] <- df[is.na(df$Suggested.Symbol), "x"]
    }
    if (sum(df$Approved) != nrow(df)) 
        warning("x contains non-approved gene symbols")

    return(df)
### The function will return a data.frame of the same number of rows as the input,
### with corrections possible from hgnc.table. 
}
