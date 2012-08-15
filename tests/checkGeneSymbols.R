library(HGNChelper)

x = c("FN1", "TP53", "UNKNOWNGENE","7-Sep", "1-Mar") 

res <- checkGeneSymbols(x)

stopifnot(sum(res[,1] != x) == 0)
stopifnot(sum(res[,2] != c(TRUE,TRUE,FALSE,FALSE,FALSE)) == 0)
stopifnot( all.equal(as.character(res[,3]), c("FN1","TP53",NA,"SEPT7","MARCH1 /// MARC1")) 
)


