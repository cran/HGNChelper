library(HGNChelper)

x = c("C21orf49", "c21ORF49",
      "MORF4L1P7", "Morf4L1P7",
      "Fn1", "FN1",
      "TP53",
      "UNKNOWNGENE",
      "7-Sep", "1-Mar", "1-MAR")

res <- checkGeneSymbols(x)

stopifnot(sum(res[,1] != x) == 0)
stopifnot(sum(res[,2] != c(TRUE, FALSE,
                           TRUE, FALSE,
                           FALSE, TRUE,
                           TRUE,
                           FALSE,
                           FALSE, FALSE, FALSE)) == 0)
stopifnot( all.equal(as.character(res[,3]), c("C21orf49", "C21orf49",
                                              "MORF4L1P7", "MORF4L1P7",
                                              "FN1", "FN1",
                                              "TP53",
                                              NA,
                                              "SEPT7", "MARCH1 /// MARC1", "MARCH1 /// MARC1")) 
)
