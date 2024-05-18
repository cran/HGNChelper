## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  checkGeneSymbols(x, unmapped.as.na = TRUE, map = NULL, species = "human")

## -----------------------------------------------------------------------------
library(HGNChelper)
human = c("FN1", "tp53", "UNKNOWNGENE","7-Sep", "9/7", "1-Mar", "Oct4", "4-Oct",
      "OCT4-PG4", "C19ORF71", "C19orf71")
checkGeneSymbols(human)

## -----------------------------------------------------------------------------
checkGeneSymbols(c("1-Feb", "Pzp", "A2m", "E430008G22Rik"), species="mouse")

## -----------------------------------------------------------------------------
dim(mouse.table)
dim(hgnc.table)

## -----------------------------------------------------------------------------
checkGeneSymbols("AAVS1", expand.ambiguous = FALSE)
checkGeneSymbols("AAVS1", expand.ambiguous = TRUE)

checkGeneSymbols(c("Cpamd8", "Mug2"), species = "mouse", expand.ambiguous = FALSE)
checkGeneSymbols(c("Cpamd8", "Mug2"), species = "mouse", expand.ambiguous = TRUE)

