#' @name getCurrentMaps
#' @title Get the current maps for correcting gene symbols
#' @aliases getCurrentHumanMap getCurrentMouseMap
#' @description Valid human and mouse gene symbols can be updated frequently. Use
#' these functions to get the most current lists of valid symbols, which you can 
#' then use as an input to the \code{map} argument of \code{checkGeneSymbols}. 
#' Make sure to change the default \code{species="human"} argument to \code{checkGeneSymbols} 
#' if you are doing this for mouse. Use \code{getCurrentHumanMap} for HGNC human gene 
#' symbols from \url{https://www.genenames.org/} and \code{getCurrentMouseMap} for MGI mouse gene 
#' symbols from \url{https://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt}.
#' @return A \code{data.frame} that can be used for \code{map} argument of \code{checkGeneSymbols} function
#' @export getCurrentMouseMap getCurrentHumanMap
#' @usage 
#' getCurrentHumanMap()
#' getCurrentMouseMap()
#' @importFrom stats complete.cases
#' @importFrom utils read.delim
#' @importFrom splitstackshape cSplit
#' @importFrom utils download.file
#' @examples
#' \dontrun{
#' ## human
#' new.hgnc.table <- getCurrentHumanMap()
#' checkGeneSymbols(c("3-Oct", "10-3", "tp53"), map=new.hgnc.table)
#' 
#' ## mouse
#' new.mouse.table <- getCurrentMouseMap()
#' ## Set species to NULL or "mouse" 
#' checkGeneSymbols(c("Gm46568", "1-Feb"), map=new.mouse.table, species="mouse")
#' }
#' 
getCurrentHumanMap <- function(){
  .fixttable <- function(hgnc.table) {
    ## remove withdrawn symbols with known new name
    hgnc.table <- hgnc.table[!(duplicated(hgnc.table$Symbol) 
                               & is.na(hgnc.table$Approved.Symbol)), ]
    hgnc.table <- hgnc.table[order(hgnc.table$Symbol), ]
    
    hgnc.table$Symbol <- as.character(hgnc.table$Symbol)
    hgnc.table$Approved.Symbol <- as.character(hgnc.table$Approved.Symbol)
    rownames(hgnc.table) <- NULL
    
    ## In the un-approved column, convert everything but orfs to upper-case
    hgnc.table$Symbol <- toupper(hgnc.table$Symbol)
    hgnc.table$Symbol <-
      sub("(.*C[0-9XY]+)ORF(.+)", "\\1orf\\2", hgnc.table$Symbol)
    
    hgnc.table <- unique(hgnc.table)
    hgnc.table <- hgnc.table[complete.cases(hgnc.table), ]
    is.ascii <-
      iconv(hgnc.table[, 1], to = "ASCII", sub = ".") == hgnc.table[, 1]
    hgnc.table <- hgnc.table[is.ascii, ]
    return(hgnc.table)
  }
  url <- paste0(
    "https://storage.googleapis.com/public-download-files/",
    "hgnc/tsv/tsv/hgnc_complete_set.txt"
  )
  if(!file.exists("hgnc_complete_set.txt")){
    message(paste("Fetching human gene symbols from", url))
    download.file(url, destfile = "hgnc_complete_set.txt")    
  }else{
    message("Using the already downloaded hgnc_complete_set.txt file.")
  }
  map <- read.delim("hgnc_complete_set.txt", as.is=TRUE)
  if(nrow(map) < 40000){
    warning("The hgnc_complete_set.txt file has fewer than 40,000 rows.",
            "As of May 17 2014, it already had 43,842 rows, so the download was probably incomplete.", 
            "Please check the download.")
  }  
  # 2 = symbol, 3=alias_symbol, 4=prev_symbol, 7=location
  # corrections
  has.corrections <- nchar(map$alias_symbol)>0 | nchar(map$prev_symbol)>0
  M  <- do.call(rbind, apply(map[has.corrections & map$status == "Approved", ], 1,
                             function(x) {
                               y <- strsplit(paste(x[c("alias_symbol", "prev_symbol")], collapse="|"), "|", fixed = TRUE)[[1]]
                               cbind(y[y!=""], x["symbol"])
                             }))
  rownames(M) <- NULL
  # valid symbols
  N <- map[, c("symbol", "symbol", "status")]
  N[N$status=="Entry Withdrawn", 2] <- NA
  N <- N[, 1:2]

  O = read.csv(system.file("extdata/mog_map.csv", package = "HGNChelper"),
               as.is=TRUE)[, 2:1]
  O$mogrified <- toupper(O$mogrified)
  
  colnames(M) <- c("Symbol", "Approved.Symbol")
  colnames(N) <- c("Symbol", "Approved.Symbol")
  colnames(O) <- c("Symbol", "Approved.Symbol")
  
  external.table <- rbind(M, N)
  ## Correct outdated symbols in extdata/mog_map.csv
  O_corrected <-
    suppressWarnings(checkGeneSymbols(O[, "Approved.Symbol"], map = external.table))
  O[, "Approved.Symbol"] <- O_corrected[, "Suggested.Symbol"]
  output <- rbind(external.table, O)
  output <- .fixttable(output)
  output <- merge(output, map[, c("location", "symbol")], by = 2)
  output <- output[, c(2, 1, 3)]
  output <- suppressWarnings(splitstackshape::cSplit(splitstackshape::cSplit(output, 3, "p"), 3, "q"))
  output <- output[, c(1, 2, 5)]
  colnames(output)[3] <- "chromosome"
  return(output)
}

getCurrentMouseMap <- function(){
  fname <- system.file("extdata/HGNChelper_mog_map_MGI_AMC_2016_03_30.csv", package = "HGNChelper")
  O = read.csv(fname, as.is=TRUE)[, 2:1]
  colnames(O) <- c("Symbol", "Approved.Symbol")
  url <- "https://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt"
  if(!file.exists("MGI_EntrezGene.rpt")){
    message(paste("Fetching gene symbols from", url))
    download.file(url, destfile = "MGI_EntrezGene.rpt")    
  }else{
    message("Using the already downloaded MGI_EntrezGene.rpt file.")
  }
  map <- read.delim("MGI_EntrezGene.rpt", as.is=TRUE, header = FALSE)
  if(nrow(map) < 700000){
    warning("The MGI_EntrezGene.rpt file has fewer than 700,000 rows.",
    "As of May 17 2014, it already had 752,802 rows, so the download was probably incomplete.", 
    "Please check the download.")
  }
  map <- map[, 2:4]
  
  map.ok <- map[grep("^O$", map[, 2]), c(1, 3)]
  map.ok <- data.frame(Symbol=unique(map.ok[, 1]), Approved.Symbol=unique(map.ok[, 1]), stringsAsFactors = FALSE)
  
  map.withdrawn <- map[grep("^W$", map[, 2]), c(1, 3)]
  map.withdrawn[, 2] <- sub("^withdrawn, ", "", map.withdrawn[, 2])
  map.withdrawn[, 2] <- sub("^= ", "", map.withdrawn[, 2])
  colnames(map.withdrawn) <- c("Symbol", "Approved.Symbol")
  
  mouse.table <- rbind(O, map.withdrawn, map.ok)
  mouse.table <- unique(mouse.table)
  mouse.table <- mouse.table[complete.cases(mouse.table), ]
  rownames(mouse.table) <- NULL
  mouse.table <- mouse.table[!is.na(iconv(mouse.table[, 1], "ASCII")), ]
  return(mouse.table)
}