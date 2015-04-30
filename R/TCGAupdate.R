#'@title TCGAUpdate
#'
#'@description replace the dataFolders.rda matrix with the
#'             new matrix from "The Cancer Genome Atlas" ftp
#'@author Davide
#'@seealso TCGAQuery
#'@export
TCGAUpdate <-function(){

}

tcga.check.updates <- function(){
  message("Looking for updates in TCGA database")

  table <- create.tcga.table ()

  message("Verifying if there is any update")
  compare <- sapply(table$name,function(x){is.element(x,tcga.db[,"name"])})

  # get new rows
  table <- table[as.integer(which(compare == F)),]

  # get barcode for new rows
  table <- tcga.get.barcode(table)
  # add disease and platform to table
  table <- tcga.db.addCol(table)

  # add table to tcga.db
  tcga.db <- rbind(tcga.db,table)

}
# This function created tcga.db
# warning to create the db takes almost 3 days
# to update it should take some minutes (function still needed to be done)
#' @import downloader XML plyr stringr
create.tcga.table <- function(){
  # get all compressed archives
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  regex <- '<table summary="Data Summary".*</a></td></tr></table>'
  message("Downloading TCGA database")
  tcga.query <- paste0("query=Archive[isLatest=1]")
  next.url <- paste0(tcga.root,tcga.query)
  pb <- txtProgressBar(min = 0, max = 32, style = 3)
  i <- 0
  while(!is.na(next.url)){
    downloader::download(next.url,"tcga.html",quiet=T)
    html <- readLines("tcga.html")
    table <- readHTMLTable(toString(str_match(html,regex)[6,]),
                           header = T,
                           stringsAsFactors = FALSE)$'NULL'
    colnames(table) <-table[1,]
    table <- table[-1,1:9]

    if(exists("db")) {
      db <- rbind(db,table)
    } else {
      db <- table
    }
    # get next table
    next.regex <-"http://tcga-data.nci.nih.gov/tcgadccws/GetHTML.*Next"
    next.url <- str_match(html,next.regex)[6,]
    next.url <- gsub("amp;","",gsub('\">Next',"",next.url))

    # update progress bar
    i <- i + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
  # remove protected data from the update
  db <- db[-grep("tcga4yeo",db$deployLocation),]
  db$addedDate <- as.Date(db$addedDate,"%m-%d-%Y")
  return(db)
}

# Using the name create two collumns Platform and Disease
tcga.db.addCol <- function(data){
  data$Platform <- ""
  data$Disease <- ""
  diseases <- disease.table$abbreviation
  for(i in seq_along(diseases)){
    idx <- grep(diseases[i],data$baseName)
    if(length(idx >0)){
      data[idx,]$Disease <- diseases[i]
    }
  }

  for(i in seq_along(platform.table$name)){
    idx <- grep(platform.table[i,]$name,data$baseName)
    if(length(idx)>0){
      data[idx,]$Platform <- platform.table[i,]$name
    }
  }
  return (data)
}



tcga.get.barcode <- function(data){
  # todo considere next page =/
  # comparar com sdrf
  # for each tcga.db id get barcodes
  message("Downloading TCGA barcodes")
  start.time <- Sys.time()
  pb <- txtProgressBar(min = 0, max = nrow(data), style = 3)
  for(j in 1:nrow(data)){
    #for(j in c(5)){
    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?query="
    tcga.query <- paste0("FileInfo&Archive[@id=",data[j,'id'],"]&roleName=fileCollection")
    # metadata files ignore barcodes
    if(length(grep(".*(aux|mage-tab).*",data[j,'name']))>0){
      print(paste("Continued for:",j))
      next
    }
    # search if a file with same serial index and revision have barcodes
    # if it has we should copy it!
    rev <- data[j,"revision"]
    index <- data[j,"serialIndex"]
    bname<- data[j,"baseName"]
    aux <- subset(data,
                  revision == rev & serialIndex == index & baseName==bname & deployStatus!="Available",
                  select = c(deployLocation,deployStatus)
    )

    if(nrow(aux)>0){
      print("Found barcode before")
      data[j,"deployStatus"] <- aux[1,"deployStatus"]
      next
    }
    next.url <- paste0(tcga.root,tcga.query)
    while(!is.na(next.url)){
      print(next.url)
      downloader::download(next.url,"tcga.html",quiet =T)
      regex <- '<table summary="Data Summary".*</a></td></tr></table>'
      html <- readLines("tcga.html")
      files <- readHTMLTable(toString(str_match(html,regex)[6,]),
                             header = T,
                             stringsAsFactors = FALSE)$'NULL'
      colnames(files) <- files[1,]
      files <- files[-1,1:5]
      idx <- grep("(README|CHANGES|DESCRIPTION|MANIFEST|DISCLAIMER).*",files$name)
      if(length(idx > 0)){files <- files[-idx,]}
      if(nrow(files) == 0) {
        next.url <- NA
        if(!exists("all.barcode")){
          all.barcode <- ""
        }
        next
      }

      #pat <- "*(TCGA)-([A-Z0-9]{2})-([A-Z0-9]{4})-(0[1-9]|[1-2][0-9])([A-Z])-(0[1-9]|[1-9][0-9])([DGHRTWX])-([A-Z0-9]{4})-([A-Z0-9]{2})*"
      pat <- "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{3}-[A-Z0-9]{3}-[A-Z0-9]{4}-[A-Z0-9]{2}"
      #pat <- "((((TCGA-[A-Z0-9]{2})-[A-Z0-9]{4})-[A-Z0-9]{3}-[A-Z0-9]{3})-[A-Z0-9]{4}-[A-Z0-9]{2})"

      barcode <- str_match(files$name,pat)[,1]
      print(Sys.time()-start.time)
      #message("Found in name")
      #print(barcode)
      if(all(is.na(barcode))){
        message("Not Found in name")
        for(i in 1:nrow( files )){
          tcga.query <- paste0("BiospecimenBarcode&FileInfo[@id=",files[i,"id"],"]&roleName=biospecimenBarcodeCollection")
          next.url <- paste0(tcga.root,tcga.query)
          downloader::download(next.url,"tcga.html",quiet =T)
          regex <- '<table summary="Data Summary".*</a></td></tr></table>'
          barcode.html <- readLines("tcga.html")
          table <- str_match(barcode.html,regex)
          idx <- which(!is.na(table))
          if(length(idx>0)){
            barcode.table <- readHTMLTable(toString(table[idx,]),
                                           header = T,
                                           stringsAsFactors = FALSE)$'NULL'
            colnames(barcode.table) <- barcode.table[1,]
            barcode.table <- barcode.table[-1,1:8]
            #print(barcode.table$barcode)
            files[[i,5]] <- list(barcode.table$barcode)
          }
        }
        barcode <- unique(unlist(files[,5]))
      }
      # get next table
      next.regex <-"http://tcga-data.nci.nih.gov/tcgadccws/GetHTML.*Next"
      next.url <- str_match(html,next.regex)[6,]
      next.url <- gsub("amp;","",gsub('\">Next',"",next.url))
      if(exists("all.barcode")){
        all.barcode <- c(all.barcode,barcode)
      } else {
        all.barcode <- barcode
      }
    }
    data[j,"deployStatus"] <- paste0(all.barcode, collapse = ",")
    rm(all.barcode)
    setTxtProgressBar(pb, j)
    print(Sys.time()-start.time)
  }
  close(pb)
  return (data)
}
