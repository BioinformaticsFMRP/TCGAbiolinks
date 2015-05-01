#'@title TCGAUpdate
#'
#'@description replace the dataFolders.rda matrix with the
#'             new matrix from "The Cancer Genome Atlas" ftp
#'@author Davide
#'@seealso TCGAQuery
#'@export
TCGAUpdate <-function(){

}

tcga.check.updates <- function(platform=NULL,type=NULL,disease=NULL){
  message("Looking for updates in TCGA database")

  table <- create.tcga.table (platform=platform,type=type,disease=disease)

  message("Verifying if there is any update")
  compare <- sapply(table$name,function(x){is.element(x,tcga.db[,"name"])})

  # get new rows
  table <- table[as.integer(which(compare == F)),]

  # get barcode for new rows
  table <- tcga.get.barcode(table)

  # add table to tcga.db
  tcga.db <- rbind(tcga.db,table)

  # remove old values
  # compare db name with the new ones, if not present it is old
  aux <- tcga.db
  if(!is.null(platform)){ aux <- subset(aux, Platform = platform)}
  if(!is.null(disease)){ aux <- subset(aux, Disease = disease)}
  if(!is.null(type)){ aux <- aux[ grepl(type,aux[,"name"]),]}
  compare <- sapply(aux$name,function(x){is.element(x,table[,"name"])})
  aux <- aux[as.integer(which(compare == F)),]
  for(i in seq_along(aux$name)){
    idx <- grep(aux[i,"name"], tcga.db$name)
    tcga.db <- tcga.db[-idx,]
  }
}
# This function created tcga.db
# warning to create the db takes almost 3 days
# to update it should take some minutes (function still needed to be done)
#' @import downloader XML plyr stringr
create.tcga.table <- function(disease=NULL, platform=NULL,type=NULL){
  # get all compressed archives
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  regex <- '<table summary="Data Summary".*</a></td></tr></table>'
  message("Downloading TCGA database")
  tcga.query <- paste0("query=Archive[isLatest=1]")

  extra <- ""
  pages <- 32
  if(!is.null(platform)){
    extra <- paste0(extra,"[Platform[@name=",platform,"]]")
    pages <- pages/4
  }
  if(!is.null(disease)){
    extra <- paste0(extra,"[Disease[@abbreviation=",disease,"]]")
    pages <- pages/4
  }
  if(!is.null(disease)){
    extra <- paste0(extra,"[ArchiveType[@type=",type,"]]")
    pages <- pages/4
  }
  url <- paste0(tcga.root,tcga.query,extra)
  db <- tcga.get.table(url,pages)

  # remove useless cols
  db <- db[,1:9]

  # remove protected data from the update
  idx <- grep("tcga4yeo",db$deployLocation)
  if(length(idx)>0){
    db <- db[-idx,]
  }
  # transoform date into same format
  db$addedDate <- as.Date(db$addedDate,"%m-%d-%Y")
  db <- tcga.db.addCol(db)
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

# A function to get a tcga table from api
# input: url (max for progress bar)
# return: table
tcga.get.table <- function(url,max){
  next.url <- url
  regex <- '<table summary="Data Summary".*</a></td></tr></table>'
  next.regex <-"http://tcga-data.nci.nih.gov/tcgadccws/GetHTML.*Next"

  if(max>0){
    pb <- txtProgressBar(min = 0, max = max, style = 3)
    i <- 0
  }

  while(!is.na(next.url)){
    downloader::download(next.url,"tcga.html",quiet=T)
    html <- readLines("tcga.html")
    table <- readHTMLTable(toString(str_match(html,regex)[6,]),
                           header = T,
                           stringsAsFactors = FALSE)$'NULL'
    colnames(table) <-table[1,]
    table <- table[-1,]

    if(exists("db")) {
      db <- rbind(db,table)
    } else {
      db <- table
    }
    # get next table
    next.url <- str_match(html,next.regex)[6,]
    next.url <- gsub("amp;","",gsub('\">Next',"",next.url))

    if(max>0){
      # update progress bar
      i <- i + 1
      setTxtProgressBar(pb, i)
    }
  }
  if(max>0){
    setTxtProgressBar(pb, max)
    close(pb)
  }
  return (db)
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
    if(length(grep(".*(aux|mage-tab|diagnostic_images|pathology_reports|MDA_RPPA_Core|tissue_images).*",data[j,'name']))>0){
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
        break
      }

      # barcode in only one page (it happens only for bio plat level 1)
      if(data[j,'Platform']=="bio"){
        all.barcode <- ""
        if(grepl("Level_1",data[j,"name"]))
        {
          tcga.query <- paste0("BiospecimenBarcode&Archive[@id=",data[j,'id'],"]&roleName=bcrBiospecimenBarcodeCollection")
          aux <- tcga.get.table(paste0(tcga.root,tcga.query),0)
          all.barcode <- aux[,1]
        }
        next.url <- NA
        break
      }
      #pat <- "*(TCGA)-([A-Z0-9]{2})-([A-Z0-9]{4})-(0[1-9]|[1-2][0-9])([A-Z])-(0[1-9]|[1-9][0-9])([DGHRTWX])-([A-Z0-9]{4})-([A-Z0-9]{2})*"
      pat <- "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{3}-[A-Z0-9]{3}-[A-Z0-9]{4}-[A-Z0-9]{2}"
      #pat <- "((((TCGA-[A-Z0-9]{2})-[A-Z0-9]{4})-[A-Z0-9]{3}-[A-Z0-9]{3})-[A-Z0-9]{4}-[A-Z0-9]{2})"
      barcode <- str_match(files$name,pat)[,1]
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
