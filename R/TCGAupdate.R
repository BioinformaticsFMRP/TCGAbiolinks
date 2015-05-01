#'@title TCGAUpdate
#'
#'@description replace the dataFolders.rda matrix with the
#'             new matrix from "The Cancer Genome Atlas" ftp
#'@author Davide
#'@seealso TCGAQuery
#'@export
getArchive <- function(id){
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- paste0("query=Archive&FileInfo[@id=",id,"]&roleName=archiveCollection")
  url <- paste0(tcga.root,tcga.query)
  db <- tcga.get.table(url)
  return(db)
}

get.samples.files <- function(id){
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- paste0("query=FileInfo&BiospecimenBarcode[@id=",id,"]&roleName=fileCollection")
  url <- paste0(tcga.root,tcga.query)
  db <- tcga.get.table(url)
  obsolete <- grep("-",db$md5sum)

  if(length(obsolete)>0){db <- db[-obsolete,]}
  files <- unique(db$name)
  #for each file, get the highest ID
  for(i in seq_along(files)){
    same.files <- grep(files[i],db$name)
    latest <- max(as.integer(db[same.files,]$id))
    if(!exists("archives")){
      archives  <- getArchive(latest)
      archives$file <- files[i]

    } else {
      aux  <- getArchive(latest)
      aux$file <- files[i]
      archives <- rbind(archives,aux)
    }
  }
  archives <- archives[,-c(10:15)]
  archives$addedDate <- as.Date(archives$addedDate,"%m-%d-%Y")
  archives <- tcga.db.addCol(archives)
  return(archives)
}

getBcrArchiveCollection <- function (id){
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- paste0("query=Archive&BiospecimenBarcode[@id=",id,"]&roleName=bcrArchiveCollection")
  url <- paste0(tcga.root,tcga.query)
  db <- tcga.get.table(url)
  db$addedDate <- as.Date(db$addedDate,"%m-%d-%Y")
  db <- tcga.db.addCol(db)
}
get.barcode.table <- function(barcode){
  tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"
  tcga.query <- paste0("query=BiospecimenBarcode[@barcode=",barcode,"]")
  url <- paste0(tcga.root,tcga.query)
  db <- tcga.get.table(url)
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
  if(!is.null(type)){
    extra <- paste0(extra,"[ArchiveType[@type=Level_",type,"]]")
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
tcga.get.table <- function(url,max=0){
  next.url <- url
  regex <- '<table summary="Data Summary".*</a></td></tr></table>'
  next.regex <-"http://tcga-data.nci.nih.gov/tcgadccws/GetHTML.*Next"

  if(max>0){
    pb <- txtProgressBar(min = 0, max = max, style = 3)
    i <- 0
  }
  #print(url)
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
