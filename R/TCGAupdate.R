#'@title TCGAUpdate
#'
#'@description replace the dataFolders.rda matrix with the
#'             new matrix from "The Cancer Genome Atlas" ftp
#'@author Davide
#'@seealso TCGAQuery
#'@export
TCGAUpdate <-function(){
  siteTCGA <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
  dataFolders <- NULL


  x <- DownloadHTML(siteTCGA)
  x <- GrepSite(x, "href")

  Tumor <- sub("/", "", x)

  siteTumor <- paste(siteTCGA, Tumor, "/", sep = "")
  for(tt in 1:length(siteTumor)){
    x <- DownloadHTML(siteTumor[tt])
    x <- GrepSite(x, "href")
    CenterType <- sub("/", "", x)

    siteCenterType <- paste(siteTumor[tt], CenterType, "/", sep = "")
    for(ct in 1:length(siteCenterType)){
      x <- DownloadHTML(siteCenterType[ct])
      x <- GrepSite(x, "href")
      Center <- sub("/", "", x)

      siteCenter <- paste(siteCenterType[ct], Center, "/", sep = "")
      for(cc in 1:length(siteCenter)){
        x <- DownloadHTML(siteCenter[cc])
        x <- GrepSite(x, "href")
        Platform <- sub("/", "", x)

        sitePlatform <- paste(siteCenter[cc],Platform, "/", sep = "")
        for(p in 1:length(sitePlatform)){
          x <- DownloadHTML(sitePlatform[p])
          x <- GrepSite(x, "href")
          kind <- sub("/", "", x)

          siteKind <- paste(sitePlatform[p],kind, "/", sep = "")
          for(k in 1:length(siteKind)){
            x <- DownloadHTML(siteKind[k])
            x <- GrepSite(x, "href")
            folder <- sub("/", "", x)

            siteFolder <- paste(siteKind[k],folder, "/", sep = "")

            if(length(Platform)==0) Platform = " "
            if(length(kind)==0) kind = " "
            if(length(folder)==0) folder = " "
            dataFolders <- rbind(dataFolders,
                                 cbind(
                                   Tumor = rep(Tumor[tt], length(folder)),
                                   CenterType = rep(CenterType[ct], length(folder)),
                                   Center = rep(Center[cc], length(folder)),
                                   Platform = rep(Platform[p],length(folder)),
                                   Kind = rep(kind[k],length(folder)),
                                   Folder = folder,
                                   Manifest = paste(siteFolder,"MANIFEST.txt",sep="")))
          }
        }
        print(paste(Tumor[tt], CenterType[ct], Center[cc], ":", paste(Platform, collapse = " - "), "done!", sep = " "))
      }
    }
  }

  save(dataFolders,
       file = paste(system.file("extdata", package="TCGAbiolinks"),
                    "dataFolders.rda",sep="/")
  )
}

# This function created tcga.db
# warning to create the db takes almost 3 days
# to update it should take some minutes (function still needed to be done)
#' @import downloader XML plyr stringr
create.tcga.table <- function(env=as.environment("package:TCGAbiolinks")){
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

    if(exists("tcga.db")) {
      tcga.db <- rbind(tcga.db,table)
    } else {
      tcga.db <- table
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
  env$tcga.db <- tcga.db
}


# save as: list (list(barcodes)) or collapse (paste0(barcode, collapse = ","))?
load.tcga.barcode <- function(){
  # todo considere next page =/
  # comparar com sdrf
  # for each tcga.db id get barcodes
  message("Downloading TCGA barcodes")
  start.time <- Sys.time()
  pb <- txtProgressBar(min = 0, max = nrow(tcga.db), style = 3)
  for(j in 1:nrow(tcga.db)){
    #for(j in c(5)){
    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?query="
    tcga.query <- paste0("FileInfo&Archive[@id=",tcga.db[j,'id'],"]&roleName=fileCollection")
    # metadata files ignore barcodes
    if(length(grep(".*(aux|mage-tab).*",tcga.db[j,'name']))>0){
      print(paste("Continued for:",j))
      next
    }
    # search if a file with same serial index and revision have barcodes
    # if it has we should copy it!
    rev <- tcga.db[j,"revision"]
    index <- tcga.db[j,"serialIndex"]
    bname<- tcga.db[j,"baseName"]
    aux <- subset(tcga.db,
                  revision == rev & serialIndex == index & baseName==bname & deployStatus!="Available",
                  select = c(deployLocation,deployStatus)
    )

    if(nrow(aux)>0){
      print("Found barcode before")
      tcga.db[j,"deployStatus"] <- aux[1,"deployStatus"]
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
      idx <- grep("(README|CHANGES|DESCRIPTION|MANIFEST).*",files$name)
      if(length(idx > 0)){files <- files[-idx,]}
      if(nrow(files) == 0) {
        next.url <- NA
        if(!exists("all.barcode")){
          all.barcode <- ""
        }
        next
      }
      pat <- "*(TCGA)-([A-Z0-9]{2})-([A-Z0-9]{4})-(0[1-9]|[1-2][0-9])([A-Z])-(0[1-9]|[1-9][0-9])([DGHRTWX])-([A-Z0-9]{4})-([A-Z0-9]{2})*"
      barcode <- str_match(files$name,pat)[,1]
      #message("Found in name")
      #print(barcode)
      if(all(is.na(barcode))){
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
    #tcga.db[[j,"deployStatus"]] <- list(barcode)
    tcga.db[j,"deployStatus"] <- paste0(all.barcode, collapse = ",")
    rm(all.barcode)
    setTxtProgressBar(pb, j)
    print(Sys.time()-start.time)
  }
  close(pb)
}


# Using the name create two collumns Platform and Disease
tcga.db.addCol <- function(x){
  tcga.db$Platform <- ""
  tcga.db$Disease <- ""
  diseases <- sapply(strsplit(biosample.tcga$biosample,split = " - "),function(x){x[1]})
  for(i in seq_along(diseases)){
    idx <- grep(diseases[i],tcga.db$baseName)
    tcga.db[idx,]$Disease <- diseases[i]
  }

  for(i in seq_along(platform.table$name)){
    idx <- grep(platform.table[i,]$name,tcga.db$baseName)
    if(length(idx)>0){
      tcga.db[idx,]$Platform <- platform.table[i,]$name
    }
  }
}
