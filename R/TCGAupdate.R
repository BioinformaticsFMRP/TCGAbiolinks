#'@title TCGAUpdate
#'
#'@description replace the dataFolders.rda matrix with the
#'             new matrix from "The Cancer Genome Atlas" ftp
#'@author Davide
#'@seealso TCGAQuery
#'@export
#'@import downloader RCurl XML
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
        if(length(Platform)==0) Platform = ""

        sitePlatform <- paste(siteCenter[cc],Platform, "/", sep = "")
        for(p in 1:length(sitePlatform)){
          x <- DownloadHTML(sitePlatform[p])
          x <- GrepSite(x, "href")
          kind <- sub("/", "", x)
          if(length(kind)==0) kind = ""

          siteKind <- paste(sitePlatform[p],kind, "/", sep = "")
          for(k in 1:length(siteKind)){
            x <- DownloadHTML(siteKind[k])
            x <- GrepSite(x, "href")
            folder <- sub("/", "", x)
            if(length(folder)==0) folder = ""

            siteFolder <- paste(siteKind[k],folder, "/", sep = "")




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
