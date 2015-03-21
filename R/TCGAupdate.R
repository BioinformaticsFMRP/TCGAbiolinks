#'@title TCGAUpdate
#'
#'@author Davide
#'
#' @export
TCGAUpdate <-function(siteTCGA="https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"){
  PlatformMat <- NULL

  x <- DownloadHTML(siteTCGA)
  x <- GrepSite(x, "href")

  Tumor <- sub("/", "", x)

  siteTumor <- paste0(siteTCGA, Tumor, "/")

  for(tt in 1:length(siteTumor)){
    x <- DownloadHTML(siteTumor[tt])
    x <- GrepSite(x, "href")
    CenterType <- sub("/", "", x)

    siteCenterType <- paste0(siteTumor[tt], CenterType, "/")
    for(ct in 1:length(siteCenterType)){
      x <- DownloadHTML(siteCenterType[ct])
      x <- GrepSite(x, "href")
      Center <- sub("/", "", x)

      siteCenter <- paste0(siteCenterType[ct], Center, "/")
      for(cc in 1:length(siteCenter)){
        x <- DownloadHTML(siteCenter[cc])
        x <- GrepSite(x, "href")
        Platform <- sub("/", "", x)

        sitePlatform <- paste0(siteCenter[cc],Platform, "/")
        for(p in 1:length(sitePlatform)){
          x <- DownloadHTML(sitePlatform[p])
          x <- GrepSite(x, "href")
          kind <- sub("/", "", x)

          siteKind <- paste0(sitePlatform[p],kind, "/")
          for(k in 1:length(siteKind)){
            x <- DownloadHTML(siteKind[k])
            x <- GrepSite(x, "href")
            folder <- sub("/", "", x)

            siteFolder <- paste0(siteKind[k],folder, "/")

            if(length(Platform)==0) Platform = " "
            if(length(kind)==0) kind = " "
            if(length(folder)==0) folder = " "
            PlatformMat <- rbind(PlatformMat,
                                 cbind(
                                   Tumor = rep(Tumor[tt], length(folder)),
                                   CenterType = rep(CenterType[ct], length(folder)),
                                   Center = rep(Center[cc], length(folder)),
                                   Platform = rep(Platform[p],length(folder)),
                                   Kind = rep(kind[k],length(folder)),
                                   Folder = folder,
                                   Manifest = paste0(siteFolder,"MANIFEST.txt")))
          }
        }
        print(paste(Tumor[tt], CenterType[ct], Center[cc], ":", paste(Platform, collapse = " - "), "done!", sep = " "))
      }
    }
  }
  save(PlatformMat,
       file = paste(system.file("extdata", package="TCGAbiolinks"),
                                 "PlatformMat.rda",sep="/")
       )
}
