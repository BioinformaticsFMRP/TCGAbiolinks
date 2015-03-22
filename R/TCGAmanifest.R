TCGAmanifest <- function(Tumor){
  FolderWd <- getwd()
  createDir("TCGAsdrf")
  setwd("./TCGAsdrf")
  on.exit(setwd(FolderWd))

  siteTCGA <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
  tmplf <- PlatformAndAssociatedData[PlatformAndAssociatedData$Tumor %in% tolower(Tumor),]
  for (w in 1:nrow(tmplf)){
    tmp <- tmplf[w,]
    if(tmp["FileName"] == "") next
    key1a <- paste(unique(tmp$CenterType),
                   unique(tmp$Center),
                   unique(tmp$Platform),
                   sep="/"
                   )
    Description <- paste0(siteTCGA,
                          tolower(tmp$Tumor),
                          "/",
                          key1a)
    DescriptionSite <- paste(Description,
                             tmp$Folder,
                             tmplf[w,"FileName"],
                             sep="/"
                             )

    file <- gsub("/","_",tmplf[w,"FileName"])
    # shoul see if the file was updated
    if (!file.exists(file))
      downloader::download(DescriptionSite, file )
    print(paste0(w, " of ", nrow(tmplf)," ", file))
  }
}
