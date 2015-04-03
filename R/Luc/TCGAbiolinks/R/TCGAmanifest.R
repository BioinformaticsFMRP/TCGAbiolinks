TCGAmanifest <-
function(Tumor){
  FolderWd <- getwd()
  siteTCGA <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"   
  .createDirectory("TCGAsdrf")
  setwd("./TCGAsdrf") 
  tmplf <- PlatformAndAssociatedData[PlatformAndAssociatedData$Tumor %in% tolower(Tumor),]
  token <- 1
  for ( w in token: nrow( tmplf)){
    tmp <- tmplf[w,]
    if(tmp["FileName"] == "") next
    key1a<- paste(unique(tmp$CenterType), unique(tmp$Center), unique(tmp$Platform), sep="/")
    Description <- paste(siteTCGA, tolower(tmp$Tumor), "/",key1a, sep="")
    DescriptionSite <- paste(Description, tmp$Folder, tmplf[w,"FileName"], sep="/")
    x <- .DownloadURL(DescriptionSite)
    #x <- .DownloadData_fromURL(DescriptionSite)
    writeLines(x, gsub("/","_",tmplf[w,"FileName"]))
    
    print(paste(w, " of ", nrow(tmplf ), " ", tmplf[w,"FileName"]),sep="")
  }
  setwd(FolderWd)
}
