#source('~/Desktop/gitTCGAbiolinks/TCGAbiolinks/R/TCGA_biolinks_functions_anto.R')
#wdpath <- "/Users/antoniocolaprico/Dropbox/tcgabiolinks/Package"
#setwd(wdpath)
#load("PlatformAndAssociatedData_2015_03_11.Rdata")

TCGAupdate <- function(PlatformAndAssociatedData, wdpath = getwd()){
#token means that you can restart from row where it stopped, problem with 76,78,86,94, 104, 106, 131
token <- 1
for( w in token: nrow(PlatformAndAssociatedData)){
  tmp <- PlatformAndAssociatedData[w,]
  key1a<- paste(unique(tmp$CenterType), unique(tmp$Center), unique(tmp$Platform), sep="/")
  Description <- paste(siteTCGA, tolower(tmp$Tumor), "/",key1a, sep="")
  key2a <- paste("/",tmp$Folder,"/",sep="")
  print( paste( w , " of ", nrow(PlatformAndAssociatedData)," ",Description, sep=""))
  PlatformAndAssociatedData[w,"FileName"] <-  .DownloaDmageTAB_sdrf(Description,key2a,tmp$Key1, tmp$Key2) 
  
  save(PlatformAndAssociatedData, file = "PlatformAndAssociatedData_2015_03_11.Rdata")
}

}