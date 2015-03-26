
# if(TRUE %in% grepl("TCGABiolinks", installed.packages())){
#   pkg <- c("package:TCGABiolinks")
#   lapply(pkg, detach, character.only = TRUE, unload = TRUE)
#   remove.packages("TCGABiolinks")
#   rm(pkg)
# }
lapply(c("package:TCGAbiolinks"), detach, character.only = TRUE, unload = TRUE)
remove.packages("TCGAbiolinks")
#set a working directory:
wd <-"~/Applications/GitHub/"
#path of the package:
pathPkg <- "TCGAbiolinks/"

setwd(wd)

#install.packages(c("downloader","RCurl","httr","devtools"))
install.packages(pathPkg, type = 'source', repos =NULL)

library("downloader")
library("httr")
library("devtools")
library('TCGAbiolinks')

#update the dataFolders.rda
#TCGAUpdate()
tumor = "acc"
centerType = "cgcc"
platform = "genome_wide_snp_6"
level = 3

# tumor = "lihc"
# platform = "mixed_dnaseq_curated"
# centerType = "gsc"
# level=2

# TCGAQuery(tumor = tumor,
#           centerType = centerType,
#           platform = platform,
#           level = level)
#
# earlyStop <- 20
# TCGADownload(earlyStop = earlyStop)




version <- TCGAVersion(tumor = tumor,
                       centerType = centerType,
                       platform = platform,
                       level = level,
                       barcode = T)





