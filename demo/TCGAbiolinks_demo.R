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
library("XML")

#update the dataFolders.rda
#TCGAUpdate()
tumor = "acc"
centerType = "cgcc"
platform = "genome_wide_snp_6"
level = 3
metadata = F
metadata = T

# tumor = "lihc"
# platform = "mixed_dnaseq_curated"
# centerType = "gsc"
# level=2
# metadata = T



ullala <- NULL


ullala <- TCGAQuery(tumor = tumor,
          centerType = centerType,
          platform = platform,
          level = level,
          metadata = metadata)


earlyStop <- 20
TCGADownload(earlyStop = earlyStop)




version <- TCGAVersion(tumor = tumor,
                       centerType = centerType,
                       platform = platform,
                       level = level,
                       barcode = T)





