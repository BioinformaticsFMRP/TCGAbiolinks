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
library('stringr')
library("XML")

#update the dataFolders.rda
#TCGAUpdate()



tumor = "acc"
centerType = "cgcc"
platform = "genome_wide_snp_6"
level = 3
metadata = T

# tumor = "lihc"
# platform = "mixed_dnaseq_curated"
# centerType = "gsc"
# level=2
# metadata = F

TCGAQuery(tumor = tumor,
          centerType = centerType,
          level = level,
          platform = platform,
          metadata = metadata)


earlyStop <- 20
TCGADownload(earlyStop = earlyStop)

