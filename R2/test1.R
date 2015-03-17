source("R2/TCGADownload.R")
source("R2/TCGAQuery.R")
source("R2/TCGAInternal.R")
source("R2/TCGAUpdate.R")

if(require(httr) == F) install.packages("httr")
if(require(downloader) == F) install.packages("downloader")

#TCGAUpdate()

#after having updated the platform names matrix using:
#TCGAUpdate()
#
#you are ready to go. (I already put the downloaded one)


x <- TCGAQuery(tumor = "acc", platform = "genome_wide_snp_6",centerType = "cgcc",level=3)

TCGADownload(x,20) #earlyStop provided. Something about the stability of the connection and the ongoing situation must be done

