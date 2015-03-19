
#ATTENTION. IN THE SAME FOLDER MUST BE PLACED "data/dataFolders.rda" AND TCGADownloader.tar.gz

#I don't know how to make that data internal, that's the cause of the first constriction. Any help?


#setwd("whatever works")
setwd("Documents/MasterThesis/proj1 - TCGAbiolinks/")
install.packages("TCGADownloader.tar.gz", type = 'source', repos =NULL)
library(TCGADownloader)

#after having updated the platform names matrix using: (I already put the downloaded one)
#TCGAUpdate()
#
#you are ready to go. 


x <- TCGAQuery(tumor = "acc", platform = "genome_wide_snp_6",centerType = "cgcc",level=3)

TCGADownload(x,20) #earlyStop provided
