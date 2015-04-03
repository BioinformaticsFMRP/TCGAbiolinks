detach("package:TCGAbiolinks", unload=TRUE)
remove.packages("TCGAbiolinks")
setwd("~/Desktop/package")
install.packages('TCGAbiolinks_1.0.tar.gz', type='source', repos = NULL)

library('TCGAbiolinks')
origDir <-"~/Desktop/package/"
setwd(origDir)


#Load info file
load('PlatformAndAssociatedData_WebSites3.RData')
sdrfFolder <- "~/Desktop/package/TCGAsdrf/"
#Download sdrf
TCGAmanifest('LIHC')
#Query
res <- TCGAQuery('LIHC', sdrfFolder, PlatformAndAssociatedData)
toDown <- intersect(res$humanmethylation450, res$illuminahiseq_rnaseqv2)

#Download
#test 1 sample
#~ TCGADownload(Tumor='BRCA', Type='mRNA', Species='RNASeqV2', PlatformAndAssociatedData, downloadFolder=paste(getwd(), 'Test_TCGAbiolinks/', sep='/'), PlatformType='illuminahiseq_rnaseqv2', nsample=1)
#test All sample with Expression + Methylation
Tumor <- 'BRCA'
downloadFolder <- "/Users/antoniocolaprico/Desktop/package"
sdrfFolder <- "/Users/antoniocolaprico/Desktop/package/TCGAsdrf/"
PlatformType <- 'humanmethylation450'
toDown <- res$humanmethylation450
listSample <- toDown

BRCA_rnaseqV2_common_Methylation <- TCGADownload(Tumor='BRCA', PlatformAndAssociatedData, downloadFolder="/Users/antoniocolaprico/Desktop/package", sdrfFolder = "/Users/antoniocolaprico/Desktop/package/TCGAsdrf/", PlatformType='humanmethylation450', listSample=toDown,newsample = F)

# Only to use after stopped first esecution with newsample = F
#BRCA_rnaseqV2_common_Methylation <- TCGADownload(Tumor='BRCA', PlatformAndAssociatedData, downloadFolder="/Users/antoniocolaprico/Desktop/package", sdrfFolder = "/Users/antoniocolaprico/Desktop/package/TCGAsdrf/", PlatformType='illuminahiseq_rnaseqv2', listSample=toDown,newsample = T)

BRCA_methylation450_version <- TCGAVersion(Tumor, PlatformType,PlatformAndAssociatedData)
     
