
# it downloads illuminahiseq_rnaseq data of breast cancer samples and stores in downloadFolder
\dontrun{
  Tumor <- "BRCA"
  downloadFolder <- "~/Data"
  setwd(downloadFolder)
  Type = "mRNA"
  Species <- "RNASeq"   
  PlatformType <- "illuminahiseq_rnaseq"
  BRCA_rnaseq <- TCGADownload(Tumor, 
                              Type, 
                              Species, 
                              PlatformAndAssociatedData,
                              downloadFolder,
                              PlatformType)
}

# it downloads illuminahiseq_rnaseqv2 data of breast cancer samples and stores in downloadFolder

\dontrun{
  Type = "mRNA"
  Species <- "RNASeqV2" 
  PlatformType <- "illuminahiseq_rnaseqv2" 
  BRCA_rnaseqv2 <- TCGADownload(Tumor, 
                                Type, 
                                Species, 
                                PlatformAndAssociatedData, 
                                downloadFolder,
                                PlatformType)
}
# it downloads agilentg4502a_07_3 data of breast cancer samples and stores in downloadFolder

\dontrun{
  Type = "mRNA"
  Species <- "Exp-Gene" 
  PlatformType <- "agilentg4502a_07_3"
  BRCA_agilent <- TCGADownload(Tumor, 
                               Type, 
                               Species, 
                               PlatformAndAssociatedData, 
                               downloadFolder,
                               PlatformType)
}

# it downloads humanmethylation27 data of breast cancer samples and stores in downloadFolder
\dontrun{
  Type = "Methylation"
  Species <- "Methyl"
  PlatformType <- "humanmethylation27" 
  BRCA_methylation27 <- TCGADownload(Tumor, 
                                     Type, 
                                     Species, 
                                     PlatformAndAssociatedData, 
                                     downloadFolder,
                                     PlatformType)
}

# it downloads humanmethylation450 data of breast cancer samples and stores in downloadFolder
\dontrun{
  Type = "Methylation"
  Species <- "Methyl"
  PlatformType <- "humanmethylation450" 
  BRCA_methylation450 <- TCGADownload(Tumor, 
                                      Type, 
                                      Species, 
                                      PlatformAndAssociatedData, 
                                      downloadFolder,
                                      PlatformType)
}

# it downloads illuminaga_mirnaseq data of breast cancer samples and stores in downloadFolder
\dontrun{
  Species <- "miRNASeq"
  Type = "miRNA"
  PlatformType <- "illuminaga_mirnaseq"
  BRCA_ga_mirnaseq <- TCGADownload(Tumor, 
                                   Type, 
                                   Species, 
                                   PlatformAndAssociatedData, 
                                   downloadFolder,
                                   PlatformType)
}

# it downloads illuminahiseq_mirnaseq data of breast cancer samples and stores in downloadFolder
\dontrun{
  Type = "miRNA"
  Species <- "miRNASeq"
  PlatformType <- "illuminahiseq_mirnaseq"
  BRCA_hiseq_mirnaseq <- TCGADownload(Tumor, 
                                      Type, 
                                      Species, 
                                      PlatformAndAssociatedData, 
                                      downloadFolder,
                                      PlatformType)
}

# it downloads genome_wide_snp_6 data of breast cancer samples and stores in downloadFolder
\dontrun{
  Type = "SNP"
  Species <- "CNV (SNP Array)"
  PlatformType <- "genome_wide_snp_6"  
  BRCA_genome_wide_snp_6 <- TCGADownload(Tumor, 
                                         Type, 
                                         Species, 
                                         PlatformAndAssociatedData, 
                                         downloadFolder,
                                         PlatformType)
}

# it downloads illuminahiseq_dnaseqc data of breast cancer samples and stores in downloadFolder
\dontrun{
  Species <- "CNV (Low Pass DNASeq)"
  PlatformType <- "illuminahiseq_dnaseqc"  
  BRCA_hiseq_dnaseqc <- TCGADownload(Tumor, 
                                     Type, 
                                     Species, 
                                     PlatformAndAssociatedData, 
                                     downloadFolder,
                                     PlatformType)
}

# it downloads mda_rppa_core data of breast cancer samples and stores in downloadFolder
\dontrun{
  Type = "Protein"
  Species <- "Exp-Protein"
  PlatformType <- "mda_rppa_core"
  BRCA_protein <- TCGADownload(Tumor, 
                               Type, 
                               Species, 
                               PlatformAndAssociatedData, 
                               downloadFolder,
                               PlatformType)
}

# it downloads illuminaga_dnaseq data of breast cancer samples and stores in downloadFolder
\dontrun{
  Type = "Exome"
  Species <- "Somatic Mutation"
  PlatformType <- "illuminaga_dnaseq"
  BRCA_mutation <- TCGADownload(Tumor, 
                                Type, 
                                Species, 
                                PlatformAndAssociatedData, 
                                downloadFolder,
                                PlatformType)
}
