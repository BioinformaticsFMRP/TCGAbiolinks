#' @title Query GDC data
#' @description
#'   Uses GDC API to search for search, it searches for both controlled and
#'   open-acess data.
#'   For GDC data arguments project, data.category, data.type and workflow.type should be used
#'   For the legacy data arguments project, data.category, platform and/or file.extension should be used.
#' @param project A valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param data.category A valid project (see list with TCGAbiolinks:::getProjectSummary(project))
#' @param data.type A data type to filter the files to download
#' @param sample.type A sample type to filter the files to download
#' @param barcode A list of barcodes to filter the files to download
#' @param legacy Search in the legacy repository
#' @param file.extension To be used in the legacy database for some platforms,
#' to define which file types to be used.
#' @param workflow.type GDC workflow type
#' @param platform Example:
#' \tabular{ll}{
#'CGH- 1x1M_G4447A                   \tab IlluminaGA_RNASeqV2   \cr
#'AgilentG4502A_07                  \tab IlluminaGA_mRNA_DGE    \cr
#'Human1MDuo                        \tab HumanMethylation450    \cr
#'HG-CGH-415K_G4124A                \tab IlluminaGA_miRNASeq    \cr
#'HumanHap550                       \tab IlluminaHiSeq_miRNASeq \cr
#'ABI                               \tab H-miRNA_8x15K  \cr
#'HG-CGH-244A                       \tab SOLiD_DNASeq                \cr
#'IlluminaDNAMethylation_OMA003_CPI \tab IlluminaGA_DNASeq_automated   \cr
#'IlluminaDNAMethylation_OMA002_CPI \tab HG-U133_Plus_2                 \cr
#'HuEx- 1_0-st-v2 \tab Mixed_DNASeq                  \cr
#'H-miRNA_8x15Kv2 \tab IlluminaGA_DNASeq_curated      \cr
#'MDA_RPPA_Core   \tab IlluminaHiSeq_TotalRNASeqV2    \cr
#'HT_HG-U133A     \tab IlluminaHiSeq_DNASeq_automated \cr
#'diagnostic_images                 \tab microsat_i                     \cr
#'IlluminaHiSeq_RNASeq              \tab SOLiD_DNASeq_curated           \cr
#'IlluminaHiSeq_DNASeqC             \tab Mixed_DNASeq_curated           \cr
#'IlluminaGA_RNASeq                 \tab IlluminaGA_DNASeq_Cont_automated  \cr
#'IlluminaGA_DNASeq                 \tab IlluminaHiSeq_WGBS             \cr
#'pathology_reports                 \tab IlluminaHiSeq_DNASeq_Cont_automated\cr
#'Genome_Wide_SNP_6                 \tab bio                            \cr
#'tissue_images                     \tab Mixed_DNASeq_automated         \cr
#'HumanMethylation27                \tab Mixed_DNASeq_Cont_curated      \cr
#'IlluminaHiSeq_RNASeqV2            \tab Mixed_DNASeq_Cont
#'}
#' @export
#' @return A data frame with the results and the parameters used
GDCquery <- function(project,
                     data.category,
                     data.type,
                     workflow.type,
                     legacy = FALSE,
                     platform,
                     file.extension,
                     barcode,
                     sample.type){

    # Check arguments
    checkProjectInput(project)
    checkDataCategoriesInput(project, data.category, legacy)
    if(!missing(sample.type)) checkBarcodeDefinition(sample.type)
    if(!legacy & !missing(platform)) message("Platform information is only available for legacy database. It will be ignored")

    # Get manifest using the API
    baseURL <- ifelse(legacy,"https://gdc-api.nci.nih.gov/legacy/files/?","https://gdc-api.nci.nih.gov/files/?")
    options.pretty <- "pretty=true"
    if(data.category == "Protein expression" & legacy) {
        options.expand <- "expand=cases.samples.portions,cases.project,center,analysis"
    } else {
        options.expand <- "expand=cases.samples.portions.analytes.aliquots,cases.project,center,analysis"
    }
    option.size <- paste0("size=",getNbFiles(project,data.category,legacy))
    options.filter <- paste0("filters=",
                             URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["'),
                             project,
                             URLencode('"]}},{"op":"in","content":{"field":"files.data_category","value":["'),
                             URLencode(data.category),
                             URLencode('"]}}]}'))
    print(paste0(baseURL,paste(options.pretty, options.expand, option.size, options.filter, sep = "&")))
    json <- fromJSON(paste0(baseURL,paste(options.pretty, options.expand,option.size, options.filter, sep = "&")), simplifyDataFrame = TRUE)

    results <- json$data$hits


    # get barcode of the samples
    # TARGET-20-PANLLX-09A-01R
    #print(results$cases[[1]])
    pat <- paste("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}-[:alnum:]{3}-[:alnum:]{2,3}-[:alnum:]{4}-[:alnum:]{2}",
                 "[:alnum:]{6}-[:alnum:]{2}-[:alnum:]{6}-[:alnum:]{3}-[:alnum:]{3}",sep = "|")
    barcodes <- na.omit(unlist(lapply(results$cases,function(x) str_extract(x,pat))))
    results$cases <- barcodes
    results$definition <- expandBarcodeInfo(barcodes)$definition

    if(legacy & !missing(platform)){
        if(!(platform %in% results$platform)) {
            stop("Please set a valid platform argument from the list below:\n  => ", paste(unique(results$platform), collapse = "\n  => "))
        }
        results <- results[results$platform %in% platform,]
    }

    # Filter by sample.type
    if(!missing(sample.type)) {
        results <- results[results$definition %in% sample.type,]
    }
    # Filter by barcode
    if(!missing(barcode)) {
        results <- results[substr(results$cases,1,str_length(barcode[1])) %in% barcode,]
    }

    # Filter by data.type
    if(!missing(data.type)) {
        if(!(data.type %in% results$data_type)) {
            stop("Please set a valid data.type argument from the list below:\n  => ", paste(unique(results$data_type), collapse = "\n  => "))
        }
        results <- results[results$data_type %in% data.type,]
    }

    # Filter by data.type
    if(!missing(workflow.type)) {
        if(!(workflow.type %in% results$analysis$workflow_type)) {
            stop("Please set a valid data.type argument from the list below:\n  => ", paste(unique(results$analysis$workflow_type), collapse = "\n  => "))
        }
        results <- results[results$analysis$workflow_type %in% workflow.type,]
    }

    # prepare output
    if(missing(sample.type)) sample.type <- NA
    if(missing(data.type)) data.type <- NA
    if(missing(barcode)) barcode <- NA
    if(missing(platform)) platform <- NA
    ret <- data.frame(results=I(list(results)), project = project,
                      data.category = data.category, data.type = data.type, legacy = legacy, platform = I(list(platform)),
                      sample.type = I(list(sample.type)), barcode = I(list(barcode)))
    return(ret)
}

expandBarcodeInfo <- function(barcode){
    ret <- data.frame(barcode = barcode,
                      patient = substr(barcode, 1, 12),
                      sample = substr(barcode, 1, 16),
                      code = substr(barcode, 14, 15))
    ret <- merge(ret,getBarcodeDefinition(), by = "code", sort = FALSE)
    ret <- ret[match(barcode,ret$barcode),]
    return(ret)
}


getBarcodeDefinition <- function(){
    code <- c('01','02','03','04','05','06','07','08','09','10','11',
              '12','13','14','20','40','50','60','61')
    shortLetterCode <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC",
                         "TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB",
                         "CELL","XP","XCL")

    definition <- c("Primary solid Tumor",
                    "Recurrent Solid Tumor",
                    "Primary Blood Derived Cancer - Peripheral Blood",
                    "Recurrent Blood Derived Cancer - Bone Marrow",
                    "Additional - New Primary",
                    "Metastatic",
                    "Additional Metastatic",
                    "Human Tumor Original Cells",
                    "Primary Blood Derived Cancer - Bone Marrow",
                    "Blood Derived Normal",
                    "Solid Tissue Normal",
                    "Buccal Cell Normal",
                    "EBV Immortalized Normal",
                    "Bone Marrow Normal",
                    "Control Analyte",
                    "Recurrent Blood Derived Cancer - Peripheral Blood",
                    "Cell Lines",
                    "Primary Xenograft Tissue",
                    "Cell Line Derived Xenograft Tissue")
    aux <- data.frame(code = code,shortLetterCode,definition)
    return(aux)
}

#' @title Searches TCGA open-access data providing also latest version of the files.
#' @description
#'    This function has been replace by GDCquery
#' @param tumor Disease Examples:
#' \tabular{lllll}{
#'OV   \tab BRCA \tab CESC \tab ESCA \tab PCPG\cr
#'LUSC \tab LGG  \tab SKCM \tab KICH \tab CHOL\cr
#'GBM  \tab UCEC \tab PRAD \tab PAAD \tab THYM\cr
#'KIRC \tab THCA \tab SARC \tab LAML \tab TGCT\cr
#'COAD \tab KIRP \tab HNSC \tab ACC  \tab UVM \cr
#'READ \tab BLCA \tab DLBC \tab UCS  \tab FPPP\cr
#'LUAD \tab LIHC \tab STAD \tab MESO \tab CNTL
#'}
#' @param platform Example:
#' \tabular{ll}{
#'CGH- 1x1M_G4447A                   \tab IlluminaGA_RNASeqV2   \cr
#'AgilentG4502A_07                  \tab IlluminaGA_mRNA_DGE    \cr
#'Human1MDuo                        \tab HumanMethylation450    \cr
#'HG-CGH-415K_G4124A                \tab IlluminaGA_miRNASeq    \cr
#'HumanHap550                       \tab IlluminaHiSeq_miRNASeq \cr
#'ABI                               \tab H-miRNA_8x15K  \cr
#'HG-CGH-244A                       \tab SOLiD_DNASeq                \cr
#'IlluminaDNAMethylation_OMA003_CPI \tab IlluminaGA_DNASeq_automated   \cr
#'IlluminaDNAMethylation_OMA002_CPI \tab HG-U133_Plus_2                 \cr
#'HuEx- 1_0-st-v2 \tab Mixed_DNASeq                  \cr
#'H-miRNA_8x15Kv2 \tab IlluminaGA_DNASeq_curated      \cr
#'MDA_RPPA_Core   \tab IlluminaHiSeq_TotalRNASeqV2    \cr
#'HT_HG-U133A     \tab IlluminaHiSeq_DNASeq_automated \cr
#'diagnostic_images                 \tab microsat_i                     \cr
#'IlluminaHiSeq_RNASeq              \tab SOLiD_DNASeq_curated           \cr
#'IlluminaHiSeq_DNASeqC             \tab Mixed_DNASeq_curated           \cr
#'IlluminaGA_RNASeq                 \tab IlluminaGA_DNASeq_Cont_automated  \cr
#'IlluminaGA_DNASeq                 \tab IlluminaHiSeq_WGBS             \cr
#'pathology_reports                 \tab IlluminaHiSeq_DNASeq_Cont_automated\cr
#'Genome_Wide_SNP_6                 \tab bio                            \cr
#'tissue_images                     \tab Mixed_DNASeq_automated         \cr
#'HumanMethylation27                \tab Mixed_DNASeq_Cont_curated      \cr
#'IlluminaHiSeq_RNASeqV2            \tab Mixed_DNASeq_Cont
#'}
#' @param level '1' '2' '3'
#' @param center center name
#' @param samples List of samples. Ex:c('TCGA-04-06','TCGA-61-1743-01A-01D-0649-04')
#' @param version List of vector with tumor/plaform/version to get old samples,
#' @examples
#' \dontrun{
#' query <- TCGAquery(tumor = "gbm")
#' }
#' @export
#' @return A dataframe with the results of the query
#'        (lastest version of the files)
#' @family data functions
TCGAquery <- function(tumor = NULL,
                      platform = NULL,
                      samples = NULL,
                      center = NULL,
                      level = NULL,
                      version = NULL) {
    stop("TCGA data moved from DCC server to GDC server. \n Please use the function GDCquery")
}

#' @title Find most studied TF in pubmed related to a specific cancer, disease, or tissue
#' @description Find most studied TF in pubmed related to a specific cancer,
#'  disease, or tissue
#' @param tumor is character such as cancer, disease, or tissue eg. BRCA or breast
#' @param dataDEGsFiltLevelTF is a table output from TCGAanalyze_LevelTab with only TFs
#' @param topgenes is the number of top genes (eg. 10) in the rownames(dataDEGsFiltLevelTF)
#' where find in pubmed if those genes or TFs are already related to that cancer or disease
#' @importFrom RCurl url.exists curlVersion
#' @examples
#' query <- TCGAquery(tumor = "lgg")
#' \dontrun{
#' TFs <- EAGenes[EAGenes$Family =="transcription regulator",]
#' TFs_inDEGs <- intersect(TFs$Gene, dataDEGsFiltLevel$mRNA )
#' dataDEGsFiltLevelTFs <- dataDEGsFiltLevel[TFs_inDEGs,]
#  # Order table DEGs TFs according to Delta decrease
#' dataDEGsFiltLevelTFs <- dataDEGsFiltLevelTFs[order(dataDEGsFiltLevelTFs$Delta,decreasing = TRUE),]
#' # Find Pubmed of TF studied related to cancer
#' tabDEGsTFPubmed <- TCGAquery_Investigate("breast", dataDEGsFiltLevelTFs, topgenes = 1)
#' }
#' @export
#' @return table with number of pubmed's publications related to tfs and disease selected
TCGAquery_Investigate <- function(tumor,dataDEGsFiltLevelTF,topgenes){
    site <- "http://www.ncbi.nlm.nih.gov/pubmed/?term="

    # GenesTofix <- c("JUN","HR","HOMEZ",
    #                 "ANKAR", "REST", "BATF", "MAX", "ECD", "FOS")
    # dataDEGsFiltLevelTF <- dataDEGsFiltLevelTF[
    #                        setdiff(dataDEGsFiltLevelTFs$mRNA,GenesTofix),]

    dataDEGsFiltLevelTF <- dataDEGsFiltLevelTF[1:topgenes,]
    Pubmed <- matrix(0, nrow(dataDEGsFiltLevelTF), 1)
    PMID <- matrix(0, nrow(dataDEGsFiltLevelTF), 1)


    dataDEGsFiltLevelTF <- cbind(dataDEGsFiltLevelTF,Pubmed,PMID)
    dataDEGsFiltLevelTF<-as.data.frame(dataDEGsFiltLevelTF)

    for (k in 1:nrow( dataDEGsFiltLevelTF)){
        CurrentGene <- dataDEGsFiltLevelTF$mRNA[k]
        site2 <- paste(site,CurrentGene, "+", tumor,sep = "")

        if(interactive() && ("ssl" %in% names(curlVersion()$features)) &&
           url.exists(site2)) {
            x = tryCatch(RCurl::getURL(site2), error = function(e) {
                RCurl::getURL(site2, ssl.verifypeer = FALSE) })
        }



        if ( length(grep("No items found.",x)) != 1) {

            if (length(grep("Display Settings",x)) == 1) {
                x6 <- 1
                dataDEGsFiltLevelTF[k,"PMID"] <- substr(
                    gsub("</dt> <dd>","",unlist(strsplit(x,"PMID:"))[2]),1,8
                )

            }


            if (length(grep("result_count",x))==1){
                x2a <- unlist(strsplit(x,"result_count"))[2]

                tmpPMID2 <- unlist(strsplit(x2a,"UidCheckBox"))
                tmpPMID3 <- tmpPMID2[grep("<span>",tmpPMID2)]
                dataDEGsFiltLevelTF[k,"PMID"] <- as.character(paste(
                    substr(tmpPMID3,1,8),collapse="; "))


                x3a <-  unlist(strsplit(x2a,"</h2>"))[1]

                if( length(grep("of",x3a))!=1){
                    x6 <- as.numeric(unlist(strsplit(x3a,": "))[2])
                } else { x6 <- as.numeric(unlist(strsplit(x3a,"of "))[2]) }
            }


            if (length(grep("following term was not found",x)) == 1) { x6 <- 0 }

            if (length(grep("Search instead for",x)) == 1) { x6 <- 1 }

            if (CurrentGene == "JUN")   { x6 <- 0 }
            if (CurrentGene == "HR")    { x6 <- 0 }
            if (CurrentGene == "HOMEZ") { x6 <- 0 }
            if (CurrentGene == "ANKAR") { x6 <- 0 }
            if (CurrentGene == "REST")  { x6 <- 0 }
            if (CurrentGene == "BATF")  { x6 <- 0 }
            if (CurrentGene == "MAX")   { x6 <- 0 }
            # if (CurrentGene =="FOS"){       x6 <- 0     }
            if (CurrentGene == "ECD")   { x6 <- 0 }
            # HR, HOMEZ, ANKAR, REST

            dataDEGsFiltLevelTF[k,"Pubmed"] <- x6
            print(paste("Cancer ", tumor, "with TF n. ",k, "of " ,
                        nrow( dataDEGsFiltLevelTF)," : ", CurrentGene,
                        "found n. ", x6, "pubmed."))

        }

        else{
            print(paste("Cancer ", tumor, "with TF n. ",k, "of " ,
                        nrow( dataDEGsFiltLevelTF)," : ", CurrentGene,
                        "no item found in pubmed."))
            dataDEGsFiltLevelTF[k,"Pubmed"] <- 0
        }

    }

    dataDEGsFiltLevelTF <- dataDEGsFiltLevelTF[order(dataDEGsFiltLevelTF$Pubmed,
                                                     decreasing = TRUE),]

    if( sum(dataDEGsFiltLevelTF$Pubmed == 1) != 0) {
        dataDEGsFiltLevelTF[dataDEGsFiltLevelTF$Pubmed == 1,][
            which( nchar(dataDEGsFiltLevelTF[
                dataDEGsFiltLevelTF$Pubmed == 1,]$PMID) > 8),"PMID"] <- substr(
                    dataDEGsFiltLevelTF[dataDEGsFiltLevelTF$Pubmed == 1,][
                        which( nchar(dataDEGsFiltLevelTF[
                            dataDEGsFiltLevelTF$Pubmed == 1,]$PMID) > 8),
                        "PMID"],1,8)
    }
    if ( sum(dataDEGsFiltLevelTF$Pubmed == 0) != 0) {
        dataDEGsFiltLevelTF[dataDEGsFiltLevelTF$Pubmed == 0,"PMID"] <- 0
    }

    tabDEGsTFPubmed$Tumor <- round(tabDEGsTFPubmed$Tumor)
    tabDEGsTFPubmed$Normal <- round(tabDEGsTFPubmed$Normal)
    tabDEGsTFPubmed$Delta <- round(tabDEGsTFPubmed$Delta)

    return(dataDEGsFiltLevelTF)
}



#' @title Retrieve open access maf files from GDC server
#' @description
#'   GDCquery_Maf uses the following guide to download maf files
#'   https://gdc-docs.nci.nih.gov/Data/Release_Notes/Data_Release_Notes/
#' @param tumor a valid tumor
#' @param save.csv Write maf file into a csv document
#' @export
#' @importFrom data.table fread
#' @import readr stringr
#' @importFrom downloader download
#' @importFrom R.utils gunzip
#' @importFrom tools md5sum
#' @examples
#' acc.maf <- GDCquery_Maf("ACC")
#' @return A data frame with the maf file information
GDCquery_Maf <- function(tumor, save.csv= FALSE){
    root <- "https://gdc-api.nci.nih.gov/data/"
    maf <- fread("https://gdc-docs.nci.nih.gov/Data/Release_Notes/Manifests/GDC_open_MAFs_manifest.txt",
                 data.table = FALSE, verbose = FALSE, showProgress = FALSE)
    maf$tumor <- unlist(lapply(maf$filename, function(x){unlist(str_split(x,"\\."))[2]}))

    # Check input
    if (missing(tumor)) stop(paste0("Please, set tumor argument. Possible values:\n => ",
                                    paste(sort(maf$tumor),collapse = "\n => ")))

    if (!(tumor %in%  maf$tumor)) stop(paste0("Please, set a valid tumor argument. Possible values:\n => ",
                                              paste(sort(maf$tumor),collapse = "\n => ")))

    #  Info to user
    message("============================================================================")
    message(" For more information about MAF data please read the following GDC manual:")
    message(" GDC manual: https://gdc-docs.nci.nih.gov/Data/PDF/Data_UG.pdf")
    message("============================================================================")
    selected <- maf[grepl(tumor,maf$tumor,ignore.case = TRUE),]

    # Download maf
    repeat{
        if (!file.exists(selected$filename)) download(file.path(root,selected$id),selected$filename)

        # check integrity
        if(md5sum(selected$filename) == selected$md5) break
        message("The data downloaded might be corrupted. We will download it again")
    }

    # uncompress file
    uncompressed <- gsub(".gz","",selected$filename)
    if (!file.exists(uncompressed)) gunzip(selected$filename, remove = FALSE)

    ret <- read_tsv(uncompressed,comment = "#")

    if(save.csv) write_csv(ret,gsub("txt","csv",uncompressed))

    return(ret)
}

#' @title Get last maf file for the tumor
#' @description
#'    This function has been replaced. Use GDCquery_maf
#' @param tumor tumor type to filter the search
#' @param center Center name to filter the search
#' @param archive.name Archive name to filter the search
#' @examples
#' \dontrun{
#'  query <- TCGAquery_maf(tumor = 'lgg')
#' }
#' @export
#' @return list of samples for a tumor
TCGAquery_maf <- function(tumor = NULL, center = NULL, archive.name = NULL){
    stop("TCGA data has moved from DCC server to GDC server. Please use GDCquery_maf function")
}
