#' @title Retrieve multiple tissue types not from the same patients.
#' @description
#'   TCGAquery_SampleTypes for a given list of samples and types,
#'    return the union of samples that are from theses type.
#' @param barcode is a list of samples as TCGA barcodes
#' @param typesample a character vector indicating tissue type to query.
#' Example:
#' \tabular{ll}{
#'TP \tab   PRIMARY SOLID TUMOR \cr
#'TR \tab   RECURRENT SOLID TUMOR \cr
#'TB \tab   Primary Blood Derived Cancer-Peripheral Blood \cr
#'TRBM \tab Recurrent Blood Derived Cancer-Bone Marrow \cr
#'TAP \tab  Additional-New Primary \cr
#'TM \tab   Metastatic \cr
#'TAM \tab  Additional Metastatic \cr
#'THOC \tab Human Tumor Original Cells \cr
#'TBM \tab  Primary Blood Derived Cancer-Bone Marrow \cr
#'NB \tab   Blood Derived Normal \cr
#'NT \tab   Solid Tissue Normal \cr
#'NBC \tab  Buccal Cell Normal \cr
#'NEBV \tab EBV Immortalized Normal \cr
#'NBM \tab  Bone Marrow Normal \cr
#'}
#' @export
#' @examples
#' # selection of normal samples "NT"
#' barcode <- c("TCGA-B0-4698-01Z-00-DX1","TCGA-CZ-4863-02Z-00-DX1")
#' # Returns the second barcode
#'  TCGAquery_SampleTypes(barcode,"TR")
#'  # Returns both barcode
#'  TCGAquery_SampleTypes(barcode,c("TR","TP"))
#' @return a list of samples / barcode filtered by type sample selected
TCGAquery_SampleTypes <- function(barcode,typesample){
    # Tumor AND Solid Tissue Normal NOT FROM THE SAME PATIENTS
    table.code <- c('01','02','03','04','05','06','07','08','09','10',
                    '11','12','13','14','20','40','50','60','61')
    names(table.code) <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC",
                           "TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB",
                           "CELL","XP","XCL")

    if (sum(is.element(typesample,names(table.code))) == length(typesample)) {

        string <- substr(barcode, 14, 15)
        barcode.all <- NULL
        for (sample.i in typesample) {
            barcode.all <- union(barcode.all,
                                 barcode[grep(table.code[sample.i], string)])
        }
        return(barcode.all)
    }else{
        return("Error message: one or more sample types do not exist")
    }
}

#' @title Retrieve multiple tissue types from the same patients.
#' @description
#'   TCGAquery_MatchedCoupledSampleTypes
#' @param barcode barcode
#' @param typesample typesample
#' @examples
#'  TCGAquery_MatchedCoupledSampleTypes(c("TCGA-B0-4698-01Z-00-DX1",
#'                              "TCGA-B0-4698-02Z-00-DX1"),
#'                              c("TP","TR"))
#' @export
#' @return a list of samples / barcode filtered by type sample selected
TCGAquery_MatchedCoupledSampleTypes <- function(barcode,typesample){
    # Tumor AND Solid Tissue Normal FROM THE SAME PATIENTS
    table.code <- c('01','02','03','04','05','06','07','08','09','10',
                    '11','12','13','14','20','40','50','60','61')
    names(table.code) <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC",
                           "TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB",
                           "CELL","XP","XCL")
    if(length(typesample)!=2){
        return("Error message: exactly two types need to be provided")
    }

    if(sum(is.element(typesample,names(table.code))) == length(typesample)) {

        string <- substr(barcode, 14, 15)
        barcode.1 <- barcode[grep(table.code[typesample[1]], string)]
        barcode.2 <- barcode[grep(table.code[typesample[2]], string)]
        barcode.common <- intersect(substr(barcode.1,1,12), substr(barcode.2,1,12))
        if(length(barcode.common) > 0){
            idx1 <- unlist(lapply(barcode.common, function(x) grep(x,barcode.1)))
            idx2 <- unlist(lapply(barcode.common, function(x) grep(x,barcode.2)))
            return(union(barcode.1[idx1], barcode.2[idx2]))
        } else {
            return("Error message: there exist no matched samples")
        }
    } else {
        return("Error message: one or more sample types do not exist")
    }
}

#' @title stage_BRCA
#' @description
#'   stage_BRCA
#' @param barcode barcode
#' @param stage stage
#' @param clinical_patient_data clinical_patient_data
#' @keywords internal
#' @return stage_BRCA
# @examples
# # clin <- TCGAquery_clinic("BRCA","clinical_patient")
# clin <- clinBRCA
# stage_BRCA(c("TCGA-3C-AALK","TCGA-A2-A04Q","TCGA-A4-A04Q"),"stage_IX",clin)
stage_BRCA <- function(barcode, stage, clinical_patient_data){
    table.stages <- c("Stage I$|Stage IA$|Stage IB$", "Stage I$", "Stage IA$",
                      "Stage IB$", "Stage II$|Stage IIA$|Stage IIB$",
                      "Stage II$", "Stage IIA$", "Stage IIB$",
                      "Stage III$|Stage IIIA$|Stage IIIB$|Stage IIIC$",
                      "Stage III$", "Stage IIIA$", "Stage IIIB$",
                      "Stage IIIC$", "Stage IV$")
    names(table.stages) <- c("stage_IX", "stage_I", "stage_IA", "stage_IB",
                             "stage_IIX", "stage_IIA", "stage_IIB",
                             "stage_IIIX","stage_IIIA", "stage_IIIB",
                             "stage_IIIC", "stage_IV")

    if (is.element(stage, names(table.stages))) {
        clinical_patient_data <- as.data.frame(clinical_patient_data)
        print(table.stages[stage])
        stage.i <- clinical_patient_data[
            grep(table.stages[stage],
                 clinical_patient_data$pathologic_stage), ]
        stage.i <- stage.i[,"bcr_patient_barcode"]
        samples <- substr(barcode, 1, 12)
        barcode <- intersect(samples,stage.i)
    }else{
        return("Error message: stage does not exist")
    }
    return(barcode)
}

#' @title gender_BRCA
#' @description
#'   gender_BRCA
#' @param barcode barcode
#' @param gender gender
#' @param clinical_patient_data clinical_patient_data
#' @keywords internal
#' @return stage_BRCA
# @examples
# # clin <- TCGAquery_clinic("BRCA","clinical_patient")
# clin <- clinBRCA
# gender_BRCA (c("TCGA-3C-AALK","TCGA-A2-A04Q","TCGA-A4-A04Q"),"FEMALE",clin)
gender_BRCA <- function(barcode, gender, clinical_patient_data){

    if (is.element(gender,c("MALE", "FEMALE"))) {
        clinical_patient_data <- as.data.frame(clinical_patient_data)
        s.gender <- as.data.frame(clinical_patient_data)[
            grep(paste0("^", gender,"$"), clinical_patient_data$gender),
            ][,"bcr_patient_barcode"]
        samples <- substr(barcode, 1, 12)
        #find common patients between FEMALE and barcode data
        barcode <- intersect(samples,s.gender)
    }else{
        return("Error message gender doesn't exist")
    }

    return(barcode)
}

#' @title ER_status_BRCA
#' @description
#'   ER_status_BRCA
#' @param barcode barcode
#' @param ER ER
#' @param clinical_patient_data clinical_patient_data
#' @keywords internal
#' @return ER_status_BRCA
# @examples
# # clin <- TCGAquery_clinic("BRCA","clinical_patient")
# clin <- clinBRCA
# ER_status_BRCA(c("TCGA-3C-AALK","TCGA-A2-A04Q","TCGA-A4-A04Q"),
# "Positive",clin)
ER_status_BRCA <- function(barcode,ER, clinical_patient_data){
    ## ER should be "Positive" or "Negative"
    # consider only barcode and ER status
    idx <- grep("estrogen_receptor_status",colnames(clinical_patient_data))
    if(length(idx) > 0) idx <- idx[1]
    if (is.element(ER, c("Positive", "Negative"))) {
        status <- as.data.frame(clinical_patient_data)[
            grep(paste0("^",ER,"$"),
                 clinical_patient_data[,idx]), ][,"bcr_patient_barcode"]
        samples <- substr(barcode, 1, 12)
        #find common patients between ER status and barcode data
        barcode <- intersect(samples,status)
        return(barcode)
    }else{
        return("Error message: ER status does not exist")
    }
}

#' @title PR_status_BRCA
#' @description
#'   PR_status_BRCA
#' @param barcode barcode
#' @param PR PR
#' @param clinical_patient_data clinical_patient_data
#' @keywords internal
#' @return PR_status_BRCA
# @examples
# # clin <- TCGAquery_clinic("BRCA","clinical_patient")
# clin <- clinBRCA
# PR_status_BRCA(c("TCGA-3C-AALK","TCGA-A2-A04Q","TCGA-A4-A04Q"),
# "Positive",clin)
PR_status_BRCA  <- function(barcode,PR, clinical_patient_data){
    ## PR should be "Positive" or "Negative"

    if(is.element(PR, c("Positive", "Negative"))){
        #for breast cancer
        status <- as.data.frame(clinical_patient_data)[
            grep(paste0("^", PR, "$"),
                 clinical_patient_data$pr_status_by_ihc), ][,"bcr_patient_barcode"]

        samples <- substr(barcode, 1, 12)
        #find common patients between PR status and barcode data
        barcode <- intersect(samples,status)
    }else{
        return("Error message: PR status does not exist")
    }

    return(barcode)

}

#' @title HER_status_BRCA
#' @description
#'   HER_status_BRCA
#' @param barcode barcode
#' @param HER HER
#' @param clinical_patient_data clinical_patient_data
#' @keywords internal
#' @return HER_status_BRCA
# @examples
# # clin <- TCGAquery_clinic("BRCA","clinical_patient")
# clin <- clinBRCA
# HER_status_BRCA(c("TCGA-3C-AALK","TCGA-A2-A04Q","TCGA-A4-A04Q"),
# "Positive",clin)
HER_status_BRCA  <- function(barcode, HER, clinical_patient_data){
    if (is.element(HER, c("Positive", "Negative"))) {
        clinical_patient_data <- as.data.frame(clinical_patient_data)
        idx <- grep("immunohistochemistry_receptor_",colnames(clinical_patient_data))
        if(length(idx) > 0) idx <- idx[1]
        #for breast cancer HER+
        status <- as.data.frame(clinical_patient_data)[
            grep(paste0("^",HER,"$"),
                 clinical_patient_data[,idx],), ][,"bcr_patient_barcode"]
        samples <- substr(barcode, 1, 12)
        #find common patients between HER+ e barcode data
        barcode <- intersect(samples,status)
        return(barcode)
    }else{
        return("Error message: HER status does not exist")
    }


}

#' @title clinical_data_site_cancer
#' @description
#'   clinical_data_site_cancer
#' @param cancer cancer
#' @keywords internal
# @examples clinical_data_site_cancer("gbm")
#' @return clinical_data_site_cancer
clinical_data_site_cancer <- function(cancer){
    return(paste0("https://tcga-data.nci.nih.gov/tcgafiles/",
                  "ftp_auth/distro_ftpusers/anonymous/tumor/",
                  cancer,"/bcr/biotab/clin/"))
}



#' @title Get DDC clinical data
#' @description
#'   GDCquery_clinic will download all clinical information from the API
#'   as the one with using the button from each project
#' @param project A valid project (see list with getGDCprojects()$project_id)]
#' @param type A valid type. Options "clinical", "Biospecimen"  (see list with getGDCprojects()$project_id)]
#' @param save.csv Write clinical information into a csv document
#' @export
#' @importFrom data.table rbindlist
#' @importFrom jsonlite fromJSON
#' @examples
#' clin <- GDCquery_clinic("TCGA-ACC", type = "clinical", save.csv = TRUE)
#' clin <- GDCquery_clinic("TCGA-ACC", type = "biospecimen", save.csv = TRUE)
#' @return A data frame with the clinical information
GDCquery_clinic <- function(project, type = "clinical", save.csv = FALSE){
    checkProjectInput(project)
    baseURL <- "https://gdc-api.nci.nih.gov/cases/?"
    options.pretty <- "pretty=true"
    if(type == "clinical"){
        options.expand <- "expand=diagnoses,demographic,exposures"
        option.size <- paste0("size=",getNbCases(project,"Clinical"))
    } else {
        options.expand <- "expand=samples,samples.portions,samples.portions.analytes,samples.portions.analytes.aliquots"
        option.size <- paste0("size=",getNbCases(project,"Biospecimen"))
    }
    options.filter <- paste0("filters=",
                             URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["'),
                             project,
                             URLencode('"]}},{"op":"in","content":{"field":"files.data_category","value":["Biospecimen"]}}]}'))
    json <- fromJSON(paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&")), simplifyDataFrame = TRUE)
    results <- json$data$hits

    if(type == "clinical"){
        diagnoses <- rbindlist(results$diagnoses, fill = TRUE)
        exposures <- rbindlist(results$exposures, fill = TRUE)
        df <- cbind(diagnoses,exposures,results$demographic)
        df$bcr_patient_barcode <- gsub("_diagnosis", "", df$submitter_id)
        df$disease <- gsub("TCGA-|TARGET-", "", project)

    } else {
        df <- rbindlist(results$samples,fill = T)
    }

    #y <- data.frame(diagnosis=I(results$diagnoses), demographic=results$demographic,exposures=I(results$exposures))
    if(save.csv){
        if(type == "biospecimen") df$portions <- NULL
        write_csv(df,paste0(project,"_",type,".csv"))
    }
    setDF(df)
    return(df)
}

#' @title Get the clinical information
#' @description
#'   This function has been replaced. Use GDCquery_clinic
#' @param tumor a character vector indicating cancer type Examples:
#' \tabular{lllll}{
#'OV   \tab BRCA \tab CESC \tab ESCA \tab PCPG\cr
#'LUSC \tab LGG  \tab SKCM \tab KICH \tab CHOL\cr
#'GBM  \tab UCEC \tab PRAD \tab PAAD \tab THYM\cr
#'KIRC \tab THCA \tab SARC \tab LAML \tab TGCT\cr
#'COAD \tab KIRP \tab HNSC \tab ACC  \tab UVM \cr
#'READ \tab BLCA \tab DLBC \tab UCS  \tab FPPP\cr
#'LUAD \tab LIHC \tab STAD \tab MESO \tab CNTL
#'}
#'
#' For information about cancer types: https://tcga-data.nci.nih.gov/tcga/
#' @param clinical_data_type a character vector indicating the types of
#' clinical data. Besides TCGA data, we created the clinical_patient_updated,
#' which is the clinical_patient file with the last follow up information from the last
#' follow up file.
#'  Example:
#' \tabular{ll}{
#' biospecimen_aliquot \tab biospecimen_analyte \cr
#' biospecimen_cqcf \tab biospecimen_diagnostic_slides \cr
#' biospecimen_normal_control \tab biospecimen_portion \cr
#' biospecimen_protocol \tab biospecimen_sample \cr
#' biospecimen_shipment_portion \tab biospecimen_slide \cr
#' biospecimen_tumor_sample \tab clinical_cqcf \cr
#' clinical_follow_up_v1.0 \tab clinical_follow_up_v1.5 \cr
#' clinical_follow_up_v2.0 \tab clinical_follow_up_v2.1 \cr
#' clinical_follow_up_v4.0 \tab clinical_follow_up_v4.0_nte \cr
#' clinical_nte \tab  clinical_omf_v4.0 \cr
#' clinical_patient \tab  clinical_radiation \cr
#' ssf_normal_controls  \tab  ssf_tumor_samples \cr
#' clinical_follow_up_v1.0_nte \cr clinical_patient_updated (TCGAbiolinks only)
#'}
#' @param samples List of barcodes to get the clinical data
#' @param path Directory to save the downloaded data default getwd()
#' @export
#' @return clinic data
#' @examples
#' \dontrun{
#' data <- TCGAquery_clinic("LGG","clinical_drug")
#' }
TCGAquery_clinic <- function(tumor, clinical_data_type, samples, path = getwd()){
    stop("TCGA data has moved from DCC server to GDC server. Please use GDCquery_clinic function")
}


update.clinical.with.last.followup <- function(clin){

    for(disease in unique(clin$disease)){
        aux <- clinical.table[,disease]
        files <- rownames(clinical.table[which(aux==1),])
        # get last follow up files not nte
        files <- sort(files[grepl("follow",files) & !grepl("nte",files)],decreasing = T)[1]

        follow <- TCGAquery_clinic(disease,files)

        colnames(follow) [colnames(follow) %in% colnames(clin)]
        aux <- plyr::ddply(follow, .(bcr_patient_barcode), function(x) x[c(nrow(x)), ])
        aux <- aux[aux$bcr_patient_barcode %in% clin$bcr_patient_barcode,colnames(aux) %in% colnames(clin)]
        clin[na.omit(match(aux$bcr_patient_barcode,clin$bcr_patient_barcode)),match(colnames(aux),colnames(clin))] <- aux
    }
    return(clin)
}

#' @title Filter samples using clinical data
#' @description
#'   This function will return the samples that matches all filters.
#'   Filters available: HER, ER,gender,PR, stage.
#' @param barcode List of barcodes
#' @param clinical_patient_data clinical_patient_data obtained with clinic function
#' Ex: clinical_patient_data <- TCGAquery_clinic("LGG","clinical_patient")
#' @param HER  her2 neu immunohistochemistry receptor status: "Positive" or "Negative"
#' @param gender "MALE" or "FEMALE"
#' @param PR  Progesterone receptor status: "Positive" or "Negative"
#' @param stage Pathologic Stage: "stage_IX", "stage_I", "stage_IA", "stage_IB", "stage_IIX",
#' "stage_IIA", "stage_IIB", "stage_IIIX","stage_IIIA", "stage_IIIB",
#' "stage_IIIC", "stage_IV" -
#' @param ER Estrogen receptor status: "Positive" or "Negative"
#' @export
#' @return List of samples that matches the filters
#' @examples
#' # clin <- TCGAquery_clinic("BRCA","clinical_patient")
#' clin <- clinBRCA
#' bar <- c("TCGA-G9-6378-02A-11R-1789-07", "TCGA-CH-5767-04A-11R-1789-07",
#'         "TCGA-G9-6332-60A-11R-1789-07", "TCGA-G9-6336-01A-11R-1789-07",
#'         "TCGA-G9-6336-11A-11R-1789-07", "TCGA-G9-7336-11A-11R-1789-07",
#'         "TCGA-G9-7336-04A-11R-1789-07", "TCGA-G9-7336-14A-11R-1789-07",
#'         "TCGA-G9-7036-04A-11R-1789-07", "TCGA-G9-7036-02A-11R-1789-07",
#'         "TCGA-G9-7036-11A-11R-1789-07", "TCGA-G9-7036-03A-11R-1789-07",
#'         "TCGA-G9-7036-10A-11R-1789-07", "TCGA-BH-A1ES-10A-11R-1789-07",
#'         "TCGA-BH-A1F0-10A-11R-1789-07", "TCGA-BH-A0BZ-02A-11R-1789-07",
#'         "TCGA-B6-A0WY-04A-11R-1789-07", "TCGA-BH-A1FG-04A-11R-1789-08",
#'         "TCGA-D8-A1JS-04A-11R-2089-08", "TCGA-AN-A0FN-11A-11R-8789-08",
#'         "TCGA-AR-A2LQ-12A-11R-8799-08", "TCGA-AR-A2LH-03A-11R-1789-07",
#'         "TCGA-BH-A1F8-04A-11R-5789-07", "TCGA-AR-A24T-04A-55R-1789-07",
#'         "TCGA-AO-A0J5-05A-11R-1789-07", "TCGA-BH-A0B4-11A-12R-1789-07",
#'         "TCGA-B6-A1KN-60A-13R-1789-07", "TCGA-AO-A0J5-01A-11R-1789-07",
#'         "TCGA-AO-A0J5-01A-11R-1789-07", "TCGA-G9-6336-11A-11R-1789-07",
#'         "TCGA-G9-6380-11A-11R-1789-07", "TCGA-G9-6380-01A-11R-1789-07",
#'         "TCGA-G9-6340-01A-11R-1789-07","TCGA-G9-6340-11A-11R-1789-07")
#'
#' TCGAquery_clinicFilt(c("TCGA-3C-AALK","TCGA-A2-A04Q","TCGA-A4-A04Q"),clin,
#' HER="Positive", gender="FEMALE",ER="Positive")
TCGAquery_clinicFilt <- function(barcode,
                                 clinical_patient_data,
                                 HER=NULL,
                                 ER=NULL,
                                 gender=NULL,
                                 PR=NULL,
                                 stage=NULL){

    x <- NULL

    if (!is.null(PR)) {
        res.pr <- PR_status_BRCA(barcode,PR, clinical_patient_data)
        message(paste0("PR ",PR," Samples:"))
        message(paste(paste("\t",res.pr,"\n")))

        if (is.null(x)) x <- res.pr
    }

    if (!is.null(ER)) {
        res.er <- ER_status_BRCA(barcode,ER, clinical_patient_data)
        message(paste0("ER ",ER," Samples:"))
        message(paste(paste("\t",res.er,"\n")))

        if (is.null(x)) {
            x <- res.er
        } else {
            x <- intersect(x,res.er)
        }
    }

    if (!is.null(HER)) {
        res.her <- HER_status_BRCA(barcode,HER, clinical_patient_data)
        message(paste0("HER ",HER," Samples:"))
        message(paste(paste("\t",res.her,"\n")))
        if (is.null(x)) {
            x <- res.her
        } else {
            x <- intersect(x,res.her)
        }

    }
    if (!is.null(stage)) {
        res.stage <- stage_BRCA(barcode,stage, clinical_patient_data)
        message(paste0("Stage ",stage," Samples:"))
        message(paste(paste("\t",res.stage,"\n")))
        if (is.null(x)) {
            x <- res.stage
        } else {
            x <- intersect(x,res.stage)
        }

    }
    if (!is.null(gender)) {
        res.gender <- gender_BRCA(barcode,gender, clinical_patient_data)
        message(paste0("GENDER ",gender," Samples:"))
        message(paste(paste("\t",res.gender,"\n")))
        if (is.null(x)) {
            x <- res.gender
        } else {
            x <- intersect(x,res.gender)
        }

    }

    return(x)
}



# This function will preprare the colData for TCGAPrepare
# The idea is to add usefull information to the object
# that will help the users to understand their samples
# ref: TCGA codeTablesReport - Table: Sample type
#' @importFrom S4Vectors DataFrame
#' @importFrom stringr str_match

colDataPrepare <- function(barcode){

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
    aux <- DataFrame(code = code,shortLetterCode,definition)

    # in case multiple equal barcode
    regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                    "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
    samples <- str_match(barcode,regex)[,1]

    ret <- DataFrame(barcode = barcode,
                     patient = substr(barcode, 1, 12),
                     sample = substr(barcode, 1, 16),
                     code = substr(barcode, 14, 15))
    ret <- merge(ret,aux, by = "code", sort = FALSE)
    ret <- ret[match(barcode,ret$barcode),]
    rownames(ret) <- gsub("\\.","-",make.names(ret$barcode,unique=TRUE))
    ret$code <- NULL
    return(DataFrame(ret))
}


load.maf <- function(){
    if (requireNamespace("xml2", quietly = TRUE) & requireNamespace("rvest", quietly = TRUE) ) {
        tables <- xml2::read_html("https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files")
        tables <-  rvest::html_table(tables)
        # Table one is junk
        tables[[1]] <- NULL

        # get which tables are from the tumor
        all.df <- data.frame()
        for(tumor in unique(TCGAbiolinks::TCGAquery()$Disease)) {

            idx <- which(mapply(function(x) {
                any(grepl(tumor,(x[,1]), ignore.case = TRUE))
            },tables) == TRUE)
            df <- lapply(idx,function(x) tables[x])

            if(length(df) == 0) next
            # merge the data frame in the lists
            if(length(idx) > 1) {
                df <- Reduce(function(...) merge(..., all=TRUE), df)
            }  else if(length(idx) == 1) {
                df <- Reduce(function(...) merge(..., all=TRUE), df)
                df <- df[[1]]
                colnames(df) <- gsub(" ",".", colnames(df))
                colnames(df) <- gsub(":",".", colnames(df))
            }

            # Remove obsolete/protected
            df <- subset(df, df$Deploy.Status == "Available")
            df <- subset(df, df$Protection.Status == "Public")

            if(nrow(df) == 0) next

            df$Tumor <- tumor
            all.df <- rbind(all.df,df)
        }

        all.df[,"Deploy.Location"] <- gsub("/dccfiles_prod/tcgafiles/",
                                           "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/",
                                           all.df[,"Deploy.Location"] )
    }
    return(all.df)
}

get.mutation.matrix <- function(barcode,query){
    maf <- load.maf()
    maf <- maf[order(as.Date(maf$Deploy.Date, "%d-%b-%Y")),]
    aux <- plyr::ddply(maf, .(Tumor), function(x) x[c(nrow(x)), ])

    ret <- NULL
    for(disease in unique(query$Disease)){
        df <- aux[aux$Tumor == disease,]
        message(paste0("Downloading maf file",df$MAF.File.Name))
        if(is.windows()) mode <- "wb" else  mode <- "w"
        if (!file.exists(df$MAF.File.Name))
            downloader::download(df$Deploy.Location,df$MAF.File.Name, quiet = FALSE,mode = mode)

        mutation.matrix <- read.table(df$MAF.File.Name, fill = TRUE,
                                      comment.char = "#", header = TRUE, sep = "\t", quote='')
        mutation.matrix <- mutation.matrix[substr(mutation.matrix$Tumor_Sample_Barcode,1,15) %in% substr(barcode,1,15),]
        if(nrow(mutation.matrix) == 0) {
            next
        }
        # Fazer um subset de acordo com as amostras que eu tenho
        mutation.matrix <- data.table::setDT(mutation.matrix)
        mutation.matrix <- reshape2::acast(mutation.matrix, Tumor_Sample_Barcode~Hugo_Symbol, value.var="Variant_Type")
        mutation.matrix[which(mutation.matrix > 0)] <- 1
        if(is.null(ret)) {
            ret <- mutation.matrix
        }  else{
            ret <- plyr::rbind.fill(as.data.frame(ret),as.data.frame(mutation.matrix))
            ret[is.na(ret)] <- 0
        }
    }
    return(ret)
}


#' @title Retrieve molecular subtypes for a given tumor
#' @description
#'   TCGAquery_subtype Retrieve molecular subtypes for a given tumor
#' @param tumor is a cancer Examples:
#' \tabular{lllll}{
#' lgg   \tab gbm \tab luad \tab stad \tab brca\cr
#' coad \tab read \tab  \tab  \tab
#'}
#' @export
#' @examples
#' dataSubt <- TCGAquery_subtype(tumor = "lgg")
#' @return a data.frame with barcode and molecular subtypes
TCGAquery_subtype <- function(tumor){
    if (grepl("acc|lgg|gbm|luad|stad|brca|coad|read|skcm|hnsc|kich|lusc|ucec|pancan|thca|prad|kirp|kirc|all",
              tumor,ignore.case = TRUE)) {

        if(tolower(tumor) == "all") {
            all.tumor <- c("acc","lgg", "gbm", "luad", "stad", "brca", "coad",
                           "skcm", "hnsc", 'kich', "lusc", "ucec", "pancan", "thca",
                           "prad","kirp","kirc")
            doi <- c("acc"="doi:10.1016/j.ccell.2016.04.002",
                     "aml"="doi:10.1056/NEJMoa1301689",
                     "blca"="doi:10.1038/nature12965",
                     "brca"="doi:10.1038/nature11412",
                     "coad"="doi:10.1038/nature11252",
                     "gbm"="doi:10.1016/j.cell.2015.12.028",
                     "lgg"="doi:10.1016/j.cell.2015.12.028",
                     "hnsc"="doi:10.1038/nature14129",
                     "kich"="doi:10.1016/j.ccr.2014.07.014",
                     "kirc"="doi:10.1038/nature12222",
                     "kirp"="doi:10.1056/NEJMoa1505917",
                     "lihc"="",
                     "luad"="doi:10.1038/nature13385",
                     "lusc"="doi:10.1038/nature11404",
                     "ovca"= "doi:10.1038/nature10166",
                     "pancan"="doi:10.1016/j.cell.2014.06.049",
                     "prad"="doi:10.1016/j.cell.2015.10.025",
                     "skcm"="doi:10.1016/j.cell.2015.05.044",
                     "stad"="doi:10.1038/nature13480",
                     "thca"="doi:10.1016/j.cell.2014.09.050",
                     "ucec"="doi:10.1038/nature12113",
                     "ucs"="")

            all <- NULL
            for(i in all.tumor){
                try({
                    aux <- TCGAquery_subtype(i)
                    aux$Disease <- toupper(i)
                    aux$doi <- doi[i]
                    if(is.null(all)) all <- aux
                    else all <- plyr::rbind.fill(all,aux)
                })
            }
            idx <- which(colnames(all) == "doi")
            all <- all[, c(idx, (1:ncol(all))[-idx])]

            idx <- which(colnames(all) == "Disease")
            all <- all[, c(idx, (1:ncol(all))[-idx])]

            return(all)
        }

        # COAD and READ are in the same object
        if(tolower(tumor) == "read") tumor <- "coad"

        # The object with the gbm and lgg classification are the same
        # source: http://dx.doi.org/10.1016/j.cell.2015.12.028
        if(tolower(tumor) %in% c("lgg","gbm")) {
            aux <- get("lgg.gbm.subtype")
            if(tolower(tumor) == "gbm"){
                aux <- subset(aux,aux$Study == "Glioblastoma multiforme")
            } else {
                aux <- subset(aux,aux$Study != "Glioblastoma multiforme")
            }
            return (aux)
        }
        return(get(paste0(tolower(tumor),".subtype")))
    } else {
        stop("For the moment we have only subtype for: acc, brca, coad, gbm, hnsc, kich, kirp, kirc, lgg, luad, lusc, prad, pancan, read, skcm, stad, thca and ucec")
    }
}
