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
#' GDCquery_clinic will download all clinical information from the API
#' as the one with using the button from each project
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
    if(!grepl("clinical|Biospecimen",type,ignore.case = TRUE)) stop("Type must be clinical or biospecemen")
    baseURL <- "https://gdc-api.nci.nih.gov/cases/?"
    options.pretty <- "pretty=true"
    if(grepl("clinical",type,ignore.case = TRUE)) {
        options.expand <- "expand=diagnoses,diagnoses.treatments,annotations,family_histories,demographic,exposures"
        option.size <- paste0("size=",getNbCases(project,"Clinical"))
        files.data_category <- "Clinical"
    } else {
        options.expand <- "expand=samples,samples.portions,samples.portions.analytes,samples.portions.analytes.aliquots"
        option.size <- paste0("size=",getNbCases(project,"Biospecimen"))
        files.data_category <- "Biospecimen"
    }
    options.filter <- paste0("filters=",
                             URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["'),
                             project,
                             URLencode('"]}},{"op":"in","content":{"field":"files.data_category","value":["'),
                             files.data_category,
                             URLencode('"]}}]}'))
    url <- paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter,"format=json", sep = "&"))
    json  <- tryCatch({
        json <- fromJSON(url, simplifyDataFrame = TRUE)
    }, error = function(e) {
        json <- fromJSON(content(GET(url), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
    })

    #message(paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&")))
    results <- json$data$hits
    if(grepl("clinical",type,ignore.case = TRUE)) {
        diagnoses <- rbindlist(results$diagnoses, fill = TRUE)
        diagnoses$submitter_id <- gsub("_diagnosis","", diagnoses$submitter_id)
        exposures <- rbindlist(results$exposures, fill = TRUE)
        exposures$submitter_id <- gsub("_exposure","", exposures$submitter_id)
        results$demographic$submitter_id <- gsub("_demographic","", results$demographic$submitter_id)
        df <- merge(diagnoses,exposures, by="submitter_id", all = TRUE)
        df <- merge(df,results$demographic, by="submitter_id", all = TRUE)
        treatments <- rbindlist(df$treatments,fill = TRUE)
        df[,treatments:=NULL]
        treatments$submitter_id <- gsub("_treatment","", treatments$submitter_id)
        df <- merge(df,treatments, by="submitter_id", all = TRUE)
        df$bcr_patient_barcode <- df$submitter_id
        df$disease <- gsub("TCGA-|TARGET-", "", project)
    } else {
        df <- rbindlist(results$samples,fill = TRUE)
    }

    #y <- data.frame(diagnosis=I(results$diagnoses), demographic=results$demographic,exposures=I(results$exposures))
    if(save.csv){
        if(grepl("biospecimen",type))  df[,portions:=NULL]
        write_csv(df,paste0(project,"_",type,".csv"))
    }
    setDF(df)
    return(df)
}

#' @title Parsing clinical xml files
#' @description
#' This function receives the query argument and parses the clinical xml files
#' based on the desired information
#' @param query Result from GDCquery, with data.category set to Clinical
#' @param clinical.info Which information should be retrieved.
#' Options Clinical: drug, admin, follow_up,radiation, patient, stage_event or new_tumor_event
#' Options Biospecimen: protocol, admin, aliquot, analyte, bio_patient, sample, portion, slide
#' @param directory Directory/Folder where the data was downloaded. Default: GDCdata
#' @importFrom xml2 read_xml xml_ns
#' @importFrom XML xmlParse getNodeSet xmlToDataFrame
#' @importFrom plyr rbind.fill
#' @importFrom dplyr mutate_each funs
#' @export
#' @examples
#' query <- GDCquery(project = "TCGA-COAD",
#'                   data.category = "Clinical",
#'                   barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))
#' GDCdownload(query)
#' clinical <- GDCprepare_clinic(query,"patient")
#' clinical.drug <- GDCprepare_clinic(query,"drug")
#' clinical.radiation <- GDCprepare_clinic(query,"radiation")
#' clinical.admin <- GDCprepare_clinic(query,"admin")
#' query <- GDCquery(project = "TCGA-COAD",
#'                   data.category = "Biospecimen",
#'                   barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))
#' GDCdownload(query)
#' clinical <- GDCprepare_clinic(query,"admin")
#' clinical.drug <- GDCprepare_clinic(query,"sample")
#' clinical.radiation <- GDCprepare_clinic(query,"portion")
#' clinical.admin <- GDCprepare_clinic(query,"slide")
GDCprepare_clinic <- function(query, clinical.info, directory = "GDCdata"){
    if(unique(query$results[[1]]$data_category) != "Biospecimen") {
        valid.clinical.info <- c("drug","admin","follow_up","radiation","patient","stage_event","new_tumor_event")
    } else  if(unique(query$results[[1]]$data_category) != "Clinical") {
        valid.clinical.info <- c("protocol","admin","aliquot","analyte","bio_patient","sample", "portion", "slide")
    } else {
        stop("Data category should be Clinical or Biospecimen")
    }
    if(missing(clinical.info)) stop(paste0("Please set clinical.info argument:\n=> ",paste(valid.clinical.info,collapse = "\n=> ")))
    if(!(clinical.info %in% valid.clinical.info)) stop(paste0("Please set a valid clinical.info argument:\n=> ",paste(valid.clinical.info,collapse = "\n=> ")))

    # Get all the clincal xml files
    source <- ifelse(query$legacy,"legacy","harmonized")
    files <- file.path(query$project, source,
                       gsub(" ","_",query$results[[1]]$data_category),
                       gsub(" ","_",query$results[[1]]$data_type),
                       gsub(" ","_",query$results[[1]]$file_id),
                       gsub(" ","_",query$results[[1]]$file_name))
    files <- file.path(directory, files)
    if(!all(file.exists(files))) stop(paste0("I couldn't find all the files from the query.",
                                             "Please check directory parameter right"))
    xpath <- NULL

    disease <- tolower(gsub("TCGA-","",query$project))
    if(tolower(clinical.info) == "drug")      xpath <- "//rx:drug"
    else if(tolower(clinical.info) == "admin")     xpath <- "//admin:admin"
    else if(tolower(clinical.info) == "radiation") xpath <- "//rad:radiation"
    else if(tolower(clinical.info) == "patient")   xpath <- paste0("//",disease,":patient")
    else if(tolower(clinical.info) == "stage_event")     xpath <- "//shared_stage:stage_event"
    else if(tolower(clinical.info) == "new_tumor_event") xpath <- paste0("//",disease,"_nte:new_tumor_event")
    # biospecimen xpaths
    else if(tolower(clinical.info) == "sample")      xpath <- "//bio:sample"
    else if(tolower(clinical.info) == "bio_patient")      xpath <- "//bio:patient"
    else if(tolower(clinical.info) == "analyte")      xpath <- "//bio:analyte"
    else if(tolower(clinical.info) == "aliquot")      xpath <- "//bio:aliquot"
    else if(tolower(clinical.info) == "protocol")      xpath <- "//bio:protocol"
    else if(tolower(clinical.info) == "portion")      xpath <- "//bio:portion"
    else if(tolower(clinical.info) == "slide")  xpath <- "//bio:slide"

    clin <- parseXML(files,xpath,clinical.info)

    if(tolower(clinical.info) == "patient") {
        composed.cols <- c("new_tumor_events","drugs","follow_ups","radiations")
        composed.cols <- composed.cols[composed.cols %in% colnames(clin)]
        message("To get the following information please change the clinical.info argument")
        message("=> new_tumor_events: new_tumor_event \n=> drugs: drug \n=> follow_ups: follow_up \n=> radiations: radiation")

        for(i in composed.cols){
            clin[,i] <- as.character(clin[,i])
            clin[which(clin[,i] != ""),i] <- "YES"
            clin[which(clin[,i] == ""),i] <- "NO"
            colnames(clin)[which(colnames(clin) == i)] <- paste0("has_",i,"_information")
        }
        if("stage_event" %in% colnames(clin)) {
            message("Adding stage event information")
            aux <- parseXML(files,"//shared_stage:stage_event","stage_event")
            colnames(aux)[1:(length(colnames(aux)) - 1)] <- paste0("stage_event_",colnames(aux)[1:(length(colnames(aux)) - 1)])
            clin <- merge(clin, aux, by = "bcr_patient_barcode" , sort = FALSE, all.x = TRUE)
            clin$stage_event <- NULL
        }
        if("primary_pathology" %in% colnames(clin)) {
            message("Adding primary pathology information")
            aux <- parseXML(files,paste0("//",disease,":primary_pathology"),"primary_pathology")
            # Last column is the idx to merge
            colnames(aux)[1:(length(colnames(aux)) - 1)] <- paste0("primary_pathology_",colnames(aux)[1:(length(colnames(aux)) - 1)])
            clin <- merge(clin, aux, by = "bcr_patient_barcode" , sort = FALSE, all.x = TRUE)
            clin$primary_pathology <- NULL
        }

    }
    if(tolower(clinical.info) == "samples") clin$samples <- NULL
    if(tolower(clinical.info) == "portion") {
        for(i in c("slides","analytes")){
            clin[,i] <- as.character(clin[,i])
            clin[which(clin[,i] != ""),i] <- "YES"
            clin[which(clin[,i] == ""),i] <- "NO"
            colnames(clin)[which(colnames(clin) == i)] <- paste0("has_",i,"_information")
        }
    }
    if(is.null(clin)) {
        message("No information found")
        return(NULL)
    }
    # Converting factor to numeric and double
    out <- clin %>%
        mutate_each(
            funs(
                type.convert(as.character(.), as.is = TRUE, numerals = "warn.loss")
            )
        )

    # Change columns back to factor
    for(i in colnames(out)[!grepl("has_",colnames(out))]) {
        if(class(out[,i]) == "character" &  length(unique(out[,i])) < nrow(out)/2)
            out[,i] <-  as.factor(out[,i])
    }
    return(out)
}


parseXML <- function(files, xpath, clinical.info ){
    clin <- NULL
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    for(i in seq_along(files)){
        xmlfile <- files[i]
        xml <- read_xml(xmlfile)
        doc = xmlParse(xmlfile)

        if(tolower(clinical.info) == "follow_up" & is.null(xpath)){
            follow_up_version <-  names(xml_ns(xml))[grepl("follow_up",names(xml_ns(xml)))]
            if(length(follow_up_version) > 1) {
                follow_up_version <- follow_up_version[length(follow_up_version)]
            }
            if(length(follow_up_version) == 1)  xpath <- paste0("//", follow_up_version, ":follow_up")
            else next
        }

        patient <- str_extract(xmlfile,"[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}")
        # Test if this xpath exists before parsing it
        if(gsub("\\/\\/","", unlist(stringr::str_split(xpath,":"))[1]) %in% names(xml_ns(xml))){
            nodes <- getNodeSet(doc,xpath)
            if(length(nodes) == 0) next;
            df <- NULL
            for(j in 1:length(nodes)) {
                df.aux <- xmlToDataFrame(nodes = nodes[j])
                if(NA %in% colnames(df.aux)) df.aux <- df.aux[,!is.na(colnames(df.aux))]
                if(nrow(df.aux) == 0) next

                if(j == 1) {
                    df <- df.aux
                } else {
                    df <- rbind.fill(df,df.aux)
                }
            }

            df$bcr_patient_barcode <- patient
            if(i == 1) {
                clin <- df
            } else {
                clin <- rbind.fill(clin,df)
            }
            setTxtProgressBar(pb, i)
        }
    }
    close(pb)
    return(clin)
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
                 "read"="doi:10.1038/nature11252",
                 "skcm"="doi:10.1016/j.cell.2015.05.044",
                 "stad"="doi:10.1038/nature13480",
                 "thca"="doi:10.1016/j.cell.2014.09.050",
                 "ucec"="doi:10.1038/nature12113",
                 "ucs"="")

        if(tolower(tumor) != "all") message(paste0("Subtype information from:", doi[tolower(tumor)]))
        if(tolower(tumor) == "all") {
            all.tumor <- c("acc","lgg", "gbm", "luad", "stad", "brca", "coad",
                           "skcm", "hnsc", 'kich', "lusc", "ucec", "pancan", "thca",
                           "prad","kirp","kirc")
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
