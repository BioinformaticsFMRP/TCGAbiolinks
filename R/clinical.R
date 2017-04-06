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


#' @title Get GDC clinical data
#' @description
#' GDCquery_clinic will download all clinical information from the API
#' as the one with using the button from each project
#' @param project A valid project (see list with getGDCprojects()$project_id)]
#' @param type A valid type. Options "clinical", "Biospecimen"  (see list with getGDCprojects()$project_id)]
#' @param save.csv Write clinical information into a csv document
#' @export
#' @importFrom data.table rbindlist as.data.table
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
    json  <- tryCatch(
        getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
        error = function(e) {
            fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
        }
    )

    #message(paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&")))
    results <- json$data$hits

    if(grepl("clinical",type,ignore.case = TRUE)) {
        if(grepl("TCGA",project)) {
            df <- rbindlist(results$diagnoses, fill = TRUE)
            df$submitter_id <- gsub("_diagnosis","", df$submitter_id)
            if("exposures" %in% colnames(results)){
                exposures <- rbindlist(results$exposures, fill = TRUE)
                exposures <- exposures[,-c("updated_datetime","state","created_datetime")]
                exposures$submitter_id <- gsub("_exposure","", exposures$submitter_id)
                df <- merge(df,exposures, by="submitter_id", all = TRUE)
            }
            if("demographic" %in% colnames(results)){
                results$demographic$submitter_id <- gsub("_demographic","", results$demographic$submitter_id)
                df <- merge(df,as.data.table(results$demographic)[,-c("updated_datetime","state","created_datetime")], by="submitter_id", all = TRUE)
            }
            if("treatments" %in% colnames(df)){
                treatments <- rbindlist(df$treatments,fill = TRUE)
                df[,treatments:=NULL]
                treatments$submitter_id <- gsub("_treatment","", treatments$submitter_id)
                df <- merge(df,as.data.table(treatments)[,-c("updated_datetime","state","created_datetime")], by="submitter_id", all = TRUE)
            }
            df$bcr_patient_barcode <- df$submitter_id
            df$disease <- gsub("TCGA-|TARGET-", "", project)
        } else {
            df <- rbindlist(results$diagnoses, fill = TRUE)
            if("exposures" %in% colnames(results)){
                exposures <- rbindlist(results$exposures, fill = TRUE)
                df <- cbind(df,exposures)
            }
            if("treatments" %in% colnames(results)){
                treatments <- rbindlist(results$treatments, fill = TRUE)
                df[,treatments:=NULL]
                df <- cbind(df,treatments)
            }
            if("demographic" %in% colnames(results)){
                df <- cbind(df,results$demographic)
            }
            df$disease <- gsub("TCGA-|TARGET-", "", project)
        }
    } else {
        df <- rbindlist(results$samples,fill = TRUE)
    }

    if(save.csv){
        if(grepl("biospecimen",type))  {
            df[,portions:=NULL]
            message("Portion column is a list, it will be removed. Please check object with save.csv argument as FALSE")
        }
        if(grepl("clinical",type))  {
            df[,treatments:=NULL]
            message("Treatments column is a list, it will be removed. Please check object with save.csv argument as FALSE")
        }
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
    if(unique(query$results[[1]]$data_category) == "Clinical") {
        valid.clinical.info <- c("drug","admin","follow_up","radiation","patient","stage_event","new_tumor_event")
    } else  if(unique(query$results[[1]]$data_category) == "Biospecimen" ) {
        valid.clinical.info <- c("protocol","admin","aliquot","analyte","bio_patient","sample", "portion", "slide")
    } else  if(unique(query$results[[1]]$data_category) == "Other") {
        valid.clinical.info <- c("admin","msi")
    } else {
        stop("Data category should be Clinical or Biospecimen or Other (auxiliary files)")
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
    else if(tolower(clinical.info)  == "msi")  xpath <- "//auxiliary:microsatellite_instability_test_result"
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
    if (grepl("acc|lgg|gbm|luad|stad|brca|coad|read|skcm|hnsc|kich|lusc|ucec|pancan|thca|prad|pcpg|kirp|kirc|all",
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
                 "pcpg"="http://dx.doi.org/10.1016/j.ccell.2017.01.001",
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
