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
#'  barcode <- c("TARGET-20-PANSBH-14A-02D","TARGET-20-PANSBH-01A-02D",
#'               "TCGA-B0-4698-01Z-00-DX1","TCGA-CZ-4863-02Z-00-DX1")
#'   TCGAquery_SampleTypes(barcode,c("TR","TP"))
#' @return a list of samples / barcode filtered by type sample selected
TCGAquery_SampleTypes <- function(barcode,typesample){
    # Tumor AND Solid Tissue Normal NOT FROM THE SAME PATIENTS
    table.code <- c('01','02','03','04','05','06','07','08','09','10',
                    '11','12','13','14','20','40','50','60','61')
    names(table.code) <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC",
                           "TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB",
                           "CELL","XP","XCL")

    if (sum(is.element(typesample,names(table.code))) == length(typesample)) {
        string <- sapply(barcode,function(x) substr(unlist(stringr::str_split(x,"-"))[4],1,2))

        barcode.all <- NULL
        for (sample.i in typesample) {
            barcode.all <- union(barcode.all,
                                 barcode[grep(table.code[sample.i], string)])
        }
        return(barcode.all)
    } else {
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
#'  barcode <- c("TARGET-20-PANSBH-02A-02D","TARGET-20-PANSBH-01A-02D",
#'               "TCGA-B0-4698-01Z-00-DX1","TCGA-CZ-4863-02Z-00-DX1",
#'               "TARGET-20-PANSZZ-02A-02D","TARGET-20-PANSZZ-11A-02D",
#'               "TCGA-B0-4699-01Z-00-DX1","TCGA-B0-4699-02Z-00-DX1"
#'               )
#'   TCGAquery_MatchedCoupledSampleTypes(barcode,c("TR","TP"))
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

        string <- sapply(barcode,function(x) substr(unlist(stringr::str_split(x,"-"))[4],1,2))
        barcode.1 <- barcode[grep(table.code[typesample[1]], string)]
        barcode.2 <- barcode[grep(table.code[typesample[2]], string)]
        barcode.common <- intersect(sapply(barcode.1,function(x) paste(unlist(stringr::str_split(x,"-"))[1:3],collapse = "-")),
                                    sapply(barcode.2,function(x) paste(unlist(stringr::str_split(x,"-"))[1:3],collapse = "-")))
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
#' \itemize{
#' \item{ BEATAML1.0-COHORT }
#' \item{ BEATAML1.0-CRENOLANIB }
#' \item{ CGCI-BLGSP }
#' \item{ CPTAC-2 }
#' \item{ CPTAC-3 }
#' \item{ CTSP-DLBCL1 }
#' \item{ FM-AD }
#' \item{ HCMI-CMDC }
#' \item{ MMRF-COMMPASS }
#' \item{ NCICCR-DLBCL }
#' \item{ OHSU-CNL }
#' \item{ ORGANOID-PANCREATIC }
#' \item{ TARGET-ALL-P1 }
#' \item{ TARGET-ALL-P2 }
#' \item{ TARGET-ALL-P3 }
#' \item{ TARGET-AML }
#' \item{ TARGET-CCSK }
#' \item{ TARGET-NBL }
#' \item{ TARGET-OS }
#' \item{ TARGET-RT }
#' \item{ TARGET-WT }
#' \item{ TCGA-ACC }
#' \item{ TCGA-BLCA }
#' \item{ TCGA-BRCA }
#' \item{ TCGA-CESC }
#' \item{ TCGA-CHOL }
#' \item{ TCGA-COAD }
#' \item{ TCGA-DLBC }
#' \item{ TCGA-ESCA }
#' \item{ TCGA-GBM }
#' \item{ TCGA-HNSC }
#' \item{ TCGA-KICH }
#' \item{ TCGA-KIRC }
#' \item{ TCGA-KIRP }
#' \item{ TCGA-LAML }
#' \item{ TCGA-LGG }
#' \item{ TCGA-LIHC }
#' \item{ TCGA-LUAD }
#' \item{ TCGA-LUSC }
#' \item{ TCGA-MESO }
#' \item{ TCGA-OV }
#' \item{ TCGA-PAAD }
#' \item{ TCGA-PCPG }
#' \item{ TCGA-PRAD }
#' \item{ TCGA-READ }
#' \item{ TCGA-SARC }
#' \item{ TCGA-SKCM }
#' \item{ TCGA-STAD }
#' \item{ TCGA-TGCT }
#' \item{ TCGA-THCA }
#' \item{ TCGA-THYM }
#' \item{ TCGA-UCEC }
#' \item{ TCGA-UCS }
#' \item{ TCGA-UVM }
#' \item{ VAREPOP-APOLLO }
#' }
#' @param type A valid type. Options "clinical", "Biospecimen"  (see list with getGDCprojects()$project_id)]
#' @param save.csv Write clinical information into a csv document
#' @export
#' @importFrom data.table rbindlist as.data.table
#' @importFrom jsonlite fromJSON
#' @examples
#' clinical <- GDCquery_clinic(
#'    project = "TCGA-ACC",
#'    type = "clinical",
#'    save.csv = FALSE
#'  )
#' clinical <- GDCquery_clinic(
#'    project = "TCGA-ACC",
#'    type = "biospecimen",
#'    save.csv = FALSE
#' )
#' \dontrun{
#' clinical_cptac_3 <- GDCquery_clinic(
#'    project = "CPTAC-3",
#'    type = "clinical"
#' )
#' clinical_cptac_2 <- GDCquery_clinic(
#'    project = "CPTAC-2",
#'    type = "clinical"
#' )
#' clinical_HCMI_CMDC <- GDCquery_clinic(
#'    project = "HCMI-CMDC",
#'    type = "clinical"
#' )
#' clinical_GCI_HTMCP_CC <- GDCquery_clinic(
#'    project = "CGCI-HTMCP-CC",
#'    type = "clinical"
#' )
#' clinical <- GDCquery_clinic(
#'    project = "NCICCR-DLBCL",
#'    type = "clinical"
#' )
#' clinical <- GDCquery_clinic(
#'    project = "ORGANOID-PANCREATIC",
#'    type = "clinical"
#' )
#' }
#' @return A data frame with the clinical information
#' @author Tiago Chedraoui Silva
GDCquery_clinic <- function(
        project,
        type = "clinical",
        save.csv = FALSE
){
    checkProjectInput(project)

    if (length(project) > 1) {
        stop("Please, project should be only one valid project")
    }

    if (!grepl("clinical|Biospecimen",type,ignore.case = TRUE)) {
        stop("Type must be clinical or Biospecimen")
    }

    baseURL <- "https://api.gdc.cancer.gov/cases/?"
    options.pretty <- "pretty=true"

    if (grepl("clinical",type,ignore.case = TRUE)) {
        options.expand <- "expand=diagnoses,diagnoses.treatments,annotations,family_histories,demographic,exposures"
        option.size <- paste0("size=",getNbCases(project,"Clinical"))
        files.data_category <- "Clinical"
    } else {
        options.expand <- "expand=samples,samples.portions,samples.portions.analytes,samples.portions.analytes.aliquots"
        option.size <- paste0("size=",getNbCases(project,"Biospecimen"))
        files.data_category <- "Biospecimen"
    }

    if (grepl("TCGA",project)){
        options.filter <- paste0(
            "filters=",
            URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["'),
            project,
            URLencode('"]}},{"op":"in","content":{"field":"files.data_category","value":["'),
            files.data_category,
            URLencode('"]}}]}')
        )
    } else {
        options.filter <- paste0(
            "filters=",
            URLencode('{"op":"in","content":{"field":"cases.project.project_id","value":["'),
            project,
            URLencode('"]}}')
        )
    }
    url <- paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter,"format=json", sep = "&"))
    json  <- tryCatch(
        getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
        error = function(e) {
            fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
        }
    )


    #message(paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&")))
    results <- json$data$hits

    if (grepl("clinical",type,ignore.case = TRUE)) {
        if (grepl("TCGA",project)) {
            df <- data.frame("submitter_id" = results$submitter_id)
            if ("diagnoses" %in% colnames(results)){
                diagnoses <- rbindlist(lapply(results$diagnoses, function(x) if(is.null(x)) data.frame(NA) else x),fill = T)
                diagnoses$submitter_id <- gsub("_diagnosis","", df$submitter_id)
                df <- merge(df,diagnoses, by="submitter_id", all = TRUE, sort = FALSE)
            }
            if ("exposures" %in% colnames(results)){
                exposures <- rbindlist(results$exposures, fill = TRUE)
                exposures <- exposures[,-c("updated_datetime","state","created_datetime")]
                exposures$submitter_id <- gsub("_exposure","", exposures$submitter_id)
                df <- merge(df,exposures, by="submitter_id", all = TRUE, sort = FALSE)
            }
            if ("demographic" %in% colnames(results)){
                results$demographic$submitter_id <- gsub("_demographic","", results$demographic$submitter_id)
                demographic <- results$demographic[!is.na(results$demographic$submitter_id),]
                df <- merge(
                    df,
                    as.data.table(demographic)[,-c("updated_datetime","state","created_datetime")],
                    by = "submitter_id",
                    all = TRUE,
                    sort = FALSE
                )
            }

            if("treatments" %in% colnames(df)){
                treatments <- rbindlist(df$treatments,fill = TRUE)
                df$treatments <- NULL
                treatments$submitter_id <- gsub("_treatment(_[0-9])?","", treatments$submitter_id)
                treatments <- treatments[,-c("updated_datetime", "state", "created_datetime")]

                # we have now two types of treatment
                treatments.pharmaceutical <- treatments[grep("Pharmaceutical",treatments$treatment_type,ignore.case = TRUE),]
                treatments.radiation <- treatments[grep("radiation",treatments$treatment_type,ignore.case = TRUE),]

                # Adding a prefix
                colnames(treatments.pharmaceutical) <- paste0("treatments_pharmaceutical_",colnames(treatments.pharmaceutical))
                colnames(treatments.radiation) <- paste0("treatments_radiation_",colnames(treatments.radiation))
                colnames(treatments.radiation)[grep("submitter",colnames(treatments.radiation))] <- "submitter_id"
                colnames(treatments.pharmaceutical)[grep("submitter",colnames(treatments.pharmaceutical))] <- "submitter_id"

                df <- merge(df, as.data.table(treatments.pharmaceutical), by = "submitter_id",  all = TRUE, sort = FALSE)
                df <- merge(df, as.data.table(treatments.radiation), by = "submitter_id",  all = TRUE, sort = FALSE)
            }

            df$bcr_patient_barcode <- df$submitter_id
            df$project <- project
            df <- df %>% dplyr::relocate(project)

        } else {

            # Although for TCGA and TARGET IDs from diagnosis, treatments, exposures etc are the same
            # for the other projects this might not be true!
            # example: ORGANOID-PANCREATIC
            # https://api.gdc.cancer.gov/cases/?pretty=true&expand=diagnoses,demographic&size=1&filters=%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22cases.project.project_id%22,%22value%22:[%22ORGANOID-PANCREATIC%22]%7D%7D&format=json
            # DEMOGRAPHIC 48, while everything else is 42
            df <- data.frame("submitter_id" = results$submitter_id)

            if("disease_type" %in% colnames(results)){
                disease_type <- data.frame("disease_type" = results$disease_type)
                df <- cbind(df,disease_type)
            }
            if("primary_site" %in% colnames(results)){
                primary_site <- data.frame("primary_site" = results$primary_site)
                df <- cbind(df,primary_site)
            }

            if("diagnoses" %in% colnames(results)){
                diagnoses <- rbindlist(
                    lapply(
                        results$diagnoses,
                        function(x) {
                            if(is.null(x)) {
                                data.frame(NA)
                            } else {
                                # HTMCP-03-06-02061 has two diagnosis
                                x$submitter_id <- gsub("_diagnosis.*","",x$submitter_id)
                                aux <- x %>% dplyr::group_by(submitter_id) %>%
                                    dplyr::summarise_each(funs(paste(unique(.), collapse = ";")))
                                aux$treatments <- list(dplyr::bind_rows(x$treatments))
                                aux
                            }
                        }
                    ),fill = T
                )
                #df$submitter_id <- gsub("^d|_diagnosis|diag-|-DX|-DIAG|-diagnosis","", df$submitter_id)
                # ^d ORGANOID-PANCREATIC
                # -DX CPTAC-2
                # -DIAG CPTAC-3
                # -diagnosis NCICCR-DLBCL
                # _diagnosis HCMI-CMDC
                df <- cbind(df,diagnoses)
            }
            if("exposures" %in% colnames(results)){
                exposures <- rbindlist(lapply(results$exposures, function(x) if(is.null(x)) data.frame(NA) else x),fill = T)
                exposures <- exposures[,-c("updated_datetime","state","created_datetime")]
                df <- cbind(df,exposures)
            }
            if("treatments" %in% colnames(results)){
                treatments <- rbindlist(results$treatments, fill = TRUE)
                df$treatments <- NULL
                df <- cbind(df,treatments)
            }
            if("submitter_sample_ids" %in% colnames(results)){
                submitter_sample_ids <- lapply(
                    results$submitter_sample_ids,
                    function(x) {
                        paste(x,collapse = ",")
                    }
                ) %>% unlist
                df <- cbind(df,submitter_sample_ids)
            }

            # We do a merge because some cases might not have demographic information
            # DEMOGRAPHIC ORGANOID-PANCREATIC
            # -demographic NCICCR-DLBCL
            # _demographic HCMI-CMDC
            # -DG and -DM CPTAC-2
            # -DEMO CPTAC-3
            if("demographic" %in% colnames(results)){
                results$demographic$submitter_id <- gsub("_demographic|DEMOGRAPHIC|_demographic|demo-|-DG|-DM|-DEMO","", results$demographic$submitter_id)
                demographic <-  results$demographic %>% dplyr::select(-c("updated_datetime","state","created_datetime"))
                df <- cbind(df,demographic)
            }

            df$project <- project
            df <- df %>% dplyr::relocate(project)
        }
        if(nrow(results) != nrow(df)){
            stop("Error: API returned more information")
        }

    } else {
        df <- rbindlist(results$samples,fill = TRUE)
    }

    if (save.csv) {
        if(grepl("biospecimen",type))  {
            df$portions <- NULL
            message("Portion column is a list, it will be removed. Please check object with save.csv argument as FALSE")
        }
        if (grepl("clinical",type))  {
            df$treatments <- NULL
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
#' @importFrom dplyr mutate_all funs
#' @return A data frame with the parsed values from the XML
#' @export
#' @examples
#' query <- GDCquery(
#'   project = "TCGA-COAD",
#'   data.category = "Clinical",
#'   file.type = "xml",
#'   barcode = c("TCGA-RU-A8FL","TCGA-AA-3972")
#' )
#' GDCdownload(query)
#' clinical <- GDCprepare_clinic(query,"patient")
#' clinical.drug <- GDCprepare_clinic(query,"drug")
#' clinical.radiation <- GDCprepare_clinic(query,"radiation")
#' clinical.admin <- GDCprepare_clinic(query,"admin")
#' \dontrun{
#' query <- GDCquery(
#'    project = "TCGA-COAD",
#'    data.category = "Biospecimen",
#'    file.type = "xml",
#'    data.type = "Biospecimen Supplement",
#'    barcode = c("TCGA-RU-A8FL","TCGA-AA-3972")
#' )
#' GDCdownload(query)
#' clinical <- GDCprepare_clinic(query,"admin")
#' clinical.drug <- GDCprepare_clinic(query,"sample")
#' clinical.radiation <- GDCprepare_clinic(query,"portion")
#' clinical.admin <- GDCprepare_clinic(query,"slide")
#' }
GDCprepare_clinic <- function(
        query,
        clinical.info,
        directory = "GDCdata"
){

    if (unique(query$results[[1]]$data_category) == "Clinical") {
        valid.clinical.info <- c(
            "drug",
            "admin",
            "follow_up",
            "radiation",
            "patient",
            "stage_event",
            "new_tumor_event"
        )
    } else  if (unique(query$results[[1]]$data_category) == "Biospecimen" ) {
        valid.clinical.info <- c(
            "protocol",
            "admin",
            "aliquot",
            "analyte",
            "bio_patient",
            "sample",
            "portion",
            "slide"
        )
    } else  if (unique(query$results[[1]]$data_category) == "Other") {
        valid.clinical.info <- c("admin","msi")
    } else {
        stop("Data category should be Clinical or Biospecimen or Other (auxiliary files)")
    }

    if (missing(clinical.info) || !(clinical.info %in% valid.clinical.info)) {
        stop(
            "Please set clinical.info argument:\n=> ",
            paste(valid.clinical.info, collapse = "\n=> ")
        )
    }

    # Get all the clincal xml files
    source <- ifelse(query$legacy,"legacy","harmonized")
    files <- file.path(
        query$results[[1]]$project, source,
        gsub(" ","_",query$results[[1]]$data_category),
        gsub(" ","_",query$results[[1]]$data_type),
        gsub(" ","_",query$results[[1]]$file_id),
        gsub(" ","_",query$results[[1]]$file_name)
    )
    files <- file.path(directory, files)
    if(!all(file.exists(files))) {
        stop(
            "I couldn't find all the files from the query.",
            "Please check directory parameter right"
        )
    }

    xpath <- NULL

    disease <- tolower(gsub("TCGA-","",unlist(query$project)))
    if(tolower(clinical.info) == "drug") xpath <- "//rx:drug"
    else if(tolower(clinical.info) == "admin") xpath <- "//admin:admin"
    else if(tolower(clinical.info) == "radiation") xpath <- "//rad:radiation"
    else if(tolower(clinical.info) == "patient") xpath <- paste0("//",disease,":patient")
    else if(tolower(clinical.info) == "stage_event") xpath <- "//shared_stage:stage_event"
    else if(tolower(clinical.info) == "new_tumor_event") xpath <- paste0("//",disease,"_nte:new_tumor_event")
    # biospecimen xpaths
    else if(tolower(clinical.info) == "sample") xpath <- "//bio:sample"
    else if(tolower(clinical.info) == "bio_patient")  xpath <- "//bio:patient"
    else if(tolower(clinical.info) == "analyte") xpath <- "//bio:analyte"
    else if(tolower(clinical.info) == "aliquot") xpath <- "//bio:aliquot"
    else if(tolower(clinical.info) == "protocol") xpath <- "//bio:protocol"
    else if(tolower(clinical.info) == "portion") xpath <- "//bio:portion"
    else if(tolower(clinical.info) == "slide") xpath <- "//bio:slide"
    else if(tolower(clinical.info) == "msi") xpath <- "//auxiliary:microsatellite_instability_test_result"

    if (tolower(clinical.info) == "follow_up") {
        clin <- parseFollowup(files, xpath, clinical.info)
    } else {
        clin <- parseXML(files,xpath,clinical.info)
        if (is.null(clin)) return(NULL)
        clin <- merge(clin,getResults(query)[,c("project","cases")],by.x = c("bcr_patient_barcode"), by.y = "cases",all.x = TRUE)
    }

    if (tolower(clinical.info) == "patient") {
        composed.cols <- c("new_tumor_events","drugs","follow_ups","radiations")
        composed.cols <- composed.cols[composed.cols %in% colnames(clin)]
        message("To get the following information please change the clinical.info argument")
        message("=> new_tumor_events: new_tumor_event \n=> drugs: drug \n=> follow_ups: follow_up \n=> radiations: radiation")

        for (i in composed.cols) {
            clin[,i] <- as.character(clin[,i])
            clin[!clin[,i] %in%  c("","NO"),i] <- "YES"
            clin[clin[,i] %in% c("","NO"),i] <- "NO"

            if (i == "new_tumor_events") {
                followup <- parseFollowup(files,xpath,clinical.info)
                barcode <- followup$bcr_patient_barcode[!followup$new_tumor_events %in% c("","NO")]
                clin[clin$bcr_patient_barcode %in% barcode,i] <- "YES"
            }
            colnames(clin)[which(colnames(clin) == i)] <- paste0("has_",i,"_information")
        }

        if ("stage_event" %in% colnames(clin)) {
            message("Adding stage event information")
            aux <- parseXML(files,"//shared_stage:stage_event","stage_event")
            colnames(aux)[grep("bcr_patient_barcode",colnames(aux),invert = TRUE)] <- paste0("stage_event_",grep("bcr_patient_barcode",colnames(aux),invert = TRUE,value = TRUE))
            clin <- merge(clin, aux, by = "bcr_patient_barcode" , sort = FALSE, all.x = TRUE)
            clin$stage_event <- NULL
        }

        if ("primary_pathology" %in% colnames(clin)) {
            message("Adding primary pathology information")
            aux <- parseXML(files,paste0("//",disease,":primary_pathology"),"primary_pathology")
            # Last column is the idx to merge
            colnames(aux)[grep("bcr_patient_barcode",colnames(aux),invert = TRUE)] <- paste0("primary_pathology_",grep("bcr_patient_barcode",colnames(aux),invert = TRUE,value = TRUE))
            clin <- merge(clin, aux, by = "bcr_patient_barcode" , sort = FALSE, all.x = TRUE)
            clin$primary_pathology <- NULL
        }

        # Update clinical data with follow-ups: days_to_last_followup and vital_status
        message("Updating days_to_last_followup and vital_status from follow_up information using last entry")
        followup <- parseFollowup(files,xpath,clinical.info)

        followup_last <- followup %>% dplyr::group_by(bcr_patient_barcode) %>% dplyr::summarise(
            days_to_last_followup = max(as.numeric(days_to_last_followup),na.rm = TRUE),
            vital_status = vital_status[
                ifelse(
                    any(followup$days_to_last_followup %in% ""),
                    which(followup$days_to_last_followup %in% ""),
                    which.max(days_to_last_followup)
                )
            ]
        )
        clin$days_to_last_followup <- followup_last$days_to_last_followup[match(clin$bcr_patient_barcode,followup_last$bcr_patient_barcode)]
        clin$vital_status <- followup_last$vital_status[match(clin$bcr_patient_barcode,followup_last$bcr_patient_barcode)]
    }

    if (tolower(clinical.info) == "sample") {
        clin$portions <- NULL
        clin$samples <- NULL
    }

    if (tolower(clinical.info) == "bio_patient") {
        clin$samples <- NULL
    }

    if (tolower(clinical.info) == "portion") {
        for (i in c("slides","analytes")) {
            clin[,i] <- as.character(clin[,i])
            clin[which(clin[,i] != ""),i] <- "YES"
            clin[which(clin[,i] == ""),i] <- "NO"
            colnames(clin)[which(colnames(clin) == i)] <- paste0("has_",i,"_information")
        }
    }

    if (is.null(clin)) {
        message("No information found")
        return(NULL)
    }

    # Converting factor to numeric and double
    out <- clin |> type.convert(as.is = TRUE, numerals = "warn.loss")

    # Change columns back to factor
    for (i in colnames(out)[!grepl("has_",colnames(out))]) {
        if (is(out[,i],"character") &  length(unique(out[,i])) < nrow(out)/2)
            out[,i] <-  as.factor(out[,i])
    }
    return(out)
}

parseFollowup <- function(files, xpath, clinical.info){
    follow_up_version <- files %>%
        sapply(function(x) {
            name <- names(xml_ns(read_xml(x)))
            name[grepl("follow_up",name)]
        }) %>% unlist() %>% unique() %>% sort()

    if (length(follow_up_version) > 1) {
        message("We found more than one follow up version!")
        message("We will parse all and add a collumn (follow_up_version) to identify each version")
    }

    clin <- plyr::adply(
        .data = follow_up_version,
        .margins =  1,
        .fun = function(x){
            message("Parsing follow up version: ", x)
            xpath <- paste0("//", x, ":follow_up")
            clin <- parseXML(files,xpath,clinical.info)
            clin$follow_up_version <- x
            return(clin)
        }, .id = "follow_up_version"
    )

    col_idx <- grep("follow_up_version", names(clin))
    clin <- clin[, c(col_idx, (1:ncol(clin))[-col_idx])]
}

parseXML <- function(files, xpath, clinical.info ){
    clin <- NULL

    pb <- txtProgressBar(min = 0, max = length(files), style = 3)

    for (i in seq_along(files)) {

        xmlfile <- files[i]
        xml <- read_xml(xmlfile)
        doc = xmlParse(xmlfile)

        patient <- str_extract(xmlfile,"[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}")
        # Test if this xpath exists before parsing it
        for (xpath.it in xpath) {
            if (gsub("\\/\\/","", unlist(stringr::str_split(xpath.it,":"))[1]) %in% names(xml_ns(xml))){
                nodes <- getNodeSet(doc,xpath.it)
                if (length(nodes) == 0) next;
                df <- NULL
                for (j in 1:length(nodes)) {
                    df.aux <- xmlToDataFrame(nodes = nodes[j])
                    if (NA %in% colnames(df.aux)) df.aux <- df.aux[,!is.na(colnames(df.aux))]
                    if (nrow(df.aux) == 0) next

                    if (j == 1) {

                        df <- df.aux
                    } else {
                        df <- rbind.fill(df,df.aux)
                    }
                }

                df$bcr_patient_barcode <- patient
                if (i == 1) {
                    clin <- df
                } else {
                    clin <- rbind.fill(clin,df)
                }
                setTxtProgressBar(pb, i)
            }
        }
    }
    close(pb)
    return(clin)
}

#' @title Retrieve table with TCGA molecular subtypes
#' @description
#'   PanCancerAtlas_subtypes is a curated table with molecular subtypes for 24 TCGA cancer types
#' @export
#' @examples
#' molecular.subtypes <- PanCancerAtlas_subtypes()
#' @return a data.frame with barcode and molecular subtypes for 24 cancer types
PanCancerAtlas_subtypes <- function(){
    return(pancan2018)
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
    available <- c("ACC",
                   "BRCA",
                   "BLCA",
                   "CESC",
                   "CHOL",
                   "COAD",
                   "ESCA",
                   "GBM",
                   "HNSC",
                   "KICH",
                   "KIRC",
                   "KIRP",
                   "LGG",
                   "LIHC",
                   "LUAD",
                   "LUSC",
                   "PAAD",
                   "PCPG",
                   "PRAD",
                   "READ",
                   "SKCM",
                   "SARC",
                   "STAD",
                   "THCA",
                   "UCEC",
                   "UCS",
                   "UVM"
    )
    if (
        grepl(
            paste(c(available,"all"),collapse = "|"), tumor,
            ignore.case = TRUE
        )
    ) {

        doi <- c(
            "acc"  = "doi:10.1016/j.ccell.2016.04.002",
            "aml"  = "doi:10.1056/NEJMoa1301689",
            "blca" = "doi:10.1016/j.cell.2017.09.007",
            "brca" = "doi.org/10.1016/j.ccell.2018.03.014",
            "cesc" = "doi:10.1038/nature21386",
            "chol" = "doi:10.1016/j.celrep.2017.02.033",
            "coad" = "doi:10.1038/nature11252",
            "esca" = "doi:10.1038/nature20805",
            "gbm"  = "doi:10.1016/j.cell.2015.12.028",
            "lgg"  = "doi:10.1016/j.cell.2015.12.028",
            "lihc"  = "https://doi.org/10.1016/j.cell.2017.05.046",
            "hnsc" = "doi:10.1038/nature14129",
            "kich" = "doi:10.1016/j.ccr.2014.07.014",
            "kirc" = "doi:10.1038/nature12222",
            "kirp" = "doi:10.1056/NEJMoa1505917",
            "lihc" = "doi:10.1016/j.cell.2017.05.046",
            "luad" = "doi:10.1038/nature13385",
            "lusc" = "doi:10.1038/nature11404",
            "paad"  = "doi:10.1016/j.ccell.2017.07.007",
            "ovca" = "doi:10.1038/nature10166",
            "pancan"="doi:10.1016/j.cell.2014.06.049",
            "pcpg" = "doi:10.1016/j.ccell.2017.01.001",
            "prad" = "doi:10.1016/j.cell.2015.10.025",
            "read" = "doi:10.1038/nature11252",
            "sarc" = "doi:10.1016/j.cell.2017.10.014",
            "skcm" = "doi:10.1016/j.cell.2015.05.044",
            "stad" = "doi:10.1038/nature13480",
            "thca" = "doi:10.1016/j.cell.2014.09.050",
            "ucec" = "doi:10.1038/nature12113",
            "ucs"  = "doi:10.1016/j.ccell.2017.02.010",
            "uvm"  = "doi:10.1016/j.ccell.2017.07.003"
        )
        if(tolower(tumor) != "all") {
            message(
                paste0(
                    tolower(tumor),
                    " subtype information from:",
                    doi[tolower(tumor)]
                )
            )
        }
        if(tolower(tumor) == "all") {
            all <- NULL
            for(i in available){
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
        stop("For the moment we have only subtype for:\no ", paste(available,collapse = "\no "))
    }
}


#' @title Retrieve molecular subtypes for given TCGA barcodes
#' @description
#'   TCGA_MolecularSubtype Retrieve molecular subtypes from TCGA consortium for a given set of barcodes
#' @param barcodes is a vector of TCGA barcodes
#' @export
#' @examples
#' TCGA_MolecularSubtype("TCGA-60-2721-01A-01R-0851-07")
#' @return List with $subtypes attribute as a dataframe with barcodes,
#' samples, subtypes, and colors. The $filtered attribute is returned as filtered samples with no subtype info
TCGA_MolecularSubtype <- function(barcodes){
    tabPanCancer <- as.data.frame(PanCancerAtlas_subtypes())
    tabSubtypeMergedNew <- subset(tabPanCancer, select = c("pan.samplesID",
                                                           "Subtype_Selected"))
    colnames(tabSubtypeMergedNew) <- c("samples", "subtype")
    rownames(tabSubtypeMergedNew) <- tabSubtypeMergedNew$samples
    tabSubtypeMergedNew <- cbind(tabSubtypeMergedNew, color = rep(0,
                                                                  nrow(tabSubtypeMergedNew)))
    TabSubtypesCol_merged <- TabSubtypesCol_merged[!duplicated(TabSubtypesCol_merged$samples),
    ]
    commonSamples <- intersect(tabSubtypeMergedNew$samples, TabSubtypesCol_merged$samples)
    rownames(TabSubtypesCol_merged) <- TabSubtypesCol_merged$samples
    tabSubtypeMergedNew[commonSamples, "color"] <- TabSubtypesCol_merged[commonSamples,
                                                                         "color"]
    dataTableSubt <- tabSubtypeMergedNew
    dataTableSubt$samples <- as.character(dataTableSubt$samples)
    dataTableSubt.samples.stripped <- sapply(as.character(dataTableSubt$samples), function(x) paste(unlist(stringr::str_split(x,
                                                                                                                              "-"))[1:3], collapse = "-"))
    dataTableSubt$patients<- dataTableSubt.samples.stripped

    barcodes <- as.character(barcodes)
    patients <- sapply(barcodes, function(x) paste(unlist(stringr::str_split(x,
                                                                             "-"))[1:3], collapse = "-"))

    #barcodes and patients with no subtype info
    filt.p <- c()
    filt.b <- c()

    df.barcodes_patID<-data.frame("barcodes"= barcodes, "patID"=unlist(patients),
                                  row.names = 1:length(barcodes), stringsAsFactors = FALSE)


    #View(df.barcodes_patID)

    for (p in patients) {

        if (p %in% unlist(dataTableSubt.samples.stripped) == FALSE){
            filt.p <- c(filt.p, p)
            #print(p)

            #filt.p is a vector containing patients with no subtype info
        }
    }

    for (b in barcodes) {

        if (b %in% unlist(dataTableSubt$samples) == FALSE){
            filt.b <- c(filt.b, b)
            #print(b)

            #filt.b is a vector containing samples with no subtype info
        }
    }

    if (length(filt.p) > 0) {
        message("the following TCGA barcodes/patients with no subtypes were filtered:")

        ###Keeping only patients/barcodes with available subtype info

        patients.with.sub <- unlist(patients[patients %in% filt.p ==
                                                 FALSE])

        idx.barcodes<-which(df.barcodes_patID$barcodes%in%filt.b ==FALSE)

        barcodes.with.sub<- df.barcodes_patID[idx.barcodes,]$barcodes

        ###indices of barcodes with available molecular subtype info
        ###if a patient Id is in filt.p and in filt.b
        idx.patient <- which(df.barcodes_patID$patID %in% filt.p ==FALSE)
        idx.barcodes <- which(df.barcodes_patID$barcodes %in% filt.b ==FALSE)
        idx.df<-intersect(idx.barcodes, idx.patient)

        ###dataframe with barcodes and patients IDs with available moolecular subtypes

        #df.barcodes_patID[idx.df,])

        idx1 <- which(dataTableSubt$samples %in% barcodes.with.sub)
        idx2 <- which(dataTableSubt$patients %in% patients.with.sub)
        idx<-intersect(idx1, idx2)

        Subtypes <- dataTableSubt[idx, ]

        #matching barcodes so they have the same order as they were provided in arguments
        Subtypes <- Subtypes[match(Subtypes$samples, df.barcodes_patID[idx.df,]$barcodes), ]

        filt <- setdiff(df.barcodes_patID$barcodes, Subtypes$samples)
        return(list(subtypes = Subtypes, filtered = filt))
    } else {
        message("All barcodes have available molecular subtype info")
        filt <- c()
        idx.patient <- which(dataTableSubt$samples %in% df.barcodes_patID$barcodes)
        Subtypes <- as.data.frame(dataTableSubt[idx.patient, ])
        Subtypes <- Subtypes[match(df.barcodes_patID$barcodes, Subtypes$samples), ]

        ##filt will be empty because all samples have molecular subtype info
        return(list(subtypes = Subtypes, filtered = filt))
    }
}

#' @title Filters TCGA barcodes according to purity parameters
#' @description
#'   TCGAtumor_purity Filters TCGA samples using 5 estimates from 5 methods as thresholds.
#' @param barcodes is a vector of TCGA barcodes
#' @param estimate uses gene expression profiles of 141 immune genes and 141 stromal genes
#' @param absolute which uses somatic copy-number data (estimations were available for only 11 cancer types)
#' @param lump (leukocytes unmethylation for purity), which averages 44 non-methylated immune-specific CpG sites
#' @param ihc as estimated by image analysis of haematoxylin and eosin stain slides produced by the Nationwide Childrens Hospital Biospecimen Core Resource
#' @param cpe CPE is a derived consensus measurement as the median purity level after normalizing levels from all methods to give them equal means and s.ds
#' @export
#' @examples
#' dataTableSubt <- TCGAtumor_purity("TCGA-60-2721-01A-01R-0851-07",
#'                           estimate = 0.6,
#'                           absolute = 0.6,
#'                           ihc = 0.8,
#'                           lump = 0.8,
#'                           cpe = 0.7)
#' @return List with $pure_barcodes attribute as a vector of pure samples and $filtered attribute as filtered samples with no purity info
TCGAtumor_purity <- function(barcodes, estimate, absolute, lump, ihc, cpe){

    Tumor.purity.L <- Tumor.purity
    barcodes <= as.character(barcodes)
    Tumor.purity.L$Sample.ID <- as.character(Tumor.purity$Sample.ID)
    Tumor.purity.L$ESTIMATE <- as.numeric(gsub(",", ".", Tumor.purity$ESTIMATE))
    Tumor.purity.L$ABSOLUTE <- as.numeric(gsub(",", ".", Tumor.purity$ABSOLUTE))
    Tumor.purity.L$LUMP <- as.numeric(gsub(",", ".", Tumor.purity$LUMP))
    Tumor.purity.L$IHC <- as.numeric(gsub(",", ".", Tumor.purity$IHC))
    Tumor.purity.L$CPE <- as.numeric(gsub(",", ".", Tumor.purity$CPE))

    #print(head(Tumor.purity.L))
    samples.id <- sapply(barcodes, function(x) paste(unlist(stringr::str_split(x, "-"))[1:4], collapse = "-"))

    df.barcodes_sampID <- data.frame(
        barcodes = barcodes,
        sampID = samples.id,
        row.names = 1:length(barcodes)
    )
    filt.s <- c()

    for(s in samples.id){
        if (s %in% Tumor.purity$Sample.ID == FALSE)
            filt.s <- c(filt.s, s)
    }

    if(length(filt.s)>0){
        message("the following TCGA barcodes do not have info on tumor purity:")
        filt <- as.character(df.barcodes_sampID[which(df.barcodes_sampID$sampID %in% filt.s),]$barcodes)
        print(filt)
    } else filt <- c()

    samples.filtered<-unlist(samples.id[samples.id %in% filt.s == FALSE])

    idx.samples <- which(
        Tumor.purity.L$Sample.ID %in% samples.filtered
        & (Tumor.purity.L$ESTIMATE >= estimate | Tumor.purity.L$ESTIMATE == 'NaN')
        & (Tumor.purity.L$ABSOLUTE >= absolute| Tumor.purity.L$ABSOLUTE == 'NaN')
        & (Tumor.purity.L$IHC >= ihc | Tumor.purity.L$IHC == 'NaN')
        & (Tumor.purity.L$LUMP >= lump | Tumor.purity.L$LUMP == 'NaN')
        & (Tumor.purity.L$CPE >= cpe | Tumor.purity.L$CPE == 'NaN')
    )


    df.purity <- Tumor.purity.L[idx.samples,]

    idx <- which(df.barcodes_sampID$sampID%in%df.purity$Sample.ID)

    filtered.barcodes <- as.character(df.barcodes_sampID[idx,]$barcodes)

    return(list(pure_barcodes = filtered.barcodes, filtered = filt))
}
