
#' @title Retrieve summary of files per sample in a project
#' @description
#'  Retrieve the numner of files under each
#'   data_category + data_type + experimental_strategy + platform
#'   Almost like https://portal.gdc.cancer.gov/exploration
#' @param project A GDC project
#' @param legacy Access legacy database ? Deafult: FALSE
#' @param files.access Filter by file access ("open" or "controlled").
#' Default: no filter
#' @export
#' @examples
#'    summary <- getSampleFilesSummary("TCGA-UCS")
#' \dontrun{
#'    summary <- getSampleFilesSummary(c("TCGA-OV","TCGA-ACC"))
#' }
#' @return A data frame with the maf file information
#' @importFrom tidyr spread unite
#' @importFrom plyr ldply count
getSampleFilesSummary <- function(project, legacy = FALSE, files.access = NA) {
    out <- NULL

    for(proj in project){
        checkProjectInput(proj)
        message("Accessing information for project: ", proj)
        url <- getSampleSummaryUrl(proj,legacy)
        x <- getURL(url,fromJSON,simplifyDataFrame = TRUE)
        y <- x$data$hits$files
        names(y) <- x$data$hits$submitter_id
        df <- ldply(y, data.frame)

        df <- df %>%
            unite("type",
                  "data_category",
                  "data_type",
                  "experimental_strategy",
                  "platform",
                  na.rm = TRUE) %>%
            plyr::count(c(".id","type")) %>%
            tidyr::spread("type", "freq")
        df$project <- proj
        df[is.na(df)] <- 0
        out <- rbind.fill(out,df)
    }
    return(out)
}

getSampleSummaryUrl <- function(project,legacy = FALSE, files.access = NA){
    # Get manifest using the API
    baseURL <- ifelse(legacy,"https://api.gdc.cancer.gov/legacy/cases/?","https://api.gdc.cancer.gov/cases/?")

    options.pretty <- "pretty=true"
    options.expand <- "expand=summary,summary.data_categories,files"
    #option.size <- paste0("size=",getNbFiles(project,data.category,legacy))
    option.size <- paste0("size=",1000)
    option.format <- paste0("format=JSON")

    options.filter <- paste0("filters=",
                             URLencode('{"op":"and","content":['),  # Start json request
                             URLencode('{"op":"in","content":{"field":"cases.project.project_id","value":["'),
                             project,
                             URLencode('"]}}'))

    if(!any(is.na(files.access))) {
        options.filter <- paste0(options.filter,addFilter("files.access", files.access))
    }
    # Close json request
    options.filter <- paste0(options.filter, URLencode(']}'))
    url <- paste0(baseURL,paste(options.pretty,
                                options.expand,
                                option.size,
                                options.filter,
                                option.format,
                                sep = "&"))
    return(url)
}



getSubmitterIDUrl <- function(project,legacy = FALSE, files.access = NA){
    # Get manifest using the API
    baseURL <- ifelse(legacy,"https://api.gdc.cancer.gov/legacy/cases/?","https://api.gdc.cancer.gov/cases/?")

    options.pretty <- "pretty=true"
    options.expand <- "expand=files.access"
    #option.size <- paste0("size=",getNbFiles(project,data.category,legacy))
    option.fields = "fields=submitter_id"
    option.size <- paste0("size=",1000)
    option.format <- paste0("format=JSON")

    options.filter <- paste0("filters=",
                             URLencode('{"op":"and","content":['),  # Start json request
                             URLencode('{"op":"in","content":{"field":"cases.project.project_id","value":["'),
                             project,
                             URLencode('"]}}'))

    if(!any(is.na(files.access))) {
        options.filter <- paste0(options.filter,addFilter("files.access", files.access))
    }
    # Close json request
    options.filter <- paste0(options.filter, URLencode(']}'))
    url <- paste0(baseURL,paste(options.pretty,
                                options.expand,
                                option.fields,
                                option.size,
                                options.filter,
                                option.format,
                                sep = "&"))
    return(url)
}

# getSubmitterID("TCGA-BRCA")
# getSubmitterID("MMRF-COMPASS")
getSubmitterID <- function(project,legacy = FALSE, files.access = NA){

    url <- getSubmitterIDUrl(project,legacy,files.access)

    json  <- tryCatch(
        getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
        error = function(e) {
            message(paste("Error: ", e, sep = " "))
            message("We will retry to access GDC!")
            fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
        }
    )
    return(unique(json$data$hits$submitter_id))



}

# getBarcodefromAliquot(c("4e06e279-5f0d-4bf5-8659-67b8069050b8","bb6e1801-b08a-49b1-bc4b-205fdefb035b"))
#' @importFrom dplyr bind_rows
getBarcodefromAliquot <- function(aliquot){
    baseURL <- "https://api.gdc.cancer.gov/cases/?"
    options.fields <- "fields=samples.portions.analytes.aliquots.aliquot_id,samples.portions.analytes.aliquots.submitter_id"
    options.pretty <- "pretty=true"
    option.size <- paste0("size=",length(aliquot))
    #message(paste(barcode,collapse = '","'))
    #message(paste0('"',paste(barcode,collapse = '","')))
    options.filter <- paste0("filters=",
                             URLencode('{"op":"and","content":[{"op":"in","content":{"field":"samples.portions.analytes.aliquots.aliquot_id","value":['),
                             paste0('"',paste(aliquot,collapse = '","')),
                             URLencode('"]}}]}'))
    #message(paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&")))
    url <- paste0(baseURL,paste(options.pretty,options.fields, option.size, options.filter, sep = "&"))
    #message(url)
    json  <- tryCatch(
        getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
        error = function(e) {
            message(paste("Error: ", e, sep = " "))
            message("We will retry to access GDC again! URL:")
            #message(url)
            fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
        }
    )
    results <- json$data$hits
    if(length(results) == 0){
        message("aliquot_id not found")
        return(NULL)
    }

    results <- plyr::ldply(results$samples,.fun = function(x){
        plyr::ldply(x$portions,.fun = function(y){
            plyr::ldply(y$analytes,.fun = function(z){
                z$aliquots %>% bind_rows()
            })
        })
    })

    idx <- is.na(match(aliquot,results$aliquot_id))
    if(any(idx)){
        message("Some aliquot IDs were not found:",aliquot[idx])
    }
    results <- results[match(aliquot,results$aliquot_id),] %>% na.omit()

    return(results)
}

#' internal function to break a huge API call into smaller ones
#' so it repects the max character limit of a string
#' @param items vector of items to be using within the function
#' (list of barcodes, aliquot ids, etc)
#' @param step How many items to be evaluated per API call
#' @param FUN function that calls the API
splitAPICall <- function(FUN, step = 20, items){
   info <- NULL
   info <- tryCatch({
        for(i in 0:(ceiling(length(items)/step) - 1)){
            start <- 1 + step * i
            end <- ifelse(((i + 1) * step) > length(items), length(items),((i + 1) * step))
            if(is.null(info)) {
                info <- FUN(items[start:end])
            } else {
                info <- plyr::rbind.fill(info, FUN(items[start:end]))
            }
        }
       info
    }, error = function(e) {
        step <- 2
        for(i in 0:(ceiling(length(items)/step) - 1)){
            start <- 1 + step * i
            end <- ifelse(((i + 1) * step) > length(items), length(items),((i + 1) * step))
            if(is.null(info)) {
                info <- FUN(items[start:end])
            } else {
                info <- plyr::rbind.fill(info, FUN(items[start:end]))
            }
        }
    })
   info
}


#' Create a Summary table for each sample in a project saying if it contains
#' or not files for a certain data category
#' @description
#' Create a Summary table for each sample in a project saying if it contains
#' or not files for a certain data category
#' @param project A GDC project
#' @param legacy Access legacy (hg19) or harmonized database (hg38).
#' @return A data frame
#' @export
#' @importFrom stats xtabs
#' @examples
#' summary <- getDataCategorySummary("TCGA-ACC", legacy = TRUE)
getDataCategorySummary <- function(project, legacy = FALSE){
    baseURL <- ifelse(legacy,"https://api.gdc.cancer.gov/legacy/files/?","https://api.gdc.cancer.gov/files/?")
    url <- paste0(baseURL,"&expand=cases&size=100000&fields=cases.submitter_id,data_category&filters=",
                  URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["'),
                  URLencode(project),
                  URLencode('"]}}]}'))

    json  <- tryCatch(
        getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
        error = function(e) {
            fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
        }
    )
    json <- json$data$hits
    json$submitter_id <- unlist(lapply(json$cases, function(x) paste0(x$submitter_id,collapse = ",")))
    json$cases <- NULL
    json <- json[!duplicated(json),]
    json <- json[stringr::str_length(json$submitter_id) == 12,]
    ret <- as.data.frame.matrix(xtabs(~ submitter_id + data_category , json))
    return(ret)
}

#' @title Get Project Summary from GDC
#' @param project A  GDC project
#' @param legacy Select between Harmonized or Legacy database
#' @examples
#' \dontrun{
#' getProjectSummary("TCGA-ACC")
#' getProjectSummary("CPTAC-2")
#' }
getProjectSummary <- function(project, legacy = FALSE){
    baseURL <- ifelse(legacy,"https://api.gdc.cancer.gov/legacy/projects/","https://api.gdc.cancer.gov/projects/")
    url <- paste0(baseURL, project,"?expand=summary,summary.data_categories&pretty=true")
    return(fromJSON(url,simplifyDataFrame = TRUE)$data$summary)
}

#' @title Get Number of cases in GDC for a project
#' @param project A  GDC project
#' @param data.category A  GDC project data category
#' @param legacy Select between Harmonized or Legacy database
#' @examples
#' \dontrun{
#' getNbCases("TCGA-ACC","Clinical")
#' getNbCases("CPTAC-2","Clinical")
#' }
getNbCases <- function(project, data.category, legacy = FALSE){
    summary <- getProjectSummary(project, legacy)
    if(data.category %in% summary$data_categories$data_category){
        summary <- getProjectSummary(project, legacy)$data_categories
        nb <- summary[summary$data_category == data.category,"case_count"]
    } else {
        nb <- summary$case_count
    }
    return(nb)
}

#' @title Get Number of files in GDC for a project
#' @param project A  GDC project
#' @param data.category A  GDC project data category
#' @param legacy Select between Harmonized or Legacy database
#' @examples
#' \dontrun{
#' getNbFiles("TCGA-ACC","Clinical")
#' getNbFiles("CPTAC-2","Clinical")
#' }
getNbFiles <- function(project, data.category, legacy = FALSE){
    summary <- getProjectSummary(project, legacy)
    if(data.category %in% summary$data_categories$data_category){
        summary <- getProjectSummary(project, legacy)$data_categories
        nb <- summary[summary$data_category == data.category,"file_count"]
    } else {
        nb <- summary$file_count
    }
    return(nb)
}

