
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
#'    summary <- getSampleFilesSummary("TCGA-LUAD")
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
            unite(type, data_category, data_type, experimental_strategy, platform, na.rm = TRUE) %>%
            plyr::count(c(".id","type")) %>%
            tidyr::spread(type, freq)
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
        x$portions[[1]]$analytes[[1]]$aliquots %>% bind_rows()
    })

    results <- results$submitter_id[match(aliquot,results$aliquot_id)]

    return(results)
}


