#' @title Download GDC data
#' @description
#'   Uses GDC API or GDC transfer tool to download gdc data
#'   The user can use query argument
#'   The data from query will be save in a folder: project/data.category
#' @param query A query for GDCquery function
#' @param token.file Token file to download controled data (only for method = "client)
#' @param method Use api method or gdc client tool (API is faster)
#' @importFrom tools md5sum
#' @importFrom utils untar
#' @import httr
#' @export
#' @return Shows the output from the GDC transfer tools
GDCdownload <- function(query, token.file, method = "api") {

    if(missing(query)) stop("Please set query argument")
    if(!(method %in% c("api","client"))) stop("method arguments possible values are: 'api' or 'client'")
    manifest <- query$results[[1]][,c("file_id","file_name","md5sum","file_size","state")]
    colnames(manifest) <- c("id","filename","md5","size","state")
    source <- ifelse(query$legacy,"legacy","harmonized")
    path <- unique(file.path(query$project, source,
                             gsub(" ","_", query$results[[1]]$data_category),
                             gsub(" ","_",query$results[[1]]$data_type)))
    # Check if the files were already downloaded by this package
    files2Download <- sapply(manifest$id, function(x) !file.exists(file.path(path,x)))
    manifest <- manifest[files2Download,]
    if(nrow(manifest) != 0 & method == "client") {
        # There exists two options to download the data, using the query or using a manifest file
        # The second option was created to let users use legacy data or the API to search

        # This will find gdc clinet, if not installed it will install it
        gdc.client.bin <- GDCclientInstall()

        # Using the query argument we will organize the files to the user
        # Creates a file with the gdc manifest format
        write_delim(manifest,"gdc_manifest.txt",delim = "\t")

        cmd <- paste0(gdc.client.bin, " download -m gdc_manifest.txt")

        if(!missing(token.file)) cmd <- paste0(cmd," -t ", token.file)

        # Download all the files in the manifest using gdc client
        message(paste0("GDCdownload will download: ",
                       humanReadableByteCount(sum(manifest$size))))
        message(paste0("Executing GDC client with the following command:\n",cmd))
        system(cmd)

        # moving the file to make it more organized
        for(i in manifest$id) move(i,file.path(path,i))

    } else if (nrow(manifest) != 0 & method =="api"){
        if(nrow(manifest) > 1) {
            name <- paste0(gsub(" |:","_",date()),".tar.gz")
            unlink(name)
            message(paste0("GDCdownload will download: ",
                           humanReadableByteCount(sum(manifest$size)),
                           " compressed in a tar.gz file"))
        } else {
            # case with one file only. This is not at tar.gz
            name <- query$results[[1]]$file_name
            message(paste0("GDCdownload will download: ",
                           humanReadableByteCount(sum(manifest$size))))
        }
        message(paste0("Downloading as: ", name))

        # Is there a better way to do it using rcurl library?
        server <- ifelse(query$legacy,"https://gdc-api.nci.nih.gov/legacy/data/", "https://gdc-api.nci.nih.gov/data/")
        body <- list(ids=list(manifest$id))
        bin <- POST(server,
                    body = body,
                    encode = "json", progress())
        writeBin(content(bin,"raw",encoding = "UTF-8"), name)

        if(nrow(manifest) > 1) {
            success <- untar(name)

            if(success != 0){
                print(success)
                stop("There was an error in the download process, please execute it again")
            }
        }
        # moving to project/data_category/data_type/file_id
        for(i in seq_along(manifest$filename)) {
            file <- manifest$filename[i]
            id <- manifest$id[i]
            # Check status
            if(!(md5sum(file) == manifest$md5[i])){
                message(paste0("File corrupted:", file))
                message("Run GDCdownload again to download it")
                unlink(file)
                next
            }
            if(file.exists(file)) move(file,file.path(path,id,file))
        }
    } else {
        message("All samples have been already downloded")
    }
}


humanReadableByteCount <- function(bytes) {
    unit <- 1000
    if (bytes < unit) return (paste0(bytes + " B"))
    exp <- floor(log(bytes) / log(unit))
    pre <- paste0(substr("KMGTPE",exp,exp))
    pre <- paste0(pre,"B")
    nb <- bytes / (unit ^ exp)
    return (paste(nb, pre))
}
GDCclientPath <- function(){
    global <- Sys.which("gdc-client")
    if(global != "") return(global)
    local <- dir(pattern = "gdc-client*[^zip]$")
    if(length(local) > 0) return(dir(pattern = "gdc-client*[^zip]$",full.names = TRUE))
    return("")
}

GDCclientExists <- function(){
    return(Sys.which("gdc-client.exe") != "" || Sys.which("gdc-client") != "" || length(dir(pattern = "gdc-client*[^zip]$") > 0))
}
#' @importFrom xml2 read_html
#' @importFrom downloader download
#' @importFrom rvest html_nodes html_attr
GDCclientInstall <- function(){
    if(GDCclientExists()) return(GDCclientPath())

    links <- read_html("https://gdc.nci.nih.gov/access-data/gdc-data-transfer-tool")  %>% html_nodes("a") %>% html_attr("href")
    bin <- links[grep("zip",links)]
    if(is.windows()) bin <- bin[grep("windows", bin)]
    if(is.mac()) bin <- bin[grep("OSX", bin)]
    if(is.linux()) bin <- bin[grep("Ubuntu", bin)]
    if(is.windows()) mode <- "wb" else  mode <- "w"
    download(paste0("https://gdc.nci.nih.gov/",bin), basename(bin), mode = mode)
    unzip(basename(bin))
    Sys.chmod("gdc-client")
    return(GDCclientPath())
}
#' @title Download the data from TCGA using as reference the output from TCGAquery
#' @description
#'      The TCGAdownload function will download the data using as reference
#'      the the lines of the TCGAquery search result.
#'
#'      There is an option to download the entire tar.gz folder or download
#'      specific files using the \emph{type} parameter or the \emph{samples}
#'      parameter
#'
#'      The outpufiles will be saved into the path parameters. If this path does
#'      not exists the package will try to create the directories.
#'
#'      By default, if a sample was already downloaded the function will not
#'      download again, unless the force parameter is set to \code{TRUE}
#'
#' @param data The TCGAquery output
#' @param path Directory to save the downloaded data
#' @param type Filter the files that will be downloaded by
#'  type. Example:"rsem.genes.results"
#' @param samples List of samples to download data
#' @param force Download files even if it was already downladed?
#' Default: \code{FALSE}
#' @examples
#' \dontrun{
#' samples <- c("TCGA-26-1442-01A-01R-1850-01")
#' query <- TCGAquery(tumor = "gbm",
#'                    platform = "IlluminaHiSeq_RNASeqV2",
#'                    level = "3",
#'                    samples = samples)
#' TCGAdownload(query,path = "RNA",
#'              samples = samples,
#'              type ="rsem.genes.results")
#' }
#' @export
#' @return Download TCGA data into the given path
#' @family data functions
TCGAdownload <- function(data = NULL, path = ".", type = NULL, samples = NULL, force = FALSE) {
    stop("TCGA data has moved from DCC server to GDC server. Please use GDCdownload function")
}

