#' @title Download GDC data
#' @description
#'   Uses GDC transfer tool to download gdc data
#'   The user can use query or manifest.file argument (not both at the same time)
#'   The data from query will be save in a folder: project/data.category
#'   The data from manifest.file  will be save in a folder: manifest.file.without.extension/
#' @param query A query for GDCquery function
#' @param manifest.file A manifest file from the GDC data portal
#' @export
#' @return Shows the output from the GDC transfer tools
GDCDownload <- function(query, token.file) {

    # There exists two options to download the data, using the query or using a manifest file
    # The second option was created to let users use legacy data or the API to search
    if(missing(query)) stop("Please set queryargument")

    # This will find gdc clinet, if not installed it will install it
    gdc.client.bin <- GDCclientInstall()

    # Using the query argument we will organize the files to the user

    # Creates a file with the gdc manifest format
    manifest <- query$results[[1]][,c("file_id","file_name","md5sum","file_size","state")]
    colnames(manifest) <- c("id","filename","md5","size","state")
    path <- unique(file.path(query$project,
                             gsub(" ","_", query$results[[1]]$data_category),
                             gsub(" ","_",query$results[[1]]$data_type)))

    # Check if the files were already downloaded by this package
    files2Download <- sapply(manifest$id, function(x) !file.exists(file.path(path,x)))
    manifest <- manifest[files2Download,]
    write_delim(manifest,"gdc_manifest.txt",delim = "\t")

    cmd <- paste0(gdc.client.bin, " download -m gdc_manifest.txt")

    if(!missing(token.file)) cmd <- paste0(cmd," -t ", token.file)

    # Download all the files in the manifest using gdc client
    message(paste0("Executing GDC client with the following command:\n",cmd))
    system(cmd)

    # moving the file to make it more organized
    for(i in manifest$id) move(i,file.path(path,i))
}


GDCclientPath <- function(){
    global <- Sys.which("gdc-client")
    if(global != "") return(global)
    local <- dir(pattern = "gdc-client*[^zip]$")
    if(local == "gdc-client") return(dir(pattern = "gdc-client*[^zip]$",full.names = TRUE))
    return("")
}

GDCclientExists <- function(){
    return(Sys.which("gdc-client") != "" | dir(pattern = "gdc-client*[^zip]$") == "gdc-client")
}
#' @importFrom xml2 read_html
#' @importFrom downloader download
#' @importFrom rvest html_nodes html_attr
GDCclientInstall <- function(){
    if(length(GDCclientExists())) return(GDCclientPath())

    links <- read_html("https://gdc.nci.nih.gov/access-data/gdc-data-transfer-tool")  %>% html_nodes("a") %>% html_attr("href")
    bin <- links[grep("zip",links)]
    if(is.windows()) bin <- bin[grep("windows", bin)]
    if(is.mac()) bin <- bin[grep("OSX", bin)]
    if(is.linux()) bin <- bin[grep("Ubuntu", bin)]

    download(paste0("https://gdc.nci.nih.gov/",bin), basename(bin))
    unzip(basename(bin))
    Sys.chmod("gdc-client")
    return(dir(pattern = "gdc-client*[^zip]$",full.names = TRUE))
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
#' @seealso \code{\link{TCGAquery}} for searching the data to download
#'
#' \code{\link{TCGAprepare}} for preparing the data for the user into
#' a Summarized experiment object, or a matrix.
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
TCGAdownload <- function(data = NULL, path = ".", type = NULL, samples = NULL,
                         force = FALSE) {
    stop("TCGA data has moved from DCC server to GDC server. Please use GDCdownload function")
}

# Filter files by barcode
#' @importFrom stringr str_extract
filterFiles <- function(data,samples,files){

    # If it if maf it is better to let download all files
    #maf <- paste("SOLiD_DNASeq_curated",
    #             "SOLiD_DNASeq", # partial barcode 2
    #             "illuminaga_dnaseq",
    #             sep = "|"
    #)

    barcodeName <- paste("IlluminaHiSeq_RNASeq",
                         #"humanmethylation",
                         "H-miRNA_8x15K",
                         "images",
                         "SOLiD_DNASeq",
                         "pathology_reports",
                         "IlluminaDNAMethylation",
                         #"HG-CGH-244A", # exception depends on the center (harvard - barcode name)
                         # mskcc.org in the mage
                         "HG-CGH-415K_G4124A",
                         "HG-U133_Plus_2",
                         "IlluminaGA_DNASeq_automated",
                         "IlluminaGA_miRNASeq",
                         "IlluminaGA_mRNA_DGE",
                         "IlluminaGA_RNASeq",
                         "IlluminaHiSeq_DNASeqC",
                         "IlluminaHiSeq_miRNASeq",
                         "IlluminaHiSeq_RNASeq", sep = "|")

    uuidName <- paste("RNASeqV2",
                      "MDA_RPPA_Core",
                      sep = "|")

    mageName <-  paste("AgilentG4502A",
                       "CGH-1x1M_G4447A",
                       "Genome_Wide_SNP_6",
                       "HT_HG-U133A",
                       "IlluminaHiSeq_WGBS",
                       sep = "|")


    level <- as.numeric(unique(substr(str_extract(data$name,"Level_[1-3]"),7,7)))

    if (grepl(uuidName,data$Platform, ignore.case = TRUE)) {
        # case uuid in name file
        regex <- paste0("[[:alnum:]]{8}-[[:alnum:]]{4}",
                        "-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}")
        uuid <- str_match(files,regex)[,1]
        map <- mapuuidbarcode(unique(na.omit(uuid)))
        idx <- unique(unlist(lapply(samples, function(x) grep(x,map$barcode))))
        idx <- which(tolower(uuid) %in% map[idx,]$uuid)
        files <- files[idx]
    } else if (grepl("IlluminaGA_DNASeq_curated|illuminaga_dnaseq_automated",data$Platform) & data$Center == "broad.mit.edu") {
        # Exception - two uuids in the name
        # case uuid in name file
        regex <- paste0("[[:alnum:]]{8}-[[:alnum:]]{4}",
                        "-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}")
        files <- files[grep(regex,files)]
        uuid <- unlist(str_match_all(files,regex))
        map <- mapuuidbarcode(unique(na.omit(uuid)))
        idx <- unique(unlist(lapply(samples, function(x) grep(x,map$barcode))))
        idx <- which(uuid %in% map[idx,]$uuid)
        files <- files[ceiling(idx / 2)]
    } else if(grepl(barcodeName, data$Platform, ignore.case = TRUE)
              | (grepl("humanmethylation", data$Platform, ignore.case = TRUE) & level != 1 )
              | (grepl("HG-CGH-244A", data$Platform, ignore.case = TRUE) & data$Center == "hms.harvard.edu" )) {
        idx <- unique(unlist(lapply(samples, function(x) grep(x,files))))
        files <- files[idx]
    } else if(grepl(mageName, data$Platform, ignore.case = TRUE)
              | (grepl("humanmethylation", data$Platform, ignore.case = TRUE) & level == 1 )
              | (grepl("HG-CGH-244A", data$Platform, ignore.case = TRUE) & data$Center == "mskcc.org" )) {
        mage <- getMage(data)
        idx <- unlist(lapply(samples,
                             function(x) grep(x,mage$Comment..TCGA.Barcode.)))
        idx <- unique(idx)
        mage <- mage[idx,]
        idx <- grep("Derived.Array.Data.Matrix.File|Array.Data.File|Derived.Data.File",
                    colnames(mage))
        names <- unique(unlist(mage[,idx]))
        idx <- unique(unlist(lapply(names, function(x) grep(x,files))))
        files <- files[idx]
    }
    return(files)
}
