#' @title Download GDC data
#' @description
#'   Uses GDC API or GDC transfer tool to download gdc data
#'   The user can use query argument
#'   The data from query will be save in a folder: project/data.category
#' @param query A query for GDCquery function
#' @param token.file Token file to download controlled data (only for method = "client")
#' @param method Uses the API (POST method) or gdc client tool. Options "api", "client".
#' API is faster, but the data might get corrupted in the download, and it might need to be executed again
#' @param directory Directory/Folder where the data was downloaded. Default: GDCdata
#' @param files.per.chunk This will make the API method only download n (files.per.chunk) files at a time.
#' This may reduce the download problems when the data size is too large. Expected a integer number (example files.per.chunk = 6)
#' @importFrom tools md5sum
#' @importFrom utils untar
#' @import httr
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' # Download clinical data from XML
#' query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical")
#' GDCdownload(query, files.per.chunk = 200)
#' query <- GDCquery(
#'   project = "TARGET-AML",
#'   data.category = "Transcriptome Profiling",
#'   data.type = "miRNA Expression Quantification",
#'   workflow.type = "BCGSC miRNA Profiling",
#'   barcode = c("TARGET-20-PARUDL-03A-01R", "TARGET-20-PASRRB-03A-01R")
#' )
#' # data will be saved in:
#' # example_data_dir/TARGET-AML/harmonized/Transcriptome_Profiling/miRNA_Expression_Quantification
#' GDCdownload(query, method = "client", directory = "example_data_dir")
#' query_acc_gbm <- GDCquery(
#'   project = c("TCGA-ACC", "TCGA-GBM"),
#'   data.category = "Transcriptome Profiling",
#'   data.type = "Gene Expression Quantification",
#'   workflow.type = "STAR - Counts"
#' )
#' GDCdownload(
#'   query = query_acc_gbm,
#'   method = "api",
#'   directory = "example",
#'   files.per.chunk = 50
#' )
#' }
#' @return Shows the output from the GDC transfer tools
#' @author Tiago Chedraoui Silva
GDCdownload <- function(
        query,
        token.file,
        method = "api",
        directory = "GDCdata",
        files.per.chunk = NULL
) {
    isServeOK()
    if (missing(query)) stop("Please set query argument")

    if (!(method %in% c("api", "client"))) {
        stop("method arguments possible values are: 'api' or 'client'")
    }

    if (length(unique(getResults(query)$data_type)) > 1) {
        print(knitr::kable(sort(unique(getResults(query)$data_type)), col.names = "data_type in query"))
        stop("We can only download one data type. Please use data.type argument in GDCquery to filter results.")
    }

    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
    for (proj in unique(unlist(query$project))) {
        message("Downloading data for project ", proj)
        query.aux <- query
        results <- getResults(query.aux)[getResults(query.aux)$project == proj, ]
        query.aux$results[[1]] <- results

        manifest <- getManifest(query.aux)

        path <- unique(
            file.path(
                proj,
                gsub(" ", "_", results$data_category),
                gsub(" ", "_", results$data_type)
            )
        )
        path <- file.path(directory, path)

        # Check if the files were already downloaded by this package
        manifest <- checkAlreadyDownloaded(path, manifest)

        # There is a bug in the API, if the files has the same name it will not download correctly
        # so method should be set to client if there are files with duplicated names
        # However for clinical XML recurrent and primary are the same file. So we will ignore that case
        if (nrow(manifest) > length(unique(manifest$filename))) method <- "client"

        if (nrow(manifest) != 0 & method == "client") {
            # There exists two options to download the data, using the query or using a manifest file
            # The second option was created to let users use legacy data or the API to search

            # This will find gdc clinet, if not installed it will install it
            gdc.client.bin <- GDCclientInstall()

            # Using the query argument we will organize the files to the user
            # Creates a file with the gdc manifest format
            readr::write_delim(manifest, "gdc_manifest.txt", delim = "\t")

            readr::write_delim(manifest, "gdc_client_configuration.dtt", delim = "\t")
            readr::write_lines(
                c("[download]", "retry_amount = 6", paste0("dir =", path)),
                file = "gdc_client_configuration.dtt"
            )
            cmd <- paste0(gdc.client.bin, " download -m gdc_manifest.txt --config gdc_client_configuration.dtt")

            dir.create(path, recursive = TRUE, showWarnings = FALSE)
            if (!missing(token.file)) cmd <- paste0(cmd, " -t ", token.file)

            # Download all the files in the manifest using gdc client
            message(paste0(
                "GDCdownload will download: ",
                humanReadableByteCount(sum(as.numeric(manifest$size)))
            ))
            message(paste0("Executing GDC client with the following command:\n", cmd))
            result <- tryCatch(
                {
                    system(cmd)
                },
                warning = function(w) {
                },
                error = function(e) {
                }
            )
        } else if (nrow(manifest) != 0 & method == "api") {
            if (nrow(manifest) > 1) {
                name <- paste0(gsub(" |:", "_", date()), ".tar.gz")
                unlink(name)
                message(
                    paste0(
                        "GDCdownload will download ", nrow(manifest), " files. A total of ",
                        humanReadableByteCount(sum(as.numeric(manifest$size)))
                    )
                )
            } else {
                # case with one file only. This is not at tar.gz
                name <- manifest$filename
                message(
                    paste0(
                        "GDCdownload will download: ",
                        humanReadableByteCount(sum(as.numeric(manifest$size)))
                    )
                )
            }

            server <- "https://api.gdc.cancer.gov/data/"

            if (is.null(files.per.chunk) & sum(as.numeric(manifest$size)) > 10^9) {
                message("The total size of files is big. We will download files in chunks")
                files.per.chunk <- floor(10^9 / mean(as.numeric(manifest$size)))
            }

            if (is.null(files.per.chunk)) {
                message(paste0("Downloading as: ", name))
                tryCatch(
                    {
                        GDCdownload.aux(server, manifest, name, path)
                    },
                    error = function(e) {
                        message("Download failed. We will retry with smaller chunks")
                        # split in groups of 100 MB
                        manifest <- checkAlreadyDownloaded(path, manifest)
                        step <- ceiling(100000000 / manifest$size[1])
                        if (step == 0) step <- 1
                        GDCdownload.by.chunk(server, manifest, name, path, step)
                    }
                )
            } else {
                step <- files.per.chunk
                # If error we will try another time.
                tryCatch(
                    {
                        GDCdownload.by.chunk(server, manifest, name, path, step)
                    },
                    error = function(e) {
                        message("At least one of the chunks download was not correct. We will retry")
                        manifest <- checkAlreadyDownloaded(path, manifest)
                        GDCdownload.by.chunk(server, manifest, name, path, step)
                    }
                )
            }
        } else {
            message("All samples have been already downloaded")
        }
    }
}

#' @title Get a Manifest from GDCquery output that can be used with GDC-client
#' @description
#' Get a Manifest from GDCquery output that can be used with GDC-client
#' @param query A query for GDCquery function
#' @param save Write Manifest to a txt file (tab separated)
#' @examples
#' query <- GDCquery(
#'   project = "TARGET-AML",
#'   data.category = "Transcriptome Profiling",
#'   data.type = "Gene Expression Quantification",
#'   workflow.type = "STAR - Counts",
#'   barcode = c("TARGET-20-PADZCG-04A-01R", "TARGET-20-PARJCR-09A-01R")
#' )
#' getManifest(query)
#' @export
getManifest <- function(query, save = FALSE) {
    manifest <- query$results[[1]][, c("file_id", "file_name", "md5sum", "file_size", "state")]
    colnames(manifest) <- c("id", "filename", "md5", "size", "state")
    if (save) {
        fname <- "gdc_manifest.txt"
        readr::write_delim(manifest, fname, delim = "\t")
        file <- file.path(getwd(), fname)
        message("Manifest saved as: ", file)
    }
    return(manifest)
}

GDCdownload.by.chunk <- function(
        server = "https://api.gdc.cancer.gov/data/",
        manifest,
        name = "TCGAbiolinks_download",
        path = ".",
        step = 1
) {
    for (idx in 0:ceiling(nrow(manifest) / step - 1)) {
        end <- ifelse(((idx + 1) * step) > nrow(manifest), nrow(manifest), ((idx + 1) * step))
        manifest.aux <- manifest[((idx * step) + 1):end, ]
        size <- humanReadableByteCount(sum(as.numeric(manifest.aux$size)))

        if (nrow(manifest.aux) > 1) {
            name.aux <- gsub("\\.tar", paste0("_", idx, ".tar"), name)
        } else {
            name.aux <- manifest.aux$filename
        }
        message(
            paste0(
                "Downloading chunk ", idx + 1, " of ", ceiling(nrow(manifest) / step),
                " (", nrow(manifest.aux), " files, size = ", size, ") ",
                "as ", name.aux
            )
        )
        repeat {
            ret <- GDCdownload.aux(server, manifest.aux, name.aux, path)
            if (ret == 1) break
        }
    }
}

GDCdownload.aux <- function(
        server = "https://api.gdc.cancer.gov/data/",
        manifest,
        name = "TCGAbiolinks_download",
        path = "."
) {
    result <- tryCatch(
        {
            bin <- getURL(
                server,
                POST,
                body = list(ids = list(manifest$id)),
                encode = "json",
                progress()
            )
            if (bin[[2]] == "405") {
                message("ERROR accessing GDC. Trying again...")
                bin <- getURL(
                    "https://api.gdc.cancer.gov/data/",
                    POST,
                    body = list(ids = list(manifest$id)),
                    encode = "json",
                    progress()
                )
            }
            writeBin(getURL(bin, content, as = "raw", encoding = "UTF-8"), name)

            if (nrow(manifest) > 1) {
                success <- untar(name)
                unlink(name) # remove tar
                if (success != 0) {
                    stop("There was an error in the download process, please execute it again")
                    return(-1)
                }
            }
            # moving to project/source/data_category/data_type/file_id
            for (i in seq_along(manifest$filename)) {
                if (nrow(manifest) > 1) file <- file.path(manifest$id[i], manifest$filename[i])
                if (nrow(manifest) == 1) file <- file.path(manifest$filename[i])
                id <- manifest$id[i]

                # Check status
                if (!(md5sum(file) == manifest$md5[i])) {
                    message(paste0("File corrupted:", file))
                    message("Run GDCdownload again to download it")
                    unlink(file)
                    next
                }
                if (nrow(manifest) > 1) {
                    move(file, file.path(path, file))
                }
                if (nrow(manifest) == 1) {
                    move(file, file.path(path, id, file))
                }
            }
            return(1)
        },
        warning = function(w) {
            return(1)
        },
        error = function(e) {
            unlink(name) # remove tar
            return(-1)
        }
    )

    if (result == -1) {
        stop(
            paste0(
                "There was an error in the download process (we might had a connection problem with GDC server).",
                "\nPlease run this function it again.",
                "\nTry using method = `client` or setting files.per.chunk to a small number."
            )
        )
    }
    message("Download completed")
}

humanReadableByteCount <- function(bytes) {
    unit <- 1000
    if (bytes < unit) {
        return(paste0(bytes + " B"))
    }
    exp <- floor(log(bytes) / log(unit))
    pre <- paste0(substr("KMGTPE", exp, exp))
    pre <- paste0(pre, "B")
    nb <- bytes / (unit^exp)
    return(paste(nb, pre))
}
GDCclientPath <- function() {
    global <- Sys.which("gdc-client")
    if (global != "") {
        return(global)
    }
    local <- dir(pattern = "gdc-client*[^zip]$")
    if (length(local) > 0) {
        return(dir(pattern = "gdc-client*[^zip]$", full.names = TRUE))
    }
    return("")
}

GDCclientExists <- function() {
    return(Sys.which("gdc-client.exe") != "" || Sys.which("gdc-client") != "" || length(dir(pattern = "gdc-client*[^zip]$") > 0))
}
#' @importFrom xml2 read_html
#' @importFrom downloader download
#' @importFrom rvest html_nodes html_attr %>%
GDCclientInstall <- function() {
    if (GDCclientExists()) {
        return(GDCclientPath())
    }

    links <- tryCatch(
        {
            read_html("https://gdc.cancer.gov/access-data/gdc-data-transfer-tool") %>%
                html_nodes("a") %>%
                html_attr("href")
        },
        error = function(e) {
            c(
                "https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_Ubuntu_x64-py3.8-ubuntu-20.04.zip",
                "https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_Windows_x64-py3.8-windows-2019.zip",
                "https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_OSX_x64-py3.8-macos-14.zip"
            )
        }
    )
    bin <- links[grep("public.*zip", links)]
    if (is.windows()) bin <- bin[grep("client*.*windows", bin, ignore.case = TRUE)]
    if (is.mac()) {
        if(Sys.info()["machine"] == "arm64"){
            bin <- bin[grep("client*.*OSX.*14", bin)]
        } else {
            bin <- bin[grep("client*.*OSX.*12", bin)]
        }
    }
    if (is.linux()) {
        if (grepl("ubuntu", Sys.info()["version"], ignore.case = TRUE)) {
            bin <- bin[grep("client*.*Ubuntu", bin)]
        } else {
            bin <- bin[grep("client*.*Cent", bin)]
        }
    }
    if (is.windows()) mode <- "wb" else mode <- "w"
    download(bin, basename(bin), mode = mode)
    unzip(basename(bin))
    if(is.mac()){
        unzip(gsub("-py3.8-macos-14|-py3.8-macos-12","",basename(bin)))
    }
    Sys.chmod("gdc-client")
    return(GDCclientPath())
}

checkAlreadyDownloaded <- function(path, manifest) {
    files2Download <- !(file.exists(file.path(path, manifest$id, manifest$filename)) | file.exists(file.path(path, manifest$filename)))
    if (any(files2Download == FALSE)) {
        message(
            "Of the ", nrow(manifest), " files for download ",
            table(files2Download)["FALSE"], " already exist."
        )
        if (any(files2Download == TRUE)) {
            message("We will download only those that are missing ones.")
        }
    }
    return(manifest[files2Download, ])
}
