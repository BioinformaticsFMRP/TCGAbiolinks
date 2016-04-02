#' @title Searches TCGA open-access data providing also latest version of the files.
#' @description
#'    The TCGAquery searches TCGA open-access data through a preprocessed table
#'    with the latest version of the files.
#'
#'    Some parameters as tumor, platorm, level, center, barcode can be used
#'    to search the data.
#'
#'    TCGAquery will return a matrix with the results found,
#'    that will be used in the other function TCGAdownload, TCGAprepare.
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
#' query <- TCGAquery(tumor = "gbm")
#'
#' query <- TCGAquery(tumor = c("gbm","lgg"),
#'                    platform = c("HumanMethylation450","HumanMethylation27"))
#'
#' query <- TCGAquery(tumor = "gbm",
#'                    platform = "HumanMethylation450",
#'                    level = "3")
#'
#' query <- TCGAquery(samples = "TCGA-61-1743-01A-01D-0649-04",
#'                    tumor = "OV",
#'                    platform = "CGH-1x1M_G4447A",
#'                    level = 3)
#'
#' # Get all LGG IlluminaHiSeq_RNASeqV2 data, but change with data version 11
#'
#' query <- TCGAquery(tumor = "LGG", platform = "IlluminaHiSeq_RNASeqV2", level = "3",
#'                    version = list(c("IlluminaHiSeq_RNASeqV2","LGG",11)))
#' @export
#' @importFrom downloader download
#' @importFrom stringr str_sub str_locate str_sub<-
# @importFrom knitr kable
#' @seealso
#'  \code{\link{TCGAdownload}} for downloading the data from the
#' search
#'
#' \code{\link{TCGAprepare}} for preparing the data for the user into
#' a Summarized experiment object, or a matrix.
#' @return A dataframe with the results of the query
#'        (lastest version of the files)
#' @family data functions
TCGAquery <- function(tumor = NULL,
                      platform = NULL,
                      samples = NULL,
                      center = NULL,
                      level = NULL,
                      version = NULL) {

    db <- tcga.db
    if (!is.null(tumor)) {
        for (j in seq_along(tumor)) {
            if (!(is.element(tolower(tumor[j]),
                             tolower(disease.table$abbreviation)))) {
                suppressWarnings(
                    df <- as.data.frame(matrix(
                        sort(unique(disease.table$abbreviation)),
                        ncol = 8))
                )
                print(knitr::kable(df, col.names = NULL, format = "pandoc",
                                   caption = "TCGA tumors"))
                cat("=======================================================\n")
                cat("ERROR: Disease not found. Select from the table above.\n")
                cat("=======================================================\n")
                return(NULL)
            }
        }
    }

    if (!is.null(platform)) {
        for (j in seq_along(platform)) {
            if (!(is.element(tolower(platform[j]),
                             tolower(platform.table$name)))) {
                suppressWarnings(
                    df <- as.data.frame(matrix(
                        sort(unique(platform.table$name)),
                        ncol = 3))
                )
                print(knitr::kable(df, col.names = NULL, format = "pandoc",
                                   caption = "TCGA Platforms"))
                cat("=======================================================\n")
                cat("ERROR: Platform not found. Select from the table above.\n")
                cat("=======================================================\n")
                return(NULL)
            }
        }
    }

    if (!is.null(center)) {
        if (!(is.element(tolower(center), tolower(center.table$name)))) {
            suppressWarnings(
                df <- as.data.frame(matrix(sort(unique(center.table$name)),
                                           ncol = 3))
            )
            print(knitr::kable(df, col.names = NULL, format = "pandoc",
                               caption = "TCGA Centers"))
            cat("=======================================================\n")
            cat("ERROR: Center not found. Select from the table above.\n")
            cat("=======================================================\n")
            return(NULL)
        }
    }
    if (!is.null(level)) {
        if (!(is.element(level, c("1", "2", "3","mage-tab")))) {
            message("Level not found. Chosse between:'1', '2','3','mage-tab'")
            return(NULL)
        }
    }

    if (!is.null(tumor)) {
        id <- sapply(tumor, function(x){
            grepl(x, db$Disease, ignore.case = TRUE)
        })
        id <- apply(id, 1,any)
        db <-  db[id,]
    }
    if (!is.null(platform)) {
        id <- sapply(platform, function(x){
            tolower(x) == tolower(db$Platform)
        })
        id <- apply(id, 1,any)
        db <-  db[id,]
    }
    if (!is.null(center)) {
        id <- sapply(center, function(x){
            grepl(x, db$Center, ignore.case = TRUE)
        })
        id <- apply(id, 1,any)
        db <-  db[id,]
    }
    if (!is.null(level)) {
        if (nchar(level) == 1) {
            id <- grep(paste0("Level_", level), db$name)
        } else {
            id <- grep("mage-tab", db$name)
        }

        if (length(id) > 0) {
            db <-  db[id,]
        } else {
            message("Sorry! No files with this level were found")
            return(NULL)
        }
    }

    # This is a workaround for working with old levels
    # The user should specify the tumor and disease and we will change
    # the path to get old version of that tumor/platform
    if( !is.null(version)) {
        message("Retrieving old version information, please wait...")

        root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?query=Archive"

        # This function will get the version using the api and add the barcode
        for(i in 1:length(version)){
            idx <- intersect(grep(version[[i]][1],db$baseName,ignore.case = TRUE),
                             grep(version[[i]][2],db$baseName,ignore.case = TRUE))

            if (length(idx > 0 )) {
                db <- db[-idx,]
            }

            extra <- paste0("[revision=",version[[i]][3],"]",
                            "[baseName=*",version[[i]][2],"*]",
                            "[deployLocation=*",version[[i]][1],"*]")

            if(!is.null(level)){
                extra <- paste0(extra, "[name=*Level_",level,"*]")
            }
            url <- paste0(root, extra)

            new <- tcgaGetTable(url)
            new <- new[, 1:9]
            new$addedDate <- as.Date(new$addedDate, "%m-%d-%Y")
            new <- tcgaDbAddCol(new)
            new <- new[order(new$Disease,new$Platform,new$Center),]

            if(any(grepl("Obsolete",new$deployStatus))){
                message("========== WARNING ==================")
                message("Some data is obsolete and it is no more available")
                message("For all version information go to the link below")
                message(url)
                message("=====================================")
                new <- new[grep("Available",new$deployStatus),]
            }
            colnames(new)[4] <- "barcode"

            for (j in 1:nrow(new)){
                message(paste0("Updating barcode for: ",new[j,]$name))
                new[j,"barcode"] <- updatebarcode(new[j,])
            }
            db <- rbind(db,new)
        }
    }

    # to be improved
    idx <- c()
    if(!is.null(samples)){
        for(i in seq_along(samples)){
            aux <- grep(samples[i],db$barcode)
            idx <- union(idx, aux)
        }
        db <- db[idx,]
    }

    return(db)
}



# Filter files by barcode
updatebarcode <- function(data){

    root <- "https://tcga-data.nci.nih.gov"
    url <- gsub(".tar.gz","",data$deployLocation)
    # maybe the folder does not exists, so this should be removed
    if (!url.exists(paste0(root,url))) return("ERROR")
    files <- getFileNames(paste0(root,url))
    idx <- grep("MANIFEST|README|CHANGES|DESCRIPTION|DATA_USE|Name|Size|Parent|Last",files)
    files <- files[-idx]

    barcodeName <- paste("IlluminaHiSeq_RNASeq",
                         "humanmethylation",
                         "H-miRNA_8x15K",
                         "images",
                         "SOLiD_DNASeq",
                         "pathology_reports",
                         "IlluminaDNAMethylation",
                         "HG-CGH-244A",
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


    if (grepl(uuidName,data$Platform, ignore.case = TRUE)) {
        # case uuid in name file
        regex <- paste0("[[:alnum:]]{8}-[[:alnum:]]{4}",
                        "-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}")
        uuid <- str_match(files,regex)[,1]
        map <- mapuuidbarcode(unique(na.omit(uuid)))
        ret <- paste0(map$barcode,collapse = ",")
    } else if (grepl("IlluminaGA_DNASeq_curated|illuminaga_dnaseq_automated",data$Platform) & data$Center == "broad.mit.edu") {
        # Exception - two uuids in the name
        # case uuid in name file
        regex <- paste0("[[:alnum:]]{8}-[[:alnum:]]{4}",
                        "-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}")
        files <- files[grep(regex,files)]
        uuid <- unlist(str_match_all(files,regex))
        map <- mapuuidbarcode(unique(na.omit(uuid)))
        ret <- paste0(map$barcode,collapse = ",")
    } else if(grepl(barcodeName, data$Platform, ignore.case = TRUE)) {
        regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                        "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
        barcode <- unlist(str_match_all(files,regex))
        ret <- paste0(barcode,collapse = ",")
    } else if(grepl(mageName, data$Platform, ignore.case = TRUE)) {
        ## TO BE IMPROVED
        mage <- getMage(data)
        ret <- paste0(mage$Comment..TCGA.Barcode.,collapse = ",")
    }
    return(ret)
}


#' @import utils
getBarcode <- function(table){

    root <- "https://tcga-data.nci.nih.gov"
    allFiles <- c()
    mages <-  tcga.db[grep("mage-tab",tcga.db$name),]
    max <- length(table[,1])
    pb <- txtProgressBar(min = 0, max = max , style = 3)

    for (i in seq_along(table[,1])) {
        message(i," - ",table[i,]$deployLocation)
        if (table[i,]$deployStatus == "Available") {

            # get the mage file for the entry
            mage <- subset(mages, mages$Disease == table[i,]$Disease &
                               mages$Platform == table[i,]$Platform &
                               mages$Center == table[i,]$Center)

            # if no mage file: not found (case PLatform:ABI)
            if (dim(mage)[1] == 0) {
                # For IlluminaGA and ABI platforms
                # Barcodes are in maf files
                if ((grepl("DNASeq",table[i,]$Platform) &
                     grepl("Level",table[i,]$name))
                ) {
                    barcodes <- tcga.get.barcode(table[i,])

                    if (is.null(barcodes)){
                        table[i,]$deployStatus <- "Not found"
                        next
                    }
                    table[i,]$deployStatus <- barcodes
                    setTxtProgressBar(pb, i)
                } else if (table[i,]$Platform == "diagnostic_images" |
                           table[i,]$Platform == "pathology_reports" |
                           table[i,]$Platform == "tissue_images"
                ) {
                    print("image or report")

                    folder <- gsub(".tar.gz","",table[i,]$deployLocation)
                    files <- getFileNames(paste0(root,folder))
                    if (is.null(files))
                        next
                    files <- files[grepl("TCGA",files)]
                    if (table[i,]$Platform == "pathology_reports") {
                        barcode <- substr(files,1,12)
                    } else {
                        barcode <- substr(files,1,23)
                    }
                    table[i,]$deployStatus <- paste0(unique(barcode),
                                                     collapse = ",")
                } else if (table[i,]$Platform == "ABI") {
                    print("ABI")
                    folder <- gsub(".tar.gz","",table[i,]$deployLocation)
                    files <- getFileNames(paste0(root,folder))
                    if (is.null(files))
                        next

                    maf <- files[grep("somatic.maf",files)]
                    print(maf)
                    if(length(maf) == 0){
                        table[i,]$deployStatus <- "Not found"
                        next
                    }
                    if ( !file.exists(maf)) {
                        download(paste0(root,folder,"/",maf), maf, quiet = TRUE)
                    }
                    df <- read.delim(file = maf,
                                     sep = "\t", comment.char = "#",
                                     stringsAsFactors = FALSE,
                                     fileEncoding = "latin1",
                                     row.names = NULL)
                    barcode <- union(df$TUMOR_SAMPLE_ID,
                                     df$MATCH_NORM_SAMPLE_ID)
                    if (is.null(barcode)) {
                        barcode <- union(df$Tumor_Sample_Barcode,
                                         df$Matched_Norm_Sample_Barcode)
                    }
                    if (!is.null(barcode)) {
                        table[i,]$deployStatus <- paste0(barcode,
                                                         collapse = ",")
                        unlink(maf)
                    } else {
                        table[i,]$deployStatus <- "Not found"
                    }


                } else if (table[i,]$Platform == "bio") {
                    print("bio")
                    if(grepl("Level_1",table[i,]$name)) {
                        folder <- gsub(".tar.gz","",table[i,]$deployLocation)
                        files <- getFileNames(paste0(root,folder))
                        if (is.null(files))
                            next
                        files <- files[grepl("TCGA",files)]

                        barcode <- unlist(strsplit(files,"\\."))
                        barcode <-  barcode[grep("TCGA",barcode)]
                        table[i,]$deployStatus <- paste0(unique(barcode),
                                                         collapse = ",")
                    } else {
                        folder <- gsub(".tar.gz","",table[i,]$deployLocation)
                        files <- getFileNames(paste0(root,folder))
                        if (is.null(files))
                            next
                        pat <- paste0("tumor_sample_|org_control|",
                                      "clinical_patient_|cqcf_")
                        idx <- grep(pat,files)
                        if(length(idx)>1){
                            files <- files[idx[1]]
                        } else {
                            files <- files[idx]
                        }
                        if ( !file.exists(files)) {
                            download(paste0(root,folder,"/",files),
                                     files, quiet = TRUE)
                        }
                        df <- read.delim(file = files,
                                         sep = "\t", comment.char = "#",
                                         stringsAsFactors = FALSE,
                                         fileEncoding = "latin1",
                                         row.names = NULL)
                        # first line is not barcode
                        df <- df[-1,]
                        barcode <- unique(df$bcr_patient_barcode)
                        table[i,]$deployStatus <- paste0(barcode,
                                                         collapse = ",")
                        unlink(files)
                    }
                }
                next
            }
            if(grepl("mage-tab|aux",table[i,]$name)){
                table[i,]$deployStatus <- ""
                next
            }
            # In case we have two files
            # This should not happen after filtering by center
            # probably this code can be remove until next comment
            file <- basename(mage$deployLocation)
            x <- stringr::str_replace_all(file, "[^[:alnum:]]", "")
            y <- stringr::str_replace_all(table[i,]$baseName,
                                          "[^[:alnum:]]", "")
            idx <- grep(y,x)

            if (length(idx) > 0) {
                file <- file[idx]
                mage <- mage[idx,]
            }
            # get name of files to be deleted
            allFiles <- c(allFiles, file)
            print(file)
            if ( !file.exists(file)) {
                download(paste0(root,mage$deployLocation), file, quiet = TRUE)
            }
            folder <- gsub(".tar.gz","",file)
            if ( !file.exists(folder)) {
                untar(file)
            }
            files <- list.files(folder)
            print("has mage-tab")

            # For this platform the barcodes are in array_design files
            if (table[i,]$Platform == "MDA_RPPA_Core") {
                sdrf <- files[grep("array_design",files)]
                # case with 2 array_design BLCA
                if (length(sdrf) > 1) {
                    sdrf <- sdrf[1]
                }
                df <- read.delim(file = file.path(folder,sdrf),
                                 sep = "\t",
                                 stringsAsFactors = FALSE,
                                 fileEncoding = "latin1")

                if (is.element("Sample.description",colnames(df))) {
                    barcode <- unique(as.character(df$Sample.description))
                } else if (is.element("Sample.barcode",colnames(df))){
                    barcode <- unique(as.character(df$Sample.barcode))
                } else {
                    barcode <- unique(as.character(df$Biospecimen.Barcode))
                }
                table[i,]$deployStatus <- paste0(barcode, collapse = ",")

            } else {
                # Platform is not MDA_RPPA_Core
                sdrf <- files[grep("sdrf",files)]
                df <- read.delim(file = file.path(folder,sdrf),
                                 sep = "\t",
                                 stringsAsFactors = FALSE,
                                 fileEncoding = "latin1")

                # The mage file could is common to the same
                # disease, center and platform.
                # So we're going to fill them so we don't need
                # To re-read the barcodes
                for (j in seq_along(table[,1])) {
                    if (table[j,]$Disease == table[i,]$Disease &&
                        table[j,]$Platform == table[i,]$Platform &&
                        table[j,]$Center == table[i,]$Center) {
                        aux <- grep("Comment..TCGA.Archive.Name",colnames(df))

                        barcode <- data.frame()
                        x <- table[j,]$name

                        for (z in seq_along(aux)){
                            barcode <- rbind(barcode,subset(df,
                                                            x == df[,aux[z]]))
                        }

                        barcode <- unique(as.character(
                            barcode$Comment..TCGA.Barcode.))
                        barcode <- barcode[barcode != "->"]
                        table[j,]$deployStatus <- paste0(barcode,
                                                         collapse = ",")
                    }
                }
            }
        }
        setTxtProgressBar(pb, i)
    }
    setTxtProgressBar(pb, max)
    close(pb)
    colnames(table)[4] <- "barcode"
    # Removing the files mess
    # porbably this can be used up in the code
    unlink(allFiles)
    folders <- gsub(".tar.gz","",allFiles)
    unlink(folders, recursive = TRUE)

    return(table)
}

#' @title Finds the number of downloads of a package on CRAN or BIOC and
#' find questions in website ("bioconductor.org", "biostars.org", "stackoverflow).
#' @description Finds the number of downloads of a package on CRAN or BIOC
#' @param siteToFind website ("bioconductor.org", "biostars.org", "stackoverflow)
#' related to TCGA or a package
#' @param listPackage list of package to find in bioconductor
#' @param KeyInfo is a key to find in biostars related to TGGA or package
#' @export
#' @return table with number of downloads about a package
#' @examples
#' TCGAquery_Social("bioconductor.org","BiocCheck")
TCGAquery_Social <- function(siteToFind=NULL, listPackage=NULL,KeyInfo=NULL){


    # Find all packages in bioconductor
    site3 <- "http://www.bioconductor.org/packages/3.1/bioc/"
    tmp <- .DownloadURL(site3)
    tmpPack <-  tmp[grep(".html",tmp)]
    tmpPackMatrix <- as.matrix(tmpPack)
    posA <- grep("bioconductor",tolower(tmpPackMatrix))[1]+1
    posB <- grep("bioconductor",tolower(tmpPackMatrix))[2]-1
    tmpPackMatrixNew <- tmpPackMatrix[posA:posB,]
    tmpPackMatrixNew <- gsub("<td><a","",tmpPackMatrixNew)
    tmpPck2 <- lapply(tmpPackMatrixNew, function(x) gsub("</a></td>","", as.matrix(unlist(strsplit(x,">")))[2]))
    tmpPck3 <- as.matrix(unlist(tmpPck2))
    tmpPck2a <- lapply(tmpPck2, function(x) gsub("</a","", x))
    BiocPackageList <- as.matrix(unlist(tmpPck2a))

    if( siteToFind == "bioconductor.org"){

        if(is.null(listPackage)) {
            msg  <- paste0(
                "\nPlease, provide a listofPackage argument\n",
                "Example:  TCGAquery_Social('bioconductor.org','BiocCheck')")
            stop(msg)
        }

        siteBioC <- "http://www.bioconductor.org/packages/stats/index.html"

        TablePackage <- matrix(0, length(listPackage), 2)
        rownames(TablePackage)<-listPackage
        colnames(TablePackage) <- c("Package","NumberDownload")
        TablePackage <- as.data.frame(TablePackage)
        TablePackage$Package <- listPackage

        tmp <- .DownloadURL(siteBioC)

        for ( i in 1:nrow(TablePackage)){
            packagetofind <- listPackage[i]
            pos <- grep(tolower(packagetofind), tolower(tmp ))
            pos1 <- pos[1]
            tmp3 <- tmp[pos1]
            tmp3a <- gsub("</A></TD></TR>","", as.matrix(unlist(strsplit(tmp3,"&nbsp;")))[2])
            tmp4 <- as.numeric(substr(as.character(tmp3a),2,nchar(tmp3a)-1))
            TablePackage[i,"NumberDownload"] <- tmp4
        }

        TablePackage <- TablePackage[order(TablePackage$NumberDownload,decreasing=TRUE),]
    }

    if( siteToFind == "biostars.org"){

        if(is.null(KeyInfo)) {
            msg  <- paste0(
                "\nPlease, provide a KeyInfo argument\n",
                "Example: TCGAquery_Social('biostars.org', KeyInfo='methylation')")
            stop(msg)
        }


        siteBioStar <- "https://www.biostars.org/local/search/page/?page=1&sort=update&limit=all%50time&q="
        siteQuestions <- "https://www.biostars.org/p"
        siteBioStarKey <- paste(siteBioStar,KeyInfo,sep="")

        tmp <- .DownloadURL(siteBioStarKey)
        tmp2 <- tmp[ grep("<h4>",tmp)]

        TableQuestions <- matrix(0, length(tmp2), 3)
        rownames(TableQuestions)<-paste("q",c(1:length(tmp2)),sep="")
        colnames(TableQuestions) <- c("question","BiostarsSite","PackageSuggested")
        TableQuestions <- as.data.frame(TableQuestions)

        for ( i in 1:nrow(TableQuestions)){
            #print(i)
            questiontofind <- tmp2[i]
            questiontofind <- gsub("<h4>","", questiontofind)
            qst_find_site <-gsub("<a href=","", as.matrix(unlist(strsplit(questiontofind,">")))[1])
            qst_find_site2 <-gsub("<a href=","", as.matrix(unlist(strsplit(qst_find_site,"p")))[2])

            qst_find_site2_sub <- substr(qst_find_site2,1, 7)
            TableQuestions[i,"BiostarsSite"] <- qst_find_site2_sub
            newsite_tofind <- paste(siteQuestions,qst_find_site2_sub,sep="")

            tmpPack <- .DownloadURL(newsite_tofind)

            if( length(grep("package",tolower(tmpPack)))!=0){
                pos <- grep("package",tolower(tmpPack))
                if( length(pos)!=1){
                    pos <- pos[1]
                }
                tmpPackage <-  tmpPack[pos]

                PackMat <- sapply(BiocPackageList, grepl, tmpPackage, ignore.case=TRUE)
                if(sum(PackMat)>=1){
                    print(which(PackMat == TRUE))
                    PackageSuggested <- paste(names(PackMat[which(PackMat == TRUE)]),collapse=";")
                    #TableQuestions[i,"PackageSuggested"] <- PackageSuggested
                    TableQuestions[i,"PackageSuggested"] <- substr(PackageSuggested,1, 64)
                }
                # print(PackageSuggested)
            }

            tmp3a <- gsub("</a","", as.matrix(unlist(strsplit(questiontofind,">")))[2])
            tmp3a <- gsub("&#39;", ".", tmp3a)
            TableQuestions[i,"question"] <- tmp3a
        }

        TableQuestions <- TableQuestions[order(TableQuestions$PackageSuggested,decreasing=TRUE),]
        TablePackage <- TableQuestions
    }

    if( siteToFind == "support.bioconductor.org"){

        if(is.null(KeyInfo)) {
            msg  <- paste0(
                "\nPlease, provide a KeyInfo argument\n",
                "Example: TCGAquery_Social('support.bioconductor.org', KeyInfo='TCGA')")
            stop(msg)
        }

        sitesupportBioc_part1 <- "https://support.bioconductor.org/local/search/page/?page="
        sitesupportBioc_part2 <-"&sort=New%20answers&limit=All%20time&q="
        TablePackage <-NULL


        for( pg in 1:2){
            message(paste("pag",pg),sep="")
            #sitesupportBioc <- "https://support.bioconductor.org/local/search/page/?q="
            sitesupportBioc<- paste(sitesupportBioc_part1,pg,sitesupportBioc_part2,sep="")
            siteQuestions <- "https://support.bioconductor.org/p"
            sitesupportBiocKey <- paste(sitesupportBioc,tolower(KeyInfo),sep="")

            tmp <- .DownloadURL(sitesupportBiocKey)
            tmp2 <- tmp[ grep("<h4>",tmp)]

            TableQuestions <- matrix(0, length(tmp2), 3)
            rownames(TableQuestions)<-paste("q",c(1:length(tmp2)),sep="")
            colnames(TableQuestions) <- c("question","supportBiocsSite","PackageSuggested")
            TableQuestions <- as.data.frame(TableQuestions)

            for ( tbi in 1:nrow(TableQuestions)){

                questiontofind <- tmp2[tbi]

                if(length(grep("MISSING",questiontofind))!=1){

                    questiontofind <- gsub("<h4>","", questiontofind)
                    qst_find_site <-gsub("<a href=","", as.matrix(unlist(strsplit(questiontofind,">")))[1])
                    qst_find_site2 <-gsub("<a href=","", as.matrix(unlist(strsplit(qst_find_site,"p")))[2])

                    qst_find_site2_sub <- substr(qst_find_site2,1, 7)
                    TableQuestions[tbi,"supportBiocsSite"] <- qst_find_site2_sub
                    newsite_tofind <- paste(siteQuestions,qst_find_site2_sub,sep="")

                    tmpPack <- .DownloadURL(newsite_tofind)
                    # starting from Question
                    # and removing similar post inside webpage
                    tmpPack<-tmpPack[grep("Question",tmpPack)[3]:grep("Similar",tmpPack)]


                    if( length(grep("package",tolower(tmpPack)))!=0){
                        pos <- grep("package",tolower(tmpPack))

                        if( length(pos)==1){
                            pos <- pos[1]
                            tmpPackage <-  tmpPack[pos]
                            PackMat <- sapply(BiocPackageList, grepl, tmpPackage, ignore.case=TRUE)
                        }

                        if( length(pos)!=1){
                            PackMatNew <- NULL
                            for( ip in 1: length(pos)){
                                tmpPackage <-  tmpPack[pos][ip]
                                PackMat <- sapply(BiocPackageList, grepl, tmpPackage, ignore.case=TRUE)
                                PackMatNew <- c(PackMatNew,PackMat)
                            }
                            PackMat<-PackMatNew
                        }

                        if(sum(PackMat)>=1){
                            print(which(PackMat == TRUE))
                            PackageSuggested <- paste(names(PackMat[which(PackMat == TRUE)]),collapse=";")
                            #TableQuestions[tbi,"PackageSuggested"] <- PackageSuggested
                            # print(PackageSuggested)
                            TableQuestions[tbi,"PackageSuggested"] <- substr(PackageSuggested,1, 64)

                        }
                    }

                    tmp3a <- gsub("</a","", as.matrix(unlist(strsplit(questiontofind,">")))[2])
                    tmp3a <- gsub("&#39;", ".", tmp3a)
                    TableQuestions[tbi,"question"] <- tmp3a
                }

                TableQuestions <- TableQuestions[order(TableQuestions$PackageSuggested,decreasing=TRUE),]
                TablePackage <- rbind(TablePackage,TableQuestions)
            }
        }
        TablePackage <- TablePackage[order(TablePackage$PackageSuggested,decreasing=TRUE),]
        TablePackage <- TablePackage[!duplicated(TablePackage$supportBiocsSite),]
    }

    return(TablePackage)

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

#' @title Filtering common samples among platforms from TCGAquery for the same tumor
#' @description In order to help the user to have an overview of the number of
#' samples in common we created the function `TCGAquery_integrate` that will receive the
#' data frame returned from `TCGAquery` and produce a matrix n platforms x n platforms
#' with the values of samples in commun.
#' @param query is the output of TCGAquery
#' @export
#' @return table with common samples among platforms from TCGAquery
#' @examples
#' query <- TCGAquery(tumor = 'brca',level = 3)
#' matSamples <- TCGAquery_integrate(query)
TCGAquery_integrate <- function(query) {

    querySamples <- TCGAquery_samplesfilter(query)

    matSamples <- matrix(0, length(querySamples), length(querySamples))
    colnames(matSamples) <- names(querySamples)
    rownames(matSamples) <- names(querySamples)

    for (i in 1:nrow(matSamples)) {
        for (j in 1:nrow(matSamples)) {
            if (i != j) {
                x <- as.vector(substr(querySamples[[i]], 1, 12))
                y <- as.vector(substr(querySamples[[j]], 1, 12))
                matSamples[i, j] <- length(intersect(x, y))
            } else {
                matSamples[i, j] <- length(querySamples[[i]])
            }
        }
    }
    return(matSamples)
}

#' @title Filtering sample output from TCGAquery
#' @description
#'    Filtering sample output from TCGAquery
#' @param query metaData output from TCGAquery
#' @examples query <- TCGAquery(tumor = 'brca',level = 3)
#' querySamples <- TCGAquery_samplesfilter(query)
#' @export
#' @return list of samples for a tumor
TCGAquery_samplesfilter <- function(query) {

    # Find unique platforms
    plat <- sort(unique(query$Platform))
    TumorDataList <- vector("list", length(plat))
    names(TumorDataList) <- plat

    for (idx in 1:length(TumorDataList)) {
        currPlatform <- query[query$Platform == names(TumorDataList)[idx], ]
        TumorDataList[[idx]] <- as.matrix(
            unlist(strsplit(currPlatform$barcode, ","))
        )
    }
    return(TumorDataList)
}

#' @title Get last maf file for the tumor
#' @description
#'    Filtering sample output from TCGAquery
#' @param tumor tumor type to filter the search
#' @param center Center name to filter the search
#' @param archive.name Archive name to filter the search
#' @importFrom rvest html_table
#' @importFrom xml2 read_html
#' @examples
#' \dontrun{
#'  query <- TCGAquery(tumor = 'lgg')
#' }
#' @export
#' @return list of samples for a tumor
TCGAquery_maf <- function(tumor = NULL, center = NULL, archive.name = NULL){
    message("Getting maf tables")
    message("Source: https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files")

    tables <- read_html("https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files")
    tables <-  html_table(tables)

    # Table one is junk
    tables[[1]] <- NULL

    if(!is.null(tumor)){
       # get which tables are from the tumor
       idx <- which(mapply(function(x) {
          any(grepl(tumor,(x[,1]), ignore.case = TRUE))
       },tables) == TRUE)
       df <- lapply(idx,function(x) tables[x])
    }
    # merge the data frame in the lists
    if(length(idx) > 1) {
        df <- Reduce(function(...) merge(..., all=TRUE), df)
     }  else if(length(idx) == 1) {
            df <- Reduce(function(...) merge(..., all=TRUE), df)
            df <- df[[1]]
            colnames(df) <- gsub(" ",".", colnames(df))
            colnames(df) <- gsub(":",".", colnames(df))
    } else {
	message("Sorry, no maf found")
        return (NULL)
    }

    # Remove obsolete/protected
    df <- subset(df, df$Deploy.Status == "Available")
    df <- subset(df, df$Protection.Status == "Public")

    if (!is.null(center)) df <- df[grepl(center,df[,"Archive.Name"],
                                         ignore.case = TRUE),]
    if (!is.null(archive.name)) df <- df[grepl(archive.name,df[,"Archive.Name"],
                                               ignore.case = TRUE),]

    message("We found these maf below:")
    print(df[,c(1,5,7)])

    if(nrow(df) > 1){
        x <- readline("Please, select the line that you want to download: ")
        df <- df[x,]
        if(nrow(df) > 1){
            message("Sorry, we have more than 1 maf file, please filter by the name")
            return (NULL)
        }
    }
    if(nrow(df) == 0){
        message("Sorry, no maf found")
        return (NULL)
    }


    # change the path to be downloaded
    df[,"Deploy.Location"] <- gsub("/dccfiles_prod/tcgafiles/",
                                   "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/",
                                   df[,"Deploy.Location"] )

    # We will order to get the file with more samples
    # we are considering this is the last. Needs to be checked
    #nb <- sapply(strsplit(df$Tumor.Samples.Normal.Samples,":"), function(x) x[[1]])
    #df -> df[order(nb,decreasing = F),]

    message("Downloading maf file")
    if (!file.exists(basename(df[1,]$Deploy.Location)))
        download(df[1,]$Deploy.Location,basename(df[1,]$Deploy.Location))

    suppressWarnings({
        ret <- read.table(basename(df[1,]$Deploy.Location), fill = TRUE,
                          comment.char = "#", header = TRUE, sep = "\t", quote='')
    })
    ret$bcr_patient_barcode <- substr(ret$Tumor_Sample_Barcode,1,12)

    return(ret)
}
