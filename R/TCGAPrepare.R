#' @title TCGAPrepare
#' @description
#'  This function will read the data from level 3 the experiments and prepare it
#'  for downstream analysis into a SummarizedExperiment object.
#'
#'  The samples are always refered by their barcode.
#'
#'  If you want to save the data into an rda file, please use the \emph{save}
#'  parameter that will save an rda object with the  \emph{filename} parameter.
#'  If no filename was set, the filename will be the concatenation of platform and
#'  Sys.time.
#'
#' List of accepted platforms:
#'\itemize{
#' \item AgilentG4502A_07_1/AgilentG4502A_07_2/AgilentG4502A_07_3
#' \item Genome_Wide_SNP_6
#' \item H-miRNA_8x15K/H-miRNA_8x15Kv2
#' \item HG-U133_Plus_2
#' \item HT_HG-U133A
#' \item HumanMethylation27
#' \item HumanMethylation450
#' \item IlluminaDNAMethylation_OMA002_CPI
#' \item IlluminaDNAMethylation_OMA003_CPI
#' \item IlluminaGA_RNASeq
#' \item IlluminaGA_RNASeqV2
#' \item IlluminaHiSeq_RNASeq
#' \item IlluminaHiSeq_RNASeqV2
#' \item IlluminaHiSeq_TotalRNASeqV2
#' }
#'  \strong{Return}The default output is a SummarizedExperiment object.
#'
#' @return A SummarizedExperiment object (If SummarizedExperiment = \code{FALSE},
#' a data.frame)
#' @param query TCGAQuery output
#' @param dir Directory with the files downloaded by TCGADownload
#' @param samples List of samples to prepare the data
#' @param type Filter the files to prepare.
#' @param save Save a rda object with the prepared object?
#'  Default: \code{FALSE}
#' @param filename Name of the saved file
#' @param toPackage For whihc package are you preparing the data?
#' Default: NULL. Options: "ELMER"
#' @param summarizedExperiment Output as SummarizedExperiment?
#' Default: \code{FALSE}
#' @examples
#' sample <- "TCGA-06-0939-01A-01D-1228-05"
#' query <- TCGAQuery(tumor = "GBM",samples = sample, level = 3)
#' TCGADownload(query,path = "exampleData",samples = sample, quiet = TRUE)
#' data <- TCGAPrepare(query, dir="exampleData")
#' @export
#' @importFrom stringr str_match str_trim str_detect str_match_all
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom limma alias2SymbolTable
#' @importFrom GenomicFeatures microRNAs
#' @importFrom BiocGenerics as.data.frame
#' @importFrom GenomicRanges GRanges nearest
#' @importFrom IRanges IRanges
#' @importFrom R.oo abort
#' @import utils data.table TxDb.Hsapiens.UCSC.hg19.knownGene
#' @seealso  \code{\link{TCGAQuery}} for searching the data to download
#'
#'  \code{\link{TCGADownload}} for downloading the data from the
#' search
#' @family data functions
TCGAPrepare <- function(query,
                        dir = NULL,
                        samples = NULL,
                        type = NULL,
                        save = FALSE,
                        filename = NULL,
                        toPackage = NULL,
                        summarizedExperiment = TRUE){

    if (is.null(dir)) {
        message("Argument dir is NULL. Plese provide the directory
                with the folders to be prepared. ")
        return(NULL)
    }

    if (length(unique(query$Platform)) > 1 |
        length(unique(query$Center)) > 2) {
        message("Sorry! But, for the moment, we can only prepare on type of
                platform per call")
        return(NULL)
    } else {
        platform <- unique(query$Platform)
    }

    # Get all files from directory except MANIFEST, README, etc
    files <- NULL
    dirs <- gsub(".tar.gz","",basename(query$deployLocation))
    for (i in seq_along(dirs)) {
        aux <- list.files(file.path(dir,dirs[i]), full.names = TRUE,
                          recursive = TRUE)
        files <- c(files, aux )
    }
    idx <- grep("MANIFEST|README|CHANGES|DESCRIPTION|DATA_USE",files)
    if (length(idx) > 0) {
        files <- files[-idx]
    }

    # Filter by type
    if (!is.null(type)) {
        files <- files[grep(type,files)]
        if(length(files) == 0){
            message("No files of that type found")
            return(NULL)
        }
    }

    # Filter by samples
    if (!is.null(samples)) {
        files <- filterFiles(query[i,],samples,files)
        if(length(files) == 0){
            message("No files for that samples found")
            return(NULL)
        }
    }


    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    df <- NULL
    sset <- NULL

    if (grepl("humanmethylation",tolower(platform))) {

        regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                        "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
        barcode <- str_match(files,regex)

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE,skip = 1,
                          colClasses=c("character", # Composite Element REF
                                       "numeric",   # beta value
                                       "character", # Gene symbol
                                       "character",   # Chromosome
                                       "integer"))  # Genomic coordinate
            setnames(data,gsub(" ", "\\.", colnames(data)))
            setnames(data,2,barcode[i])

            if (i == 1) {
                setcolorder(data,c(1, 3:5, 2))
                df <- data
            } else {
                data <- subset(data,select = c(1,2))
                df <- merge(df, data, by = "Composite.Element.REF")
            }

            setTxtProgressBar(pb, i)
        }

        if (summarizedExperiment) {

            gene.GR <- GRanges(seqnames = paste0("chr", gene.location$chromosome_name),
                               ranges = IRanges(start = gene.location$start_position,
                                                end = gene.location$end_position),
                               strand = gene.location$strand,
                               symbol = gene.location$external_gene_name,
                               EntrezID = gene.location$entrezgene)

            rowRanges <- GRanges(seqnames = paste0("chr", df$Chromosome),
                                 ranges = IRanges(start = df$Genomic_Coordinate,
                                                  end = df$Genomic_Coordinate),
                                 probeID = df$Composite.Element.REF)

            # closest gene to each 450k probe ##your data
            neargene <- as.data.frame(nearest(rowRanges, gene.GR))
            neargene <- gene.location[neargene[,1],]
            rowRanges$Gene_Symbol <-  neargene$external_gene_name

            names(rowRanges) <- as.character(df$Composite.Element.REF)
            colData <-  colDataPrepare(colnames(df)[5:ncol(df)])
            assay <- data.matrix(subset(df,select = c(5:ncol(df))))

            sset <- SummarizedExperiment(assays = assay,
                                         rowRanges = rowRanges,
                                         colData = colData)

        } else {
            setDF(df)
            rownames(df) <- df$Composite.Element.REF
            df$Composite.Element.REF <- NULL
            df[,3:ncol(df)] <- sapply(df[,3:ncol(df)], as.numeric)
        }
    }

    if (grepl("mda_rppa_core",tolower(platform))) {
        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE)
            sample <- gsub("\\.", "-", colnames(data)[2])
            colnames(data) <- data[1,]
            data <- data[-1,] # removing Composite Element REF
            colnames(data)[2] <- sample

            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data,by = "Composite Element REF")
            }
            setTxtProgressBar(pb, i)
        }
        rownames(df) <- df[,1]
        df[,1] <- NULL
        # get array_design.txt from mage folder
        # and change uuid by Barcode
        uuid <- colnames(df)
        idx <- grep("Sample|Control",uuid)
        if(length(idx) > 0){
            uuid <- uuid[-idx]
        }
        map <- mapuuidbarcode(uuid)
        idx <- which(colnames(df) %in% map$uuid)
        colnames(df)[idx] <- as.character(map$barcode)
    }


    if (grepl("illuminadnamethylation_oma",
              platform, ignore.case = TRUE)) {

        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE, check.names = FALSE,
                               comment.char = "#",fill = TRUE)
            data <- data[-1,] # removing Composite Element REF
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data,by = "Hybridization REF")
            }
            setTxtProgressBar(pb, i)
        }

        rownames(df) <- df[,1]
        df[,1] <- NULL
    }

    if (tolower(platform) == "illuminaga_rnaseq" ||
        tolower(platform) == "illuminahiseq_rnaseq") {
        # Barcode in the name
        regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                        "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
        barcode <- str_match(files,regex)

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)

            setnames(data,colnames(data)[2:ncol(data)],
                     paste0(colnames(data)[2:ncol(data)],"_",barcode[i]))
            # removing duplicated rows
            data <- subset(data,subset=(!duplicated(data[,1,with=FALSE])))
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data, colnames(data)[1])
            }
            setTxtProgressBar(pb, i)
        }
        if(summarizedExperiment){
            if(grepl("gene",colnames(df)[1])){

                aux <- strsplit(df$gene,"\\|")
                GeneID <- unlist(lapply(aux,function(x) x[2]))
                df$entrezgene <- as.numeric(GeneID)
                merged <- merge(df,gene.location,by="entrezgene")
                rowRanges <- GRanges(seqnames = paste0("chr", merged$chromosome_name),
                                     ranges = IRanges(start = merged$start_position,
                                                      end = merged$end_position),
                                     strand=merged$strand,
                                     gene_id = merged$external_gene_name,
                                     entrezgene = merged$entrezgene)
                names(rowRanges) <- as.character(merged$gene)
                assays <- SimpleList(
                    raw_counts=data.matrix(subset(merged,select=seq(3,ncol(df),3))),
                    median_length_normalized=data.matrix(subset(merged,select=seq(4,ncol(df),3))),
                    RPKM=data.matrix(subset(merged,select=seq(5,ncol(df),3))))

            }  else if(grepl("junction",colnames(df)[1])){

                # junction example: chr1:12227:+,chr1:12595:+
                aux    <- strsplit(df$junction,":")
                name   <- unlist(lapply(aux,function(x) x[1]))
                x <- as.numeric(unlist(lapply(aux,function(x) x[2])))
                y    <- as.numeric(unlist(lapply(aux,function(x) x[4])))
                start <- apply(data.frame(x,y),1,min)
                end <- apply(data.frame(x,y),1,max)

                strand <- unlist(lapply(aux,function(x) x[5]))

                rowRanges <- GRanges(seqnames = name,
                                     ranges = IRanges(start = start, end = end),
                                     strand = strand)
                names(rowRanges) <- as.character(df$junction)
                assays <- SimpleList(raw_counts=data.matrix(subset(df,select=2:ncol(df))))

            } else if(grepl("exon",colnames(df)[1])){
                # exon chr1:11874-12227:+
                aux       <- strsplit(df$exon,":")
                name      <- unlist(lapply(aux,function(x) x[1]))
                start.end <- strsplit(unlist(lapply(aux,function(x) x[2])),"-")
                x       <- as.numeric(unlist(lapply(start.end,function(x) x[1])))
                y     <- as.numeric(unlist(lapply(start.end,function(x) x[2])))
                start <- apply(data.frame(x,y),1,min)
                end <- apply(data.frame(x,y),1,max)
                strand    <- unlist(lapply(aux,function(x) x[3]))

                rowRanges <- GRanges(seqnames = name,
                                     ranges = IRanges(start = start, end = end),
                                     strand = strand)
                names(rowRanges) <- as.character(df$exon)
                assays <- SimpleList(
                    raw_counts=data.matrix(subset(df,select=seq(2,ncol(df),3))),
                    median_length_normalized=data.matrix(subset(df,select=seq(3,ncol(df),3))),
                    RPKM=data.matrix(subset(df,select=seq(4,ncol(df),3))))

            }
            colData <- colDataPrepare(as.character(barcode))
            sset <- SummarizedExperiment(assays=assays,
                                         rowRanges=rowRanges,
                                         colData=colData)
        }else {
            setDF(df)
            rownames(df) <- df[,1]
            df[,1] <- NULL
        }
    }

    if (tolower(platform) == tolower("HT_HG-U133A")) {
        # Barcode in the mage file
        if(is.vector(query)){
            mage <- getMage(query)
        } else {
            mage <- getMage(query[1,])
        }

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t", skip= 1,
                          stringsAsFactors = FALSE)

            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data, colnames(data)[1])
            }
            setTxtProgressBar(pb, i)
        }
        names <- sapply(files,
                        function(x) {
                            y <- fread(x, header = FALSE,
                                       stringsAsFactors = FALSE,
                                       nrows=1)$V2
                            idx <- grep(y,mage$Hybridization.Name)
                            mage[idx,]$Comment..TCGA.Barcode.
                        })
        setnames(df,2:ncol(df),names)

        if(summarizedExperiment){
            # TODO create GRanges
            df$external_gene_name <-  alias2SymbolTable(df$`Composite Element REF`)
            merged <- merge(df,gene.location,by="external_gene_name")
            rowRanges <- GRanges(seqnames = paste0("chr", merged$chromosome_name),
                                 ranges = IRanges(start = merged$start_position,
                                                  end = merged$end_position),
                                 strand=merged$strand,
                                 gene_id = merged$external_gene_name,
                                 entrezgene = merged$entrezgene,
                                 alias = merged$`Composite Element REF`)
            names(rowRanges) <- as.character(merged$`Composite Element REF`)

            regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                            "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
            barcode <- unique(unlist(str_match_all(colnames(merged),regex)))
            colData <- colDataPrepare(barcode)

            assays <- SimpleList(raw_counts=data.matrix(
                subset(merged,select=seq(3,2+length(barcode)))
            )
            )

            sset <- SummarizedExperiment(assays=assays,
                                         rowRanges=rowRanges,
                                         colData=colData)
        }
    }

    if (tolower(platform) == tolower("HG-U133_Plus_2") ||
        grepl("H-miRNA_8x15K|agilent",platform, ignore.case = TRUE)) {

        for (i in seq_along(files)) {
            header <- fread(files[i], sep = "\t",
                            stringsAsFactors = FALSE, nrows =  1, header = FALSE,
                            colClasses=c("character",
                                         "character"))
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, skip = 1,
                          colClasses=c("character",
                                       "numeric"))
            setnames(data,1:2,as.character(header))
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data, colnames(data)[1])
            }
            #setTxtProgressBar(pb, i)
        }
        df <- df[-1,]

        if(summarizedExperiment){
            if(!grepl("H-miRNA_8x15K|agilent",platform, ignore.case = TRUE)){
                # TODO create GRanges
                df$external_gene_name <-  alias2SymbolTable(df$`Hybridization REF`)
                merged <- merge(df,gene.location,by="external_gene_name")
                rowRanges <- GRanges(seqnames = paste0("chr", merged$chromosome_name),
                                     ranges = IRanges(start = merged$start_position,
                                                      end = merged$end_position),
                                     strand=merged$strand,
                                     gene_id = merged$external_gene_name,
                                     entrezgene = merged$entrezgene,
                                     alias = merged$`Hybridization REF`)
                names(rowRanges) <- as.character(merged$`Hybridization REF`)
            } else {
                microRNA <- as.data.frame(microRNAs(TxDb.Hsapiens.UCSC.hg19.knownGene))
                df$mirna_id <- tolower(df$`Hybridization REF`)
                merged <- merge(df,microRNA, by="mirna_id")
                rowRanges <- GRanges(seqnames = merged$seqnames,
                                     ranges = IRanges(start = merged$start,
                                                      end = merged$end),
                                     strand=merged$strand,
                                     mirna_id = merged$mirna_id)
                names(rowRanges) <- as.character(merged$mirna_id)
            }
            regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                            "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
            barcode <- unique(unlist(str_match_all(colnames(merged),regex)))
            colData <- colDataPrepare(barcode)

            assays <- SimpleList(raw_counts=data.matrix(
                subset(merged,select=seq(3,2+length(barcode)))
            )
            )

            sset <- SummarizedExperiment(assays=assays,
                                         rowRanges=rowRanges,
                                         colData=colData)
        }

    }
    if (grepl("rnaseqv2",platform, ignore.case = TRUE)) {

        if(is.null(type) || (type != "rsem.genes.results" &&
                             type != "rsem.isoforms.results" &&
                             type != "rsem.genes.normalized_results" &&
                             type != "rsem.isoforms.normalized_results" &&
                             type != "bt.exon_quantification" &&
                             type != "junction_quantification")
        ){
            msg <- paste0("Plase select a type. \n Possibilities:\n",
                          " = rsem.genes.results\n",
                          " = rsem.isoforms.results\n",
                          " = rsem.genes.normalized_results\n",
                          " = rsem.isoforms.normalized_results\n",
                          " = bt.exon_quantification\n",
                          " = junction_quantification"
            )
            message(msg)
            return()
        }

        if(type == "rsem.genes.results")               pat <- "rsem.genes.results"
        if(type == "rsem.isoforms.results")            pat <- "rsem.isoforms.results"
        if(type == "rsem.genes.normalized_results")    pat <- "rsem.genes.normalized_results"
        if(type == "rsem.isoforms.normalized_results") pat <- "rsem.isoforms.normalized_results"
        if(type == "junction_quantification")          pat <- "junction_quantification"
        if(type == "bt.exon_quantification")           pat <- "bt.exon_quantification"



        files <- files[grep(pat,files, perl = TRUE)]

        regex <- paste0("[:alnum:]{8}-[:alnum:]{4}",
                        "-[:alnum:]{4}-[:alnum:]{4}-[:alnum:]{12}")
        uuid <- str_match(files,regex)
        map <- mapuuidbarcode(uuid)

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
            x <- subset(map, uuid == uuid[i])

            if(summarizedExperiment){
                setnames(data,colnames(data)[2:ncol(data)],
                         paste0(colnames(data)[2:ncol(data)],"_",x$barcode))
            } else {
                setnames(data,2, as.character(x$barcode))
            }
            # removing duplicated rows
            data <- subset(data,subset=(!duplicated(data[,1,with=FALSE])))
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data,by = colnames(df)[1])
            }
            setTxtProgressBar(pb, i)
        }

        if(summarizedExperiment){
            if(grepl("gene_id",colnames(df)[1])){
                aux <- strsplit(df$gene_id,"\\|")
                GeneID<-unlist(lapply(aux,function(x) x[2]))
                df$entrezgene <- as.numeric(GeneID)
                merged <- merge(df,gene.location,by="entrezgene")


                rowRanges <- GRanges(seqnames = paste0("chr", merged$chromosome_name),
                                     ranges = IRanges(start = merged$start_position,
                                                      end = merged$end_position),
                                     strand=merged$strand,
                                     gene_id = merged$external_gene_name,
                                     entrezgene = merged$entrezgene )
                names(rowRanges) <- as.character(merged$gene_id)

                if(length(colnames(data))>2){
                    assays <- SimpleList(
                        raw_counts=data.matrix(subset(merged,select=seq(3,ncol(df),3))),
                        scaled_estimate=data.matrix(subset(merged,select=seq(4,ncol(df),3))),
                        transcript_id=data.matrix(subset(merged,select=seq(5,ncol(df),3))))
                } else {
                    assays <- SimpleList(
                        raw_counts=data.matrix(subset(merged,select=seq(3,ncol(df)))))
                }
            } else if(grepl("junction",colnames(df)[1])){
                aux    <- strsplit(df$junction,":")
                name   <- unlist(lapply(aux,function(x) x[1]))
                x <- as.numeric(unlist(lapply(aux,function(x) x[2])))
                y    <- as.numeric(unlist(lapply(aux,function(x) x[4])))
                start <- apply(data.frame(x,y),1,min)
                end <- apply(data.frame(x,y),1,max)

                strand <- unlist(lapply(aux,function(x) x[5]))

                rowRanges <- GRanges(seqnames = name,
                                     ranges = IRanges(start = start, end = end),
                                     strand = strand)
                names(rowRanges) <- as.character(df$junction)
                assays <- SimpleList(raw_counts=data.matrix(subset(df,select=2:ncol(df))))
            } else if(grepl("exon",colnames(df)[1])){
                # exon chr1:11874-12227:+
                aux       <- strsplit(df$exon,":")
                name      <- unlist(lapply(aux,function(x) x[1]))
                start.end <- strsplit(unlist(lapply(aux,function(x) x[2])),"-")
                x       <- as.numeric(unlist(lapply(start.end,function(x) x[1])))
                y     <- as.numeric(unlist(lapply(start.end,function(x) x[2])))
                start <- apply(data.frame(x,y),1,min)
                end <- apply(data.frame(x,y),1,max)
                strand    <- unlist(lapply(aux,function(x) x[3]))

                rowRanges <- GRanges(seqnames = name,
                                     ranges = IRanges(start = start, end = end),
                                     strand = strand)
                names(rowRanges) <- as.character(df$exon)
                assays <- SimpleList(
                    raw_counts=data.matrix(subset(df,select=seq(2,ncol(df),3))),
                    median_length_normalized=data.matrix(subset(df,select=seq(3,ncol(df),3))),
                    RPKM=data.matrix(subset(df,select=seq(4,ncol(df),3))))
            } else if(grepl("isoform",colnames(df)[1])){
                message("TBD")
            }
            regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                            "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
            barcode <- unique(unlist(str_match_all(colnames(merged),regex)))
            colData <- colDataPrepare(barcode)

            sset <- SummarizedExperiment(assays=assays,
                                         rowRanges=rowRanges,
                                         colData=colData)
        } else {
            setDF(df)
            rownames(df) <- df[,1]
            df[,1] <- NULL
        }
    }

    if (grepl("illuminahiseq_mirnaseq",platform, ignore.case = TRUE)) {

        if(is.null(type) || (type != "hg19.mirna" && type != "mirna")){
            msg <- paste0("Plase select a type. \n Possibilities:\n",
                          " = hg19.mirna\n = mirna")
            message(msg)
            return()
        }

        if(type == "hg19.mirna")   pat <- "(hg19.)mirna"
        if(type == "mirna")        pat <- "(?<!hg19\\.)mirna"

        files <- files[grep(pat,files, perl = TRUE)]

        if(length(files) == 0){
            message("No mirna files of that type found")
            return(NULL)
        }

        regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                        "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
        barcode <- str_match(files,regex)

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
            data <- subset(data,select=c(1:3))
            setnames(data,2:ncol(data),
                     paste0(as.character(barcode[i]),".",colnames(data)[2:ncol(data)]))
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data, by=colnames(data)[1])
            }
        }
        setDF(df)
        rownames(df) <- df[,1]
        df <- df[,-1]
    }

    if (grepl("bio",platform,ignore.case = TRUE)) {
        if (!is.null(type)) {
            files <- files[grep(type,files)]
        }
        if (length(files) == 1) {
            df <- read.table(files, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE,
                             comment.char = "#",fill = TRUE)
            regex <- paste0("[[:alnum:]]{8}-[[:alnum:]]{4}",
                            "-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}")
            idx <- grepl(regex,df$bcr_patient_uuid)
            df <- df[idx,]
        } else {
            message("We're preaparing for the moment only one clinical file")
            return(NULL)
        }
    }
    if (grepl("genome_wide_snp_6",tolower(platform))) {

        close(pb)
        message("Preparing h19 files...")
        idx <- grep("nocnv|hg18", files)
        if(length(idx)>0){
            files <- files[-idx]
        }

        pb <- txtProgressBar(min = 0, max = length(files), style = 3)

        if(is.vector(query)){
            mage <- getMage(query)
        } else {
            mage <- getMage(query[1,])
        }

        genes <- sort(unique(gene.location$external_gene_name))

        df <- matrix(0, nrow = length(genes), ncol = length(files))
        rownames(df) <- genes

        colNames <- rep("", ncol(df)) #check barcode
        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
            ####Check barcodes
            colNames[i] <- data$Sample[1]

            for(j in 1:nrow(data)){
                gg <- sort(unique(gene.location[
                    gene.location$start_position >= data$Start[j]
                    & gene.location$end_position <= data$End[j],
                    "external_gene_name"]))
                df[gg, i] <- df[gg, i] + data$Segment_Mean[j]
            }
            setTxtProgressBar(pb, i)
        }
        id <- data.frame(id=colNames)
        names <- merge(id,mage,by.x="id",by.y="Hybridization.Name")
        colnames(df) <- names$Comment..TCGA.Barcode.

        if (summarizedExperiment){
            message("Preparing summarizedExperiment")
            df <- DataFrame(df)
            df$external_gene_name <- rownames(df)
            merged <- merge(df,gene.location,by="external_gene_name")
            rowRanges <- GRanges(seqnames = paste0("chr", merged$chromosome_name),
                                 ranges = IRanges(start = merged$start_position,
                                                  end = merged$end_position),
                                 strand=merged$strand,
                                 gene_id = merged$external_gene_name,
                                 entrezgene = merged$entrezgene)
            names(rowRanges) <- as.character(merged$external_gene_name)

            idx  <- grep("TCGA",colnames(merged))
            colData <- colDataPrepare(colnames(merged)[idx])
            assays <- SimpleList(sum_Segment_Mean=data.matrix(
                subset(merged,select = idx))
            )

            sset <- SummarizedExperiment(assays = assays,
                                         rowRanges = rowRanges,
                                         colData = colData)

        }
    }
    close(pb)

    if (save) {
        if (is.null(filename)) {
            filename <- paste0(platform,"_",gsub(" ","_",Sys.time()),".rda")
        }
        if (!is.null(sset)) {
            save(sset,file = filename)
        } else {
            save(df,file = filename)
        }
    }

    if (!is.null(toPackage) && !summarizedExperiment) {
        df <- prepareToPackage(df, platform,toPackage)
    }
    if (!is.null(sset)) {
        return(sset)
    }
    return(df)
}

# This function will help the user to prepare the data to an specific package
prepareToPackage <- function(df, platform, toPackage){

    if(grepl("elmer", toPackage, ignore.case = TRUE)){

        if (grepl("illuminahiseq_rnaseqv2|illuminahiseq_totalrnaseqv2",
                  platform, ignore.case = TRUE)){
            message("============ Pre-pocessing expression data =============")
            message(paste0("1 - expression = log2(expression + 1): ",
                           "To linearize \n    relation between ",
                           "methylation and expression"))
            df <- log2(df+1)
            message("2 - rownames  (gene|loci) => ('ID'loci) ")
            aux <- strsplit(rownames(df),"\\|")
            GeneID <- unlist(lapply(aux,function(x) x[2]))
            row.names(df) <- paste0("ID",GeneID)
            df <- data.matrix(df)
        }

        if (grepl("humanmethylation", platform, ignore.case = TRUE)) {
            message("============ Pre-pocessing methylation data =============")
            msg <- paste0("1 - Removing Columns: \n  * Gene_Symbol  \n",
                          "  * Chromosome  \n  * Genomic_Coordinate")
            message(msg)
            df <- subset(df,select = 4:ncol(df))
            msg <- paste0("2 - Removing probes with ",
                          "NA values in more than 0.80% samples")
            message(msg)
            df <- df[rowMeans(is.na(df))<0.2,]
        }
    }

    if(grepl("limma", toPackage, ignore.case = TRUE)){
        message("TBD")
    }

    if(grepl("methylmix", toPackage, ignore.case = TRUE)){
        message("TBD")
    }

    if(grepl("biomics", toPackage, ignore.case = TRUE)){
        message("TBD")
    }
    return(df)
}

#  Internal function
#  Get a list of barcode from a list of uuid
#  example mapuuidbarcode(c("011bb13f-e0e8-4f4b-b7a5-4867bbe3b30a",
#                           "048615c7-c08c-4199-b394-c59160337d67"))
#' @importFrom rjson fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom RCurl postForm
mapuuidbarcode <- function(uuid){
    # Using tcga api: https://goo.gl/M1uQLR
    serv <- "https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/uuid/batch"
    header <- c('Content-Type' = 'text/plain')
    uuids <- paste0(uuid, collapse = ",")
    ans <- fromJSON(postForm(serv,
                             .opts = list(postfields = uuids,
                                          httpheader = header,
                                          ssl.verifypeer = FALSE)))
    if(length(uuid) == 1){
        x <- data.frame(ans$uuidMapping$barcode,ans$uuidMapping$uuid,
                        stringsAsFactors = FALSE)
        colnames(x) <- c("barcode","uuid")
    } else {
        # Extract patient barcode from sample barcode.
        x <- (do.call("rbind.fill", lapply(ans$uuidMapping, as.data.frame)))
    }
    return(x)
}

#  Internal function
#  Get a list of barcode from a list of uuid
#  example mapuuidbarcode(c("011bb13f-e0e8-4f4b-b7a5-4867bbe3b30a",
#                           "048615c7-c08c-4199-b394-c59160337d67"))
#' @importFrom rjson fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom RCurl postForm
mapbarcodeuuid <- function(barcode){
    # Using tcga api: https://goo.gl/M1uQLR
    barcodes <- paste0(barcode, collapse = ",")
    serv <- paste0("https://tcga-data.nci.nih.gov/",
                   "uuid/uuidws/mapping/json/barcode/batch")
    header <- c('Content-Type' = 'text/plain')
    ans <- fromJSON(postForm(serv,
                             .opts = list(postfields = barcodes,
                                          httpheader = header,
                                          ssl.verifypeer = FALSE)))

    if(length(barcode) == 1){
        x <- data.frame(ans$uuidMapping$barcode,ans$uuidMapping$uuid,
                        stringsAsFactors = FALSE)
        colnames(x) <- c("barcode","uuid")
    } else {
        # Extract patient barcode from sample barcode.
        x <- (do.call("rbind.fill", lapply(ans$uuidMapping, as.data.frame)))
    }
    return(x)
}

# Get sdrf file/array_design of a line
# example
# query <- TCGAQuery(tumor = "BRCA")
# getMagecontent(query[1,])
# Obs: delete the file after reading
#      is it better to save it?
getMage <- function(line){

    root <- "https://tcga-data.nci.nih.gov"
    path <- "mages"
    dir.create(path,showWarnings = FALSE)
    mages <-  tcga.db[grep("mage-tab",tcga.db$name),]
    # get the mage file for the entry
    mage <- subset(mages, mages$Disease == line$Disease &
                       mages$Platform == line$Platform &
                       mages$Center == line$Center)

    if (dim(mage)[1] != 0) {
        file <- file.path(path,basename(mage$deployLocation))
        if ( !file.exists(file)) {
            if(!is.windows()){
                download(paste0(root,mage$deployLocation), file,
                         quiet = TRUE)
            } else {
                download(paste0(root,mage$deployLocation), file,
                         quiet = TRUE,method="auto")
            }
        }
        folder <- gsub(".tar.gz","",file)
        if ( !file.exists(folder)) {
            untar(file,exdir = "mages")
        }
        files <- list.files(folder)
        if (line$Platform == "MDA_RPPA_Core") {
            sdrf <- files[grep("array_design",files)]
        } else {
            # Platform is not MDA_RPPA_Core
            sdrf <- files[grep("sdrf",files)]
        }
        if (length(sdrf) > 1) {
            sdrf <- sdrf[1]
        }

        df <- read.delim(file = file.path(folder,sdrf),
                         sep = "\t",
                         stringsAsFactors = FALSE,
                         fileEncoding = "latin1")
    }
    return(df)
}

