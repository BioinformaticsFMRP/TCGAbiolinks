
#' @title Prepare GDC data
#' @description
#'   Reads the data downloaded and prepare it into an R object
#' @param query A query for GDCquery function
#' @param save Save result as RData object?
#' @param save.filename Name of the file to be save if empty an automatic will be created
#' @export
#' @return A summarizedExperiment or a data.frame
GDCPrepare <- function(query, save = FALSE, save.filename, summarizedExperiment = TRUE){

    if(missing(query)) stop("Please set query parameter")

    # We save the files in project/data.category/data.type/file_id/file_name
    files <- file.path(query$project,
                       gsub(" ","_",query$results[[1]]$data_category),
                       gsub(" ","_",query$results[[1]]$data_type),
                       gsub(" ","_",query$results[[1]]$file_id),
                       gsub(" ","_",query$results[[1]]$file_name))

    if(query$data.category == "Transcriptome Profiling"){
        data <- readTranscriptomeProfiling(files,unique(query$results[[1]]$analysis$workflow_type),query$results[[1]]$cases)
    } else if(query$data.category == "Copy Number Variation") {
        data <- readCopyNumberVariantion(files, query$results[[1]]$cases)
    }  else if(query$data.category == "DNA methylation") {
        data <- readDNAmethylation(files, query$results[[1]]$cases, summarizedExperiment, unique(query$platform))
    }  else if(query$data.category == "Protein expression") {
        data <- readProteinExpression(files, query$results[[1]]$cases)
    } else if (query$data.category == "Gene expression") {
        if(query$data.type == "Gene expression quantification")
            data <- readGeneExpressionQuantification(files,
                                                     query$results[[1]]$cases,
                                                     summarizedExperiment,
                                                     unique(query$platform))
        #if(query$data.type == "Gene expression quantification") data <- readGeneExpression()

    }

    if(save){
        if(missing(save.filename) & !missing(query)) save.filename <- paste0(query$project,gsub(" ","_", query$data.category),gsub(" ","_",date()),".RData")
        message(paste0("Saving file:",save.filename))
        save(data, file = save.filename)
        message("File saved")
    }
    return(data)
}

readGeneExpressionQuantification <- function(files, cases, summarizedExperiment = TRUE, platform){
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    skip = 0
    if(grepl("Agilent", platform)) {
        skip = 1
        summarizedExperiment = FALSE
    }
    for (i in seq_along(files)) {
        data <- fread(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE,skip = skip)

        if(!missing(cases)) {
            assay.list <- colnames(data)[2:ncol(data)]
            setnames(data,colnames(data)[2:ncol(data)],
                     paste0(colnames(data)[2:ncol(data)],"_",cases[i]))
        }
        if (i == 1) {
            df <- data
        } else {
            df <- merge(df, data, by=colnames(data)[1], all = TRUE)
        }
        setTxtProgressBar(pb, i)
    }
    setDF(df)

    if (summarizedExperiment) df <- makeSEfromGeneExpressionQuantification(df,assay.list)
    return(df)
}
makeSEfromGeneExpressionQuantification <- function(df, assay.list, genome="hg19"){
    gene.location <- get.GRCh.bioMart(genome)
    if(all(grepl("\\|",df[,1]))){
        aux <- strsplit(df$gene_id,"\\|")
        GeneID <- unlist(lapply(aux,function(x) x[2]))
        df$entrezid <- as.numeric(GeneID)
        GeneSymbol <- unlist(lapply(aux,function(x) x[1]))
        df$external_gene_id <- as.character(GeneSymbol)
    }
    df <- merge(df, gene.location, by="external_gene_id")
    print(assay.list)
    print(head(df))
    if("transcript_id" %in% assay.list){
        rowRanges <- GRanges(seqnames = paste0("chr", df$chromosome_name),
                             ranges = IRanges(start = df$start_position,
                                              end = df$end_position),
                             strand = df$strand,
                             gene_id = df$external_gene_id,
                             entrezgene = df$entrezid,
                             transcript_id = subset(df, select = 5))
        names(rowRanges) <- as.character(df$gene_id)
        assay.list <- assay.list[which(assay.list != "transcript_id")]
    } else {
        rowRanges <- GRanges(seqnames = paste0("chr", df$chromosome_name),
                             ranges = IRanges(start = df$start_position,
                                              end = df$end_position),
                             strand = df$strand,
                             gene_id = df$external_gene_id,
                             entrezgene = df$entrezid)
        names(rowRanges) <- as.character(df$gene_id)
    }
    assays <- lapply(assay.list, function (x) {
        return(data.matrix(subset(df, select = grep(x,colnames(df)))))
    })
    names(assays) <- assay.list
    regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                    "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
    samples <- na.omit(unique(str_match(colnames(df),regex)[,1]))
    colData <-  colDataPrepare(samples)

    assays <- lapply(assays, function(x){
        colnames(x) <- NULL
        rownames(x) <- NULL
        return(x)
    })
    save(assays,rowRanges,colData,file = "test2.rda")
    rse <- SummarizedExperiment(assays=assays,
                                rowRanges=rowRanges,
                                colData=colData)
    return(rse)
}

makeSEfromDNAmethylation <- function(df, genome="hg19"){
    gene.location <- get.GRCh.bioMart(genome)
    gene.GR <- GRanges(seqnames = paste0("chr", gene.location$chromosome_name),
                       ranges = IRanges(start = gene.location$start_position,
                                        end = gene.location$end_position),
                       strand = gene.location$strand,
                       symbol = gene.location$external_gene_id,
                       EntrezID = gene.location$entrezgene)

    rowRanges <- GRanges(seqnames = paste0("chr", df$Chromosome),
                         ranges = IRanges(start = df$Genomic_Coordinate,
                                          end = df$Genomic_Coordinate),
                         probeID = df$Composite.Element.REF,
                         Gene_Symbol = df$Gene_Symbol)

    names(rowRanges) <- as.character(df$Composite.Element.REF)
    colData <-  colDataPrepare(colnames(df)[5:ncol(df)])
    assay <- data.matrix(subset(df,select = c(5:ncol(df))))
    colnames(assay) <- rownames(colData)
    rownames(assay) <- as.character(df$Composite.Element.REF)

    rse <- SummarizedExperiment(assays = assay, rowRanges = rowRanges, colData = colData)

}

# We will try to make this function easier to use this function in its own data
# In case it is not TCGA I should not consider that there is a barcode in the header
# Instead the user should be able to add the names to his data
# The only problem is that the data from the user will not have all the columns
# TODO: Improve this function to be more generic as possible
readDNAmethylation <- function(files, cases, summarizedExperiment = TRUE, platform){
    if (grepl("OMA00",platform)){
        pb <- txtProgressBar(min = 0, max = length(files), style = 3)
        for (i in seq_along(files)) {
            print(files[i])
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE,skip = 1,
                          na.strings="N/A",
                          colClasses=c("character", # Composite Element REF
                                       "numeric"))   # beta value
            setnames(data,gsub(" ", "\\.", colnames(data)))
            if(!missing(cases)) setnames(data,2,cases[i])
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data, by = "Composite.Element.REF")
            }
            setTxtProgressBar(pb, i)
        }
        setDF(df)
        rownames(df) <- df$Composite.Element.REF
        df$Composite.Element.REF <- NULL
    } else {
        pb <- txtProgressBar(min = 0, max = length(files), style = 3)
        for (i in seq_along(files)) {
            print(files[i])
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE,skip = 1,
                          colClasses=c("character", # Composite Element REF
                                       "numeric",   # beta value
                                       "character", # Gene symbol
                                       "character", # Chromosome
                                       "integer"))  # Genomic coordinate
            setnames(data,gsub(" ", "\\.", colnames(data)))
            if(!missing(cases)) setnames(data,2,cases[i])
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
            df <- makeSEfromDNAmethylation(df)
        } else {
            setDF(df)
            rownames(df) <- df$Composite.Element.REF
            df$Composite.Element.REF <- NULL
        }
    }
    return(df)
}

colDataPrepare <- function(barcode){

    code <- c('01','02','03','04','05','06','07','08','09','10','11',
              '12','13','14','20','40','50','60','61')
    shortLetterCode <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC",
                         "TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB",
                         "CELL","XP","XCL")

    definition <- c("Primary solid Tumor",
                    "Recurrent Solid Tumor",
                    "Primary Blood Derived Cancer - Peripheral Blood",
                    "Recurrent Blood Derived Cancer - Bone Marrow",
                    "Additional - New Primary",
                    "Metastatic",
                    "Additional Metastatic",
                    "Human Tumor Original Cells",
                    "Primary Blood Derived Cancer - Bone Marrow",
                    "Blood Derived Normal",
                    "Solid Tissue Normal",
                    "Buccal Cell Normal",
                    "EBV Immortalized Normal",
                    "Bone Marrow Normal",
                    "Control Analyte",
                    "Recurrent Blood Derived Cancer - Peripheral Blood",
                    "Cell Lines",
                    "Primary Xenograft Tissue",
                    "Cell Line Derived Xenograft Tissue")
    aux <- DataFrame(code = code,shortLetterCode,definition)

    # in case multiple equal barcode
    regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                    "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
    samples <- str_match(barcode,regex)[,1]

    ret <- DataFrame(barcode = barcode,
                     patient = substr(barcode, 1, 12),
                     sample = substr(barcode, 1, 16),
                     code = substr(barcode, 14, 15))
    ret <- merge(ret,aux, by = "code", sort = FALSE)
    ret <- ret[match(barcode,ret$barcode),]
    rownames(ret) <- gsub("\\.","-",make.names(ret$barcode,unique=TRUE))
    ret$code <- NULL
    return(DataFrame(ret))
}

#' @importFrom biomaRt getBM useMart
get.GRCh.bioMart <- function(genome="hg19") {
    message(paste0("Downloading genome information: ",genome))
    if (genome == "hg19"){
        # for hg19
        ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           host = "feb2014.archive.ensembl.org",
                           path = "/biomart/martservice" ,
                           dataset = "hsapiens_gene_ensembl")
    } else {
        # for hg38
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    }

    chrom <- c(1:22, "X", "Y")
    gene.location <- getBM(attributes = c("chromosome_name",
                                          "start_position",
                                          "end_position", "strand",
                                          "external_gene_id",
                                          "entrezgene"),
                           filters = c("chromosome_name"),
                           values = list(chrom), mart = ensembl)

    return(gene.location)
}


readProteinExpression <- function(files,cases) {
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    for (i in seq_along(files)) {
        data <- read_tsv(file = files[i], col_names = TRUE)
        data <- data[-1,]
        if(i == 1) df <- data
        if(i != 1) df <- merge(df, data,all=TRUE, by="Sample REF")
        setTxtProgressBar(pb, i)
    }
    close(pb)
    if(!missing(cases))  colnames(df)[2:length(colnames(df))] <- cases

    return(df)
}

readTranscriptomeProfiling <- function(files, workflow.type, cases) {
    # Status working for:
    #  - htseq
    #  - FPKM
    #  - FPKM-UQ
    if(grepl("HTSeq",workflow.type)){
        pb <- txtProgressBar(min = 0, max = length(files), style = 3)
        for (i in seq_along(files)) {
            data <- read_tsv(file = files[i], col_names = FALSE)
            if(!missing(cases))  colnames(data)[2] <- cases[i]
            if(i == 1) df <- data
            if(i != 1) df <- merge(df, data, by=colnames(df)[1],all = TRUE)
            setTxtProgressBar(pb, i)
        }
        close(pb)
    }
    return(df)
}

# Reads Copy Number Variantion files to a data frame, basically it will rbind it
readCopyNumberVariantion <- function(files, cases){

    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    for (i in seq_along(files)) {
        data <- read_tsv(file = files[i])
        aux <- query$results[[1]]
        if(!missing(cases)) data$Barcode <- cases[i]
        if(i == 1) df <- data
        if(i != 1) df <- rbind(df, data, make.row.names = FALSE)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    return(df)
}

# Source: https://stackoverflow.com/questions/10266963/moving-files-between-folders
move <- function(from, to) {
    todir <- dirname(to)
    if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE,showWarnings = FALSE)
    file.rename(from = from,  to = to)
}

#' @title Read the data from level 3 the experiments and prepare it
#'  for downstream analysis into a SummarizedExperiment object.
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
#' @param query TCGAquery output
#' @param dir Directory with the files downloaded by TCGAdownload
#' @param samples List of samples to prepare the data
#' @param type Filter the files to prepare.
#' @param save Save a rda object with the prepared object?
#'  Default: \code{FALSE}
#' @param filename Name of the saved file
#' @param add.mutation.genes Integrate information about genes mutation? DEFAULT: FALSE
#' @param reannotate Reannotate genes?  Source http://grch37.ensembl.org/.
#' DEFAULT: FALSE. (For the moment only working for methylation data)
#' @param summarizedExperiment Output as SummarizedExperiment?
#' Default: \code{FALSE}
#' @param add.subtype Add subtype information from tcgaquery_subtype? Default: \code{FALSE}
#' @param add.clinical Add clinical information from TCGAquery_clinic?
#' (The information file add will be: clinical_patient) Default: \code{FALSE}
#' @examples
#' \dontrun{
#' sample <- "TCGA-06-0939-01A-01D-1228-05"
#' query <- TCGAquery(tumor = "GBM",samples = sample, level = 3)
#' TCGAdownload(query,path = "exampleData",samples = sample)
#' data <- TCGAprepare(query, dir="exampleData")
#' }
#' @export
#' @importFrom stringr str_match str_trim str_detect str_match_all
#' @importFrom SummarizedExperiment SummarizedExperiment metadata<- colData<-
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom limma alias2SymbolTable
#' @importFrom GenomicFeatures microRNAs
#' @importFrom BiocGenerics as.data.frame
#' @importFrom GenomicRanges GRanges distanceToNearest
#' @importFrom data.table data.table
#' @importFrom IRanges IRanges
#' @import utils TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom data.table fread setnames setcolorder setDF data.table
TCGAprepare <- function(query,
                        dir = NULL,
                        samples = NULL,
                        type = NULL,
                        save = FALSE,
                        filename = NULL,
                        add.mutation.genes = FALSE,
                        add.clinical = TRUE,
                        reannotate = FALSE,
                        summarizedExperiment = TRUE,
                        add.subtype = FALSE){
    stop("TCGA data has moved from DCC server to GDC server. Please use GDCprepare function")
    if(add.subtype){
        for (i in unique(query$Disease)) {
            if (!grepl("lgg|gbm|luad|stad|brca|coad|read|skcm|hnsc|kich|lusc|ucec|pancan|thca|prad|kirp|kirc|all",
                       i,ignore.case = TRUE)) {
                message(paste0("Sorry we don't have subtypes for: ",i))
                message("add.subtype set to FALSE")
                add.subtype <- FALSE
            }
        }
    }

    if (is.null(dir)) {
        message("Argument dir is NULL. Plese provide the directory
                with the folders to be prepared. ")
        return(NULL)
    }

    if (length(unique(query$Platform)) > 1 |
        length(unique(query$Center)) > 2) {
        # This case (27k/450k)accepts two platforms
        if (all(grepl("HumanMethylation[0-9]{2,3}",unique(query$Platform)))){
            platform <- "humanMethylation"
        } else {
            message("Sorry! But, for the moment, we can only prepare on type of
                platform per call")
            return(NULL)
        }
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
    rse <- NULL

    if (grepl("illuminaga_dnaseq",tolower(platform))) {

        message(
            paste("Sorry, but for this platform we haven't prepared",
                  "the data into a summarizedExperiment object.",
                  "\nBut we will do it soon! The return is a data frame of the maf file")
        )

        idx <- grep("maf",files)
        if (length(idx) == 0) {
            message("Sorry, we are preparing only maf files")
            return(NULL)
        }
        files <- files[idx]
        for (i in seq_along(files)) {
            data <- read.table(files[i], fill = TRUE,
                               comment.char = "#", header = TRUE, sep = "\t", quote="")

            # some center has different ways to name the collums
            idx.pos <- grep("Start_position", colnames(data))
            colnames(data)[idx.pos] <- "Start_Position"
            idx.pos <- grep("End_position", colnames(data))
            colnames(data)[idx.pos] <- "End_Position"

            if (i == 1) {
                df <- data
            } else {
                df <-  plyr::rbind.fill(df, data)
            }
            setTxtProgressBar(pb, i)
        }
    }

    if (tolower(platform) == "illuminaga_rnaseq" ||
        tolower(platform) == "illuminahiseq_rnaseq") {


        if(is.null(type) || (type != "exon.quantification" &&
                             type != "spljxn.quantification" &&
                             type != "gene.quantification")
        ){
            msg <- paste0("Plase select a type. \n Possibilities:\n",
                          " = gene.quantification\n",
                          " = spljxn.quantification\n",
                          " = exon.quantification\n"
            )
            message(msg)
            return()
        }

        files <- files[grep(type,files)]
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
                                     gene_id = merged$external_gene_id,
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
            colData <- colDataPrepare(as.character(barcode),
                                      query,add.subtype = add.subtype,add.clinical = add.clinical, add.mutation.genes = add.mutation.genes)
            rse <- SummarizedExperiment(assays=assays,
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
            data <- fread(files[i], header = TRUE, sep = "\t", skip = 1,
                          stringsAsFactors = FALSE, data.table = FALSE)

            if (i == 1) {
                df <- data
            } else {
                df <- cbind(df, data[,2])
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
            #df$external_gene_id <-  alias2SymbolTable(df$`Composite Element REF`)
            df$external_gene_id <-  df$`Composite Element REF`
            merged <- merge(df,gene.location,by="external_gene_id")
            rowRanges <- GRanges(seqnames = paste0("chr", merged$chromosome_name),
                                 ranges = IRanges(start = merged$start_position,
                                                  end = merged$end_position),
                                 strand=merged$strand,
                                 gene_id = merged$external_gene_id,
                                 entrezgene = merged$entrezgene,
                                 alias = merged$`Composite Element REF`)
            names(rowRanges) <- as.character(merged$`Composite Element REF`)

            regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                            "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
            barcode <- unique(unlist(str_match_all(colnames(merged),regex)))
            colData <- colDataPrepare(barcode,query,
                                      add.subtype = add.subtype, add.clinical = add.clinical, add.mutation.genes = add.mutation.genes)

            assays <- SimpleList(raw_counts=data.matrix(
                subset(merged,select=seq(3,2+length(barcode)))
            )
            )

            rse <- SummarizedExperiment(assays=assays,
                                        rowRanges=rowRanges,
                                        colData=colData)
        }
    }

    if (tolower(platform) == tolower("HG-U133_Plus_2") ||
        grepl("H-miRNA_8x15K|agilent",platform, ignore.case = TRUE)) {

        for (i in seq_along(files)) {

            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
            if (i == 1) {
                df <- data
            } else {
                df <- merge(df, data, colnames(data)[1])
            }
            #setTxtProgressBar(pb, i)
        }
        df <- df[-1,]

        if(summarizedExperiment){
            message("===================================================================================")
            message(" As we can't map some miRNA to genomic positions this step might loose some rows .")
            message(" Please, for all rows run TCGAprepare with summarizedExperiment=F")
            message("====================================================================================")

            if(grepl("HG-U133_Plus_2|agilent",platform, ignore.case = TRUE)){
                suppressWarnings(
                    df$external_gene_id <-  alias2SymbolTable(df$`Hybridization REF`)
                )
                merged <- merge(df,gene.location,by="external_gene_id")
                rowRanges <- GRanges(seqnames = paste0("chr", merged$chromosome_name),
                                     ranges = IRanges(start = merged$start_position,
                                                      end = merged$end_position),
                                     strand=merged$strand,
                                     gene_id = merged$external_gene_id,
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
            colData <- colDataPrepare(barcode, query,
                                      add.subtype = add.subtype, add.clinical = add.clinical, add.mutation.genes = add.mutation.genes)

            suppressWarnings(
                assays <- SimpleList(raw_counts=data.matrix(
                    subset(merged,select=seq(3,2+length(barcode)))))
            )

            rse <- SummarizedExperiment(assays=assays,
                                        rowRanges=rowRanges,
                                        colData=colData)
        } else {
            df <- as.data.frame(df)
            rownames(df) <- df[,1]
            df[,1] <- NULL
            df <- data.matrix(df)
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

        files <- files[grep(type,files, perl = TRUE)]

        regex <- paste0("[:alnum:]{8}-[:alnum:]{4}",
                        "-[:alnum:]{4}-[:alnum:]{4}-[:alnum:]{12}")
        uuid <- str_match(files,regex)
        map <- mapuuidbarcode(uuid)

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
            x <- subset(map, uuid == uuid[i])

            if( length(x$barcode)!=0){

                if (summarizedExperiment) {
                    setnames(data,colnames(data)[2:ncol(data)],
                             paste0(colnames(data)[2:ncol(data)],"_",x$barcode))
                } else {
                    setnames(data,2, as.character(x$barcode))
                }

            } else {
                next
            }

            if (i == 1) {
                df <- data
            } else {
                df <- cbind(df, data[,2:ncol(data),with = FALSE])
            }
            setTxtProgressBar(pb, i)
        }

        if (summarizedExperiment){
            if (grepl("gene_id",colnames(df)[1])) {
                aux <- strsplit(df$gene_id,"\\|")
                GeneID <- unlist(lapply(aux,function(x) x[2]))
                df$entrezid <- as.numeric(GeneID)
                GeneSymbol <- unlist(lapply(aux,function(x) x[1]))
                df$external_gene_id <- as.character(GeneSymbol)
                df <- merge(df,gene.location,by="external_gene_id")

                rowRanges <- GRanges(seqnames = paste0("chr", df$chromosome_name),
                                     ranges = IRanges(start = df$start_position,
                                                      end = df$end_position),
                                     strand = df$strand,
                                     gene_id = df$external_gene_id,
                                     entrezgene = df$entrezid,
                                     transcript_id = subset(df, select = 5))
                names(rowRanges) <- as.character(df$gene_id)

                if (length(colnames(data)) > 2) {
                    assays <- SimpleList(
                        raw_counts = data.matrix(
                            subset(df,
                                   select = grep("raw_count",colnames(df)))
                        ),
                        scaled_estimate = data.matrix(
                            subset(df,
                                   select = grep("scaled_estimate",colnames(df)))
                        )
                    )
                } else {
                    normalized_count <- data.matrix(subset(df, select = grep("normalized_count",colnames(df))))
                    assays <- list(normalized_count = normalized_count)
                }
            } else if (grepl("junction",colnames(df)[1])){
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
                assays <- SimpleList(
                    raw_counts = data.matrix(subset(df,select=2:ncol(df)))
                )
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
                raw_counts = data.matrix(
                    df[,grep("raw_count",colnames(df)),with = FALSE]
                )
                median_length_normalized=data.matrix(
                    df[,grep("median_length",colnames(df)),with = FALSE]
                )
                RPKM=data.matrix(
                    df[,grep("RPKM",colnames(df)),with = FALSE]
                )
                assays <- list( raw_counts = raw_counts, median_length_normalized=median_length_normalized, RPKM=RPKM)
            } else if(grepl("isoform",colnames(df)[1])){
                message("TBD")
                return (NULL)
            }
            regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                            "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")

            barcode <- unique(unlist(str_match_all(colnames(df),regex)))
            colData <- colDataPrepare(barcode,query,
                                      add.subtype = add.subtype,add.clinical = add.clinical, add.mutation.genes = add.mutation.genes)

            assays <- lapply(assays, function(x){
                colnames(x) <- barcode
                return(x)
            })

            rse <- SummarizedExperiment(assays=assays,
                                        rowRanges=rowRanges,
                                        colData=colData)
        } else {
            setDF(df)
            rownames(df) <- df[,1]
            df[,1] <- NULL
        }
    }

    if (grepl("illuminahiseq_mirnaseq",platform, ignore.case = TRUE) ||
        grepl("illuminaga_mirnaseq",platform, ignore.case = TRUE)) {

        if (is.null(type) || (type != "isoform.quantification" &&
                              type != "hg19.mirbase20.isoform.quantification" &&
                              type != "hg19.mirbase20.mirna.quantification" &&
                              type != "mirna.quantification")){
            msg <- paste0("Plase select a type. \n Possibilities:\n",
                          "\n = hg19.mirbase20.mirna.quantification",
                          "\n = mirna.quantification",
                          "\n = hg19.mirbase20.isoform.quantification",
                          "\n = isoform.quantification")
            message(msg)
            return()
        }

        if(type == "hg19.mirbase20.mirna.quantification") regex <- "hg19.mirbase20.mirna.quantification"
        if(type == "hg19.mirbase20.isoform.quantification") regex <- "hg19.mirbase20.isoform.quantification"
        if(type == "isoform.quantification" ) regex <- "[^(hg19.mirbase20)].isoform.quantification"
        if(type == "mirna.quantification" ) regex <- "[^(hg19.mirbase20)].mirna.quantification"

        files <- files[grep(regex,files)]

        if(length(files) == 0){
            message("No mirna files of that type found")
            return(NULL)
        }

        if(summarizedExperiment){
            message(
                paste("Sorry, but for this platform we haven't prepared",
                      "the data into a summarizedExperiment object.",
                      "\nBut we will do it soon! The return is a data frame")
            )
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
            setTxtProgressBar(pb, i)
        }
        setDF(df)
        rownames(df) <- df[,1]
        df <- df[,-1]
    }

    if (grepl("bio",platform,ignore.case = TRUE)) {
        if (!is.null(type)) {
            if(grepl("clinical_follow_up",type))  type <- paste0(type,"_[^(nte)]")
            files <- files[grep(type,files)]
        }
        if (length(files) == 1) {
            df <- read.table(files, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, check.names = FALSE,
                             comment.char = "#",fill = TRUE,quote="")
            regex <- paste0("[[:alnum:]]{8}-[[:alnum:]]{4}",
                            "-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}")
            if (grepl("clinical",type)) colnames(df) <- df[1,]
            idx <- grepl(regex,df$bcr_patient_uuid)
            df <- df[idx,]
        } else {
            message("We're preaparing for the moment only one clinical file")
            return(NULL)
        }
    }

    if (grepl("HuEx-1_0-st-v2",tolower(platform),ignore.case = TRUE)) {
        print("HuEx-1_0-st-v2")
        if(!(type %in% c("FIRMA","gene"))){
            stop("Please, set type argument to FIRMA or gene (type = 'gene' or type = 'FIRMA'")
        }

        files <- files[grep(type,files)]

        if(length(files) == 0){
            message("No files of that type found")
            return(NULL)
        }

        if(is.vector(query)){
            mage <- getMage(query)
        } else {
            mage <- getMage(query[1,])
        }

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, data.table = FALSE)
            if(i == 1) df <- data
            if(i != 1) df <- merge(df,data,by="Hybridization REF")
            setTxtProgressBar(pb, i)
        }


        # The df object is:
        #
        # Hibridation ref 123 345
        # gene 1
        # gene 2
        #
        # The mage has the map btw  Hibridation ref and barcode
        codes <- colnames(df)
        codes <- codes[-1] # remove the Hybridization REF
        barcodes <- mage[match(as.integer(codes), as.integer(mage$Hybridization.Name)),"Comment..TCGA.Barcode."]
        not.found <- which(barcodes == "->")
        barcodes[not.found] <- codes[not.found] # some does not have a barcode
        colnames(df)[2:length(colnames(df))] <- barcodes

        if (summarizedExperiment & type == "FIRMA" ) {
            message("Sorry we can prepare FIRMA type into summarizedExperiment please set summarizedExperiment = FALSE")
            summarizedExperiment <- FALSE
        }
        if (summarizedExperiment){
            df$external_gene_id <- df[,1]

            df <- merge(df,gene.location,by="external_gene_id")

            rowRanges <- GRanges(seqnames = paste0("chr", df$chromosome_name),
                                 ranges = IRanges(start = df$start_position,
                                                  end = df$end_position),
                                 strand = df$strand,
                                 gene_id = df$external_gene_id,
                                 entrezgene = df$entrezgene)
            names(rowRanges) <- as.character(df$external_gene_id)
            assays <- SimpleList(
                signal = data.matrix(df[,3:(length(barcodes)+2)],rownames.force = T)
            )

            colData <- colDataPrepare(barcodes,query,
                                      add.subtype = add.subtype, add.clinical = add.clinical, add.mutation.genes = add.mutation.genes)

            rownames(colData) <- barcodes
            assays <- lapply(assays, function(x){
                colnames(x) <- barcodes
                rownames(x) <- as.character(df$external_gene_id)
                return(x)
            })
            rse <- SummarizedExperiment(assays=assays,
                                        rowRanges=rowRanges,
                                        colData=colData)
        }

    }

    if (grepl("CGH-1x1M_G4447A",tolower(platform),ignore.case = TRUE)) {

        if(summarizedExperiment){
            message(
                paste("Sorry, but for this platform we haven't prepared",
                      "the data into a summarizedExperiment object.",
                      "\n The return is a data frame")
            )
            summarizedExperiment = FALSE
        }

        if(is.vector(query)){
            mage <- getMage(query)
        } else {
            mage <- getMage(query[1,])
        }

        for (i in seq_along(files)) {
            data <- read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)

            if(i == 1) df <- data
            if(i != 1) df <- rbind(df, data, make.row.names = FALSE)

            setTxtProgressBar(pb, i)
        }
        mage <- mage[,c("Comment..TCGA.Barcode.","Hybridization.Name")]
        mage <- mage[!mage$Comment..TCGA.Barcode. == "->",]
        print(head(mage))
        print(head(df))
        df$sample <- gsub("\\.","-",df$sample)
        df <- merge(df,mage,
                    by.x="sample",by.y="Hybridization.Name", sort=FALSE)
        print(head(df))
        df[,1] <- df[,ncol(df)]
        df[,ncol(df)] <- NULL

    }

    if (grepl("HG-CGH-244A|HG-CGH-415K_G4124A",tolower(platform),ignore.case = TRUE)) {

        if(summarizedExperiment){
            message(
                paste("Sorry, but for this platform we haven't prepared",
                      "the data into a summarizedExperiment object.",
                      "\n The return is a data frame")
            )
            summarizedExperiment = FALSE
        }

        if(is.vector(query)){
            mage <- getMage(query)
            center <- query$Center
        } else {
            mage <- getMage(query[1,])
            center <- query[1,]$Center
        }
        if(query[1,]$Center == "hms.harvard.edu"){
            regex <- paste0("[:alnum:]{4}-[:alnum:]{2}-[:alnum:]{4}",
                            "-[:alnum:]{3}-[:alnum:]{3}-[:alnum:]{4}-[:alnum:]{2}")
            barcode <- str_match(files,regex)
            for (i in seq_along(files)) {
                data <- read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
                data$sample <- barcode[i]
                if(i == 1) df <- data
                if(i != 1) df <- rbind(df, data, make.row.names = FALSE)

                setTxtProgressBar(pb, i)
            }
            df <- df[,c(ncol(df),1:(ncol(df)-1))]
        } else {
            for (i in seq_along(files)) {
                data <- read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
                if(i == 1) df <- data
                if(i != 1) df <- rbind(df, data, make.row.names = FALSE)
                setTxtProgressBar(pb, i)
            }
            mage <- mage[,c("Comment..TCGA.Barcode.","Hybridization.Name")]
            mage <- mage[!mage$Comment..TCGA.Barcode. == "->",]
            df$sample <- gsub("\\.","-",df$sample)
            df <- merge(df,mage,
                        by.x="sample",by.y="Hybridization.Name", sort=FALSE)
            df[,1] <- df[,ncol(df)]
            df[,ncol(df)] <- NULL
        }
    }

    if (grepl("genome_wide_snp_6",tolower(platform))) {

        while(!(type %in% c("nocnv_hg18","nocnv_hg19","cnv_hg18","cnv_hg19",
                            "nocnv_hg18.seg","hg18.seg","hg19.seg","nocnv_hg19.seg"))){
            type <- readline(
                paste("Which type do you want?",
                      "(Options: nocnv_hg19,nocnv_hg18,cnv_hg18,cnv_hg19, cancel)  ")
            )
            if(type == "cancel") return(NULL)
        }

        if(type == "nocnv_hg18" | type == "nocnv_hg18.seg") regex <- "nocnv_hg18"
        if(type == "cnv_hg18" | type == "hg18.seg") regex <- "[^nocnv_]hg18.seg"
        if(type == "nocnv_hg19" | type == "nocnv_hg19.seg") regex <- "nocnv_hg19"
        if(type == "cnv_hg19" | type == "hg19.seg") regex <- "[^nocnv_]hg19.seg"

        files <- files[grep(regex,files)]

        if(length(files) == 0){
            message("No files of that type found")
            return(NULL)
        }

        idx <- grep(regex, files)
        if (length(idx) > 0){
            files <- files[idx]
        } else {
            message("No files of that type found")
            return (NULL)
        }

        if(is.vector(query)){
            mage <- getMage(query)
        } else {
            mage <- getMage(query[1,])
        }

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, data.table = FALSE,
                          colClasses=c("character", # ID
                                       "character",   # chrom
                                       "numeric", # start
                                       "numeric",   # end
                                       "integer", # num_probes
                                       "numeric"))  # seg mean
            if(i == 1) df <- data
            if(i != 1) df <- rbind(df, data, make.row.names = FALSE)

            setTxtProgressBar(pb, i)
        }
        mage <- mage[,c("Comment..TCGA.Barcode.","Hybridization.Name")]
        df <- merge(df,mage,
                    by.x="Sample",by.y="Hybridization.Name", sort=FALSE)

        df[,1] <- df[,7]
        df[,7] <- NULL

    }
    close(pb)

    if (!is.null(rse)) {
        message("Adding metadata to the rse object...")

        finf <- c()
        finf <- file.info(files)
        rownames(finf) <- basename(rownames(finf))
        finf <- finf[,c("mtime","ctime")]

        metadata(rse) <- list("Query:"=list(query),
                              "TCGAprepareParameters"=c("dir"=dir,
                                                        "samples"=samples,
                                                        "type"=type,
                                                        "save"=save,
                                                        "add.mutation.genes"=add.mutation.genes,
                                                        "filename"=filename),
                              "FilesInfo:"=list(finf))
    }

    if (save) {
        message("Saving the data...")

        if (is.null(filename)) {
            filename <- paste0(platform,"_",gsub(" ","_",Sys.time()),".rda")
        }
        if (!is.null(rse)) {
            save(rse,file = filename)
        } else {
            save(df,file = filename)
        }
        message(paste("Data saved in:",filename))
    }

    if (!is.null(rse)) {
        return(rse)
    }
    return(df)
}

#' @title Prepare the data for ELEMR package
#' @description Prepare the data for ELEMR package
#' @return Matrix prepared for fetch.mee function
#' @param data A data frame or summarized experiment from TCGAPrepare
#' @param platform platform of the data. Example: "HumanMethylation450", "IlluminaHiSeq_RNASeqV2"
#' @param met.na.cut Define the percentage of NA that the line should have to
#'  remove the probes for humanmethylation platforms.
#' @param save Save object? Default: FALSE.
#' Names of the files will be: "Exp_elmer.rda" (object Exp) and "Met_elmer.rda" (object Met)
#' @export
#' @examples
#' df <- data.frame(runif(200, 1e5, 1e6),runif(200, 1e5, 1e6))
#' rownames(df) <- sprintf("?|%03d", 1:200)
#' df <- TCGAprepare_elmer(df,platform="IlluminaHiSeq_RNASeqV2")
TCGAprepare_elmer <- function(data,
                              platform,
                              met.na.cut = 0.2,
                              save = FALSE){
    # parameters veryfication

    if (missing(data))  stop("Please set the data parameter")
    if (missing(platform))  stop("Please set the platform parameter")

    if (grepl("illuminahiseq_rnaseqv2|illuminahiseq_totalrnaseqv2",
              platform, ignore.case = TRUE)) {
        message("============ Pre-pocessing expression data =============")
        message(paste0("1 - expression = log2(expression + 1): ",
                       "To linearize \n    relation between ",
                       "methylation and expression"))
        if(typeof(data) == typeof(SummarizedExperiment())){
            data <- assay(data)
        }

        data <- log2(data+1)
        message("2 - rownames  (gene|loci) => ('ID'loci) ")
        aux <- strsplit(rownames(data),"\\|")
        GeneID <- unlist(lapply(aux,function(x) x[2]))
        row.names(data) <- paste0("ID",GeneID)
        Exp <- data.matrix(data)

        if (save)  save(Exp,file = "Exp_elmer.rda")
        return(Exp)
    }

    if (grepl("humanmethylation", platform, ignore.case = TRUE)) {
        message("============ Pre-pocessing methylation data =============")
        if (class(data) == class(data.frame())){
            msg <- paste0("1 - Removing Columns: \n  * Gene_Symbol  \n",
                          "  * Chromosome  \n  * Genomic_Coordinate")
            message(msg)
            data <- subset(data,select = 4:ncol(data))
        }
        if(typeof(data) == typeof(SummarizedExperiment())){
            data <- assay(data)
        }
        msg <- paste0("2 - Removing probes with ",
                      "NA values in more than 20% samples")
        message(msg)
        data <- data[rowMeans(is.na(data)) < met.na.cut,]
        Met <- data.matrix(data)
        if (save)  save(Met,file = "Met_elmer.rda")
        return (Met)
    }
}

#' @title Prepare CEL files into an AffyBatch.
#' @description Prepare CEL files into an AffyBatch.
#' @param ClinData write
#' @param PathFolder write
#' @param TabCel write
#' @importFrom affy ReadAffy
#' @importFrom affy rma
#' @importFrom Biobase exprs
#' @examples
#' \dontrun{
#' to add example
#' }
#' @export
#' @return Normalizd Expression data from Affy eSets
TCGAprepare_Affy <- function(ClinData, PathFolder, TabCel){

    affy_batch <- ReadAffy(filenames=as.character(paste(TabCel$samples, ".CEL", sep="")))

    eset <- rma(affy_batch)

    mat <- exprs(eset)

    return(mat)

}



mutation.genes <- function(tumor = NULL, data=NULL){
    df <- TCGAquery_maf(tumor)
    DT <- data.table(df)
    mutated.genes <- with(DT, {
        mutated.genes <- DataFrame(
            DT[, list(genes = list(as.character(Hugo_Symbol))), by = "bcr_patient_barcode"]
        )
    })
    colnames(mutated.genes)[1] <- "patient"

    if(!is.null(data)){
        df <- merge(data,mutated.genes,all.x = TRUE, sort = FALSE, all.y= FALSE)
        df <- df[match(data$patient,df$patient),]
    } else {
        df <- mutated.genes
    }
    return(df)
}
