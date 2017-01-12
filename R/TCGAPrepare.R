
#' @title Prepare GDC data
#' @description
#'   Reads the data downloaded and prepare it into an R object
#' @param query A query for GDCquery function
#' @param save Save result as RData object?
#' @param save.filename Name of the file to be save if empty an automatic will be created
#' @param directory Directory/Folder where the data was downloaded. Default: GDCdata
#' @param summarizedExperiment Create a summarizedExperiment? Default TRUE (if possible)
#' @param remove.files.prepared Remove the files read? Default: FALSE
#' This argument will be considered only if save argument is set to true
#' @param add.gistic2.mut If a list of genes (gene symbol) is given columns with gistic2 results from GDAC firehose and
#' a column indicating if there is or not mutation in that gene (TRUE or FALSE - use the maf file for more information)
#' will be added to the sample matrix in the summarized Experiment object
#' @param mut.pipeline If add.gistic2.mut is not NULL this field will be taken in consideration.
#' Four separate variant calling pipelines are implemented for GDC data harmonization.
#' Options: muse, varscan2, somaticsniper, MuTect2. For more information:
#' https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/
#' @export
#' @examples
#' query <- GDCquery(project = "TCGA-KIRP",
#'                   data.category = "Simple Nucleotide Variation",
#'                   data.type = "Masked Somatic Mutation",
#'                   workflow.type = "MuSE Variant Aggregation and Masking")
#' GDCdownload(query, method = "api", directory = "maf")
#' maf <- GDCprepare(query, directory = "maf")
#'
#' query <- GDCquery(project = "TCGA-ACC",
#'                    data.category =  "Copy number variation",
#'                    legacy = TRUE,
#'                    file.type = "hg19.seg",
#'                    barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01", "TCGA-OR-A5LJ-10A-01D-A29K-01"))
#' # data will be saved in  GDCdata/TCGA-ACC/legacy/Copy_number_variation/Copy_number_segmentation
#' GDCdownload(query, method = "api")
#' acc.cnv <- GDCprepare(query)
#' \dontrun{
#'  query <- GDCquery(project = "TCGA-GBM",
#'                    legacy = TRUE,
#'                    data.category = "Gene expression",
#'                    data.type = "Gene expression quantification",
#'                    platform = "Illumina HiSeq",
#'                    file.type = "normalized_results",
#'                    experimental.strategy = "RNA-Seq")
#'  GDCdownload(query, method = "api")
#'  data <- GDCprepare(query,add.gistic2.mut = c("PTEN","FOXJ1"))
#' }
#' @return A summarizedExperiment or a data.frame
#' @importFrom  S4Vectors DataFrame
#' @importFrom data.table setcolorder setnames
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
GDCprepare <- function(query,
                       save = FALSE,
                       save.filename,
                       directory = "GDCdata",
                       summarizedExperiment = TRUE,
                       remove.files.prepared = FALSE,
                       add.gistic2.mut = NULL,
                       mut.pipeline = "mutect2"){

    if(missing(query)) stop("Please set query parameter")
    if(any(duplicated(query$results[[1]]$cases))) {
        dup <- query$results[[1]]$cases[duplicated(query$results[[1]]$cases)]
        dup <- query$results[[1]][sapply(dup, function(x) grep(x,query$results[[1]]$cases)),c(2,4,9,13,15,18)]
        print(knitr::kable(dup))
        stop("There are samples duplicated. We will not be able to preapre it")
    }
    if(!save & remove.files.prepared) {
        stop("To remove the files, please set save to TRUE. Otherwise, the data will be lost")
    }
    # We save the files in project/source/data.category/data.type/file_id/file_name
    source <- ifelse(query$legacy,"legacy","harmonized")
    files <- file.path(query$project, source,
                       gsub(" ","_",query$results[[1]]$data_category),
                       gsub(" ","_",query$results[[1]]$data_type),
                       gsub(" ","_",query$results[[1]]$file_id),
                       gsub(" ","_",query$results[[1]]$file_name))

    files <- file.path(directory, files)

    if(!all(file.exists(files))) stop(paste0("I couldn't find all the files from the query. ",
                                             "Please check if the directory parameter right or GDCdownload downloaded the samples."))

    if(grepl("Transcriptome Profiling", query$data.category, ignore.case = TRUE)){
        data <- readTranscriptomeProfiling(files = files,
                                           data.type = ifelse(!is.na(query$data.type),  as.character(query$data.type),  unique(query$results[[1]]$data_type)),
                                           workflow.type = unique(query$results[[1]]$analysis$workflow_type),
                                           cases = query$results[[1]]$cases,
                                           summarizedExperiment)
    } else if(grepl("Copy Number Variation",query$data.category,ignore.case = TRUE)) {
        data <- readCopyNumberVariation(files, query$results[[1]]$cases)
    }  else if(grepl("DNA methylation",query$data.category, ignore.case = TRUE)) {
        data <- readDNAmethylation(files, query$results[[1]]$cases, summarizedExperiment, unique(query$platform))
    }  else if(grepl("Protein expression",query$data.category,ignore.case = TRUE)) {
        data <- readProteinExpression(files, query$results[[1]]$cases)
    }  else if(grepl("Simple Nucleotide Variation",query$data.category,ignore.case = TRUE)) {
        if(grepl("Masked Somatic Mutation",query$data.type,ignore.case = TRUE) | source == "legacy")
            suppressWarnings(data <- readSimpleNucleotideVariationMaf(files))
    }  else if(grepl("Clinical|Biospecimen", query$data.category, ignore.case = TRUE)){
        message("Mot working yet")
        # data <- readClinical(files, query$results[[1]]$cases)
    } else if (grepl("Gene expression",query$data.category,ignore.case = TRUE)) {
        if(query$data.type == "Gene expression quantification")
            data <- readGeneExpressionQuantification(files = files,
                                                     cases = query$results[[1]]$cases,
                                                     summarizedExperiment = summarizedExperiment,
                                                     experimental.strategy = unique(query$results[[1]]$experimental_strategy))

        if(query$data.type == "miRNA gene quantification")
            data <- readGeneExpressionQuantification(files = files,
                                                     cases = query$results[[1]]$cases,
                                                     summarizedExperiment = F,
                                                     experimental.strategy = unique(query$results[[1]]$experimental_strategy))
        if(query$data.type == "miRNA isoform quantification")
            data <- readmiRNAIsoformQuantification(files = files,
                                                   cases = query$results[[1]]$cases)

        if(query$data.type == "Isoform expression quantification")
            data <- readIsoformExpressionQuantification(files = files, cases = query$results[[1]]$cases)

    }

    if((!is.null(add.gistic2.mut)) & summarizedExperiment) {
        message("=> Adding GISTIC2 and mutation information....")
        genes <- tolower(levels(EAGenes$Gene))
        if(!all(tolower(add.gistic2.mut) %in% genes)) message(paste("These genes were not found:\n",
                                                                    paste(add.gistic2.mut[! tolower(add.gistic2.mut) %in% genes],collapse = "\n=> ")))
        add.gistic2.mut <- add.gistic2.mut[tolower(add.gistic2.mut) %in% tolower(genes)]
        if(length(add.gistic2.mut) > 0){
            info <- get.mut.gistc.information(colData(data),query$project, add.gistic2.mut, mut.pipeline = mut.pipeline)
            colData(data) <- info
        }
    }

    if(save){
        if(missing(save.filename) & !missing(query)) save.filename <- paste0(query$project,gsub(" ","_", query$data.category),gsub(" ","_",date()),".RData")
        message(paste0("Saving file:",save.filename))
        save(data, file = save.filename)
        message("File saved")

        # save is true, due to the check in the beggining of the code
        if(remove.files.prepared){
            # removes files and empty directories
            remove.files.recursively(files)
        }
    }

    return(data)
}

remove.files.recursively <- function(files){
    files2rm <- dirname(files)
    unlink(files2rm,recursive = TRUE)
    files2rm <- dirname(files2rm) # data category
    if(length(list.files(files2rm)) == 0) remove.files.recursively(files2rm)
}


readIsoformExpressionQuantification <- function (files, cases){
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)

    for (i in seq_along(files)) {
        data <- fread(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)

        if(!missing(cases)) {
            assay.list <- gsub(" |\\(|\\)|\\/","_",colnames(data)[2:ncol(data)])
            # We will use this because there might be more than one col for each samples
            setnames(data,colnames(data)[2:ncol(data)],
                     paste0(gsub(" |\\(|\\)|\\/","_",colnames(data)[2:ncol(data)]),"_",cases[i]))
        }
        if (i == 1) {
            df <- data
        } else {
            df <- merge(df, data, by=colnames(data)[1], all = TRUE)
        }
        setTxtProgressBar(pb, i)
    }
    setDF(df)
    rownames(df) <- df[,1]
    df[,1] <- NULL
    return(df)

}
readmiRNAIsoformQuantification <- function (files, cases){
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)

    for (i in seq_along(files)) {
        data <- fread(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        data$barcode <- cases[i]
        if (i == 1) {
            df <- data
        } else {
            df <- rbind(df, data)
        }
        setTxtProgressBar(pb, i)
    }
    setDF(df)

}
readSimpleNucleotideVariationMaf <- function(files){
    ret <- read_tsv(files,
                    comment = "#",
                    col_types = cols(
                        Entrez_Gene_Id = col_integer(),
                        Start_Position = col_integer(),
                        End_Position = col_integer(),
                        t_depth = col_integer(),
                        t_ref_count = col_integer(),
                        t_alt_count = col_integer(),
                        n_depth = col_integer(),
                        ALLELE_NUM = col_integer(),
                        TRANSCRIPT_STRAND = col_integer(),
                        PICK = col_integer(),
                        TSL = col_integer(),
                        HGVS_OFFSET = col_integer(),
                        MINIMISED = col_integer()),
                    progress = TRUE)
    if(ncol(ret) == 1) ret <- read_csv(files,
                                       comment = "#",
                                       col_types = cols(
                                           Entrez_Gene_Id = col_integer(),
                                           Start_Position = col_integer(),
                                           End_Position = col_integer(),
                                           t_depth = col_integer(),
                                           t_ref_count = col_integer(),
                                           t_alt_count = col_integer(),
                                           n_depth = col_integer(),
                                           ALLELE_NUM = col_integer(),
                                           TRANSCRIPT_STRAND = col_integer(),
                                           PICK = col_integer(),
                                           TSL = col_integer(),
                                           HGVS_OFFSET = col_integer(),
                                           MINIMISED = col_integer()),
                                       progress = TRUE)
    return(ret)
}

readGeneExpressionQuantification <- function(files, cases, summarizedExperiment = TRUE, experimental.strategy, platform){
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)

    skip <- unique((ifelse(experimental.strategy == "Gene expression array",1,0)))

    if(length(skip) > 1) stop("It is not possible to handle this different platforms together")

    for (i in seq_along(files)) {
        data <- fread(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE,skip = skip)

        if(!missing(cases)) {
            assay.list <- gsub(" |\\(|\\)|\\/","_",colnames(data)[2:ncol(data)])
            # We will use this because there might be more than one col for each samples
            setnames(data,colnames(data)[2:ncol(data)],
                     paste0(gsub(" |\\(|\\)|\\/","_",colnames(data)[2:ncol(data)]),"_",cases[i]))
        }
        if (i == 1) {
            df <- data
        } else {
            df <- merge(df, data, by=colnames(data)[1], all = TRUE)
        }
        setTxtProgressBar(pb, i)
    }
    setDF(df)

    if (summarizedExperiment) {
        df <- makeSEfromGeneExpressionQuantification(df,assay.list)
    } else {
        rownames(df) <- df$gene_id
        df$gene_id <- NULL
    }
    return(df)
}
makeSEfromGeneExpressionQuantification <- function(df, assay.list, genome="hg19"){
    gene.location <- get.GRCh.bioMart(genome)
    if(all(grepl("\\|",df[,1]))){
        aux <- strsplit(df$gene_id,"\\|")
        GeneID <- unlist(lapply(aux,function(x) x[2]))
        df$entrezgene <- as.numeric(GeneID)
    } else {
        df$external_gene_id <- as.character(df[,1])
    }

    df <- merge(df, gene.location, by="entrezgene")

    if("transcript_id" %in% assay.list){
        rowRanges <- GRanges(seqnames = paste0("chr", df$chromosome_name),
                             ranges = IRanges(start = df$start_position,
                                              end = df$end_position),
                             strand = df$strand,
                             gene_id = df$external_gene_id,
                             entrezgene = df$entrezgene,
                             ensembl_gene_id = df$ensembl_gene_id,
                             transcript_id = subset(df, select = 5))
        names(rowRanges) <- as.character(df$gene_id)
        assay.list <- assay.list[which(assay.list != "transcript_id")]
    } else {
        rowRanges <- GRanges(seqnames = paste0("chr", df$chromosome_name),
                             ranges = IRanges(start = df$start_position,
                                              end = df$end_position),
                             strand = df$strand,
                             gene_id = df$external_gene_id,
                             entrezgene = df$entrezgene,
                             ensembl_gene_id = df$ensembl_gene_id)
        names(rowRanges) <- as.character(df$external_gene_id)
    }
    assays <- lapply(assay.list, function (x) {
        return(data.matrix(subset(df, select = grep(x,colnames(df),ignore.case = TRUE))))
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
    rse <- SummarizedExperiment(assays=assays,
                                rowRanges=rowRanges,
                                colData=colData)
    return(rse)
}

makeSEfromDNAmethylation <- function(df, probeInfo=NULL){
    if(is.null(probeInfo)) {
        gene.location <- get.GRCh.bioMart()
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
    } else {
        rowRanges <- makeGRangesFromDataFrame(probeInfo, keep.extra.columns = T)
        colData <-  colDataPrepare(colnames(df)[(ncol(probeInfo) + 1):ncol(df)])
        assay <- data.matrix(subset(df,select = c((ncol(probeInfo) + 1):ncol(df))))
    }
    colnames(assay) <- rownames(colData)
    rownames(assay) <- as.character(df$Composite.Element.REF)
    rse <- SummarizedExperiment(assays = assay, rowRanges = rowRanges, colData = colData)
}

# We will try to make this function easier to use this function in its own data
# In case it is not TCGA I should not consider that there is a barcode in the header
# Instead the user should be able to add the names to his data
# The only problem is that the data from the user will not have all the columns
# TODO: Improve this function to be more generic as possible
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tibble as_data_frame
readDNAmethylation <- function(files, cases, summarizedExperiment = TRUE, platform){
    if (grepl("OMA00",platform)){
        pb <- txtProgressBar(min = 0, max = length(files), style = 3)
        for (i in seq_along(files)) {
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
        skip <- ifelse(all(grepl("hg38",files)), 0,1)
        colClasses <- NULL
        if(!all(grepl("hg38",files))) colClasses <- c("character", # Composite Element REF
                                                      "numeric",   # beta value
                                                      "character", # Gene symbol
                                                      "character", # Chromosome
                                                      "integer")

        for (i in seq_along(files)) {
            data <- fread(files[i], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE,skip = skip, colClasses = colClasses)
            setnames(data,gsub(" ", "\\.", colnames(data)))
            if(!missing(cases)) setnames(data,2,cases[i])
            if (i == 1) {
                setcolorder(data,c(1, 3:ncol(data), 2))
                df <- data
            } else {
                data <- subset(data,select = c(1,2))
                df <- merge(df, data, by = "Composite.Element.REF")
            }
            setTxtProgressBar(pb, i)
        }
        if (summarizedExperiment) {
            if(skip == 0) {
                df <- makeSEfromDNAmethylation(df, probeInfo = as_data_frame(df)[,grep("TCGA",colnames(df),invert = TRUE)])
            } else {
                df <- makeSEfromDNAmethylation(df)
            }
        } else {
            setDF(df)
            rownames(df) <- df$Composite.Element.REF
            df$Composite.Element.REF <- NULL
        }
    }
    return(df)
}

colDataPrepareTARGET <- function(barcode){
    message("Adding description to TARGET samples")
    tissue.code <- c('01','02','03','04','05','06','07','08','09','10','11',
                     '12','13','14','15','16','17','20','40','41','42','50','60','61','99')

    tissue.definition <- c("Primary solid Tumor", # 01
                           "Recurrent Solid Tumor", # 02
                           "Primary Blood Derived Cancer - Peripheral Blood", # 03
                           "Recurrent Blood Derived Cancer - Bone Marrow", # 04
                           "Additional - New Primary", # 05
                           "Metastatic", # 06
                           "Additional Metastatic", # 07
                           "Tissue disease-specific post-adjuvant therapy", # 08
                           "Primary Blood Derived Cancer - Bone Marrow", # 09
                           "Blood Derived Normal", # 10
                           "Solid Tissue Normal",  # 11
                           "Buccal Cell Normal",   # 12
                           "EBV Immortalized Normal", # 13
                           "Bone Marrow Normal", # 14
                           "Fibroblasts from Bone Marrow Normal", # 15
                           "Mononuclear Cells from Bone Marrow Normal", # 16
                           "Lymphatic Tissue Normal (including centroblasts)", # 17
                           "Control Analyte", # 20
                           "Recurrent Blood Derived Cancer - Peripheral Blood", # 40
                           "Blood Derived Cancer- Bone Marrow, Post-treatment", # 41
                           "Blood Derived Cancer- Peripheral Blood, Post-treatment", # 42
                           "Cell line from patient tumor", # 50
                           "Xenograft from patient not grown as intermediate on plastic tissue culture dish", # 60
                           "Xenograft grown in mice from established cell lines", #61
                           "Granulocytes after a Ficoll separation") # 99
    aux <- DataFrame(tissue.code = tissue.code,tissue.definition)

    # in case multiple equal barcode
    regex <- paste0("[:alnum:]{5}-[:alnum:]{2}-[:alnum:]{6}",
                    "-[:alnum:]{3}-[:alnum:]{3}")
    samples <- str_match(barcode,regex)[,1]

    ret <- DataFrame(barcode = barcode,
                     patient = substr(barcode, 1, 16),
                     tumor.code = substr(barcode, 8, 9),
                     case.unique.id = substr(barcode, 11, 16),
                     tissue.code = substr(barcode, 18, 19),
                     nucleic.acid.code = substr(barcode, 24, 24))

    ret <- merge(ret,aux, by = "tissue.code", sort = FALSE)

    tumor.code <- c('00','01','02','03','04','10','15','20','21','30','40',
                    '41','50','51','52','60','61','62','63','64','65','70','71','80','81')

    tumor.definition <- c("Non-cancerous tissue", # 00
                          "Diffuse Large B-Cell Lymphoma (DLBCL)", # 01
                          "Lung Cancer (all types)", # 02
                          "Cervical Cancer (all types)", # 03
                          "Anal Cancer (all types)", # 04
                          "Acute lymphoblastic leukemia (ALL)", # 10
                          "Mixed phenotype acute leukemia (MPAL)", # 15
                          "Acute myeloid leukemia (AML)", # 20
                          "Induction Failure AML (AML-IF)", # 21
                          "Neuroblastoma (NBL)", # 30
                          "Osteosarcoma (OS)",  # 40
                          "Ewing sarcoma",   # 41
                          "Wilms tumor (WT)", # 50
                          "Clear cell sarcoma of the kidney (CCSK)", # 51
                          "Rhabdoid tumor (kidney) (RT)", # 52
                          "CNS, ependymoma", # 60
                          "CNS, glioblastoma (GBM)", # 61
                          "CNS, rhabdoid tumor", # 62
                          "CNS, low grade glioma (LGG)", # 63
                          "CNS, medulloblastoma", # 64
                          "CNS, other", # 65
                          "NHL, anaplastic large cell lymphoma", # 70
                          "NHL, Burkitt lymphoma (BL)", # 71
                          "Rhabdomyosarcoma", #80
                          "Soft tissue sarcoma, non-rhabdomyosarcoma") # 81
    aux <- DataFrame(tumor.code = tumor.code,tumor.definition)
    ret <- merge(ret,aux, by = "tumor.code", sort = FALSE)

    nucleic.acid.code <- c('D','E','W','X','Y','R','S')
    nucleic.acid.description <-  c("DNA, unamplified, from the first isolation of a tissue",
                                   "DNA, unamplified, from the first isolation of a tissue embedded in FFPE",
                                   "DNA, whole genome amplified by Qiagen (one independent reaction)",
                                   "DNA, whole genome amplified by Qiagen (a second, separate independent reaction)",
                                   "DNA, whole genome amplified by Qiagen (pool of 'W' and 'X' aliquots)",
                                   "RNA, from the first isolation of a tissue",
                                   "RNA, from the first isolation of a tissue embedded in FFPE")
    aux <- DataFrame(nucleic.acid.code = nucleic.acid.code,nucleic.acid.description)
    ret <- merge(ret,aux, by = "nucleic.acid.code", sort = FALSE)


    ret <- ret[match(barcode,ret$barcode),]
    rownames(ret) <- gsub("\\.","-",make.names(ret$barcode,unique=TRUE))
    ret$code <- NULL
    return(DataFrame(ret))
}

colDataPrepareTCGA <- function(barcode){
    # For the moment this will work only for TCGA Data
    # We should search what TARGET data means

    code <- c('01','02','03','04','05','06','07','08','09','10','11',
              '12','13','14','20','40','50','60','61')
    shortLetterCode <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC",
                         "TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB",
                         "CELL","XP","XCL")

    definition <- c("Primary solid Tumor", # 01
                    "Recurrent Solid Tumor", # 02
                    "Primary Blood Derived Cancer - Peripheral Blood", # 03
                    "Recurrent Blood Derived Cancer - Bone Marrow", # 04
                    "Additional - New Primary", # 05
                    "Metastatic", # 06
                    "Additional Metastatic", # 07
                    "Human Tumor Original Cells", # 08
                    "Primary Blood Derived Cancer - Bone Marrow", # 09
                    "Blood Derived Normal", # 10
                    "Solid Tissue Normal",  # 11
                    "Buccal Cell Normal",   # 12
                    "EBV Immortalized Normal", # 13
                    "Bone Marrow Normal", # 14
                    "Control Analyte", # 20
                    "Recurrent Blood Derived Cancer - Peripheral Blood", # 40
                    "Cell Lines", # 50
                    "Primary Xenograft Tissue", # 60
                    "Cell Line Derived Xenograft Tissue") # 61
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

colDataPrepare <- function(barcode){

    # For the moment this will work only for TCGA Data
    # We should search what TARGET data means
    message("Starting to add information to samples")

    if(all(grepl("TARGET",barcode))) ret <- colDataPrepareTARGET(barcode)
    if(all(grepl("TCGA",barcode))) ret <- colDataPrepareTCGA(barcode)
    message(" => Add clinical information to samples")
    # There is a limitation on the size of the string, so this step will be splited in cases of 100
    patient.info <- NULL
    step <- 50 # more than 100 gives a bug =/
    for(i in 0:(ceiling(length(ret$patient)/step) - 1)){
        start <- 1 + step * i
        end <- ifelse(((i + 1) * step) > length(ret$patient), length(ret$patient),((i + 1) * step))
        if(is.null(patient.info)) {
            patient.info <- getBarcodeInfo(ret$patient[start:end])
        } else {
            patient.info <- rbind(patient.info,getBarcodeInfo(ret$patient[start:end]))
        }
    }
    ret <- merge(ret,patient.info, by.x = "patient", by.y = "submitter_id", all.x = TRUE )

    if(!"project_id" %in% colnames(ret)) {
        aux <- getGDCprojects()[,5:6]
        aux <- aux[aux$disease_type == unique(ret$disease_type),2]
        ret$project_id <- as.character(aux)
    }
    # na.omit should not be here, exceptional case

    if(grepl("TCGA",na.omit(unique(ret$project_id)))) {
        message(" => Adding subtype information to samples")
        tumor <- gsub("TCGA-","",na.omit(unique(ret$project_id)))
        if (grepl("acc|lgg|gbm|luad|stad|brca|coad|read|skcm|hnsc|kich|lusc|ucec|pancan|thca|prad|kirp|kirc|all",
                  tumor,ignore.case = TRUE)) {
            subtype <- TCGAquery_subtype(tumor)
            colnames(subtype) <- paste0("subtype_", colnames(subtype))

            if(all(str_length(subtype$subtype_patient) == 12)){
                # Subtype information were to primary tumor in priority
                subtype$sample <- paste0(subtype$subtype_patient,"-01A")
            }

            ret <- merge(ret,subtype, by = "sample", all.x = TRUE)

        }
    }

    # Add purity information from http://www.nature.com/articles/ncomms9971
    # purity  <- getPurityinfo()
    # ret <- merge(ret, purity, by = "sample", all.x = TRUE, sort = FALSE)

    # Put data in the right order
    ret <- ret[match(barcode,ret$barcode),]
    rownames(ret) <- ret$barcode
    return(ret)
}

#' @importFrom biomaRt getBM useMart listDatasets
get.GRCh.bioMart <- function(genome="hg19") {
    tries <- 0L
    msg <- character()
    while (tries < 3L) {
        gene.location <- tryCatch({
            if (genome == "hg19"){
                # for hg19
                ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                   host = "feb2014.archive.ensembl.org",
                                   path = "/biomart/martservice" ,
                                   dataset = "hsapiens_gene_ensembl")
                attributes <- c("chromosome_name",
                                "start_position",
                                "end_position", "strand",
                                "ensembl_gene_id", "entrezgene",
                                "external_gene_id")
            } else {
                # for hg38
                ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
                attributes <- c("chromosome_name",
                                "start_position",
                                "end_position", "strand",
                                "ensembl_gene_id",
                                "external_gene_name")
            }
            description <- listDatasets(ensembl)[listDatasets(ensembl)$dataset=="hsapiens_gene_ensembl",]$description
            message(paste0("Downloading genome information (try:", tries,") Using: ", description))

            filename <-  paste0(gsub("[[:punct:]]| ", "_",description),".rda")
            if(!file.exists(filename)) {
                chrom <- c(1:22, "X", "Y")
                gene.location <- getBM(attributes = attributes,
                                       filters = c("chromosome_name"),
                                       values = list(chrom), mart = ensembl)
                save(gene.location, file = filename)
            } else {
                message("Loading from disk")
                gene.location <- get(load(filename))
            }
            gene.location
        }, error = function(e) {
            msg <<- conditionMessage(e)
            tries <<- tries + 1L
        })
        if(!is.null(gene.location)) break
    }
    if (tries == 3L) stop("failed to get URL after 3 tries:", "\n  error: ", msg)

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

makeSEfromTranscriptomeProfiling <- function(data, cases, assay.list){
    # How many cases do we have?
    # We wil consider col 1 is the ensemble gene id, other ones are data
    size <- ncol(data)

    # Prepare Patient table
    colData <-  colDataPrepare(cases)

    gene.location <- get.GRCh.bioMart("hg38")
    aux <- strsplit(data$X1,"\\.")
    data$ensembl_gene_id <- as.character(unlist(lapply(aux,function(x) x[1])))

    found.genes <- table(data$ensembl_gene_id %in% gene.location$ensembl_gene_id)
    if("FALSE" %in% names(found.genes))
        message(paste0("From the ", nrow(data), " genes we couldn't map ", found.genes[["FALSE"]]))

    data <- merge(data, gene.location, by="ensembl_gene_id")

    # Prepare data table
    # Remove the version from the ensembl gene id
    assays <- list(data.matrix(data[,2:size+1]))
    names(assays) <- assay.list
    assays <- lapply(assays, function(x){
        colnames(x) <- NULL
        rownames(x) <- NULL
        return(x)
    })

    # Prepare rowRanges
    rowRanges <- GRanges(seqnames = paste0("chr", data$chromosome_name),
                         ranges = IRanges(start = data$start_position,
                                          end = data$end_position),
                         strand = data$strand,
                         ensembl_gene_id = data$ensembl_gene_id,
                         external_gene_name = data$external_gene_name,
                         original_ensembl_gene_id = data$X1)
    names(rowRanges) <- as.character(data$ensembl_gene_id)
    rse <- SummarizedExperiment(assays=assays,
                                rowRanges=rowRanges,
                                colData=colData)

    return(rse)
}

readTranscriptomeProfiling <- function(files, data.type, workflow.type, cases,summarizedExperiment) {
    if(grepl("Gene Expression Quantification", data.type, ignore.case = TRUE)){

        # Status working for:
        #  - htseq
        #  - FPKM
        #  - FPKM-UQ
        if(grepl("HTSeq",workflow.type)){
            pb <- txtProgressBar(min = 0, max = length(files), style = 3)
            for (i in seq_along(files)) {
                data <- read_tsv(file = files[i],
                                 col_names = FALSE,
                                 col_types = cols(
                                     X1 = col_character(),
                                     X2 = col_double()
                                 ))
                if(!missing(cases))  colnames(data)[2] <- cases[i]
                if(i == 1) df <- data
                if(i != 1) df <- merge(df, data, by=colnames(df)[1],all = TRUE)
                setTxtProgressBar(pb, i)
            }
            close(pb)
            if(summarizedExperiment) df <- makeSEfromTranscriptomeProfiling(df,cases,workflow.type)
        }
    } else if(grepl("miRNA", workflow.type, ignore.case = TRUE) & grepl("miRNA", data.type, ignore.case)) {
        pb <- txtProgressBar(min = 0, max = length(files), style = 3)
        for (i in seq_along(files)) {
            data <- read_tsv(file = files[i], col_names = TRUE,col_types = "cidc")
            if(!missing(cases))
                colnames(data)[2:ncol(data)] <- paste0(colnames(data)[2:ncol(data)],"_",cases[i])

            if(i == 1) df <- data
            if(i != 1) df <- merge(df, data, by=colnames(df)[1],all = TRUE)
            setTxtProgressBar(pb, i)
        }
        close(pb)

    } else if(grepl("Isoform Expression Quantification", data.type, ignore.case = TRUE)){
        pb <- txtProgressBar(min = 0, max = length(files), style = 3)
        for (i in seq_along(files)) {
            data <- read_tsv(file = files[i], col_names = TRUE, col_types = c("ccidcc"))
            if(!missing(cases)) data$barcode <- cases[i] else data$file <- i
            if(i == 1) df <- data
            if(i != 1) df <- rbind(df,data)
            setTxtProgressBar(pb, i)
        }
        close(pb)
    }
    return(df)
}

# Reads Copy Number Variation files to a data frame, basically it will rbind it
readCopyNumberVariation <- function(files, cases){
    message("Reading a copy  number variation")
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    for (i in seq_along(files)) {
        data <- read_tsv(file = files[i], col_names = TRUE, col_types = "ccnnnd")
        if(!missing(cases)) data$Sample <- cases[i]
        if(i == 1) df <- data
        if(i != 1) df <- rbind(df, data, make.row.names = FALSE)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    return(df)
}

getBarcodeInfo <- function(barcode) {
    baseURL <- "https://gdc-api.nci.nih.gov/cases/?"
    options.pretty <- "pretty=true"
    options.expand <- "expand=project,diagnoses,diagnoses.treatments,annotations,family_histories,demographic,exposures"
    option.size <- paste0("size=",length(barcode))
    #message(paste(barcode,collapse = '","'))
    #message(paste0('"',paste(barcode,collapse = '","')))
    options.filter <- paste0("filters=",
                             URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.submitter_id","value":['),
                             paste0('"',paste(barcode,collapse = '","')),
                             URLencode('"]}}]}'))
    #message(paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&")))
    url <- paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&"))
    json  <- tryCatch(
        fromJSON(url, simplifyDataFrame = TRUE),
        error = function(e) {
            fromJSON(content(GET(url), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
        }
    )

    results <- json$data$hits
    submitter_id <- results$submitter_id

    # We dont have the same cols for TCGA and TARGET so we need to check them
    if(!is.null(results$diagnoses)) {
        diagnoses <- rbindlist(results$diagnoses, fill = TRUE)
        if(any(grepl("submitter_id", colnames(diagnoses)))) {
            diagnoses$submitter_id <- gsub("_diagnosis","", diagnoses$submitter_id)
        }  else {
            diagnoses$submitter_id <- submitter_id
        }
        df <- diagnoses
    } else {
        df <- as.data.frame(submitter_id)

    }

    if(!is.null(results$exposures) > 0) {
        exposures <- rbindlist(results$exposures, fill = TRUE)
        if(any(grepl("submitter_id", colnames(exposures)))) {
            exposures$submitter_id <- gsub("_exposure","", exposures$submitter_id)
        }  else {
            exposures$submitter_id <- submitter_id
        }
        df <- merge(df,exposures, by="submitter_id", all = TRUE,sort = FALSE)
    }


    if(!is.null(results$demographic)) {
        demographic <- results$demographic
        if(any(grepl("submitter_id", colnames(demographic)))) {
            demographic$submitter_id <- gsub("_demographic","", results$demographic$submitter_id)
        } else {
            demographic$submitter_id <-submitter_id
        }
        df <- merge(df,demographic, by="submitter_id", all = TRUE,sort = FALSE)
    }

    treatments <- rbindlist(results$treatments,fill = TRUE)
    if(nrow(treatments) > 0) {
        df[,treatments:=NULL]

        if(any(grepl("submitter_id", colnames(treatments)))) {
            treatments$submitter_id <- gsub("_treatment","", treatments$submitter_id)
        } else {
            treatments$submitter_id <-submitter_id
        }
        df <- merge(df,treatments, by="submitter_id", all = TRUE,sort = FALSE)
    }
    df$bcr_patient_barcode <- df$submitter_id
    df <- cbind(df,results$project)

    # Adding in the same order
    df <- df[match(barcode,df$submitter_id)]
    # This line should not exists, but some patients does not have clinical data
    # case: TCGA-R8-A6YH"
    # this has been reported to GDC, waiting answers
    # So we will remove this NA cases
    df <- df[!is.na(df$submitter_id),]
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
            row.names(data) <- paste0("ID",values(data)$entrezgene)
            data <- assay(data)
        }

        if(all(grepl("\\|",rownames(data)))){
            message("2 - rownames  (gene|loci) => ('ID'loci) ")
            aux <- strsplit(rownames(data),"\\|")
            GeneID <- unlist(lapply(aux,function(x) x[2]))
            row.names(data) <- paste0("ID",GeneID)
        }
        data <- log2(data+1)
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

# Is this
# @import From gdata read.xls
#getPurityinfo <- function(){
#    message("Adding purity information from: doi:10.1038/ncomms9971")
#    x <- read.xls("http://www.nature.com/article-assets/npg/ncomms/2015/151204/ncomms9971/extref/ncomms9971-s2.xlsx", na.strings=c("NA","#DIV/0!","NaN"), skip = 2)
#    x <- x[,c(1,3:7)]
#    colnames(x)[1] <- "sample"
#    return(x)
#}


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
