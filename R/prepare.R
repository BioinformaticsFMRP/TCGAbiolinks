
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
#' @param add.gistic2.mut If a list of genes (gene symbol) is given, columns with gistic2 results from GDAC firehose (hg19)
#' and a column indicating if there is or not mutation in that gene (hg38)
#' (TRUE or FALSE - use the MAF file for more information)
#' will be added to the sample matrix in the summarized Experiment object.
#' @param mut.pipeline If add.gistic2.mut is not NULL this field will be taken in consideration.
#' Four separate variant calling pipelines are implemented for GDC data harmonization.
#' Options: muse, varscan2, somaticsniper, MuTect2. For more information:
#' https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/
#' @param mutant_variant_classification List of mutant_variant_classification that will be
#' consider a sample mutant or not. Default: "Frame_Shift_Del", "Frame_Shift_Ins",
#' "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "In_Frame_Del",
#' "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation"
#' @export
#' @examples
#' \dontrun{
#' query <- GDCquery(project = "TCGA-KIRP",
#'                   data.category = "Simple Nucleotide Variation",
#'                   data.type = "Masked Somatic Mutation",
#'                   workflow.type = "MuSE Variant Aggregation and Masking")
#' GDCdownload(query, method = "api", directory = "maf")
#' maf <- GDCprepare(query, directory = "maf")
#'
#' # Get GISTIC values
#' gistic.query <- GDCquery(project = "TCGA-ACC",
#'                          data.category = "Copy Number Variation",
#'                          data.type = "Gene Level Copy Number Scores",
#'                          access = "open")
#' GDCdownload(gistic.query)
#' gistic <- GDCprepare(gistic.query)
#' }
#' @return A summarizedExperiment or a data.frame
#' @importFrom  S4Vectors DataFrame
#' @importFrom SummarizedExperiment metadata<-
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
                       mut.pipeline = "mutect2",
                       mutant_variant_classification = c("Frame_Shift_Del",
                                                         "Frame_Shift_Ins",
                                                         "Missense_Mutation",
                                                         "Nonsense_Mutation",
                                                         "Splice_Site",
                                                         "In_Frame_Del",
                                                         "In_Frame_Ins",
                                                         "Translation_Start_Site",
                                                         "Nonstop_Mutation")){

  isServeOK()
  if(missing(query)) stop("Please set query parameter")

  test.duplicated.cases <- (any(duplicated(query$results[[1]]$cases)) &
                              !(query$data.type %in% c("Clinical data",
                                                       "Protein expression quantification",
                                                       "Raw intensities",
                                                       "Clinical Supplement",
                                                       "Biospecimen Supplement")))

  if(test.duplicated.cases) {
    dup <- query$results[[1]]$cases[duplicated(query$results[[1]]$cases)]
    cols <- c("tags","cases","experimental_strategy","analysis_workflow_type")
    cols <- cols[cols %in% colnames(query$results[[1]])]
    dup <- query$results[[1]][query$results[[1]]$cases %in% dup,cols]
    dup <- dup[order(dup$cases),]
    print(knitr::kable(dup))
    stop("There are samples duplicated. We will not be able to prepare it")
  }
  if(!save & remove.files.prepared) {
    stop("To remove the files, please set save to TRUE. Otherwise, the data will be lost")
  }
  # We save the files in project/source/data.category/data.type/file_id/file_name
  source <- ifelse(query$legacy,"legacy","harmonized")
  files <- file.path(query$results[[1]]$project, source,
                     gsub(" ","_",query$results[[1]]$data_category),
                     gsub(" ","_",query$results[[1]]$data_type),
                     gsub(" ","_",query$results[[1]]$file_id),
                     gsub(" ","_",query$results[[1]]$file_name))

  files <- file.path(directory, files)

  if(!all(file.exists(files))) stop(paste0("I couldn't find all the files from the query. ",
                                           "Please check if the directory parameter is right or `GDCdownload` downloaded the samples."))

  cases <- ifelse(grepl("TCGA|TARGET",query$results[[1]]$project %>% unlist()),query$results[[1]]$cases,query$results[[1]]$sample.submitter_id)
  if(grepl("Transcriptome Profiling", query$data.category, ignore.case = TRUE)){
    data <- readTranscriptomeProfiling(files = files,
                                       data.type = ifelse(!is.na(query$data.type),  as.character(query$data.type),  unique(query$results[[1]]$data_type)),
                                       workflow.type = unique(query$results[[1]]$analysis_workflow_type),
                                       cases = cases,
                                       summarizedExperiment)
  } else if(grepl("Copy Number Variation",query$data.category,ignore.case = TRUE)) {
    if(unique(query$results[[1]]$data_type) == "Gene Level Copy Number Scores") {
      data <- readGISTIC(files, query$results[[1]]$cases)
    } else {
      data <- readCopyNumberVariation(files, query$results[[1]]$cases)
    }
  }  else if(grepl("DNA methylation",query$data.category, ignore.case = TRUE)) {
    data <- readDNAmethylation(files, cases = cases, summarizedExperiment, unique(query$platform))
  }  else if(grepl("Raw intensities",query$data.type, ignore.case = TRUE)) {
    # preparing IDAT files
    data <- readIDATDNAmethylation(files, barcode = cases, summarizedExperiment, unique(query$platform), query$legacy)
  }  else if(grepl("Protein expression",query$data.category,ignore.case = TRUE)) {
    data <- readProteinExpression(files, cases = cases)
  }  else if(grepl("Simple Nucleotide Variation",query$data.category,ignore.case = TRUE)) {
    if(grepl("Masked Somatic Mutation",query$results[[1]]$data_type,ignore.case = TRUE) | source == "legacy")
      suppressWarnings(data <- readSimpleNucleotideVariationMaf(files))
  }  else if(grepl("Clinical|Biospecimen", query$data.category, ignore.case = TRUE)){
    data <- readClinical(files, query$data.type, cases = cases)
    summarizedExperiment <- FALSE
  } else if (grepl("Gene expression",query$data.category,ignore.case = TRUE)) {
    if(query$data.type == "Gene expression quantification")
      data <- readGeneExpressionQuantification(files = files,
                                               cases = cases,
                                               summarizedExperiment = summarizedExperiment,
                                               genome = ifelse(query$legacy,"hg19","hg38"),
                                               experimental.strategy = unique(query$results[[1]]$experimental_strategy))

    if(query$data.type == "miRNA gene quantification")
      data <- readGeneExpressionQuantification(files = files,
                                               cases = cases,
                                               summarizedExperiment = FALSE,
                                               genome = ifelse(query$legacy,"hg19","hg38"),
                                               experimental.strategy = unique(query$results[[1]]$experimental_strategy))
    if(query$data.type == "miRNA isoform quantification")
      data <- readmiRNAIsoformQuantification(files = files,
                                             cases = query$results[[1]]$cases)

    if(query$data.type == "Isoform expression quantification")
      data <- readIsoformExpressionQuantification(files = files, cases = cases)

    if(query$data.type == "Exon quantification")
      data <- readExonQuantification(files = files,
                                     cases = cases,
                                     summarizedExperiment = summarizedExperiment)

  }
  # Add data release to object
  if(summarizedExperiment & !is.data.frame(data)){
    metadata(data) <- list("data_release" = getGDCInfo()$data_release)
  }


  if((!is.null(add.gistic2.mut)) & summarizedExperiment) {
    message("=> Adding GISTIC2 and mutation information....")
    genes <- tolower(levels(EAGenes$Gene))
    if(!all(tolower(add.gistic2.mut) %in% genes)) message(paste("These genes were not found:\n",
                                                                paste(add.gistic2.mut[! tolower(add.gistic2.mut) %in% genes],collapse = "\n=> ")))
    add.gistic2.mut <- add.gistic2.mut[tolower(add.gistic2.mut) %in% tolower(genes)]
    if(length(add.gistic2.mut) > 0){
      info <- colData(data)
      for(i in unlist(query$project)){
        info <- get.mut.gistc.information(info,
                                          i,
                                          add.gistic2.mut,
                                          mut.pipeline = mut.pipeline,
                                          mutant_variant_classification = mutant_variant_classification)
      }
      colData(data) <- info
    }
  }
  if("samples" %in% colnames(data)){
    if(any(duplicated(data$sample))) {
      message("Replicates found.")
      if(any(data$is_ffpe)) message("FFPE should be removed. You can modify the data with the following command:\ndata <- data[,!data$is_ffpe]")
      print(as.data.frame(colData(data)[data$sample %in% data$sample[duplicated(data$sample)],c("is_ffpe"),drop=F]))
    }
  }


  if(save){
    if(missing(save.filename) & !missing(query)) save.filename <- paste0(query$project,gsub(" ","_", query$data.category),gsub(" ","_",date()),".RData")
    message(paste0("=> Saving file: ",save.filename))
    save(data, file = save.filename)
    message("=> File saved")

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


readClinical <- function(files, data.type, cases){
  if(data.type == "Clinical data"){
    suppressMessages({
      ret <- plyr::alply(files,.margins = 1,readr::read_tsv, .progress = "text")
    })
    names(ret) <- gsub("nationwidechildrens.org_","",gsub(".txt","",basename(files)))
  } else if(data.type %in% c("Clinical Supplement","Biospecimen Supplement")){
    ret <- plyr::alply(files,.margins = 1,function(f) {
      readr::read_tsv(f,col_types = readr::cols())
    }, .progress = "text")
    names(ret) <- gsub("nationwidechildrens.org_","",gsub(".txt","",basename(files)))
  }
  return(ret)
}


#' @importFrom tidyr separate
readExonQuantification <- function (files, cases, summarizedExperiment = TRUE){
  pb <- txtProgressBar(min = 0, max = length(files), style = 3)
  assay.list <- NULL

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
  df <- df %>% separate(exon,into = c("seqnames","coordinates","strand"),sep = ":") %>%
    separate(coordinates,into = c("start","end"),sep = "-")
  if(summarizedExperiment) {
    suppressWarnings({
      assays <- lapply(assay.list, function (x) {
        return(data.matrix(subset(df, select = grep(x,colnames(df),ignore.case = TRUE))))
      })
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
    rowRanges <- makeGRangesFromDataFrame(df)
    rse <- SummarizedExperiment(assays=assays,
                                rowRanges=rowRanges,
                                colData=colData)
    return(rse)
  }
  return(df)

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


#' @importFrom purrrogress with_progress
#' @importFrom purrr reduce
readGeneExpressionQuantification <- function(files,
                                             cases,
                                             genome = "hg19",
                                             summarizedExperiment = TRUE,
                                             experimental.strategy,
                                             platform){
  skip <- unique((ifelse(experimental.strategy == "Gene expression array",1,0)))

  if(length(skip) > 1) stop("It is not possible to handle those different platforms together")

  print.header(paste0("Reading ", length(files)," files"),"subsection")
  ret <- plyr::alply(seq_along(files),1,.fun = function(i,cases){
    data <- fread(files[i],
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  skip = skip)
    if(!missing(cases)) {
      assay.list <<- gsub(" |\\(|\\)|\\/","_",colnames(data)[2:ncol(data)])
      # We will use this because there might be more than one col for each samples
      setnames(data,colnames(data)[2:ncol(data)],
               paste0(gsub(" |\\(|\\)|\\/","_",colnames(data)[2:ncol(data)]),"_",cases[i]))
    }
    data
  },.progress = "time",cases = cases)

  print.header(paste0("Merging ", length(files)," files"),"subsection")
  merging.col <- colnames(ret[[1]])[1]
  df <- purrr::reduce(ret,
                      purrrogress::with_progress(dplyr::full_join,
                                                 length(ret),
                                                 title = "Merging files"),
                      by = merging.col)

  if (summarizedExperiment) {
    df <- makeSEfromGeneExpressionQuantification(df,assay.list, genome = genome)
  } else {
    rownames(df) <- df$gene_id
    df$gene_id <- NULL
  }
  return(df)
}
makeSEfromGeneExpressionQuantification <- function(df, assay.list, genome = "hg19"){
  gene.location <- get.GRCh.bioMart(genome)
  if(all(grepl("\\|",df[,1]))){
    aux <- strsplit(df$gene_id,"\\|")
    GeneID <- unlist(lapply(aux,function(x) x[2]))
    df$entrezgene_id <- as.numeric(GeneID)
    gene.location <- gene.location[!duplicated(gene.location$entrezgene_id),]
    df <- merge(df, gene.location, by = "entrezgene_id")
  } else {
    df$external_gene_name <- as.character(df[,1])
    df <- merge(df, gene.location, by = "external_gene_name")
  }


  if("transcript_id" %in% assay.list){
    rowRanges <- GRanges(seqnames = paste0("chr", df$chromosome_name),
                         ranges = IRanges(start = df$start_position,
                                          end = df$end_position),
                         strand = df$strand,
                         gene_id = df$external_gene_name,
                         entrezgene = df$entrezgene_id,
                         ensembl_gene_id = df$ensembl_gene_id,
                         transcript_id = subset(df, select = 5))
    names(rowRanges) <- as.character(df$gene_id)
    assay.list <- assay.list[which(assay.list != "transcript_id")]
  } else {
    rowRanges <- GRanges(seqnames = paste0("chr", df$chromosome_name),
                         ranges = IRanges(start = df$start_position,
                                          end = df$end_position),
                         strand = df$strand,
                         gene_id = df$external_gene_name,
                         entrezgene = df$entrezgene_id,
                         ensembl_gene_id = df$ensembl_gene_id)
    names(rowRanges) <- as.character(df$external_gene_name)
  }
  suppressWarnings({
    assays <- lapply(assay.list, function (x) {
      return(data.matrix(subset(df, select = grep(x,colnames(df),ignore.case = TRUE))))
    })
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


#' @importFrom downloader download
#' @importFrom S4Vectors DataFrame
makeSEFromDNAMethylationMatrix <- function(betas, genome = "hg38", met.platform = "450K") {
  message("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  message("Creating a SummarizedExperiment from DNA methylation input")

  # Instead of looking on the size, it is better to set it as a argument as the annotation is different
  annotation <-   getInfiniumAnnotation(met.platform, genome)

  rowRanges <- annotation[names(annotation) %in% rownames(betas),,drop=FALSE]

  colData <-  DataFrame(samples = colnames(betas))
  betas <- betas[rownames(betas) %in% names(rowRanges),,drop = FALSE]
  betas <- betas[names(rowRanges),,drop = FALSE]
  assay <- data.matrix(betas)
  betas <- SummarizedExperiment(assays=assay,
                                rowRanges=rowRanges,
                                colData=colData)
  return(betas)
}


getInfiniumAnnotation <- function(platform, genome){
  base <- "http://zwdzwd.io/InfiniumAnnotation/current/"
  path <- file.path(base,platform,paste(platform,"hg19.manifest.rds", sep ="."))
  if (grepl("hg38", genome)) path <- gsub("hg19","hg38",path)
  if(platform == "EPIC") {
    annotation <- paste0(base,"EPIC/EPIC.hg19.manifest.rds")
  } else if(platform == "450K") {
    annotation <- paste0(base,"hm450/hm450.hg19.manifest.rds")
  } else {
    annotation <- paste0(base,"hm27/hm27.hg19.manifest.rds")
  }
  if(grepl("hg38", genome)) annotation <- gsub("hg19","hg38",annotation)
  if(!file.exists(basename(annotation))) {
    if(Sys.info()["sysname"] == "Windows") mode <- "wb" else  mode <- "w"
    downloader::download(annotation, basename(annotation), mode = mode)
  }
  readRDS(basename(annotation))
}


makeSEfromDNAmethylation <- function(df, probeInfo=NULL){
  if(is.null(probeInfo)) {
    rowRanges <- GRanges(seqnames = paste0("chr", df$Chromosome),
                         ranges = IRanges(start = df$Genomic_Coordinate,
                                          end = df$Genomic_Coordinate),
                         probeID = df$Composite.Element.REF,
                         Gene_Symbol = df$Gene_Symbol)

    names(rowRanges) <- as.character(df$Composite.Element.REF)
    colData <-  colDataPrepare(colnames(df)[5:ncol(df)])
    assay <- data.matrix(subset(df,select = c(5:ncol(df))))
  } else {
    rowRanges <- makeGRangesFromDataFrame(probeInfo, keep.extra.columns = TRUE)
    colData <-  colDataPrepare(colnames(df)[(ncol(probeInfo) + 1):ncol(df)])
    assay <- data.matrix(subset(df,select = c((ncol(probeInfo) + 1):ncol(df))))
  }
  colnames(assay) <- rownames(colData)
  rownames(assay) <- as.character(df$Composite.Element.REF)
  rse <- SummarizedExperiment(assays = assay, rowRanges = rowRanges, colData = colData)
}

readIDATDNAmethylation <- function(files,
                                   barcode,
                                   summarizedExperiment,
                                   platform,
                                   legacy){

  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("sesame package is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Check if moved files would be moved outside of scope folder, if so, path doesn't change
  moved.files <- sapply(files,USE.NAMES=FALSE,function(x){
    if(grepl("Raw_intensities",dirname(dirname(x)))){
      return(file.path(dirname(dirname(x)), basename(x)))
    }
    return(x)
  })

  # for each file move it to upper parent folder if necessary
  plyr::a_ply(files, 1,function(x){
    if(grepl("Raw_intensities",dirname(dirname(x)))){
      tryCatch(move(x,file.path(dirname(dirname(x)), basename(x)),keep.copy = FALSE),error = function(e){})
    }
  })

  samples <- unique(gsub("_Grn.idat|_Red.idat","",moved.files))
  message("Processing  IDATs with Sesame - http://bioconductor.org/packages/sesame/")
  message("Running opensesame - applying quality masking and nondetection masking (threshold P-value 0.05)")
  message("Please cite: doi: 10.1093/nar/gky691 and 10.1093/nar/gkt090")
  betas <- sesame::openSesame(samples)
  barcode <- unique(data.frame("file" = gsub("_Grn.idat|_Red.idat","",basename(moved.files)), "barcode" = barcode))
  colnames(betas) <- barcode$barcode[match(basename(samples),barcode$file)]

  if(summarizedExperiment){
    met.platform <- "EPIC"
    if(grepl("450",platform)) met.platform <- "450K"
    if(grepl("27",platform)) met.platform <- "27K"
    betas <- makeSEFromDNAMethylationMatrix(betas,genome = ifelse(legacy,"hg19","hg38"),met.platform = met.platform)
    colData(betas) <- DataFrame(colDataPrepare(colnames(betas)))
  }
  return(betas)

}

# We will try to make this function easier to use this function in its own data
# In case it is not TCGA I should not consider that there is a barcode in the header
# Instead the user should be able to add the names to his data
# The only problem is that the data from the user will not have all the columns
# TODO: Improve this function to be more generic as possible
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tibble as_data_frame
readDNAmethylation <- function(files, cases, summarizedExperiment = TRUE, platform){
  if(missing(cases)) cases <- NULL
  if (grepl("OMA00",platform)){
    pb <- txtProgressBar(min = 0, max = length(files), style = 3)
    for (i in seq_along(files)) {
      data <- fread(files[i], header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE,skip = 1,
                    na.strings="N/A",
                    colClasses=c("character", # Composite Element REF
                                 "numeric"))   # beta value
      setnames(data,gsub(" ", "\\.", colnames(data)))
      if(!is.null(cases)) setnames(data,2,cases[i])
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
    skip <- ifelse(all(grepl("hg38",files)), 0,1)
    colClasses <- NULL
    if(!all(grepl("hg38",files))) colClasses <- c("character", # Composite Element REF
                                                  "numeric",   # beta value
                                                  "character", # Gene symbol
                                                  "character", # Chromosome
                                                  "integer")


    x <- plyr::alply(files,1, function(f) {
      data <- fread(f, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE,skip = skip, colClasses = colClasses)
      setnames(data,gsub(" ", "\\.", colnames(data)))
      if(!is.null(cases)) setnames(data,2,cases[which(f == files)])
      setcolorder(data,c(1, 3:ncol(data), 2))
    }, .progress = "time")


    print.header(paste0("Merging ", length(files)," files"),"subsection")
    df <- purrr::reduce(x,
                        purrrogress::with_progress(dplyr::left_join,
                                                   length(x)))

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

# Barcode example MMRF_1358_1_BM_CD138pos_T1_TSMRU_L02337
colDataPrepareMMRF <- function(barcode){
  DataFrame(barcode = barcode,
            sample = barcode,
            patient = substr(barcode,1,9)
  )
}

colDataPrepareTARGET <- function(barcode){
  message("Adding description to TARGET samples")
  tissue.code <- c('01','02','03','04','05','06','07','08','09','10','11',
                   '12','13','14','15','16','17','20','40','41','42','50','60','61','99')

  definition <- c("Primary solid Tumor", # 01
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
  aux <- DataFrame(tissue.code = tissue.code,definition)

  # in case multiple equal barcode
  regex <- paste0("[:alnum:]{5}-[:alnum:]{2}-[:alnum:]{6}",
                  "-[:alnum:]{3}-[:alnum:]{3}")
  samples <- str_match(barcode,regex)[,1]

  ret <- DataFrame(barcode = barcode,
                   sample = substr(barcode, 1, 20),
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

#' @title Create samples information matrix for GDC samples
#' @description Create samples information matrix for GDC samples add subtype information
#' @param barcode TCGA or TARGET barcode
#' @importFrom plyr rbind.fill
#' @importFrom dplyr left_join
#' @examples
#'  metadata <- colDataPrepare(c("TCGA-OR-A5K3-01A","C3N-00321-01"))
#'  metadata <- colDataPrepare(c("BLGSP-71-06-00157-01A",
#'                               "BLGSP-71-22-00332-01A"))
#' @export
colDataPrepare <- function(barcode){
  # For the moment this will work only for TCGA Data
  # We should search what TARGET data means
  message("Starting to add information to samples")
  ret <- NULL

  if(all(grepl("TARGET",barcode))) ret <- colDataPrepareTARGET(barcode)
  if(all(grepl("TCGA",barcode))) ret <- colDataPrepareTCGA(barcode)
  if(all(grepl("MMRF",barcode))) ret <- colDataPrepareMMRF(barcode)
  if(is.null(ret)) ret <- data.frame(sample = barcode,stringsAsFactors = FALSE)

  message(" => Add clinical information to samples")
  # There is a limitation on the size of the string, so this step will be splited in cases of 100
  patient.info <- NULL

  patient.info <- splitAPICall(FUN = getBarcodeInfo,
                               step = 20,
                               items = ret$sample)

  if(!is.null(patient.info)) {
    ret$sample_submitter_id <- ret$sample %>% as.character()
    ret <- left_join(ret %>% as.data.frame, patient.info, by = "sample_submitter_id")
  }
  ret$bcr_patient_barcode <- ret$sample %>% as.character()
  ret$sample_submitter_id <- ret$sample %>% as.character()

  if(!"project_id" %in% colnames(ret)) {
    if("disease_type" %in% colnames(ret)){
      aux <- getGDCprojects()[,c(5,7)]
      aux <- aux[aux$disease_type == unique(ret$disease_type),2]
      ret$project_id <- as.character(aux)
    }
  }
  # There is no subtype info for target, return as it is
  if(any(grepl("TCGA",barcode))) {
    ret <- addSubtypeInfo(ret)
  }

  # na.omit should not be here, exceptional case
  if(is.null(ret)) return(data.frame(row.names = barcode, barcode,stringsAsFactors = FALSE))

  # Add purity information from http://www.nature.com/articles/ncomms9971
  # purity  <- getPurityinfo()
  # ret <- merge(ret, purity, by = "sample", all.x = TRUE, sort = FALSE)

  # Put data in the right order
  idx <- sapply(substr(barcode,1,min(str_length(ret$bcr_patient_barcode))), function(x) {
    grep(x,ret$bcr_patient_barcode)
  })

  ret <- ret[idx,]
  rownames(ret) <- barcode
  return(ret)
}

#' @title Get hg19 or hg38 information from biomaRt
#' @description Get hg19 or hg38 information from biomaRt
#' @param genome hg38 or hg19
#' @param as.granges Output as GRanges or data.frame
#' @importFrom biomaRt getBM useMart listDatasets useEnsembl
#' @export
get.GRCh.bioMart <- function(genome = "hg19", as.granges = FALSE) {
  tries <- 0L
  msg <- character()
  while (tries < 3L) {
    gene.location <- tryCatch({
      host <- ifelse(genome == "hg19", "grch37.ensembl.org",
                     "www.ensembl.org")
      mirror <- list(NULL, "useast", "uswest", "asia")[[tries + 1]]

      ensembl <- tryCatch({
        message(ifelse(is.null(mirror),
                       paste0("Accessing ", host, " to get gene information"),
                       paste0("Accessing ", host," (mirror ", mirror,")")))
        useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host = host, mirror = mirror)
      }, error = function(e) {
        message(e)
        return(NULL)
      })

      attributes <- c("chromosome_name",
                      "start_position",
                      "end_position", "strand",
                      "ensembl_gene_id",
                      "entrezgene_id",
                      "external_gene_name")

      db.datasets <- listDatasets(ensembl)
      description <- db.datasets[db.datasets$dataset == "hsapiens_gene_ensembl",]$description
      message(paste0("Downloading genome information (try:", tries,") Using: ", description))

      chrom <- c(1:22, "X", "Y")
      gene.location <- getBM(attributes = attributes,
                             filters = c("chromosome_name"),
                             values = list(chrom), mart = ensembl)
      gene.location
    }, error = function(e) {
      msg <<- conditionMessage(e)
      tries <<- tries + 1L
      NULL
    })
    if(!is.null(gene.location)) break
  }
  if (tries == 3L) stop("failed to get URL after 3 tries:", "\n  error: ", msg)

  if(as.granges) {
    gene.location$strand[gene.location$strand == 1] <- "+"
    gene.location$strand[gene.location$strand == -1] <- "-"
    gene.location$chromosome_name <- paste0("chr",gene.location$chromosome_name)
    gene.location <- makeGRangesFromDataFrame(gene.location, seqnames.field = "chromosome_name",
                                              start.field = "start_position",
                                              end.field = "end_position",
                                              keep.extra.columns = TRUE) # considering the whole gene no their promoters
  }
  return(gene.location)
}

addSubtypeInfo <- function(ret){
  out <- NULL
  message(" => Adding TCGA molecular information from marker papers")
  message(" => Information will have prefix 'paper_' ")

  for(proj in na.omit(unique(ret$project_id))){
    if(grepl("TCGA",proj,ignore.case = TRUE)) {
      # remove letter from 01A 01B etc
      ret$sample.aux <- substr(ret$sample,1,15)

      tumor <- gsub("TCGA-","",proj)
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
                     "UVM")
      if (grepl(paste(c(available,"all"),collapse = "|"),tumor,ignore.case = TRUE)) {
        subtype <- TCGAquery_subtype(tumor)
        colnames(subtype) <- paste0("paper_", colnames(subtype))

        if(all(str_length(subtype$paper_patient) == 12)){
          # Subtype information were to primary tumor in priority
          subtype$sample.aux <- paste0(subtype$paper_patient,"-01")
        }
        ret.aux <- ret[ret$sample.aux %in% subtype$sample.aux,]
        ret.aux <- merge(ret.aux,subtype, by = "sample.aux", all.x = TRUE)
        out <- rbind.fill(as.data.frame(out),as.data.frame(ret.aux))
      }
    }
  }
  if(is.null(out)) return(ret)

  # We need to put together the samples with subtypes with samples without subytpes
  ret.aux <- ret[!ret$sample %in% out$sample,]
  ret <- rbind.fill(as.data.frame(out),as.data.frame(ret.aux))
  ret$sample.aux <- NULL

  return(ret)
}

readProteinExpression <- function(files,cases) {
  pb <- txtProgressBar(min = 0, max = length(files), style = 3)
  for (i in seq_along(files)) {
    data <- read_tsv(file = files[i], col_names = TRUE,skip = 1, col_types = c("cn"))
    if(!missing(cases))  colnames(data)[2] <- cases[i]
    if(i == 1) df <- data
    if(i != 1) df <- merge(df, data,all=TRUE, by="Composite Element REF")
    setTxtProgressBar(pb, i)
  }
  close(pb)

  return(df)
}


makeSEfromTranscriptomeProfilingSTAR <- function(data, cases, assay.list){
  # How many cases do we have?
  # We wil consider col 1 is the ensemble gene id, other ones are data
  size <- ncol(data)
  # Prepare Patient table
  colData <-  colDataPrepare(cases)

  # one ensemblID can be mapped to multiple entrezgene ID
  gene.location <- get.GRCh.bioMart("hg38")
  gene.location <- gene.location[!duplicated(gene.location$ensembl_gene_id),]

  data$ensembl_gene_id <-  as.character(gsub("\\.[0-9]*","",data$`#gene`))
  metrics <- subset(data, !grepl("ENSG", data$ensembl_gene_id))
  data <- subset(data, grepl("ENSG", data$ensembl_gene_id))
  found.genes <- table(data$ensembl_gene_id %in% gene.location$ensembl_gene_id)

  if("FALSE" %in% names(found.genes))
    message(paste0("From the ", nrow(data), " genes we couldn't map ", found.genes[["FALSE"]]))

  data <- merge(data, gene.location, by = "ensembl_gene_id")

  # Prepare data table
  # Remove the version from the ensembl gene id
  assays <- list(data.matrix(data[,grep("unstranded",colnames(data))]),
                 data.matrix(data[,grep("stranded_first",colnames(data))]),
                 data.matrix(data[,grep("stranded_second",colnames(data))]))
  names(assays) <- c("unstranded","stranded_first","stranded_second")
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
                       original_ensembl_gene_id = data$`#gene`)
  names(rowRanges) <- as.character(data$ensembl_gene_id)
  rse <- SummarizedExperiment(assays = assays,
                              rowRanges = rowRanges,
                              colData = colData)
  metadata(rse) <- metrics
  return(rse)
}


makeSEfromTranscriptomeProfiling <- function(data, cases, assay.list){
  # How many cases do we have?
  # We wil consider col 1 is the ensemble gene id, other ones are data
  size <- ncol(data)
  # Prepare Patient table
  colData <-  colDataPrepare(cases)

  # one ensemblID can be mapped to multiple entrezgene ID
  gene.location <- get.GRCh.bioMart("hg38")
  gene.location <- gene.location[!duplicated(gene.location$ensembl_gene_id),]

  data$ensembl_gene_id <- as.character(gsub("\\.[0-9]*","",data$X1))
  data <- subset(data, grepl("ENSG", data$ensembl_gene_id))
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
  rse <- SummarizedExperiment(assays = assays,
                              rowRanges = rowRanges,
                              colData = colData)

  return(rse)
}


#' @importFrom purrr reduce
#' @importFrom dplyr left_join
#' @importFrom plyr alply
readTranscriptomeProfiling <- function(files, data.type, workflow.type, cases,summarizedExperiment) {

  if(grepl("Gene Expression Quantification", data.type, ignore.case = TRUE)){
    # Status working for:
    #  - htseq
    #  - FPKM
    #  - FPKM-UQ
    if(grepl("HTSeq",workflow.type)){

      x <- plyr::alply(files,1, function(f) {
        readr::read_tsv(file = f,
                        col_names = FALSE,
                        progress = FALSE,
                        col_types = c("cd"))
      }, .progress = "time")
      df <- x %>% purrr::reduce(left_join, by = "X1")
      if(!missing(cases))  colnames(df)[-1] <- cases
      if(summarizedExperiment) df <- makeSEfromTranscriptomeProfiling(df,cases,workflow.type)
    } else  if(grepl("STAR",workflow.type)){
      x <- plyr::alply(files,1, function(f) {
        readr::read_tsv(file = f,
                        col_names = TRUE,
                        progress = FALSE)
      }, .progress = "time")

      df <- x %>% purrr::reduce(left_join, by = "#gene")
      if(!missing(cases))  colnames(df)[-1] <- sapply(cases, function(x){stringr::str_c(c("unstranded_","stranded_first_", "stranded_second_"),x)}) %>% as.character()
      if(summarizedExperiment) df <- makeSEfromTranscriptomeProfilingSTAR(df,cases,workflow.type)
    }
  } else if(grepl("miRNA", workflow.type, ignore.case = TRUE) & grepl("miRNA", data.type, ignore.case = TRUE)) {
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

#' @importFrom purrr reduce
readGISTIC <- function(files, cases){
  message("Reading GISTIC file")
  gistic.df <- NULL
  gistic.list <- plyr::alply(files,1,.fun = function(file) {
    message("Reading file: ", file)
    data <- read_tsv(file = file, col_names = TRUE, progress = TRUE,col_types = readr::cols())

    aliquot <- colnames(data)[-c(1:3)]
    info <- splitAPICall(FUN = getBarcodefromAliquot,
                         step = 20,
                         items = aliquot)

    barcode <- as.character(info$submitter_id)[match(aliquot,as.character(info$aliquot_id))]
    colnames(data)[-c(1:3)] <- barcode
    return(data)
  })
  gistic.df <- gistic.list %>% purrr::reduce(dplyr::full_join, by = c("Gene Symbol","Gene ID","Cytoband"))

  return(gistic.df)
}

# Reads Copy Number Variation files to a data frame, basically it will rbind it
readCopyNumberVariation <- function(files, cases){
  message("Reading copy number variation files")
  pb <- txtProgressBar(min = 0, max = length(files), style = 3)
  for (i in seq_along(files)) {
    data <- read_tsv(file = files[i], col_names = TRUE, col_types = "ccnnnd")
    if(!missing(cases)) data$Sample <- cases[i]
    if(i == 1) df <- data
    if(i != 1) df <- rbind(df, data)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(df)
}

# getBarcodeInfo(c("TCGA-A6-6650-01B"))
addFFPE <- function(df) {
  message("Add FFPE information. More information at: \n=> https://cancergenome.nih.gov/cancersselected/biospeccriteria \n=> http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html")
  barcode <- df$barcode
  patient <- df$patient

  ffpe.info <- NULL
  ffpe.info <- splitAPICall(FUN = getFFPE,
                            step = 20,
                            items = df$patient)

  df <- merge(df, ffpe.info,by.x = "sample", by.y = "submitter_id")
  df <- df[match(barcode,df$barcode),]
  return(df)
}

# getFFPE("TCGA-A6-6650")
#' @importFrom plyr rbind.fill
#' @importFrom httr content
getFFPE <- function(patient){
  baseURL <- "https://api.gdc.cancer.gov/cases/?"
  options.pretty <- "pretty=true"
  options.expand <- "expand=samples"
  option.size <- paste0("size=",length(patient))
  options.filter <- paste0("filters=",
                           URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.submitter_id","value":['),
                           paste0('"',paste(patient,collapse = '","')),
                           URLencode('"]}}]}'))
  url <- paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&"))
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
  results <- rbind.fill(results$samples)[,c("submitter_id","is_ffpe")]
  return(results)
}

getAliquot_ids <- function(barcode){
  baseURL <- "https://api.gdc.cancer.gov/cases/?"
  options.fields <- "fields=samples.portions.analytes.aliquots.aliquot_id,samples.portions.analytes.aliquots.submitter_id"
  options.pretty <- "pretty=true"
  option.size <- paste0("size=",length(barcode))
  #message(paste(barcode,collapse = '","'))
  #message(paste0('"',paste(barcode,collapse = '","')))
  options.filter <- paste0("filters=",
                           URLencode('{"op":"and","content":[{"op":"in","content":{"field":"cases.submitter_id","value":['),
                           paste0('"',paste(barcode,collapse = '","')),
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
  results <- unlist(json$data$hits$samples)
  results.barcode <- grep("TCGA",results,value = TRUE)
  results.aliquot <- grep("TCGA",results,value = TRUE,invert = TRUE)

  df <- data.frame(results.aliquot,results.barcode)
  colnames(df) <- c("aliquot_id","barcode")
  return(df)
}


# getBarcodeInfo(c("TCGA-OR-A5K3-01A","C3N-00321-01"))
# barcode is: sample_submitter_id
#' @importFrom dplyr bind_cols
getBarcodeInfo <- function(barcode) {
  baseURL <- "https://api.gdc.cancer.gov/cases/?"
  options.pretty <- "pretty=true"
  options.expand <- "expand=samples,project,diagnoses,diagnoses.treatments,annotations,family_histories,demographic,exposures"
  option.size <- paste0("size=",length(barcode))
  options.filter <- paste0("filters=",
                           URLencode('{"op":"or","content":[{"op":"in","content":{"field":"cases.submitter_id","value":['),
                           paste0('"',paste(barcode,collapse = '","')),
                           URLencode('"]}},'),
                           URLencode('{"op":"in","content":{"field":"submitter_sample_ids","value":['),
                           paste0('"',paste(barcode,collapse = '","')),
                           URLencode('"]}},'),
                           URLencode('{"op":"in","content":{"field":"submitter_aliquot_ids","value":['),
                           paste0('"',paste(barcode,collapse = '","')),
                           URLencode('"]}}'),
                           URLencode(']}')
  )
  url <- paste0(baseURL,paste(options.pretty,options.expand, option.size, options.filter, sep = "&"))
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

  # no results
  if(length(results) == 0){
    return(data.frame(barcode,stringsAsFactors = FALSE))
  }

  submitter_id <- results$submitter_id
  submitter_aliquot_ids <- results$submitter_aliquot_ids

  if(!is.null(results$samples)) {
    samples <- rbindlist(results$samples, fill = TRUE)
    samples <- samples[match(barcode,samples$submitter_id),]
    samples$sample_submitter_id <- str_extract_all(samples$submitter_id,paste(barcode,collapse = "|")) %>%
      unlist %>% as.character

    tryCatch({
      samples$submitter_id <-
        str_extract_all(samples$submitter_id,
                        paste(submitter_id, collapse = "|"),
                        simplify = TRUE) %>% as.character
    }, error = function(e){
      samples$submitter_id <- submitter_id
    })
    df <- samples[!is.na(samples$submitter_id),]
    suppressWarnings({
      df[,c("updated_datetime","created_datetime")] <- NULL
    })
  }


  # We dont have the same cols for TCGA and TARGET so we need to check them
  if(!is.null(results$diagnoses)) {
    diagnoses <- rbindlist(lapply(results$diagnoses, function(x) if(is.null(x)) data.frame(NA) else x),fill = TRUE)
    diagnoses[,c("updated_datetime","created_datetime","state")] <- NULL
    if(any(grepl("submitter_id", colnames(diagnoses)))) {
      diagnoses$submitter_id <- gsub("_diagnosis.*|-DIAG|diag-","", diagnoses$submitter_id)
    }  else {
      diagnoses$submitter_id <- submitter_id
    }
    # this is required since the sample might not have a diagnosis
    if(!any(df$submitter_id %in% diagnoses$submitter_id)){
      diagnoses$submitter_id <- NULL
      df <- dplyr::bind_cols(df,diagnoses)
    } else {
      df <- left_join(df, diagnoses, by = "submitter_id")
    }
  }
  if(!is.null(results$exposures)) {
    exposures <- rbindlist(lapply(results$exposures, function(x) if(is.null(x)) data.frame(NA) else x),fill = TRUE)
    exposures[,c("updated_datetime","created_datetime","state")] <- NULL
    if(any(grepl("submitter_id", colnames(exposures)))) {
      exposures$submitter_id <- gsub("_exposure.*|-EXP","", exposures$submitter_id)
    }  else {
      exposures$submitter_id <- submitter_id
    }
    if(!any(df$submitter_id %in% exposures$submitter_id)){
      exposures$submitter_id <- NULL
      df <- dplyr::bind_cols(df,exposures)
    } else {
      df <- left_join(df, exposures, by = "submitter_id")
    }
  }


  if(!is.null(results$demographic)) {
    demographic <- results$demographic
    demographic[,c("updated_datetime","created_datetime","state")] <- NULL
    if(any(grepl("submitter_id", colnames(demographic)))) {
      demographic$submitter_id <- gsub("_demographic.*|-DEMO|demo-","", results$demographic$submitter_id)
    } else {
      demographic$submitter_id <- submitter_id
    }

    if(!any(df$submitter_id %in% demographic$submitter_id)){
      demographic$submitter_id <- NULL
      df <- dplyr::bind_cols(df,demographic)
    } else {
      df <- left_join(df,demographic, by = "submitter_id")
    }
  }

  df$bcr_patient_barcode <- df$submitter_id %>% as.character()
  projects.info <- results$project
  projects.info <- results$project[,grep("state",colnames(projects.info),invert = TRUE)]

  if(any(submitter_id %in% df$submitter_id)){
    projects.info <-  cbind("submitter_id" = submitter_id, projects.info)

    suppressWarnings({
      df <- left_join(df,
                      projects.info,
                      by = "submitter_id")
    })
  } else {
    df <- dplyr::bind_cols(df,projects.info)
  }

  # Adding in the same order


  if(any(substr(barcode,1,str_length(df$submitter_id)) %in% df$submitter_id)){
    df <- df[match(substr(barcode,1,str_length(df$sample_submitter_id)),df$sample_submitter_id),]
    # This line should not exists, but some patients does not have clinical data
    # case: TCGA-R8-A6YH"
    # this has been reported to GDC, waiting answers
    # So we will remove this NA cases
    df <- df[!is.na(df$submitter_id),]
  } else {
    idx <- sapply(substr(barcode,1,str_length(df$submitter_aliquot_ids) %>% max),
                  FUN = function(x){
                    grep(x,df$submitter_aliquot_ids)
                  })
    df <- df[idx,]
  }

  # remove empty columns
  df <- df %>% as.data.frame() %>% dplyr::select(which(colSums(is.na(df)) < nrow(df)))
  return(df)
}

#' @title Prepare CEL files into an AffyBatch.
#' @description Prepare CEL files into an AffyBatch.
#' @param ClinData write
#' @param PathFolder write
#' @param TabCel write
#' @examples
#' \dontrun{
#' to add example
#' }
#' @export
#' @return Normalized Expression data from Affy eSets
TCGAprepare_Affy <- function(ClinData, PathFolder, TabCel){
  if (!requireNamespace("affy", quietly = TRUE)) {
    stop("affy package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("affy package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  affy_batch <- affy::ReadAffy(filenames = as.character(paste(TabCel$samples, ".CEL", sep = "")))

  eset <- affy::rma(affy_batch)

  mat <- Biobase::exprs(eset)

  return(mat)

}
