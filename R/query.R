#' @title Query GDC data
#' @description
#'   Uses GDC API to search for search, it searches for both controlled and
#'   open-acess data.
#'   For GDC data arguments project, data.category, data.type and workflow.type should be used
#'   For the legacy data arguments project, data.category, platform and/or file.extension should be used.
#'   Please, see the vignette for a table with the possibilities.
#' @param project A list of valid project (see list with TCGAbiolinks:::getGDCprojects()$project_id)]
#' @param data.category A valid project (see list with TCGAbiolinks:::getProjectSummary(project))
#' @param data.type A data type to filter the files to download
#' @param sample.type A sample type to filter the files to download
#' @param barcode A list of barcodes to filter the files to download
#' @param legacy Search in the legacy repository
#' @param data.format Data format filter ("VCF", "TXT", "BAM","SVS","BCR XML","BCR SSF XML",
#' "TSV", "BCR Auxiliary XML", "BCR OMF XML", "BCR Biotab", "MAF", "BCR PPS XML", "XLSX")
#' @param file.type To be used in the legacy database for some platforms,
#' to define which file types to be used.
#' @param workflow.type GDC workflow type
#' @param experimental.strategy Filter to experimental stratey. Harmonized: WXS, RNA-Seq, miRNA-Seq, Genotyping Array.
#' Legacy:  WXS, RNA-Seq, miRNA-Seq, Genotyping Array,
#' DNA-Seq, Methylation array, Protein expression array, WXS,CGH array, VALIDATION, Gene expression array,WGS,
#' MSI-Mono-Dinucleotide Assay, miRNA expression array, Mixed strategies, AMPLICON, Exon array,
#' Total RNA-Seq, Capillary sequencing, Bisulfite-Seq
#' @param access Filter by access type. Possible values: controlled, open
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
#' @export
#' @examples
#' query <- GDCquery(project = "TCGA-ACC",
#'                   data.category = "Copy Number Variation",
#'                   data.type = "Copy Number Segment")
#' \dontrun{
#' query <- GDCquery(project = "TARGET-AML",
#'                   data.category = "Transcriptome Profiling",
#'                   data.type = "miRNA Expression Quantification",
#'                   workflow.type = "BCGSC miRNA Profiling",
#'                   barcode = c("TARGET-20-PARUDL-03A-01R","TARGET-20-PASRRB-03A-01R"))
#' query <- GDCquery(project = "TARGET-AML",
#'                   data.category = "Transcriptome Profiling",
#'                   data.type = "Gene Expression Quantification",
#'                   workflow.type = "HTSeq - Counts",
#'                   barcode = c("TARGET-20-PADZCG-04A-01R","TARGET-20-PARJCR-09A-01R"))
#' query <- GDCquery(project = "TCGA-ACC",
#'                   data.category =  "Copy Number Variation",
#'                   data.type = "Masked Copy Number Segment",
#'                   sample.type = c("Primary solid Tumor"))
#' query.met <- GDCquery(project = c("TCGA-GBM","TCGA-LGG"),
#'                       legacy = TRUE,
#'                       data.category = "DNA methylation",
#'                       platform = "Illumina Human Methylation 450")
#' query <- GDCquery(project = "TCGA-ACC",
#'                   data.category =  "Copy number variation",
#'                   legacy = TRUE,
#'                   file.type = "hg19.seg",
#'                   barcode = c("TCGA-OR-A5LR-01A-11D-A29H-01"))
#' }
#' @return A data frame with the results and the parameters used
#' @importFrom  jsonlite fromJSON
#' @importFrom knitr kable
#' @importFrom httr timeout
GDCquery <- function(project,
                     data.category,
                     data.type,
                     workflow.type,
                     legacy = FALSE,
                     access,
                     platform,
                     file.type,
                     barcode,
                     data.format,
                     experimental.strategy,
                     sample.type){

    isServeOK()
    suppressWarnings({
        # prepare output
        if(missing(sample.type)) {
            sample.type <- NA
        } else if(all(sample.type == FALSE)) {
            sample.type <- NA
        }
        if(missing(data.type)) {
            data.type <- NA
        } else if(data.type == FALSE) {
            data.type <- NA
        }
        if(missing(barcode)) {
            barcode <- NA
        } else if(length(barcode) == 1) {
            if(barcode == FALSE) barcode <- NA
        }
        if(missing(platform)) {
            platform <- NA
        } else if(platform == FALSE) {
            platform <- NA
        }
        if(missing(file.type)) {
            file.type <- NA
        } else if(file.type == FALSE) {
            file.type <- NA
        }
        if(missing(workflow.type)) {
            workflow.type <- NA
        } else if(workflow.type == FALSE) {
            workflow.type <- NA
        }
        if(missing(experimental.strategy)) {
            experimental.strategy <- NA
        } else if(experimental.strategy == FALSE) {
            experimental.strategy <- NA
        }
        if(missing(access)) {
            access <- NA
        } else if(access == FALSE) {
            access <- NA
        }
        if(missing(data.format)) {
            data.format <- NA
        } else if(data.format == FALSE) {
            data.format <- NA
        }
    })
    print.header("GDCquery: Searching in GDC database","section")
    message("Genome of reference: ",ifelse(legacy,"hg19","hg38"))
    # Check arguments
    checkProjectInput(project)
    checkDataCategoriesInput(project, data.category, legacy)
    if(!is.na(data.type)) checkDataTypeInput(legacy = legacy, data.type = data.type)
    if(!any(is.na(sample.type))) checkBarcodeDefinition(sample.type)

    results <- NULL
    print.header("Accessing GDC. This might take a while...","subsection")
    for(proj in project){
        url <- getGDCquery(project = proj,
                           data.category = data.category,
                           data.type = data.type,
                           legacy = legacy,
                           workflow.type = workflow.type,
                           platform = platform,
                           file.type = file.type,
                           files.access = access,
                           experimental.strategy = experimental.strategy,
                           sample.type = sample.type)
        message("ooo Project: ", proj)
        json  <- tryCatch(
            getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
            error = function(e) {
                message(paste("Error: ", e, sep = " "))
                message("We will retry to access GDC!")
                fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
            }
        )
        if(json$data$pagination$count == 0) {
            url <- getGDCquery(project = proj,
                               data.category = data.category,
                               data.type = data.type,
                               legacy = legacy,
                               workflow.type = NA,
                               platform = NA,
                               file.type = file.type,
                               experimental.strategy = experimental.strategy,
                               files.access = access,
                               sample.type = sample.type)
            json  <- tryCatch(
                getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
                error = function(e) {
                    message(paste("Error: ", e, sep = " "))
                    message("We will retry to access GDC!")
                    fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
                }
            )
        }
        json$data$hits$acl <- NULL
        json$data$hits$project <- proj
        if("archive" %in% colnames(json$data$hits)){
            if(is.data.frame(json$data$hits$archive)){
                archive <- json$data$hits$archive
                colnames(archive)[1:ncol(archive)] <- paste0("archive_", colnames(archive)[1:ncol(archive)])
                json$data$hits$archive <- NULL
                json$data$hits <- cbind(json$data$hits, archive)
            }
        }
        if("analysis" %in% colnames(json$data$hits)){
            if(is.data.frame(json$data$hits$analysis)){
                analysis <- json$data$hits$analysis
                colnames(analysis)[2:ncol(analysis)] <- paste0("analysis_", colnames(analysis)[2:ncol(analysis)])
                json$data$hits$analysis <- NULL
                json$data$hits <- cbind(json$data$hits, analysis)
            }
        }
        if("center" %in% colnames(json$data$hits)){
            if(is.data.frame(json$data$hits$center)){
                center <- json$data$hits$center
                colnames(center)[2:ncol(center)] <- paste0("center_", colnames(center)[2:ncol(center)])
                json$data$hits$center <- NULL
                json$data$hits <- cbind(json$data$hits, center)
            }
        }
        results <-  plyr::rbind.fill(as.data.frame(results),as.data.frame(json$data$hits))

    }
    if(ncol(results) == 1) {
        message("Sorry! There is no result for your query. Please check in GDC the data available or if there is no error in your query.")
        return (NULL)
    }
    print.header("Filtering results","subsection")
    if(!any(is.na(platform))){
        if(!(all(platform %in% results$platform))){
            stop("Please set a valid platform argument from the list below:\n  => ", paste(unique(results$platform), collapse = "\n  => "))
        }
        message("ooo By platform")
        results <- results[tolower(results$platform) %in% tolower(platform),]
    }

    # Filter by access
    if(!is.na(access)) {
        message("ooo By access")
        results <- results[grepl(access,results$access,ignore.case = TRUE),]
    }
    # Filter by experimental strategy
    if(!is.na(experimental.strategy)) {
        if(all(tolower(experimental.strategy) %in%  tolower(results$experimental_strategy))) {
            message("ooo By experimental.strategy")
            results <- results[tolower(results$experimental_strategy) %in% tolower(experimental.strategy),]
        } else {
            message(paste0("The argument experimental_strategy does not match any of the results.\nPossible values:",
                           paste(unique(results$experimental_strategy),collapse = "\n=>")))
        }
    }

    if(!is.na(data.format)) {
        if(all(tolower(data.format) %in%  tolower(results$data_format))) {
            message("ooo By data.format")
            results <- results[tolower(results$data_format) %in% tolower(data.format),]
        } else {
            message(paste0("The argument experimental_strategy does not match any of the results.\nPossible values:",
                           paste(unique(results$data_format),collapse = "\n=>")))
        }
    }

    # Filter by data.type
    if(!is.na(data.type)) {
        if(!(tolower(data.type) %in% tolower(results$data_type))) {
            stop("Please set a valid data.type argument from the list below:\n  => ", paste(unique(results$data_type), collapse = "\n  => "))
        }
        message("ooo By data.type")
        results <- results[tolower(results$data_type) %in% tolower(data.type),]
    }

    # Filter by workflow.type
    if(!is.na(workflow.type)) {
        if(!(workflow.type %in% results$analysis_workflow_type)) {
            stop("Please set a valid workflow.type argument from the list below:\n  => ", paste(unique(results$analysis_workflow_type), collapse = "\n  => "))
        }
        message("ooo By workflow.type")
        results <- results[results$analysis_workflow_type %in% workflow.type,]
    }


    # Filter by file.type
    if(!is.na(file.type)){
        message("ooo By file.type")
        pat <- file.type
        invert <- FALSE
        if(file.type == "normalized_results") pat <- "normalized_results"
        if(file.type == "results") pat <- "[^normalized_]results"
        if(file.type == "nocnv_hg18" | file.type == "nocnv_hg18.seg") pat <- "nocnv_hg18"
        if(file.type == "cnv_hg18" | file.type == "hg18.seg") pat <- "[^nocnv_]hg18.seg"
        if(file.type == "nocnv_hg19" | file.type == "nocnv_hg19.seg") pat <- "nocnv_hg19"
        if(file.type == "cnv_hg19" | file.type == "hg19.seg") pat <- "[^nocnv_]hg19.seg"
        if(file.type == "mirna") {
            pat <-  "hg19.*mirna"
            invert <- TRUE
        }
        # if(file.type == "hg19.mirna") pat <- "hg19.mirna"
        # if(file.type == "hg19.mirbase20.mirna") pat <- "hg19.mirbase20.mirna"
        if(file.type == "hg19.isoform") pat <- "hg19.*isoform"
        if(file.type == "isoform") {
            pat <-  "hg19.*isoform"
            invert <- TRUE
        }
        idx <- grep(pat,results$file_name,invert = invert)
        if(length(idx) == 0) {
            print(knitr::kable(sort(results$file_name)[1:10],col.names = "Files"))
            stop("We were not able to filter using this file type. Examples of available files are above. Please check the vignette for possible entries")
        }
        results <- results[idx,]
    }

    # get barcode of the samples
    # 1) Normally for each sample we will have only single information
    # however the mutation call uses both normal and tumor which are both
    # reported by the API
    if(!data.category %in% c("Clinical","Biospecimen","Other","Simple Nucleotide Variation", "Simple nucleotide variation")){

        # we also need to deal with pooled samples (mixed from different patients)
        # example CPT0000870008

        if("portions" %in% (results$cases[[1]]$samples[[1]] %>% names)) {
            aux <- plyr::laply(results$cases,
                               function(x) {
                                   summarize(x$samples[[1]],
                                             submitter_id = paste(submitter_id,collapse = ";"),
                                             is_ffpe = any(is_ffpe),
                                             sample_type =  paste(sample_type,collapse = ";"),
                                             aliquot.submiter.id = x$samples[[1]]$portions[[1]]$analytes[[1]]$aliquots[[1]]$submitter_id)
                               }) %>% as.data.frame
        } else {
            aux <- plyr::laply(results$cases,
                               function(x) {
                                   summarize(x$samples[[1]],
                                             submitter_id = paste(submitter_id,collapse = ";"),
                                             is_ffpe = any(is_ffpe),
                                             sample_type =  paste(sample_type,collapse = ";"))
                               }) %>% as.data.frame

        }
        results$sample_type <- aux$sample_type %>% as.character()
        results$is_ffpe <- aux$is_ffpe %>% as.logical

        # ORGANOID-PANCREATIC does not have aliquots
        if("aliquot.submiter.id" %in% colnames(aux)){
            results$cases <- aux$aliquot.submiter.id  %>% as.character()
            results$sample.submitter_id <- aux$submitter_id  %>% as.character()
        } else{
            results$cases <- aux$submitter_id  %>% as.character()
            results$sample.submitter_id <- aux$submitter_id  %>% as.character()
        }

    } else  if(data.category %in% c("Clinical")){
        # Clinical has another structure
        aux <- plyr::laply(results$cases,
                           function(x) {
                               unlist(x,recursive = T)[c("submitter_id")]
                           }) %>% as.data.frame
        results$cases <- aux  %>% dplyr::pull(1) %>% as.character()

    } else  if(data.category %in% c("Biospecimen")){
        # Biospecimen has another structure
        aux <- plyr::laply(results$cases,
                           function(x) {
                               paste(x$submitter_id,collapse = ",")
                           })
        results$cases <- aux
    }  else  if(data.category == "Other"){
        # Auxiliary test files does not have information linked toit.
        # get frm file names
        results$cases <- str_extract_all(results$file_name,"TCGA-[:alnum:]{2}-[:alnum:]{4}") %>% unlist
    } else if(data.category == "Simple nucleotide variation"){
        aux <- plyr::laply(results$cases,
                           function(x) {
                               lapply(x$samples,FUN = function(y)   unlist(y,recursive = T)[c("portions.analytes.aliquots.submitter_id")]) %>% unlist %>% paste(collapse = ",")
                           }) %>% as.data.frame %>% pull(1) %>% as.character()
        results$cases <- aux
    } else if(data.category == "Simple Nucleotide Variation"){

        if(data.type %in% "Masked Somatic Mutation"){
            # MAF files are one single file for all samples
            aux <- plyr::laply(results$cases[[1]]$samples,
                               function(x) {
                                   unlist(x,recursive = T)[c("portions.analytes.aliquots.submitter_id","sample_type1","sample_type2","is_ffpe1","is_ffpe2")]
                               }) %>% as.data.frame

            results$cases <- aux$portions.analytes.aliquots.submitter_id  %>% as.character() %>% paste(collapse = ",")
            if(!is.na(sample.type)) sample.type <- NA # ensure no filtering will be applied
        } else {
            # TODO: Add comnetary with case
            aux <- plyr::laply(results$cases,
                               function(x) {
                                   unlist(x$samples[[1]],recursive = T)[c("portions.analytes.aliquots.submitter_id","sample_type1","sample_type2","is_ffpe1","is_ffpe2")]
                               }) %>% as.data.frame
            results$sample_type1 <- aux$sample_type1 %>% as.character()
            results$sample_type2 <- aux$sample_type2 %>% as.character()
            results$is_ffpe1 <- aux$is_ffpe1 %>% as.logical
            results$is_ffpe2 <- aux$is_ffpe2 %>% as.logical
            results$cases <- aux$portions.analytes.aliquots.submitter_id  %>% as.character()
            if(!is.na(sample.type)) sample.type <- NA # ensure no filtering will be applied
        }
    }
    # Filter by barcode
    if(!any(is.na(barcode))) {
        message("ooo By barcode")
        idx <- unique(unlist(sapply(barcode, function(x) grep(x, results$cases,ignore.case = TRUE))))
        if(length(idx) == 0)  {
            print(knitr::kable(results$cases,col.names = "Available barcodes"))
            stop("None of the barcodes were matched. Available barcodes are above")
        }

        results <- results[idx,]
    }

    # Filter by sample.type
    if(!any(is.na(sample.type))) {
        if(!any(tolower(results$sample_type) %in% tolower(sample.type))) {
            aux <- as.data.frame(table(results$sample_type))
            aux <- aux[aux$Freq > 0,]
            print(kable(aux,row.names = FALSE,col.names = c("sample.type","Number of samples")))
            stop("Please set a valid sample.type argument from the list above.")
        }
        message("ooo By sample.type")
        results <- results[tolower(results$sample_type) %in% tolower(sample.type),]
    }
    # some how there are duplicated files in GDC we should remove them
    # Example of problematic query
    # query.exp <- GDCquery(project = "TCGA-BRCA",
    #                  legacy = TRUE,
    #                  data.category = "Gene expression",
    #                  data.type = "Gene expression quantification",
    #                  platform = "Illumina HiSeq",
    #                  file.type = "results",
    #                  experimental_strategy = "RNA-Seq",
    #                  sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
    #
    print.header("Checking data","subsection")

    message("ooo Check if there are duplicated cases")
    if(any(duplicated(results$cases))) {
        message("Warning: There are more than one file for the same case. Please verify query results. You can use the command View(getResults(query)) in rstudio")
    }

    message("ooo Check if there results for the query")
    if(nrow(results) == 0) stop("Sorry, no results were found for this query")

    print.header("Preparing output","section")
    ret <- data.frame(results = I(list(results)),
                      project = I(list(project)),
                      data.category = data.category,
                      data.type = data.type,
                      legacy = legacy,
                      access = I(list(access)),
                      experimental.strategy =  I(list(experimental.strategy)),
                      file.type = file.type,
                      platform = I(list(platform)),
                      sample.type = I(list(sample.type)),
                      barcode = I(list(barcode)),
                      workflow.type = workflow.type)
    return(ret)
}

getGDCquery <- function(project, data.category, data.type, legacy, workflow.type,platform,file.type,files.access,sample.type,experimental.strategy){
    # Get manifest using the API
    baseURL <- ifelse(legacy,"https://api.gdc.cancer.gov/legacy/files/?","https://api.gdc.cancer.gov/files/?")
    options.pretty <- "pretty=true"
    if(data.category == "Protein expression" & legacy) {
        options.expand <- "fields=archive.revision,archive.file_name,md5sum,state,data_category,file_id,platform,file_name,file_size,md5sum,submitter_id,data_type&expand=cases.samples.portions,cases.project,center,analysis"
    } else if(data.category %in% c("Clinical","Biospecimen")) {
        options.expand <- "expand=cases,cases.project,center,analysis"
    } else {
        options.expand <- "expand=cases.samples.portions.analytes.aliquots,cases.project,center,analysis,cases.samples"
    }
    option.size <- paste0("size=",getNbFiles(project,data.category,legacy))
    option.format <- paste0("format=JSON")

    options.filter <- paste0("filters=",
                             URLencode('{"op":"and","content":['),  # Start json request
                             URLencode('{"op":"in","content":{"field":"cases.project.project_id","value":["'),
                             project,
                             URLencode('"]}}'))

    if(!is.na(experimental.strategy))  options.filter <- paste0(options.filter,addFilter("files.experimental_strategy", experimental.strategy))
    if(!is.na(data.category))  options.filter <- paste0(options.filter,addFilter("files.data_category", data.category))
    if(!is.na(data.type))  options.filter <- paste0(options.filter,addFilter("files.data_type", data.type))
    if(!is.na(workflow.type))  options.filter <- paste0(options.filter,addFilter("files.analysis.workflow_type", workflow.type))
    if(!any(is.na(platform))) options.filter <- paste0(options.filter,addFilter("files.platform", platform))
    if(!any(is.na(file.type))) {
        if(file.type == "results" & legacy) options.filter <- paste0(options.filter,addFilter("files.tags", "unnormalized"))
        if(file.type == "normalized_results" & legacy) options.filter <- paste0(options.filter,addFilter("files.tags", "normalized"))
        if(file.type == "nocnv_hg19.seg" & legacy) options.filter <- paste0(options.filter,addFilter("files.tags", "nocnv"))
        if(file.type == "hg19.isoform" & legacy) options.filter <- paste0(options.filter,addFilter("files.tags", "hg19"))
    }
    if(!any(is.na(files.access))) {
        options.filter <- paste0(options.filter,addFilter("files.access", files.access))
    }
    if(!any(is.na(sample.type))) {
        if("Primary solid Tumor" %in% sample.type) sample.type[sample.type == "Primary solid Tumor"] <- "Primary Tumor"
        if("Recurrent Solid Tumor" %in% sample.type) sample.type[sample.type == "Recurrent Solid Tumor"] <- "Recurrent Tumor"
        options.filter <- paste0(options.filter,addFilter("cases.samples.sample_type", sample.type))
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

addFilter <- function(field, values){
    ret <- paste0(
        URLencode(',{"op":"in","content":{"field":"'),
        URLencode(field),
        URLencode('","value":["'),
        URLencode(paste0(values, collapse = '","')),
        URLencode('"]}}')
    )
    return(ret)
}

expandBarcodeInfo <- function(barcode){
    if(any(grepl("TARGET",barcode))) {
        ret <- DataFrame(barcode = barcode,
                         code = substr(barcode, 8, 9),
                         case.unique.id = substr(barcode, 11, 16),
                         tissue.code = substr(barcode, 18, 19),
                         nucleic.acid.code = substr(barcode, 24, 24))
        ret <- merge(ret,getBarcodeDefinition(), by = "tissue.code", sort = FALSE, all.x = TRUE)
        ret <- ret[match(barcode,ret$barcode),]
    }
    if(any(grepl("TCGA",barcode))) {
        ret <- data.frame(barcode = barcode,
                          patient = substr(barcode, 1, 12),
                          sample = substr(barcode, 1, 16),
                          tissue.code = substr(barcode, 14, 15))
        ret <- merge(ret,getBarcodeDefinition(), by = "tissue.code", sort = FALSE, all.x = TRUE)
        ret <- ret[match(barcode,ret$barcode),]
    }
    return(ret)
}


getBarcodeDefinition <- function(type = "TCGA"){
    if(type == "TCGA"){
        tissue.code <- c('01','02','03','04','05','06','07','08','09','10','11',
                         '12','13','14','20','40','50','60','61')
        shortLetterCode <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC",
                             "TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB",
                             "CELL","XP","XCL")

        tissue.definition <- c("Primary Tumor",
                               "Recurrent Tumor",
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
        aux <- data.frame(tissue.code = tissue.code,shortLetterCode,tissue.definition)
    } else {
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

    }
    return(aux)
}


#' @title Retrieve open access maf files from GDC server
#' @description
#'   GDCquery_Maf uses the following guide to download maf files
#'   https://gdc-docs.nci.nih.gov/Data/Release_Notes/Data_Release_Notes/
#' @param pipelines Four separate variant calling pipelines are implemented for GDC data harmonization.
#' Options: muse, varscan2, somaticsniper, mutect2. For more information:
#' https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/
#' @param tumor a valid tumor
#' @param save.csv Write maf file into a csv document
#' @param directory Directory/Folder where the data will downloaded. Default: GDCdata
#' @export
#' @importFrom data.table fread
#' @import readr stringr
#' @importFrom downloader download
#' @importFrom R.utils gunzip
#' @importFrom tools md5sum
#' @examples
#' \dontrun{
#'    acc.muse.maf <- GDCquery_Maf("ACC", pipelines = "muse")
#'    acc.varscan2.maf <- GDCquery_Maf("ACC", pipelines = "varscan2")
#'    acc.somaticsniper.maf <- GDCquery_Maf("ACC", pipelines = "somaticsniper")
#'    acc.mutect.maf <- GDCquery_Maf("ACC", pipelines = "mutect2")
#' }
#' @return A data frame with the maf file information
GDCquery_Maf <- function(tumor,
                         save.csv = FALSE,
                         directory = "GDCdata",
                         pipelines = NULL){

    if(is.null(pipelines)) stop("Please select the pipeline argument (muse, varscan2, somaticsniper, mutect2)")
    if(grepl("varscan",pipelines, ignore.case = TRUE)) {
        workflow.type <- "VarScan2 Variant Aggregation and Masking"
    } else if(pipelines == "muse") {
        workflow.type <- "MuSE Variant Aggregation and Masking"
    } else if(pipelines == "somaticsniper") {
        workflow.type <- "SomaticSniper Variant Aggregation and Masking"
    } else if(grepl("mutect",pipelines, ignore.case = TRUE)) {
        workflow.type <-  "MuTect2 Variant Aggregation and Masking"
    } else {
        stop("Please select the pipeline argument (muse, varscan2, somaticsniper, mutect2)")
    }

    #  Info to user
    message("============================================================================")
    message(" For more information about MAF data please read the following GDC manual and web pages:")
    message(" GDC manual: https://gdc-docs.nci.nih.gov/Data/PDF/Data_UG.pdf")
    message(" https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/")
    message(" https://gdc.cancer.gov/about-gdc/variant-calling-gdc")
    message("============================================================================")

    query <- GDCquery(paste0("TCGA-",tumor),
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      workflow.type = workflow.type,
                      access = "open")

    if(nrow(query$results[[1]]) == 0) stop("No MAF file found for this type of workflow")
    maf <- tryCatch({
        tryCatch({
            GDCdownload(query, directory = directory, method = "api")
        }, error = function(e) {
            GDCdownload(query, directory = directory, method = "client")
        })
        maf <- GDCprepare(query, directory = directory)
        maf
    }, error = function(e) {
        manifest <- getManifest(query)
        GDCdownload.aux( "https://api.gdc.cancer.gov/data/", manifest, manifest$filename, ".")
        maf <- readSimpleNucleotideVariationMaf(file.path(manifest$id,manifest$filename))
        maf
    })

    if(save.csv) {
        fout <- file.path(directory,gsub("\\.gz", "\\.csv",getResults(query)$file_name))
        write_csv(maf, fout)
        message(paste0("File created: ", fout))
    }
    return(maf)
}



#' @title Retrieve open access mc3 MAF file from GDC server
#' @description
#'   Download data from https://gdc.cancer.gov/about-data/publications/mc3-2017
#'   https://gdc-docs.nci.nih.gov/Data/Release_Notes/Data_Release_Notes/
#' @examples
#' \dontrun{
#'    maf <- getMC3MAF()
#' }
#' @return A data frame with the MAF file information from https://gdc.cancer.gov/about-data/publications/mc3-2017
#' @export
getMC3MAF <- function(){
    fout <- "mc3.v0.2.8.PUBLIC.maf.gz"
    fpath <- "https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc"
    if(is.windows()) mode <- "wb" else  mode <- "w"
    message(rep("-",100))
    message("o Starting to download Publi MAF from GDC")
    message("o More information at: https://gdc.cancer.gov/about-data/publications/mc3-2017")
    message("o Please, cite: Cell Systems. Volume 6 Issue 3: p271-281.e7, 28 March 2018 10.1016/j.cels.2018.03.002")
    if(!file.exists(gsub("\\.gz", "", fout))){
        download(fpath, fout, mode = mode)
        message("o Uncompressing file")
        gunzip(fout, remove = FALSE)
    }
    message("o Reading MAF")
    maf <- readr::read_tsv(gsub("\\.gz", "", fout),progress = TRUE, col_types = readr::cols())
    message("o Adding project_id information")
    project <- grep("TCGA",sort(getGDCprojects()$project_id),value = TRUE)
    df <- plyr::adply(project,
                .margins = 1,
                .fun = function(proj) {
                    samples <- getSubmitterID(proj)
                    return(data.frame(proj,samples))
                    }
                )
    maf$project_id <- df$proj[match(substr(maf$Tumor_Sample_Barcode,1,12),df$samples)] %>% as.character
    message(rep("-",100))
}


#' @title Query gene counts of TCGA and GTEx data from the Recount2 project
#' @description
#'   TCGArecount2_query queries and downloads data produced by the Recount2 project. User can specify which project and which tissue to query
#' @param project is a string denoting which project the user wants. Options are "tcga" and "gtex"
#' @param tissue a vector of tissue(s) to download. Options are "adipose tissue", "adrenal", "gland", "bladder","blood", "blood vessel", "bone marrow", "brain", "breast","cervix uteri", "colon", "esophagus", "fallopian tube","heart", "kidney", "liver", "lung", "muscle", "nerve", "ovary","pancreas", "pituitary", "prostate", "salivary", "gland", "skin", "small intestine", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina"
#' @export
#' @examples
#' \dontrun{
#'  brain.rec<-TCGAquery_recount2(project = "gtex", tissue = "brain")
#' }
#' @return List with $subtypes attribute as a dataframe with barcodes, samples, subtypes, and colors. The $filtered attribute is returned as filtered samples with no subtype info
TCGAquery_recount2<-function(project, tissue=c()){
    tissuesGTEx <- c(
        "adipose_tissue",
        "adrenal_gland",
        "bladder",
        "blood",
        "blood_vessel",
        "bone_marrow",
        "brain",
        "breast",
        "cervix_uteri",
        "colon",
        "esophagus",
        "fallopian_tube",
        "heart",
        "kidney",
        "liver",
        "lung",
        "muscle",
        "nerve",
        "ovary",
        "pancreas",
        "pituitary",
        "prostate",
        "salivary_gland",
        "skin",
        "small_intestine",
        "spleen",
        "stomach",
        "testis",
        "thyroid",
        "uterus",
        "vagina"
    )

    tissuesTCGA <- c(
        "adrenal_gland",
        "bile_duct",
        "bladder",
        "bone_marrow",
        "brain",
        "breast",
        "cervix",
        "colorectal",
        "esophagus",
        "eye",
        "head_and_neck",
        "kidney",
        "liver",
        "lung",
        "lymph_nodes",
        "ovary",
        "pancreas",
        "pleura",
        "prostate",
        "skin",
        "soft_tissue",
        "stomach",
        "testis",
        "thymus",
        "thyroid",
        "uterus")


    tissue<-paste(unlist(strsplit(tissue, " ")), collapse="_")
    Res<-list()

    if(tolower(project)=="gtex"){
        for(t_i in tissue){
            if(tissue%in%tissuesGTEx){
                con<-"http://duffel.rail.bio/recount/v2/SRP012682/rse_gene_"
                con<-paste0(con,tissue,".Rdata")

                message(paste0("downloading Range Summarized Experiment for: ", tissue))
                load(url(con))
                Res[[paste0(project,"_", t_i)]]<-rse_gene
            }
            else stop(paste0(tissue, " is not an available tissue on Recount2"))
        }
        return(Res)
    }
    else if(tolower(project)=="tcga"){
        for(t_i in tissue){
            if(tissue%in%tissuesTCGA){
                con<-"http://duffel.rail.bio/recount/v2/TCGA/rse_gene_"
                con<-paste0(con,tissue,".Rdata")
                message(paste0("downloading Range Summarized Experiment for: ", tissue))
                load(url(con))
                Res[[paste0(project,"_", t_i)]]<-rse_gene

            }
            else stop(paste0(tissue, " is not an available tissue on Recount2"))
        }
        return(Res)
    }
    else stop(paste0(project, " is not a valid project"))

}

#' @title Retrieve open access ATAC-seq files from GDC server
#' @description
#'   Retrieve open access ATAC-seq files from GDC server
#'  https://gdc.cancer.gov/about-data/publications/ATACseq-AWG
#'  Manifest available at: https://gdc.cancer.gov/files/public/file/ATACseq-AWG_Open_GDC-Manifest.txt
#' @param tumor a valid tumor
#' @param file.type Write maf file into a csv document
#' @export
#' @examples
#' \dontrun{
#'    query <- GDCquery_ATAC_seq(file.type = "txt")
#'    GDCdownload(query)
#'    query <- GDCquery_ATAC_seq(file.type = "bigWigs")
#'    GDCdownload(query)
#' }
#' @return A data frame with the maf file information
GDCquery_ATAC_seq <- function(tumor = NULL,
                              file.type = NULL) {
    isServeOK()
    results <- readr::read_tsv("https://gdc.cancer.gov/files/public/file/ATACseq-AWG_Open_GDC-Manifest.txt")

    if(!is.null(tumor)) results <- results[grep(tumor,results$filename,ignore.case = T),]
    if(!is.null(file.type))  results <- results[grep(file.type,results$filename,ignore.case = T),]

    colnames(results) <- c("file_id", "file_name", "md5sum", "file_size")
    results$state <- "released"
    results$data_type <- "ATAC-seq"
    results$data_category <- "ATAC-seq"
    results$project <- "ATAC-seq"
    ret <- data.frame(results=I(list(results)),
                      tumor = I(list(tumor)),
                      project = I(list("ATAC-seq")),
                      data.type = I(list("ATAC-seq")),
                      data.category = I(list("ATAC-seq")),
                      legacy = I(list(FALSE)))

    return(ret)
}


