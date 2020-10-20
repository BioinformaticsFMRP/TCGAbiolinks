#' @title Get hg19 or hg38 information from biomaRt
#' @description Get hg19 or hg38 information from biomaRt
#' @param genome hg38 or hg19
#' @param as.granges Output as GRanges or data.frame
#' @import TCGAbiolinksGUI.data
#' @export
get.GRCh.bioMart <- function(
    genome = c("hg19", "hg38"),
    as.granges = FALSE
) {

    genome <- match.arg(genome)
    # Since the amout of users complaining about the access we
    # also added the data into a data package
    check_package("TCGAbiolinksGUI.data")
    if (genome == "hg19") {
        gene.location <- get(
            data("gene.location.hg19",
                 package = "TCGAbiolinksGUI.data",
                 envir = environment())
        )
    } else {
        gene.location <- get(
            data(
                "gene.location.hg38",
                package = "TCGAbiolinksGUI.data",
                envir = environment()
            )
        )
    }

    if (as.granges) {
        gene.location$strand[gene.location$strand == 1] <- "+"
        gene.location$strand[gene.location$strand == -1] <- "-"
        gene.location$chromosome_name <- paste0("chr",gene.location$chromosome_name)
        gene.location <- makeGRangesFromDataFrame(
            gene.location, seqnames.field = "chromosome_name",
            start.field = "start_position",
            end.field = "end_position",
            keep.extra.columns = TRUE
        )
    }
    return(gene.location)
}

#' @noRd
map.ensg <- function(genome = "hg38", genes) {
    gene.location <- get.GRCh.bioMart(genome)
    gene.location <- gene.location[match(genes, gene.location$ensembl_gene_id), ]
    return(gene.location)
}



#' getTSS to fetch GENCODE gene annotation (transcripts level) from Bioconductor package biomaRt
#' If upstream and downstream are specified in TSS list, promoter regions of GENCODE gene will be generated.
#' @description
#' getTSS to fetch GENCODE gene annotation (transcripts level) from Bioconductor package biomaRt
#' If upstream and downstream are specified in TSS list, promoter regions of GENCODE gene will be generated.
#' @param TSS A list. Contains upstream and downstream like TSS=list(upstream, downstream).
#'  When upstream and downstream is specified, coordinates of promoter regions with gene annotation will be generated.
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @return GENCODE gene annotation if TSS is not specified. Coordinates of GENCODE gene promoter regions if TSS is specified.
#' @examples
#' # get GENCODE gene annotation (transcripts level)
#' \dontrun{
#'     getTSS <- getTSS()
#'     getTSS <- getTSS(genome.build = "hg38", TSS=list(upstream=1000, downstream=1000))
#' }
#' @importFrom GenomicRanges makeGRangesFromDataFrame promoters
#' @importFrom biomaRt useEnsembl listDatasets getBM
getTSS <- function(
    genome = "hg38",
    TSS = list(upstream = NULL, downstream = NULL)
) {
    host <- ifelse(genome == "hg19",  "grch37.ensembl.org", "www.ensembl.org")
    ensembl <- tryCatch({
        useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host =  host)
    },  error = function(e) {
        useEnsembl(
            "ensembl",
            dataset = "hsapiens_gene_ensembl",
            mirror = "uswest",
            host =  host
        )
    })
    attributes <- c(
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
        "transcript_start",
        "transcription_start_site",
        "transcript_end",
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "entrezgene_id",
        "external_gene_name"
    )

    chrom <- c(1:22, "X", "Y", "M", "*")
    db.datasets <- listDatasets(ensembl)
    description <- db.datasets[db.datasets$dataset == "hsapiens_gene_ensembl", ]$description
    message(paste0("Downloading transcripts information. Using: ", description))

    tss <- getBM(
        attributes = attributes,
        filters = c("chromosome_name"),
        values = list(chrom),
        mart = ensembl
    )
    tss <- tss[!duplicated(tss$ensembl_transcript_id), ]
    if (genome == "hg19") tss$external_gene_name <- tss$external_gene_id
    tss$chromosome_name <-  paste0("chr", tss$chromosome_name)
    tss$strand[tss$strand == 1] <- "+"
    tss$strand[tss$strand == -1] <- "-"

    tss <- makeGRangesFromDataFrame(
        tss,
        start.field = "transcript_start",
        end.field = "transcript_end",
        keep.extra.columns = TRUE
    )

    if (!is.null(TSS$upstream) & !is.null(TSS$downstream)) {
        tss <- promoters(
            tss,
            upstream = TSS$upstream,
            downstream = TSS$downstream
        )
    }

    return(tss)
}
