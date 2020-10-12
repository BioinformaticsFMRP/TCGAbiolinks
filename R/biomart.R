#' @title Get hg19 or hg38 information from biomaRt
#' @description Get hg19 or hg38 information from biomaRt
#' @param genome hg38 or hg19
#' @param as.granges Output as GRanges or data.frame
#' @importFrom biomaRt getBM useMart listDatasets useEnsembl
#' @export
get.GRCh.bioMart <- function(
    genome = c("hg19", "hg38"),
    as.granges = FALSE
) {

    genome <- match.arg(genome)
    msg <- character()
    gene.location <- tryCatch({
        host <- ifelse(
            genome == "hg19",
            "grch37.ensembl.org",
            "www.ensembl.org"
        )

        message("Accessing ", host, " to get gene information")
        ensembl <- useEnsembl(
            biomart = "ensembl",
            dataset = "hsapiens_gene_ensembl",
            host = host
        )

        attributes <- c(
            "chromosome_name",
            "start_position",
            "end_position",
            "strand",
            "ensembl_gene_id",
            "entrezgene_id",
            "external_gene_name"
        )

        db.datasets <- listDatasets(ensembl)
        description <- db.datasets[db.datasets$dataset == "hsapiens_gene_ensembl",]$description
        message(paste0("Downloading genome information (try:", tries,") Using: ", description))

        chrom <- c(1:22, "X", "Y")
        gene.location <- getBM(
            attributes = attributes,
            filters = c("chromosome_name"),
            values = list(chrom), mart = ensembl
        )
        gene.location
    }, error = function(e) {

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
        return(gene.location)
    })

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
