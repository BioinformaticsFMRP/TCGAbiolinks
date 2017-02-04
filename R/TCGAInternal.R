# @title .onAttach
# @description  Load required data into gloval enviroment
# @keywords internal
.onAttach <- function (libname, pkgname){

    if (!interactive() || stats::runif(1) > 0.1) return()
    welcome.message <- paste0(
        " =============================================================\n",
        " ______  ___  ____   ___                                        \n",
        "   ||   |    |      |   | |    o  __  |   o  _         __         \n",
        "   ||   |    | ___  |___| |__  | |  | |   | | | | |_/ |__         \n",
        "   ||   |___ |____| |   | |__| | |__| |__ | | |_| | \\  __|       \n",
        " ------------------------------------------------------------\n",
        " Query, download & analyze - GDC                  \n",
        " Version:",utils::packageVersion("TCGAbiolinks"),"\n",
        " ==============================================================\n"
    )
    packageStartupMessage(welcome.message)

}

#' @title Check GDC server status
#' @description
#'   Check GDC server status using the api
#'   https://gdc-api.nci.nih.gov/status
#' @export
#' @importFrom jsonlite fromJSON
#' @examples
#' status <- isServeOK()
#' @return Return true if status is ok
isServeOK <- function(){
    status <- fromJSON("https://gdc-api.nci.nih.gov/status",simplifyDataFrame = TRUE)$status
    if(status != "OK") stop("GDC server down, try to use this package later")
    return(TRUE)
}


checkProjectInput <- function(project){
    projects <- getGDCprojects()
    if(missing(project)) {
        print(knitr::kable(projects[,c(4:6,8)]))
        stop("Please set a project argument from the column project_id above")
    }
    for(proj in project)
        if(!(proj %in% projects$project_id)) {
            print(knitr::kable(projects[,c(4:6,8)]))
            stop("Please set a valid project argument from the column project_id above. Project ", proj, " was not found.")
        }
}

checkLegacyPlatform <- function(project,data.category, legacy = FALSE){
    project.summary <- getProjectSummary(project, legacy)
    if(missing(data.category)) {
        print(knitr::kable(project.summary$data_categories))
        stop("Please set a data.category argument from the column data_category above")
    }
    if(!(data.category %in% project.summary$data_categories$data_category)) {
        print(knitr::kable(project.summary$data_categories))
        stop("Please set a valid data.category argument from the column data_category above")
    }
}

checkDataTypeInput <- function(legacy, data.type){
    if(legacy){
        legacy.data.type <- c("Copy number segmentation",
                              "Raw intensities",
                              "Aligned reads",
                              "Copy number estimate",
                              "Simple nucleotide variation",
                              "Gene expression quantification",
                              "Coverage WIG",
                              "miRNA gene quantification",
                              "Genotypes",
                              "miRNA isoform quantification",
                              "Normalized copy numbers",
                              "Isoform expression quantification",
                              "Normalized intensities",
                              "Tissue slide image",
                              "Exon quantification",
                              "Exon junction quantification",
                              "Methylation beta value",
                              "Unaligned reads",
                              "Diagnostic image",
                              "CGH array QC",
                              "Biospecimen Supplement",
                              "Pathology report",
                              "Clinical Supplement",
                              "Intensities",
                              "Protein expression quantification",
                              "Microsatellite instability",
                              "Structural variation",
                              "Auxiliary test",
                              "Copy number QC metrics",
                              "Intensities Log2Ratio",
                              "Methylation array QC metrics",
                              "Clinical data",
                              "Copy number variation",
                              "ABI sequence trace",
                              "Biospecimen data",
                              "Simple somatic mutation",
                              "Bisulfite sequence alignment",
                              "Methylation percentage",
                              "Sequencing tag",
                              "Sequencing tag counts",
                              "LOH")
        if(!data.type %in% legacy.data.type) {
            print(knitr::kable(as.data.frame(sort(legacy.data.type))))
            stop("Please set a data.type argument from the column legacy.data.type above")
        }
    } else {
        harmonized.data.type <- c("Gene Expression Quantification",
                                  "Copy Number Segment",
                                  "Masked Copy Number Segment",
                                  "Isoform Expression Quantification",
                                  "miRNA Expression Quantification",
                                  "Biospecimen Supplement",
                                  "Clinical Supplement",
                                  "Masked Somatic Mutation")
        if(!data.type %in% harmonized.data.type) {
            print(knitr::kable(as.data.frame(sort(harmonized.data.type))))
            stop("Please set a data.type argument from the column harmonized.data.type above")
        }
    }
}

checkDataCategoriesInput <- function(project,data.category, legacy = FALSE){
    for(proj in project){
        project.summary <- getProjectSummary(proj, legacy)
        if(missing(data.category)) {
            print(knitr::kable(project.summary$data_categories))
            stop("Please set a data.category argument from the column data_category above")
        }
        if(!(data.category %in% project.summary$data_categories$data_category)) {
            print(knitr::kable(project.summary$data_categories))
            stop("Please set a valid data.category argument from the column data_category above. We could not validade the data.category for project ", proj)
        }
    }
}

checkBarcodeDefinition <- function(definition){
    for(i in definition){
        if(!(i %in% getBarcodeDefinition()$tissue.definition)){
            print(knitr::kable(getBarcodeDefinition()))
            stop(i, " was not found. Please select a difinition from the table above ")
        }
    }
}

#' @title Retrieve all GDC projects
#' @description
#'   getGDCprojects uses the following api to get projects
#'   https://gdc-api.nci.nih.gov/projects
#' @export
#' @import readr stringr
#' @examples
#' projects <- getGDCprojects()
#' @return A data frame with last GDC projects
getGDCprojects <- function(){

    projects  <- tryCatch({
        projects <- read_tsv("https://gdc-api.nci.nih.gov/projects?size=1000&format=tsv", col_types = "ccccccc")
        projects$tumor <- unlist(lapply(projects$project_id, function(x){unlist(str_split(x,"-"))[2]}))
        return(projects)
    }, error = function(e) {
        Sys.sleep(1)
        url <- "https://gdc-api.nci.nih.gov/projects?size=1000&format=json"
        json <- fromJSON(content(GET(url), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
        projects <- json$data$hits
        projects$tumor <- unlist(lapply(projects$project_id, function(x){unlist(str_split(x,"-"))[2]}))
        return(projects)
    })
    if(nrow(projects) == 0) stop("I couldn't access GDC API. Please, check if it is not down.")
    return(projects)
}

# Source: https://stackoverflow.com/questions/10266963/moving-files-between-folders
move <- function(from, to) {
    if(file.exists(from)){
        todir <- dirname(to)
        if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE,showWarnings = FALSE)
        file.rename(from = from,  to = to)
        if(dirname(from) != ".") unlink(dirname(from),recursive=TRUE,force = TRUE)
    }
}

getProjectSummary <- function(project, legacy = FALSE){
    baseURL <- ifelse(legacy,"https://gdc-api.nci.nih.gov/legacy/projects/","https://gdc-api.nci.nih.gov/projects/")
    url <- paste0(baseURL, project,"?expand=summary,summary.data_categories&pretty=true")
    return(fromJSON(url,simplifyDataFrame = TRUE)$data$summary)
}

getNbCases <- function(project, data.category, legacy = FALSE){
    summary <- getProjectSummary(project, legacy)$data_categories
    nb <- summary[summary$data_category == data.category,"case_count"]
    return(nb)
}
getNbFiles <- function(project, data.category, legacy = FALSE){
    summary <- getProjectSummary(project, legacy)$data_categories
    nb <- summary[summary$data_category == data.category,"file_count"]
    return(nb)
}

.quantileNormalization <-
    function(wd, distribution) {
        n <- nrow(wd)
        m <- ncol(wd)
        if(!missing(distribution))
            if(m != length(distribution))
                stop("The reference distribution has length
                 different from the col dimension of the data matrix.") else
                     distribution  <-  sort(distribution)

        o <- matrix(0, n, m)
        for(i in 1:n)
            o[i,] <- order(wd[i,])

        j <- 1
        tmp <- rep(0, n)

        while(j <= m) {
            if(missing(distribution)) {
                for(i in 1:n)
                    tmp[i] <- wd[i,o[i,j]]
                value <- mean(tmp)
            } else value  <- distribution[j]
            for(i in 1:n)
                wd[i,o[i,j]] <- value
            j <- j + 1
        }
        return(wd)
    }

is.windows <- function() {
    Sys.info()["sysname"] == "Windows"
}

is.mac <- function() {
    Sys.info()["sysname"] == "Darwin"
}

is.linux <- function() {
    Sys.info()["sysname"] == "Linux"
}

#
#  ggbiplot.r
#
#  Copyright 2011 Vincent Q. Vu.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#' Biplot for Principal Components using ggplot2
#'
#' @param pcobj           an object returned by prcomp() or princomp()
#' @param choices         which PCs to plot
#' @param scale           covariance biplot (scale = 1), form biplot (scale = 0).
#'                         When scale = 1, the inner product between the variables
#'                         approximates the covariance and the distance between
#'                         the points approximates the Mahalanobis distance.
#' @param obs.scale       scale factor to apply to observations
#' @param var.scale       scale factor to apply to variables
#' @param pc.biplot       for compatibility with biplot.princomp()
#' @param groups          optional factor variable indicating the groups that
#'                         the observations belong to. If provided the points
#'                         will be colored according to groups
#' @param ellipse         draw a normal data ellipse for each group?
#' @param ellipse.prob    size of the ellipse in Normal probability
#' @param labels          optional vector of labels for the observations
#' @param labels.size     size of the text used for the labels
#' @param alpha           alpha transparency value for the points
#'                        (0 = transparent, 1 = opaque)
#' @param circle          draw a correlation circle?
#'                        (only applies when prcomp was called with scale = TRUE
#'                        and when var.scale = 1)
#' @param circle.prob     definition of circle.prob
#' @param var.axes        draw arrows for the variables?
#' @param varname.size    size of the text for variable names
#' @param varname.adjust  adjustment factor the placement of the variable names,
#'                         >= 1 means farther from the arrow
#' @param varname.abbrev  whether or not to abbreviate the variable names
#' @import ggplot2
#' @importFrom plyr ddply
#' @importFrom scales muted
#' @import stats
#' @keywords internal
#' @return  a ggplot2 plot
#' @author Vincent Q. Vu.
# @examples
# data(iris)
# iris.pca <- prcomp(iris[,1:4], scale. = TRUE)
# print(ggbiplot(iris.pca, obs.scale = 1, var.scale = 1, groups = iris[,5],
# ellipse = TRUE, circle = TRUE))
ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
                     obs.scale = 1 - scale, var.scale = scale,
                     groups = NULL, ellipse = FALSE, ellipse.prob = 0.68,
                     labels = NULL, labels.size = 3, alpha = 1,
                     var.axes = TRUE,
                     circle = FALSE, circle.prob = 0.69,
                     varname.size = 3, varname.adjust = 1.5,
                     varname.abbrev = FALSE)
{
    xvar <- NULL
    yvar <- NULL
    varname <-  NULL
    angle <- NULL
    hjust <- NULL

    stopifnot(length(choices) == 2)

    # Recover the SVD
    if(inherits(pcobj, 'prcomp')){
        nobs.factor <- sqrt(nrow(pcobj$x) - 1)
        d <- pcobj$sdev
        u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
        v <- pcobj$rotation
    } else if(inherits(pcobj, 'princomp')) {
        nobs.factor <- sqrt(pcobj$n.obs)
        d <- pcobj$sdev
        u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
        v <- pcobj$loadings
    } else if(inherits(pcobj, 'PCA')) {
        nobs.factor <- sqrt(nrow(pcobj$call$X))
        d <- unlist(sqrt(pcobj$eig)[1])
        u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
        v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
    } else if(inherits(pcobj, "lda")) {
        nobs.factor <- sqrt(pcobj$N)
        d <- pcobj$svd
        u <- predict(pcobj)$x/nobs.factor
        v <- pcobj$scaling
        d.total <- sum(d^2)
    } else {
        stop('Expected a object of class prcomp, princomp, PCA, or lda')
    }

    # Scores
    choices <- pmin(choices, ncol(u))
    df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))

    # Directions
    v <- sweep(v, 2, d^var.scale, FUN='*')
    df.v <- as.data.frame(v[, choices])

    names(df.u) <- c('xvar', 'yvar')
    names(df.v) <- names(df.u)

    if(pc.biplot) {
        df.u <- df.u * nobs.factor
    }

    # Scale the radius of the correlation circle so that it corresponds to
    # a data ellipse for the standardized PC scores
    r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)

    # Scale directions
    v.scale <- rowSums(v^2)
    df.v <- r * df.v / sqrt(max(v.scale))

    # Change the labels for the axes
    if(obs.scale == 0) {
        u.axis.labs <- paste('standardized PC', choices, sep='')
    } else {
        u.axis.labs <- paste('PC', choices, sep='')
    }

    # Append the proportion of explained variance to the axis labels
    u.axis.labs <- paste(u.axis.labs,
                         sprintf('(%0.1f%% explained var.)',
                                 100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))

    # Score Labels
    if(!is.null(labels)) {
        df.u$labels <- labels
    }

    # Grouping variable
    if(!is.null(groups)) {
        df.u$groups <- groups
    }

    # Variable Names
    if(varname.abbrev) {
        df.v$varname <- abbreviate(rownames(v))
    } else {
        df.v$varname <- rownames(v)
    }

    # Variables for text label placement
    df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
    df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)

    # Base plot
    g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) +
        xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()

    if(var.axes) {
        # Draw circle
        if(circle)
        {
            theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
            circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
            g <- g + geom_path(data = circle, color = muted('white'),
                               size = 1/2, alpha = 1/3)
        }

        # Draw directions
        g <- g +
            geom_segment(data = df.v,
                         aes(x = 0, y = 0, xend = xvar, yend = yvar),
                         arrow = arrow(length = unit(1/2, 'picas')),
                         color = muted('red'))
    }

    # Draw either labels or points
    if(!is.null(df.u$labels)) {
        if(!is.null(df.u$groups)) {
            g <- g + geom_text(aes(label = labels, color = groups),
                               size = labels.size)
        } else {
            g <- g + geom_text(aes(label = labels), size = labels.size)
        }
    } else {
        if(!is.null(df.u$groups)) {
            g <- g + geom_point(aes(color = groups), alpha = alpha)
        } else {
            g <- g + geom_point(alpha = alpha)
        }
    }

    # Overlay a concentration ellipse if there are groups
    if(!is.null(df.u$groups) && ellipse) {
        theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
        circle <- cbind(cos(theta), sin(theta))

        ell <- ddply(df.u, 'groups', function(x) {
            if(nrow(x) <= 2) {
                return(NULL)
            }
            sigma <- var(cbind(x$xvar, x$yvar))
            mu <- c(mean(x$xvar), mean(x$yvar))
            ed <- sqrt(qchisq(ellipse.prob, df = 2))
            data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                       groups = x$groups[1])
        })
        names(ell)[1:2] <- c('xvar', 'yvar')
        g <- g + geom_path(data = ell, aes(color = groups, group = groups))
    }

    # Label the variable axes
    if(var.axes) {
        g <- g +
            geom_text(data = df.v,
                      aes(label = varname, x = xvar, y = yvar,
                          angle = angle, hjust = hjust),
                      color = 'darkred', size = varname.size)
    }
    # Change the name of the legend for groups
    # if(!is.null(groups)) {
    #   g <- g + scale_color_brewer(name = deparse(substitute(groups)),
    #                               palette = 'Dark2')
    # }

    # TODO: Add a second set of axes

    return(g)
}

#' @title GenesCutID
#' @description
#'   GenesCutID
#' @param GeneList GeneList
#' @return list of gene symbol without IDs
# @examples
# GenesCutID(c("CRKL|1399","TADA2A|6871","KRT76|51350"))
#' @keywords internal
GenesCutID <- function(GeneList){
    GeneListCutID <- as.matrix(matrix(unlist(strsplit(as.character(GeneList),
                                                      "|",fixed = TRUE)),length(GeneList),2,byrow = TRUE))[,1]
    return(as.matrix(GeneListCutID))
}

#' @title GeneSplitRegulon
#' @description
#'   GeneSplitRegulon
#' @param Genelist Genelist
#' @param Sep Sep
#' @keywords internal
#' @return GeneSplitRegulon
# @examples
# GeneSplitRegulon("CRKL;TADA2A;KRT76",Sep =";")
GeneSplitRegulon <- function(Genelist,Sep){
    RegSplitted <- as.matrix(unlist(strsplit(as.character(Genelist), Sep)))

    return(RegSplitted)
}


getGistic <- function(disease) {
    base <- paste0("http://gdac.broadinstitute.org/runs/analyses__latest/data/", disease)
    x <- tryCatch({
        read_html(base)  %>% html_nodes("a") %>% html_attr("href")
    }, error = function(e) {
        return(NULL)
    })
    if(is.null(x)) {
        message("No GISTIC file found")
        return(NULL)
    }
    base <- file.path(base,tail(x,n=1))
    x <- read_html(base)  %>% html_nodes("a") %>% html_attr("href")
    x <- x[grep("CopyNumber_Gistic2.Level_4",x)]

    if(!file.exists(x[1])) {
        if(Sys.info()["sysname"] == "Windows") mode <- "wb" else  mode <- "w"
        downloader::download(file.path(base,x[1]),x[1], mode = mode)
    }
    # Check if downlaod was not corrupted
    md5 <- readr::read_table(file.path(base,x[2]), col_names = FALSE, progress = FALSE)$X1
    if(tools::md5sum(x[1]) != md5) stop("Error while downloading CNV data")
    file <- paste0(gsub(".tar.gz","",x[1]),"/all_thresholded.by_genes.txt")
    if(!file.exists(file)) {
        compressed.files <- untar(x[1], list = TRUE)
        compressed.files <- compressed.files[grepl("all_thresholded.by_genes.txt", compressed.files)]
        untar(x[1],files = compressed.files )
    }
    ret <- tryCatch({
        fread(file, data.table = FALSE, colClasses = "character")
    }, error = function(e) {
        file <- dir(pattern = "all_thresholded.by_genes.txt", recursive = T, full.names = T)
        file <- file[grep(disease,file,ignore.case = TRUE)]
        fread(file, data.table = FALSE, colClasses = "character")
    })
    return(ret)
}
get.cnv <- function(project, genes){
    if(missing(project)) stop("Argument project is missing")
    if(missing(genes)) stop("Argument genes is missing")

    gistic <- getGistic(gsub("TCGA-","",project))
    cnv.annotation <- t(gistic[tolower(gistic[,1]) %in% tolower(genes),-c(2:3)])
    colnames(cnv.annotation) <- paste0("gistic2_",cnv.annotation[1,])
    cnv.annotation <- cnv.annotation[-1,, drop = FALSE]
    rownames(cnv.annotation) <- substr(gsub("\\.","-",rownames(cnv.annotation)),1,15)
    return(cnv.annotation)
}

get.mutation <- function(project, genes, pipeline = pipeline){
    if(missing(project)) stop("Argument project is missing")
    if(missing(genes)) stop("Argument genes is missing")

    # Get mutation annotation file
    maf <- GDCquery_Maf(gsub("TCGA-","",project),pipelines = pipeline)
    mut <- NULL
    for(i in genes) {
        if(!i %in% maf$Hugo_Symbol) next
        aux <-  data.frame(patient = substr(unique(maf[grepl(i,maf$Hugo_Symbol,ignore.case = TRUE),]$Tumor_Sample_Barcode),1,15), mut = TRUE)
        colnames(aux)[2] <- paste0("mut_",i)
        if(is.null(mut)) {
            mut <- aux
        } else {
            mut <- merge(mut, aux, by = "patient", all = TRUE)
        }
    }
    rownames(mut) <- mut$patient; mut$patient <- NULL

    # Lets replaces NA to FALSE
    # TRUE: has mutation
    # FALSE: has no mutation
    mut <- !is.na(mut)

    return(mut)
}
get.mut.gistc <- function(project, genes,mut.pipeline) {
    if(missing(project)) stop("Argument project is missing")
    if(missing(genes)) stop("Argument genes is missing")
    mut <- get.mutation(project, genes, pipeline = mut.pipeline )
    cnv <- get.cnv(project, genes)
    if(!is.null(mut) & !is.null(cnv)) {
        annotation <- merge(mut, cnv, by = 0 , sort = FALSE,all=TRUE)
        mut.idx <- grep("mut_",colnames(annotation))
        annotation[,mut.idx] <- !is.na(annotation[,mut.idx]) & annotation[,mut.idx] != FALSE
        rownames(annotation) <- annotation$Row.names
        annotation$Row.names <- NULL
        return(annotation)
    } else if(!is.null(mut) & is.null(cnv)) {
        return(mut)
    } else if(is.null(mut) & !is.null(cnv)) {
        return(cnv)
    }
    return(NULL)
}
get.mut.gistc.information <- function(df, project, genes, mut.pipeline = "mutect2") {
    order <- rownames(df)
    for(i in genes) if(!tolower(i) %in% tolower(EAGenes$Gene)) message(paste("Gene not found:", i))
    info <- get.mut.gistc(project, genes, mut.pipeline = mut.pipeline)
    if(is.null(info)) return(df)
    info$aux <- rownames(info)
    df$aux <- substr(df$barcode,1,15)
    df <- merge(df,info,by = "aux", all.x = TRUE, sort = FALSE)
    df$aux <- NULL
    mut.idx <- grep("mut_",colnames(df))
    df[,mut.idx] <- !is.na(df[,mut.idx]) & df[,mut.idx] != FALSE
    mut.idx <- grep("mut_",colnames(df))
    for(i in paste0("mut_",genes)){
        if(!i %in% colnames(df)) {
            df$aux <- FALSE
            colnames(df)[grep("aux",colnames(df))] <- i
        }
    }
    rownames(df) <- df$barcode
    df <- DataFrame(df[order,])
    return(df)
}


print.header <- function(text, type ="section"){
    message(paste(rep("-",nchar(text) + 3),collapse = ""))
    message(paste(ifelse(type=="section","o","oo"),text))
    message(paste(rep("-",nchar(text)+ 3),collapse = ""))
}

#' @title Get the results table from query
#' @description
#' Get the results table from query, it can select columns with cols argument
#' and return a number of rows using rows argument.
#' @param query A object from GDCquery
#' @param rows Rows identifiers (row numbers)
#' @param cols Columns identifiers (col names)
#' @export
#' @examples
#' query <- GDCquery(project = "TCGA-GBM",
#'                   data.category = "Transcriptome Profiling",
#'                   data.type = "Gene Expression Quantification",
#'                   workflow.type = "HTSeq - Counts",
#'                   barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"))
#' results <- getResults(query)
getResults <- function(query, rows, cols){
    if(missing(cols) & missing(rows)) return(query$results[[1]])
    if(missing(cols) & !missing(rows)) return(query$results[[1]][rows,])
    if(!missing(cols) & missing(rows)) return(query$results[[1]][,cols])
    if(!missing(cols) & !missing(rows)) return(query$results[[1]][rows,cols])
}
