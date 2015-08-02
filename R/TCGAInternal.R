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
        " Query, download & analyze - TCGA                  \n",
        " Version:",utils::packageVersion("TCGAbiolinks"),"\n",
        " ==============================================================\n"
    )
    packageStartupMessage(welcome.message)

}

# Updates tcga platform and diseases
# @param env package environment
#' @importFrom stringr str_match
#' @importFrom XML readHTMLTable
#' @importFrom downloader download
#' @keywords internal
load.tcga <- function(env) {
    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?"

    # Get platform table
    tcga.query <- "query=Platform"
    next.url <- paste0(tcga.root, tcga.query)
    platform.table <- tcgaGetTable(next.url)
    platform.table <- platform.table[, 1:4]
    platform.table <- platform.table[order(platform.table$name,
                                           decreasing = TRUE),]

    # Get disease table
    tcga.query <- "query=Disease"
    next.url <- paste0(tcga.root, tcga.query)
    disease.table <- tcgaGetTable(next.url)
    disease.table <- disease.table[, 1:3]

    # Get center table
    tcga.query <- "query=Center"
    next.url <- paste0(tcga.root, tcga.query)
    center.table  <- tcgaGetTable(next.url)
    center.table <- center.table[, 1:3]

    env$platform.table <- platform.table
    env$center.table <- center.table
    env$disease.table <- disease.table

    # Get tcga folders with barcodes without private folders
    tcga.db <- createTcgaTable()
    env$tcga.db <- tcga.db
    save(tcga.db, paste0(system.file("extdata", package = "TCGAbiolinks"),
                         "/tcgadb.rda"))
    step <- 200
    for(i in seq(1, nrow(tcga.db) , by = step)){
        print(i)
        j <- i + step
        if(j > nrow(tcga.db)){
            j <-nrow(tcga.db)
        }
        tcga.db[i:j,]$deployStatus <- getBarcode(tcga.db[i:j,])$barcode
        save(tcga.db, paste0(system.file("extdata", package = "TCGAbiolinks"),
                             "/tcgadb.rda"))
    }

    env$tcga.db <- tcga.db
    save(platform.table, disease.table, tcga.db, center.table,
         file = paste0(system.file("extdata", package = "TCGAbiolinks"),
                       "/dataFolders.rda")
    )
}

# Get all files in the ftp directory @keywords internal
getFileNames <- function(url) {

    if(is.windows()){
        suppressWarnings(
            download(url, "tmp.html",
                     quiet = TRUE,
                     mode = "wb",
                     method = "internal"
            )
        )
    } else {
        suppressWarnings(
            download(url,
                     "tmp.html",
                     mode = "wb",
                     method = "internal",
                     quiet = 1)
        )
    }
    x <- capture.output(XML::htmlTreeParse("tmp.html"))
    unlink("tmp.html")
    x <- x[grep("href", x)]
    if (is.null(x)){
        return(NULL)
    }
    x = sapply(strsplit(x, ">"), function(y) y[2])
    if (is.null(x)){
        return(NULL)
    }
    x = sapply(strsplit(x, "<"), function(y) y[1])
    return(x)
}

tcga.get.barcode <- function(data){
    # todo considere next page =/
    # comparar com sdrf
    # for each tcga.db id get barcodes
    message("Downloading TCGA barcodes")
    all.barcode <- c()

    tcga.root <- "http://tcga-data.nci.nih.gov/tcgadccws/GetHTML?query="
    tcga.query <- paste0("FileInfo&Archive[@id=",data$id,
                         "]&roleName=fileCollection")

    next.url <- paste0(tcga.root,tcga.query)
    files <- tcgaGetTable(next.url)
    #print(files)
    if (nrow(files) == 0) {
        return(NULL)
    }
    files <- files[,1:4]
    idx <- grep("somatic.maf",files$name)
    if (length(idx > 0)) {
        files <- files[idx,]
    } else {
        # no maf file
        # try in the name
        pat <- paste0("((((TCGA-[A-Z0-9]{2})-[A-Z0-9]{4})-[A-Z0-9]{3}-",
                      "[A-Z0-9]{3})-[A-Z0-9]{4}-[A-Z0-9]{2})")
        barcode <- str_match(files$name,pat)[,1]
        #message("Found in name")
        if (!all(is.na(barcode))) {
            message("Found in name")
            barcode <- barcode[!is.na(barcode)]
            barcode <- paste0(unique(barcode), collapse = ",")
            return(barcode)
        }

        return(NULL)
    }

    # maybe two maf files
    for (i in  seq_along(files$name)) {
        tcga.query <- paste0("BiospecimenBarcode&FileInfo[@id=",files[i,"id"],
                             "]&roleName=biospecimenBarcodeCollection")
        next.url <- paste0(tcga.root,tcga.query)
        print(next.url)
        barcode.table <- tcgaGetTable(next.url)
        barcode.table <- barcode.table[,1:8]
        all.barcode <- union(all.barcode, unique(barcode.table$barcode))
    }

    all.barcode <- paste0(unique(all.barcode), collapse = ",")

    return(all.barcode)
}

#' @import utils
#' @importFrom RCurl getURL
.DownloadURL <-
    function(Site){
        # setInternet2(use = TRUE)
        Site <- URLencode(Site)
        x=  getURL(Site, ssl.verifypeer = FALSE)
        x <- unlist(strsplit(x,"\n"))
        return(x)
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
# @import plyr
# @import scales
#' @import grid
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
