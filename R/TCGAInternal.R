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

    suppressWarnings(
        download(url, "tmp.html", quiet = TRUE)
    )
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

.scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
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




.heatmap.plus.sm <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,
                             distfun = dist, hclustfun = hclust,
                             reorderfun = function(d, w) reorder(d, w),
                             add.expr, symm = FALSE,
                             revC = identical(Colv,"Rowv"),
                             scale = c("row", "column", "none"), na.rm = TRUE,
                             margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 +
                                 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
                             labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE,
                             verbose = getOption("verbose"), breaks, ...)
{
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv))
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv))
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram"))
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv)
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr)))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram"))
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc)
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm)
                x
                else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv)
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc)))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow))
        if (is.null(rownames(x)))
            (1:nr)[rowInd]
    else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol))
        if (is.null(colnames(x)))
            (1:nc)[colInd]
    else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0,
              4)
    if (!missing(ColSideColors)) {
        if (!is.matrix(ColSideColors))
            stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || dim(ColSideColors)[1] !=
            nc)
            stop("'ColSideColors' dim()[2] must be of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1], 0.6, lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors))
            stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || dim(RowSideColors)[1] !=
            nr)
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1),
                                       1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei,
            "; lmat=\n")
        print(lmat)
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = RowSideColors[rowInd, ]
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
            rsc.colors[rsc.i] = rsc.name
            rsc[rsc == rsc.name] = rsc.i
            rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(colnames(RowSideColors)) > 0) {
            axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors),
                 las = 2, tick = FALSE)
        }
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, ]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
            csc.colors[csc.i] = csc.name
            csc[csc == csc.name] = csc.i
            csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 0) {
            axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), colnames(ColSideColors),
                 las = 2, tick = FALSE)
        }
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") {
        x <- t(x)
    }
    if (revC) {
        iy <- nr:1
        ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
              c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend)
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
    if (doCdend)
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main))
        frame()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    invisible(list(rowInd = rowInd, colInd = colInd,
                   Rowv = if (keep.dendro && doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}
