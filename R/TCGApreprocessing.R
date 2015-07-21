#' @title TCGApreprocessing
#' @description TCGApreprocessing perform Array Array Intensity correlation (AAIC).
#' It defines a square symmetric matrix of pearson correlation among samples.
#' According this matrix and boxplot of correlation samples by samples it is possible
#' to find samples with low correlation that can be identified as possible outliers.
#' @param object object of class eset or gene expresion
#' @param tabGroupCol table with group samples information in tabGroupCol$Color
#' @param ArrayIndex ArrayIndex
#' @param pmat pmat
#' @importFrom grDevices dev.list
#' @export
#' @return Plot with array array intensity correlation and boxplot of correlation samples by samples
TCGApreprocessing<- function(object,tabGroupCol,ArrayIndex,pmat){

    if (!(is.null(dev.list()["RStudioGD"]))){dev.off()}

    png("PreprocessingOutput.png", width = 1200, height = 1200)

    # array array IC after RMA
    #object<-eset
    #ArrayIndex = as.character(1:length(sampleNames(object)))
    #pmat <- as.matrix(pData(phenoData(object)))
    phenodepth <- min(ncol(pmat), 3)
    order <- switch(phenodepth + 1, ArrayIndex, order(pmat[, 1]), order(pmat[, 1], pmat[, 2]), order(pmat[, 1], pmat[, 2], pmat[, 3]))
    arraypos <- (1:length(ArrayIndex)) * (1/(length(ArrayIndex) - 1)) - (1/(length(ArrayIndex) - 1))
    arraypos2 = seq(1:length(ArrayIndex) - 1)
    for (i in 2:length(ArrayIndex)) { arraypos2[i - 1] <- (arraypos[i] + arraypos[i - 1])/2 }
    layout(matrix(c(1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 3, 3, 3, 4), 4, 4, byrow = TRUE))

    c <- cor(exprs(object)[, order], method = "spearman")
    colnames(c) <- gsub(".CEL","",colnames(c))
    rownames(c) <- gsub(".CEL","",rownames(c))
    image(c, xaxt = "n", yaxt = "n", xlab = "Array Samples", ylab = "Array Samples",  main = "Array-Array Intensity Correlation after RMA")
    #abline(h = arraypos2, v = arraypos2)

    for ( i in 1:length(names(table(tabGroupCol$Color)) )){
        currentCol <- names(table(tabGroupCol$Color))[i]
        pos.col <- arraypos[which(tabGroupCol$Color == currentCol)]
        lab.col <- colnames(c)[which(tabGroupCol$Color == currentCol)]
        axis(1, labels = lab.col , at = pos.col, col = currentCol,lwd = 6,las =2)
        axis(2, labels = lab.col , at = pos.col, col = currentCol,lwd = 6,las =2)
    }

    m = matrix(pretty(c, 10), nrow = 1, ncol = length(pretty(c, 10)))
    image(m, xaxt = "n", yaxt = "n", ylab = "Correlation Coefficient")
    axis(2, labels = as.list(pretty(c, 10)), at = seq(0, 1, by = (1/(length(pretty(c,  10)) - 1))))
    abline(h = seq((1/(length(pretty(c, 10)) - 1))/2, 1 - (1/(length(pretty(c, 10)) - 1)), by = (1/(length(pretty(c, 10)) - 1))))

    boxplot(c, outline = FALSE,las =2, lwd = 6,col = tabGroupCol$Color, main ="Boxplot of correlation samples by samples after RMA")

    dev.off()

    return(c)
}
