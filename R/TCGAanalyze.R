#' @title Array Array Intensity correlation (AAIC) and correlation boxplot to define outlier
#' @description TCGAanalyze_Preprocessing perform Array Array Intensity correlation (AAIC).
#' It defines a square symmetric matrix of pearson correlation among samples.
#' According this matrix and boxplot of correlation samples by samples it is possible
#' to find samples with low correlation that can be identified as possible outliers.
#' @param object of gene expression of class RangedSummarizedExperiment from TCGAprepare
#' @importFrom grDevices dev.list
#' @export
#' @return Plot with array array intensity correlation and boxplot of correlation samples by samples
TCGAanalyze_Preprocessing<- function(object){

    if (!(is.null(dev.list()["RStudioGD"]))){dev.off()}

    png("PreprocessingOutput.png", width = 1200, height = 1200)

    # array array IC after RMA
    #object <-BRCARnaseq_assay

    #object<-eset_COMBAT
    #ArrayIndex = as.character(1:length(sampleNames(object)))
    ArrayIndex = as.character(1:length( colData(object)$sample))

    pmat_new <- matrix(0, length(ArrayIndex),4)
    colnames(pmat_new) <-c("Disease","platform","SampleID","Study")
    rownames(pmat_new)<- as.character(colData(object)$sample)
    pmat_new <- as.data.frame(pmat_new)
    pmat_new$Disease <-as.character(colData(object)$shortLetterCode)
    pmat_new$platform <-"platform"
    pmat_new$SampleID <- as.character(colData(object)$sample)
    pmat_new$Study<-"study"

    tabGroupCol <-cbind(pmat_new, Color = matrix(0,nrow(pmat_new),1))
    tabGroupCol[which(tabGroupCol$Disease=="TP"),"Color"]<-"red"
    tabGroupCol[which(tabGroupCol$Disease=="NT"),"Color"]<-"blue"

    #    pmat <- as.matrix(pData(phenoData(object)))
    pmat <- pmat_new
    phenodepth <- min(ncol(pmat), 3)
    order <- switch(phenodepth + 1, ArrayIndex, order(pmat[, 1]), order(pmat[, 1], pmat[, 2]), order(pmat[, 1], pmat[, 2], pmat[, 3]))
    arraypos <- (1:length(ArrayIndex)) * (1/(length(ArrayIndex) - 1)) - (1/(length(ArrayIndex) - 1))
    arraypos2 = seq(1:length(ArrayIndex) - 1)
    for (i in 2:length(ArrayIndex)) { arraypos2[i - 1] <- (arraypos[i] + arraypos[i - 1])/2 }
    layout(matrix(c(1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 3, 3, 3, 4), 4, 4, byrow = TRUE))

    #c <- cor(exprs(object)[, order], method = "spearman")
    c <- cor(assay(object,"raw_counts")[, order], method = "spearman")

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

#' @title survival analysis (SA) univariate with Kaplan-Meier (KM) method.
#' @description TCGAanalyze_SurvivalKM perform an univariate Kaplan-Meier (KM) survival analysis (SA).
#' It performed Kaplan-Meier survival univariate using complte follow up with all days
#' taking one gene a time from Genelist of gene symbols.
#' For each gene according its level of mean expression in cancer samples,
#' defining two thresholds for quantile
#' expression of that gene in all samples (default ThreshTop=0.67,ThreshDown=0.33) it is possible
#' to define a threshold of intensity of gene expression to divide the samples in 3 groups
#' (High, intermediate, low).
#' TCGAanalyze_SurvivalKM performs SA between High and low groups using following functions
#' from survival package
#' \enumerate{
#' \item survival::Surv
#' \item survival::survdiff
#' \item survival::survfit
#' }
#' @param clinical_patient is a data.frame using function 'clinic' with information
#' related to barcode / samples such as bcr_patient_barcode, days_to_death ,
#' days_to_last_followup , vital_status, etc
#' @param dataGE is a matrix of Gene expression (genes in rows, samples in cols) from TCGAprepare
#' @param Genelist is a list of gene symbols where perform survival KM.
#' @param Survresult is a parameter (default = FALSE) if is TRUE will show KM plot and results.
#' @param ThreshTop is a quantile threshold to identify samples with high expression of a gene
#' @param ThreshDown is a quantile threshold to identify samples with low expression of a gene
#' @importFrom survival Surv survdiff survfit
#' @export
#' @return table with survival genes pvalues from KM.
#' @examples
#' \dontrun{
#' clinical_patient_Cancer <- TCGAquery_clinic("brca","clinical_patient")
#' dataBRCAcomplete <- log2(BRCA_rnaseqv2)
#' # Selecting only 10 genes for example
#' dataBRCAcomplete <- dataBRCAcomplete[1:10,]
#' tabSurvKM<-TCGAanalyze_SurvivalKM(clinical_patient_Cancer,dataBRCAcomplete,
#' Genelist = rownames(dataBRCAcomplete), Survresult = FALSE,ThreshTop=0.67,ThreshDown=0.33)
#' # Filtering by pvalue < 0.01
#' tabSurvKM <- tabSurvKM[tabSurvKM$pvalue < 0.01,]
#' tabSurvKM <- tabSurvKM[!duplicated(tabSurvKM$mRNA),]
#' rownames(tabSurvKM) <-tabSurvKM$mRNA
#' tabSurvKM <- tabSurvKM[,-1]
#' tabSurvKM <- tabSurvKM[order(tabSurvKM$pvalue, decreasing=FALSE),]
#' }
TCGAanalyze_SurvivalKM<-function(clinical_patient,dataGE,Genelist, Survresult,ThreshTop=0.67, ThreshDown=0.33){
    samplesNT <- TCGAquery_SampleTypes(colnames(dataGE), typesample = c("NT"))
    samplesTP <- TCGAquery_SampleTypes(colnames(dataGE), typesample = c("TP"))
    Genelist <- intersect(rownames(dataGE),Genelist)
    dataCancer <- dataGE[Genelist,samplesTP]
    dataNormal <- dataGE[Genelist,samplesNT]
    colnames(dataCancer)  <- substr(colnames(dataCancer),1,12)
    cfu<-clinical_patient[clinical_patient[,"bcr_patient_barcode"] %in% substr(colnames(dataCancer),1,12),]
    cfu <- as.data.frame(subset(cfu, select=c("bcr_patient_barcode","days_to_death","days_to_last_followup","vital_status"))  )
    cfu[which(cfu$vital_status=="Alive"),"days_to_death"]<-"-Inf"
    cfu[which(cfu$vital_status=="Dead"),"days_to_last_followup"]<-"-Inf"

    followUpLevel<-FALSE
    Survresult<-FALSE

    #FC_FDR_table_mRNA
    mRNAselected_surv_results_Matrix<-matrix(0,nrow(as.matrix(rownames(dataNormal))),8)
    colnames(mRNAselected_surv_results_Matrix)<-c("mRNA","pvalue","Cancer Deaths","Cancer Deaths with Top","Cancer Deaths with Down","Mean Tumor Top","Mean Tumor Down","Mean Normal")

    mRNAselected_surv_results_Matrix<-as.data.frame(mRNAselected_surv_results_Matrix)

    cfu$days_to_death<-as.numeric(as.character(cfu$days_to_death))
    cfu$days_to_last_followup<-as.numeric(as.character(cfu$days_to_last_followup))
    rownames(cfu) <- cfu[, "bcr_patient_barcode" ] #mod1
    cfu_complete<-cfu
    ngenes<-nrow(as.matrix(rownames(dataNormal)))

    for( i in 1: nrow(as.matrix(rownames(dataNormal))))  {
        #print(i)
        cat(paste( (ngenes-i),".",sep=""))

        mRNAselected<-as.matrix(rownames(dataNormal))[i]
        mRNAselected_surv_results_Matrix[i,"mRNA"]<-mRNAselected

        mRNAselected_values<-dataCancer[rownames(dataCancer) == mRNAselected,]
        mRNAselected_values_normal<-dataNormal[rownames(dataNormal) == mRNAselected,]

        mRNAselected_values_ordered<-sort(mRNAselected_values,decreasing=TRUE)
        mRNAselected_values_ordered_top<-as.numeric(quantile(mRNAselected_values_ordered,ThreshTop)[1])
        mRNAselected_values_ordered_down<-as.numeric(quantile(mRNAselected_values_ordered,ThreshDown)[1])

        mRNAselected_values_newvector<-mRNAselected_values


        if (is.na(mRNAselected_values_ordered_top)!=1){

            numberOfSamples<-nrow(as.matrix(mRNAselected_values_ordered))
            lastelementTOP<-round(numberOfSamples/3)

            firstelementDOWN<- numberOfSamples  - lastelementTOP

            samples_top_mRNA_selected<-rownames( as.matrix(mRNAselected_values_ordered[1: (lastelementTOP-1)  ] ))
            samples_down_mRNA_selected<-rownames( as.matrix(mRNAselected_values_ordered[ (firstelementDOWN+1) : numberOfSamples] ))

            samples_UNCHANGED_mRNA_selected<-rownames(as.matrix(which((mRNAselected_values_newvector) > mRNAselected_values_ordered_down & mRNAselected_values_newvector < mRNAselected_values_ordered_top )))

            cfu_onlyTOP<-cfu_complete[cfu_complete[,"bcr_patient_barcode"] %in% samples_top_mRNA_selected,]
            cfu_onlyDOWN<-cfu_complete[cfu_complete[,"bcr_patient_barcode"] %in% samples_down_mRNA_selected,]
            cfu_onlyUNCHANGED<-cfu_complete[cfu_complete[,"bcr_patient_barcode"] %in% samples_UNCHANGED_mRNA_selected,]

            #if( followUpLevel == TRUE)
            #{
            # samplesTop_over_followUplevel<- !(cfu_onlyTOP[,"days_to_death"] < 0 &   cfu_onlyTOP[,"days_to_last_followup"] < Thresh_followUP)
            #  cfu_onlyTOP<- cfu_onlyTOP[samplesTop_over_followUplevel,]
            # samplesDown_over_followUplevel<- !(cfu_onlyDOWN[,"days_to_death"] < 0 &   cfu_onlyDOWN[,"days_to_last_followup"] < Thresh_followUP)
            #cfu_onlyDOWN<- cfu_onlyDOWN[samplesDown_over_followUplevel,]
            #print(paste("Processing ... with followUP level >",Thresh_followUP," days and",nrow(cfu),"clinical samples"))
            #  }

            # else {
            #  print(paste("Processing ... without followUP level and", nrow(as.matrix(cfu)),"clinical samples"))

            #}
            cfu_ordered<-NULL
            cfu_ordered<-rbind(cfu_onlyTOP,cfu_onlyDOWN)
            cfu<-cfu_ordered

            # print(dim(cfu))

            # } #end else with all samples

            ttime <- as.numeric(cfu[, "days_to_death"])

            #ttime <- cfu[, "days_to_death"]
            sum(status <- ttime > 0) # morti
            deads_complete <- sum(status <- ttime > 0)

            ttime_only_top <- cfu_onlyTOP[, "days_to_death"]
            deads_top<- sum(ttime_only_top > 0)


            if(  dim(cfu_onlyDOWN)[1] >= 1) {
                ttime_only_down <- cfu_onlyDOWN[, "days_to_death"]
                deads_down<- sum(ttime_only_down > 0)
            }

            else {deads_down <-0 }


            #print(paste("deaths =",deads_complete))
            mRNAselected_surv_results_Matrix[i,"Cancer Deaths"]<-deads_complete
            mRNAselected_surv_results_Matrix[i,"Cancer Deaths with Top"]<- deads_top
            mRNAselected_surv_results_Matrix[i,"Cancer Deaths with Down"]<- deads_down

            mRNAselected_surv_results_Matrix[i,"Mean Normal"]<-  mean(mRNAselected_values_normal)





            dataCancer_onlyTop_sample<-dataCancer[,samples_top_mRNA_selected]
            dataCancer_onlyTop_sample_mRNASelected<- dataCancer_onlyTop_sample[rownames(dataCancer_onlyTop_sample) == mRNAselected,]


            dataCancer_onlyDown_sample<-dataCancer[,samples_down_mRNA_selected]
            dataCancer_onlyDown_sample_mRNASelected<- dataCancer_onlyDown_sample[rownames(dataCancer_onlyDown_sample) == mRNAselected,]


            mRNAselected_surv_results_Matrix[i,"Mean Tumor Top"]<- mean(dataCancer_onlyTop_sample_mRNASelected)
            mRNAselected_surv_results_Matrix[i,"Mean Tumor Down"]<- mean(dataCancer_onlyDown_sample_mRNASelected)

            ttime[!status] <- as.numeric(cfu[!status, "days_to_last_followup"])
            #ttime[!status] <- cfu[!status, "days_to_last_followup"]

            ttime[which(ttime== -Inf)]<-0


            ttime <- Surv(ttime, status)
            rownames(ttime) <- rownames(cfu)
            length(ttime)
            #plot(survfit(ttime ~ 1))

            #plot(survfit(ttime ~ c(rep("top", nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN)), rep("unchanged", nrow(cfu_onlyUNCHANGED)))), col = c("red", "green","grey"))

            #   plot(survfit(ttime ~ c(rep("top", nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN)))), col = c("red", "green"),main= mRNAselected)





            legendHigh<- paste(mRNAselected,"High")
            legendLow<- paste(mRNAselected,"Low")






            mRNAselected_surv_results<-survdiff(ttime  ~ c(rep("top", nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN)) ))
            mRNAselected_surv_results_chis<-unlist(mRNAselected_surv_results)$chisq

            mRNAselected_surv_results_pvalue <- as.numeric(1 - pchisq(abs(mRNAselected_surv_results$chisq), df = 1))
            #miRselected_surv_results_pvalue <- as.numeric(round(as.numeric(1 - pchisq(abs(miRselected_surv_results$chisq), df = 1)),6))

            mRNAselected_surv_results_Matrix[i,"pvalue"]<-mRNAselected_surv_results_pvalue


            #print(paste(i,"....",mRNAselected,"pvalue=",mRNAselected_surv_results_pvalue))

            if (Survresult ==TRUE) {
                titlePlot<- paste("Kaplan-Meier Survival analysis, pvalue=",mRNAselected_surv_results_pvalue )


                plot(survfit(ttime ~ c(rep("low", nrow(cfu_onlyTOP)), rep("high", nrow(cfu_onlyDOWN)))), col = c("green", "red"),main= titlePlot,xlab="Days",ylab="Survival")


                legend(100, 1, legend = c(legendLow,legendHigh), col = c("green", "red"), text.col = c("green", "red"), pch = 15)

                print(mRNAselected_surv_results)
            }
        } #end if

    } #end for

    mRNAselected_surv_results_Matrix[mRNAselected_surv_results_Matrix=="-Inf"]<-0

    return(mRNAselected_surv_results_Matrix)
}


#' @title Filtering mRNA transcripts and miRNA selecting a threshold.
#' @description
#'    TCGAanalyze_Filtering allows user to filter mRNA transcripts and miRNA,
#'    selecting a threshold. For istance returns all mRNA or miRNA with mean across all
#'    samples, higher than the threshold defined quantile mean across all samples.
#' @param TableRnaseq is a dataframe or numeric matrix, each row represents a gene,
#' each column represents a sample come from TCGAPrepare
#' @param QuantileThresh is threshold selected as mean for filtering
#' @export
#' @return A filtered dataframe or numeric matrix where each row represents a gene,
#' each column represents a sample
#' @examples
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
#' dataFilt <- TCGAanalyze_Filtering(dataNorm, 0.25)
TCGAanalyze_Filtering <- function(TableRnaseq,QuantileThresh ){
    GeneThresh <- as.numeric(quantile(rowMeans(TableRnaseq), QuantileThresh))
    geneFiltered <- names(which(rowMeans(TableRnaseq) > GeneThresh))
    Table_Rnaseq_Rawcount_Filt <- TableRnaseq[geneFiltered, ]
    return( Table_Rnaseq_Rawcount_Filt)
}

#' @title normalization mRNA transcripts and miRNA using EDASeq package.
#' @description
#'   TCGAanalyze_Normalization allows user to normalize mRNA transcripts and miRNA,
#'    using EDASeq package.
#'
#'    Normalization for RNA-Seq Numerical and graphical
#'     summaries of RNA-Seq read data. Within-lane normalization procedures
#'    to adjust for GC-content effect (or other gene-level effects) on read counts:
#'    loess robust local regression, global-scaling, and full-quantile normalization
#'    (Risso et al., 2011). Between-lane normalization procedures to adjust for
#'    distributional differences between lanes (e.g., sequencing depth):
#'    global-scaling and full-quantile normalization (Bullard et al., 2010).
#'
#'    For istance returns all mRNA or miRNA with mean across all
#'    samples, higher than the threshold defined quantile mean across all samples.
#'
#'    TCGAanalyze_Normalization performs normalization using following functions
#'    from EDASeq
#'    \enumerate{
#'    \item  EDASeq::newSeqExpressionSet
#'    \item  EDASeq::withinLaneNormalization
#'    \item  EDASeq::betweenLaneNormalization
#'    \item  EDASeq::counts
#'    }
#' @param TCGA_RnaseqTable Rnaseq numeric matrix, each row represents a gene,
#' each column represents a sample
#' @param geneInfo Information matrix of 20531 genes about geneLength and gcContent
#' @param method is method of normalization such as 'gcContent' or 'geneLength'
#' @importFrom EDASeq newSeqExpressionSet withinLaneNormalization
#'  betweenLaneNormalization exprs counts offst
#' @export
#' @return Rnaseq matrix normalized with counts slot holds the count data as a matrix
#' of non-negative integer count values, one row for each observational unit (gene or the like),
#' and one column for each sample.
#' @examples
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
TCGAanalyze_Normalization <- function(TCGA_RnaseqTable,geneInfo,method = "geneLength"){

    TCGA_RnaseqTable <- TCGA_RnaseqTable[ !(GenesCutID(as.matrix(rownames(TCGA_RnaseqTable))) == "?"),]
    TCGA_RnaseqTable <- TCGA_RnaseqTable[ !(GenesCutID(as.matrix(rownames(TCGA_RnaseqTable))) == "SLC35E2"),]
    rownames(TCGA_RnaseqTable) <- GenesCutID(as.matrix(rownames(TCGA_RnaseqTable)))
    TCGA_RnaseqTable <- TCGA_RnaseqTable[rownames(TCGA_RnaseqTable) != "?", ]
    TCGA_RnaseqTable <- TCGA_RnaseqTable[!duplicated(rownames(TCGA_RnaseqTable)), !duplicated(colnames(TCGA_RnaseqTable))]
    #TCGA_RnaseqTable <- TCGA_RnaseqTable[, which(substr(colnames(TCGA_RnaseqTable), 14, 15) != "02")]
    TCGA_RnaseqTable <- TCGA_RnaseqTable[rownames(TCGA_RnaseqTable) %in% rownames(geneInfo),]
    TCGA_RnaseqTable <- as.matrix(TCGA_RnaseqTable)

    if(method == "gcContent"){
        rawCounts<- TCGA_RnaseqTable
    wwhich <- which(rownames(geneInfo) %in% rownames(rawCounts))
    geneInfo <- geneInfo[wwhich,]
    geneInfo <- geneInfo[rownames(rawCounts),]

    timeEstimated <- format(ncol(TCGA_RnaseqTable)*nrow(TCGA_RnaseqTable)/80000,digits = 2)
    print(messageEstimation <- paste("I Need about ", timeEstimated,
                                     "seconds for this Complete Normalization Upper Quantile",
                                     " [Processing 80k elements /s]  "))

    ffData  <- as.data.frame(geneInfo)
    rawCounts <- floor(rawCounts)
    print("Step 1 of 4: newSeqExpressionSet ...")
    tmp <- newSeqExpressionSet(rawCounts, featureData = ffData)

    #fData(tmp)[, "gcContent"] <- as.numeric(geneInfo[, "gcContent"])

    print("Step 2 of 4: withinLaneNormalization ...")
    tmp <- withinLaneNormalization(tmp, "gcContent", which = "upper", offset = TRUE)
    print("Step 3 of 4: betweenLaneNormalization ...")
    tmp <- betweenLaneNormalization(tmp, which = "upper", offset = TRUE)
    normCounts <-  log(rawCounts + .1) + offst(tmp)
    normCounts <-  floor(exp(normCounts) - .1)
    print("Step 4 of 4: .quantileNormalization ...")
    tmp <- t(.quantileNormalization(t(normCounts)))
    TCGA_RnaseqTable_norm <- floor(tmp)
    }

    if(method == "geneLength"){
    geneInfo <- geneInfo[rownames(geneInfo) %in% rownames(TCGA_RnaseqTable), ]
    geneInfo <- geneInfo[!duplicated(rownames(geneInfo)), ]
    toKeep <- which(geneInfo[, "geneLength"] != 0)
    geneInfo <- geneInfo[toKeep, ]
    TCGA_RnaseqTable <- TCGA_RnaseqTable[toKeep, ]
    geneInfo <- as.data.frame(geneInfo)
    TCGA_RnaseqTable <- round(TCGA_RnaseqTable)
    commonGenes <- intersect(rownames(TCGA_RnaseqTable),rownames(geneInfo))

    TCGA_RnaseqTable <- TCGA_RnaseqTable[commonGenes,]
    geneInfo <- geneInfo[commonGenes,]

    timeEstimated <- format(ncol(TCGA_RnaseqTable)*nrow(TCGA_RnaseqTable)/80000,digits = 2)
    print(messageEstimation <- paste("I Need about ", timeEstimated,
                                     "seconds for this Complete Normalization Upper Quantile",
                                     " [Processing 80k elements /s]  "))

    print("Step 1 of 4: newSeqExpressionSet ...")
    system.time(TCGA_RnaseqTable_norm <- EDASeq::newSeqExpressionSet(TCGA_RnaseqTable, featureData = geneInfo))
    print("Step 2 of 4: withinLaneNormalization ...")
    system.time(TCGA_RnaseqTable_norm <- EDASeq::withinLaneNormalization(TCGA_RnaseqTable_norm, "geneLength", which = "upper", offset = FALSE))
    print("Step 3 of 4: betweenLaneNormalization ...")
    system.time(TCGA_RnaseqTable_norm <- EDASeq::betweenLaneNormalization(TCGA_RnaseqTable_norm, which = "upper", offset = FALSE))
    print("Step 4 of 4: exprs ...")

    #system.time(TCGA_RnaseqTable_norm <- EDASeq::exprs(TCGA_RnaseqTable_norm))
    system.time(TCGA_RnaseqTable_norm <- EDASeq::counts(TCGA_RnaseqTable_norm))
    }


    return(TCGA_RnaseqTable_norm)
}

#' @title Differentially expression analysis (DEA) using edgeR package.
#' @description
#'    TCGAanalyze_DEA allows user to perform Differentially expression analysis (DEA),
#'    using edgeR package to identify differentially expressed genes (DEGs).
#'     It is possible to do a two-class analysis.
#'
#'     TCGAanalyze_DEA performs DEA using following functions from edgeR:
#'     \enumerate{
#'     \item edgeR::DGEList converts the count matrix into an edgeR object.
#'     \item edgeR::estimateCommonDisp each gene gets assigned the same dispersion estimate.
#'     \item edgeR::exactTest performs pair-wise tests for differential expression between two groups.
#'     \item edgeR::topTags takes the output from exactTest(), adjusts the raw p-values using the
#'     False Discovery Rate (FDR) correction, and returns the top differentially expressed genes.
#'     }
#' @param mat1 numeric matrix, each row represents a gene,
#' each column represents a sample with Cond1type
#' @param mat2 numeric matrix, each row represents a gene,
#' each column represents a sample with Cond2type
#' @param Cond1type a string containing the class label of the samples in mat1
#'  (e.g., control group)
#' @param Cond2type a string containing the class label of the samples in mat2
#' (e.g., case group)
#' @param method is 'glmLRT' (1) or 'exactTest' (2).
#' (1) Fit a negative binomial generalized log-linear model to
#' the read counts for each gene
#' (2) Compute genewise exact tests for differences in the means between
#' two groups of negative-binomially distributed counts.
#' @importFrom edgeR DGEList estimateCommonDisp exactTest topTags estimateGLMCommonDisp
#' estimateGLMTagwiseDisp glmFit glmLRT
#' @export
#' @examples
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
#' dataFilt <- TCGAanalyze_Filtering(dataNorm, 0.25)
#' samplesNT <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
#' samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))
#' dataDEGs <- TCGAanalyze_DEA(dataFilt[,samplesNT],
#'                       dataFilt[,samplesTP],"Normal", "Tumor")
#' @return table with DEGs containing for each gene logFC, logCPM, pValue,and FDR
TCGAanalyze_DEA <- function(mat1,mat2,Cond1type,Cond2type,method = "exactTest") {

    TOC <- cbind(mat1,mat2)
    Cond1num <- ncol(mat1)
    Cond2num <- ncol(mat2)

    print(message1 <- paste( "there are Cond1 type", Cond1type ,"in ",
                             Cond1num, "samples"))
    print(message2 <- paste( "there are Cond2 type", Cond2type ,"in ",
                             Cond2num, "samples"))
    print(message3 <- paste( "there are ", nrow(TOC) ,
                             "features as miRNA or genes "))

    timeEstimated <- format(ncol(TOC)*nrow(TOC)/30000,digits = 2)
    print(messageEstimation <- paste("I Need about ", timeEstimated,
                                     "seconds for this DEA. [Processing 30k elements /s]  "))

    # Reading in the data and creating a DGEList object
    colnames(TOC) <- paste0('s',1:ncol(TOC))
    #DGE <- DGEList(TOC,group=rep(c("Normal","Tumor"),c(NormalSample,
    #TumorSample)))

    if (method == "exactTest"){
    DGE <- edgeR::DGEList(TOC,group = rep(c(Cond1type,Cond2type),
                                          c(Cond1num,Cond2num)))
    # Analysis using common dispersion
    disp <- edgeR::estimateCommonDisp(DGE) # Estimating the common dispersion
    #tested <- exactTest(disp,pair=c("Normal","Tumor")) # Testing
    tested <- edgeR::exactTest(disp,pair = c(Cond1type,Cond2type)) # Testing
    # Results visualization
    logFC_table <- tested$table
    tableDEA <- edgeR::topTags(tested,n = nrow(tested$table))$table
    }

    if (method == "glmLRT"){
        tumorType <- rep(c(Cond1type,Cond2type),
                         c(Cond1num,Cond2num))
    design <- model.matrix(~as.factor(tumorType))
    aDGEList <- edgeR::DGEList(counts = TOC, group = as.factor(tumorType))
    aDGEList <- edgeR::estimateGLMCommonDisp(aDGEList, design)
    aDGEList <- edgeR::estimateGLMTagwiseDisp(aDGEList, design)
    aGlmFit <- edgeR::glmFit(aDGEList, design, dispersion = aDGEList$tagwise.dispersion,
                             prior.count.total=0)
    aGlmLRT <- edgeR::glmLRT(aGlmFit, coef = 2)
    tableDEA <- aGlmLRT
    }

    return(tableDEA)

}

#' @title Adding information related to DEGs genes from DEA as mean values in two conditions.
#' @description
#'    TCGAanalyze_LevelTab allows user to add information related to DEGs genes from
#'    Differentially expression analysis (DEA) such as mean values and in two conditions.
#' @param FC_FDR_table_mRNA Output of dataDEGs filter by abs(LogFC) >=1
#' @param typeCond1 a string containing the class label of the samples
#'  in TableCond1  (e.g., control group)
#' @param typeCond2 a string containing the class label of the samples
#' in TableCond2  (e.g., case group)
#' @param TableCond1 numeric matrix, each row represents a gene, each column
#'  represents a sample with Cond1type
#' @param TableCond2 numeric matrix, each row represents a gene, each column
#' represents a sample with Cond2type
#' @param typeOrder typeOrder
#' @importFrom edgeR DGEList estimateCommonDisp exactTest topTags
#' @export
#' @return table with DEGs, log Fold Change (FC), false discovery rate (FDR),
#' the gene expression level
#' for samples in  Cond1type, and Cond2type, and Delta value (the difference
#' of gene expression between the two
#' conditions multiplied logFC)
#' @examples
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
#' dataFilt <- TCGAanalyze_Filtering(dataNorm, 0.25)
#' samplesNT <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
#' samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))
#' dataDEGs <- TCGAanalyze_DEA(dataFilt[,samplesNT], dataFilt[,samplesTP],
#' "Normal", "Tumor")
#' dataDEGsFilt <- dataDEGs[abs(dataDEGs$logFC) >= 1,]
#' dataTP <- dataFilt[,samplesTP]
#' dataTN <- dataFilt[,samplesNT]
#' dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGsFilt,"Tumor","Normal",
#' dataTP,dataTN)
TCGAanalyze_LevelTab <- function(FC_FDR_table_mRNA,typeCond1,typeCond2,
                                 TableCond1,TableCond2,typeOrder = TRUE) {

    TF_enriched <- as.matrix(rownames(FC_FDR_table_mRNA))
    TableLevel <- matrix(0,nrow(TF_enriched),6)
    TableLevel <- as.data.frame(TableLevel)

    colnames(TableLevel) <- c("mRNA","logFC","FDR",typeCond1,typeCond2,"Delta")


    TableLevel[,"mRNA"] <- TF_enriched
    Tabfilt <- FC_FDR_table_mRNA[which( rownames(FC_FDR_table_mRNA) %in%
                                            TF_enriched),]
    TableLevel[,"logFC"] <- as.numeric(Tabfilt[TF_enriched,][,"logFC"])
    TableLevel[,"FDR"] <- as.numeric(Tabfilt[TF_enriched,][,"FDR"])


    MeanTumor <- matrix(0,nrow(TF_enriched),1)
    MeanDiffTumorNormal <- matrix(0,nrow(TF_enriched),1)


    for (i in 1:nrow(TF_enriched)) {
        #print(paste(i, "of", nrow(TF_enriched),TF_enriched[i]))
        TableLevel[i,typeCond1] <- mean(TableCond1[rownames(TableCond1) %in%
                                                       TF_enriched[i] , ])
        TableLevel[i,typeCond2] <- mean(TableCond2[rownames(TableCond2) %in%
                                                       TF_enriched[i] , ])
    }


    TableLevel[,"Delta"] <- as.numeric(abs(TableLevel[,"logFC"]) *
                                           TableLevel[,typeCond1]  )

    TableLevel <- TableLevel[order( as.numeric(TableLevel[,"Delta"]),
                                    decreasing = typeOrder),]

    rownames(TableLevel) <-  TableLevel[,"mRNA"]
    return(TableLevel)
}

#' @title Enrichment analysis for Gene Ontology (GO) [BP,MF,CC] and Pathways
#' @description
#'   Researchers, in order to better understand the underlying biological
#'   processes, often want to retrieve a functional profile of a set of genes
#'   that might have an important role. This can be done by performing an
#'   enrichment analysis.
#'
#'We will perform an enrichment analysis on gene sets using the TCGAanalyze_EAcomplete
#'function. Given a set of genes that are
#'up-regulated under certain conditions, an enrichment analysis will find
#'identify classes of genes or proteins that are #'over-represented using
#'annotations for that gene set.
#' @param TFname is the name of the list of genes or TF's regulon.
#' @param RegulonList List of genes such as TF's regulon or DEGs where to find enrichment.
#' @export
#' @return Enrichment analysis GO[BP,MF,CC] and Pathways complete table enriched by genelist.
#' @examples
#' Genelist <- c("FN1","COL1A1")
#' ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist)
#' \dontrun{
#' Genelist <- rownames(dataDEGsFiltLevel)
#' system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))
#' }
TCGAanalyze_EAcomplete <- function(TFname, RegulonList){

    print(paste("I need about ", "1 minute to finish complete ",
                "Enrichment analysis GO[BP,MF,CC] and Pathways... "))

    ResBP <- TCGAanalyze_EA(TFname,RegulonList,DAVID_BP_matrix,
                            EAGenes,GOtype = "DavidBP")
    print("GO Enrichment Analysis BP completed....done")
    ResMF <- TCGAanalyze_EA(TFname,RegulonList,DAVID_MF_matrix,
                            EAGenes,GOtype = "DavidMF")
    print("GO Enrichment Analysis MF completed....done")
    ResCC <- TCGAanalyze_EA(TFname,RegulonList,DAVID_CC_matrix,
                            EAGenes,GOtype = "DavidCC")
    print("GO Enrichment Analysis CC completed....done")
    ResPat <- TCGAanalyze_EA(TFname,RegulonList,listEA_pathways,
                             EAGenes,GOtype = "Pathway")
    print("Pathway Enrichment Analysis completed....done")

    ans <- list(ResBP = ResBP, ResMF = ResMF, ResCC = ResCC, ResPat = ResPat)
    return(ans)
}

#' @title Enrichment analysis of a gene-set with GO [BP,MF,CC]  and pathways.
#' @description
#' The rational behind a enrichment analysis ( gene-set, pathway etc) is to compute
#' statistics of whether the overlap between the focus list (signature) and the gene-set
#' is significant. ie the confidence that overlap between the list is not due to chance.
#'  The Gene Ontology project describes genes (gene products) using terms from
#'  three structured vocabularies: biological process, cellular component and molecular function.
#'  The Gene Ontology Enrichment component, also referred to as the GO Terms" component, allows
#'  the genes in any such "changed-gene" list to be characterized using the Gene Ontology terms
#'  annotated to them. It asks, whether for any particular GO term, the fraction of genes
#'  assigned to it in the "changed-gene" list is higher than expected by chance
#'  (is over-represented), relative to the fraction of genes assigned to that term in the
#'  reference set.
#'  In statistical terms it peform the analysis tests the null hypothesis that,
#'  for any particular ontology term, there is no diffeerence in the proportion of genes
#'  annotated to it in the reference list and the proportion annotated to it in the test list.
#'  We adopted a Fisher Exact Test to perform the EA.
#' @param GeneName is the name of gene signatures list
#' @param TableEnrichment is a table related to annotations of gene symbols such as
#' GO[BP,MF,CC] and Pathways. It was created from DAVID gene ontology on-line.
#' @param RegulonList is a gene signature (lisf of genes) in which perform EA.
#' @param GOtype is type of gene ontology Biological process (BP), Molecular Function (MF),
#' Cellular componet (CC)
#' @param FDRThresh pvalue corrected (FDR) as threshold to selected significant
#' BP, MF,CC, or pathways. (default FDR < 0.01)
#' @param EAGenes is a table with informations about genes
#' such as ID, Gene, Description, Location and Family.
# @export
#' @import stats
#' @return Table with enriched GO or pathways by selected gene signature.
#' @examples
#' \dontrun{
#' EAGenes <- get("EAGenes")
#' RegulonList <- rownames(dataDEGsFiltLevel)
#' ResBP <- TCGAanalyze_EA(GeneName="DEA genes Normal Vs Tumor",
#'                            RegulonList,DAVID_BP_matrix,
#'                            EAGenes,GOtype = "DavidBP")
#'}
TCGAanalyze_EA <- function(GeneName,RegulonList,TableEnrichment,
                           EAGenes,GOtype,FDRThresh=0.01) {
    topPathways <- nrow(TableEnrichment)
    topPathways_tab <- matrix(0,1,topPathways)
    topPathways_tab <- as.matrix(topPathways_tab)
    rownames(topPathways_tab) <- GeneName

    rownames(EAGenes) <- toupper(rownames(EAGenes) )
    EAGenes <- EAGenes[!duplicated(EAGenes[,"ID"]),]
    rownames(EAGenes) <- EAGenes[,"ID"]
    allgene <- EAGenes[,"ID"]
    current_pathway_from_EA <- as.matrix(TableEnrichment[,GOtype]) # genes from EA pathways

    TableNames <- gsub("David","",paste("Top ", GOtype, " n. ", 1:topPathways,
                                        " of ", topPathways, sep = ""))
    colnames(topPathways_tab) <- TableNames
    topPathways_tab <- as.data.frame(topPathways_tab)

    table_pathway_enriched <- matrix(1, nrow(current_pathway_from_EA),7)
    colnames(table_pathway_enriched) <- c("Pathway","GenesInPathway","Pvalue",
                                          "FDR","CommonGenesPathway",
                                          "PercentPathway","PercentRegulon")
    table_pathway_enriched <- as.data.frame(table_pathway_enriched)

    for (i in 1:nrow(current_pathway_from_EA)) {
        table_pathway_enriched[i,"Pathway"] <- as.character(current_pathway_from_EA[i,])

        if (nrow(TableEnrichment) == 589) {
            genes_from_current_pathway_from_EA <- GeneSplitRegulon(TableEnrichment[ TableEnrichment[GOtype] == as.character(current_pathway_from_EA[i,]) ,][,"Molecules"], ",")
        }
        else {
            genes_from_current_pathway_from_EA <- GeneSplitRegulon(TableEnrichment[ TableEnrichment[GOtype] == as.character(current_pathway_from_EA[i,]) ,][,"Molecules"], ", ")
        }

        genes_common_pathway_TFregulon <- as.matrix(intersect(toupper(RegulonList),toupper(genes_from_current_pathway_from_EA)))



        if (length(genes_common_pathway_TFregulon) != 0) {
            current_pathway_commongenes_num <- length(genes_common_pathway_TFregulon)
            seta <-  allgene %in% RegulonList
            setb <-  allgene %in% genes_from_current_pathway_from_EA
            ft <- fisher.test(seta,setb)
            FisherpvalueTF <- ft$p.value
            table_pathway_enriched[i,"Pvalue"] <- as.numeric(FisherpvalueTF)
            if (FisherpvalueTF < 0.01) {
                current_pathway_commongenes_percent <- paste("(",format( (current_pathway_commongenes_num/length(genes_from_current_pathway_from_EA)) * 100,digits = 2),"%)")
                current_pathway_commongenes_num_with_percent <- gsub(" ","",paste(current_pathway_commongenes_num, current_pathway_commongenes_percent,"pv=",format(FisherpvalueTF,digits=2)))
                table_pathway_enriched[i,"CommonGenesPathway"] <- length(genes_common_pathway_TFregulon)
                table_pathway_enriched[i,"GenesInPathway"] <- length(genes_from_current_pathway_from_EA)
                table_pathway_enriched[i,"PercentPathway"] <-  as.numeric(table_pathway_enriched[i,"CommonGenesPathway"]) / as.numeric(table_pathway_enriched[i,"GenesInPathway"])  *100
                table_pathway_enriched[i,"PercentRegulon"] <-  as.numeric(table_pathway_enriched[i,"CommonGenesPathway"]) / length(RegulonList)  *100
            } }
    }
    table_pathway_enriched <- table_pathway_enriched[order(table_pathway_enriched[,"Pvalue"],decreasing = FALSE),]
    table_pathway_enriched <- table_pathway_enriched[table_pathway_enriched[,"Pvalue"] < 0.01 ,]
    table_pathway_enriched[,"FDR"] <- p.adjust(table_pathway_enriched[,"Pvalue"],method = "fdr")
    table_pathway_enriched <- table_pathway_enriched[table_pathway_enriched[,"FDR"] < FDRThresh ,]
    table_pathway_enriched <- table_pathway_enriched[order(table_pathway_enriched[,"FDR"],decreasing = FALSE),]

    tmp <- table_pathway_enriched[1:topPathways,]
    tmp <- paste(tmp[,"Pathway"],"; FDR= ", format(tmp[,"FDR"],digits = 3),"; (ng="   ,round(tmp[,"GenesInPathway"]),"); (ncommon=", format(tmp[,"CommonGenesPathway"],digits = 2), ")" ,sep = "")
    tmp <- as.matrix(tmp)
    topPathways_tab[1,] <- tmp
    rm(tmp)

    return(topPathways_tab)
}

