#' @title SurvivalKMtable
#' @description SurvivalKMtable perform an univariate Kaplan-Meier (KM) survival analysis (SA).
#' It performed Kaplan-Meier survival univariate using complte follow up with all days
#' taking one gene a time from Genelist of gene symbols.
#' For each gene according its level of mean expression in cancer samples, defining two thresholds for quantile
#' expression of that gene in all samples (default ThreshTop=0.67,ThreshDown=0.33) it is possible
#' to define a threshold of intensity of gene expression to divide the samples in 3 groups
#' (High, intermediate, low).
#' SurvivalKMtable performs SA between High and low groups using following functions
#' from survival package
#'    1. survival::Surv
#'    2. survival::survdiff
#'    3. survival::survfit
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
#' clinical_patient_Cancer <- clinic("brca","clinical_patient")
#' dataBRCAcomplete <- log2(BRCA_rnaseqv2)
#' # Selecting only 10 genes for example
#' dataBRCAcomplete <- dataBRCAcomplete[1:10,]
#' tabSurvKM<-SurvivalKMunivariate(clinical_patient_Cancer,dataBRCAcomplete,
#' Genelist = rownames(dataBRCAcomplete), Survresult = FALSE,ThreshTop=0.67,ThreshDown=0.33)
#' # Filtering by pvalue < 0.01
#' tabSurvKM <- tabSurvKM[tabSurvKM$pvalue < 0.01,]
#' tabSurvKM <- tabSurvKM[!duplicated(tabSurvKM$mRNA),]
#' rownames(tabSurvKM) <-tabSurvKM$mRNA
#' tabSurvKM <- tabSurvKM[,-1]
#' tabSurvKM <- tabSurvKM[order(tabSurvKM$pvalue, decreasing=FALSE),]
SurvivalKMunivariate<-function(clinical_patient,dataGE,Genelist, Survresult,ThreshTop=0.67, ThreshDown=0.33){

    samplesNT <- MultiSampleTypes(colnames(dataGE), typesample = c("NT"))
    samplesTP <- MultiSampleTypes(colnames(dataGE), typesample = c("TP"))

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


