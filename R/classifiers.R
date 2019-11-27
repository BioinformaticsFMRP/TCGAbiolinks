#' @title Gliomar classifier
#' @description Classify DNA methylation gliomas using data from
#' https://doi.org/10.1016/j.cell.2015.12.028
#' @param data DNA methylation matrix or Summarized Experiments with samples on columns
#' and probes on the rows
#' @export
#' @return A list of 3 data frames:
#' 1) Sample final classification
#' 2) Each model final classification
#' 3) Each class probability of classification
#' @examples
#' \donrun{
#' query <- GDCquery(project= "TCGA-GBM",
#'                   data.category = "DNA methylation",
#'                   barcode = c("TCGA-06-0122","TCGA-14-1456"),
#'                   platform = "Illumina Human Methylation 27",
#'                   legacy = TRUE)
#' GDCdownload(query)
#' data.hg19 <- GDCprepare(query)
#' classification <- gliomaClassifier(data.hg19)
#'
#' # Comparing reslts
#' TCGAquery_subtype("GBM") %>%
#' dplyr::filter(patient %in% c("TCGA-06-0122","TCGA-14-1456")) %>%
#' dplyr::select("patient","Supervised.DNA.Methylation.Cluster")
#' }
gliomaClassifier <- function(data){

    if (!requireNamespace("TCGAbiolinksGUI.data", quietly = TRUE)) {
        stop("TCGAbiolinksGUI.data package is needed for this function to work. Please install it.",
             call. = FALSE)
    }

    message("Plese cite https://doi.org/10.1016/j.cell.2015.12.028")
    message("Training model described at https://bioconductor.org/packages/TCGAbiolinksGUI.data/")

    if(is(data, "RangedSummarizedExperiment")) {
        met <- assay(data) %>% as.matrix  %>% t
    } else {
        met <- data %>% as.matrix  %>% t
    }


    df.all <- NULL
    models <- c("idh","gcimp","idhwt","idhmut")
    models <- paste("glioma",models,"model",sep = ".")
    data(list = models, package = "TCGAbiolinksGUI.data")
    for(i in models){
        model <- get(i)
        # If it is a Summarized Experiment object

        # keep only probes used in the model
        aux <- met[,colnames(met) %in% colnames(model$trainingData),drop=FALSE]

        # This should not happen!
        if(any(apply(aux,2,function(x) all(is.na(x))))) {
            print("NA columns")
            aux[,apply(aux,2,function(x) all(is.na(x)))] <- 0.5
        }
        if(any(apply(aux,2,function(x) any(is.na(x))))) {
            print("NA values")
            colMedians <- colMedians(aux,na.rm = T)
            x <- which(is.na(aux),arr.ind = T)
            for(l in 1:nrow(x)){
                aux[x[l,1],x[l,2]] <- colMedians[x[l,2]]
            }
        }

        pred <- predict(model, aux)
        pred.prob <- predict(model, aux, type = "prob")
        colnames(pred.prob) <- paste0(i,"_", colnames(pred.prob))
        pred.prob$samples <- rownames(pred.prob)
        df <- data.frame(samples = rownames(aux),
                         groups.classified = pred,
                         stringsAsFactors = FALSE)
        colnames(df)[2] <- paste0(i,"_groups.classified")
        if(is.null(df.all)) {
            df.all <- df
            df.prob <- pred.prob
        } else {
            df.all <- merge(df.all,df, by = "samples")
            df.prob <- merge(df.prob,pred.prob,by = "samples")
        }
    }
    colnames(df.all) <- c("samples",models)
    fctr.cols <- sapply(df.all, is.factor)
    df.all[, fctr.cols] <- sapply(df.all[, fctr.cols], as.character)
    df.all[grep("6|5|4",df.all$glioma.idh.model),c("glioma.gcimp.model","glioma.idhmut.model")]  <- NA
    df.all[grep("3|2|1",df.all$glioma.idh.model),c("glioma.idhwt.model")]  <- NA
    df.all[grep("3",df.all$glioma.idhmut.model),c("glioma.gcimp.model")]  <- "Codel"
    df.all[grep("1",df.all$glioma.idhwt.model),c("glioma.idhwt.model")]  <- "Classic-like"
    df.all[grep("2",df.all$glioma.idhwt.model),c("glioma.idhwt.model")]  <- "Mesenchymal-like"
    df.all[grep("3",df.all$glioma.idhwt.model),c("glioma.idhwt.model")]  <- "PA-like"

    # Final column with results
    df.all$glioma.DNAmethylation.subtype <- NA
    df.all$glioma.DNAmethylation.subtype <- df.all$glioma.idhwt.model
    idx <- which(is.na(df.all$glioma.DNAmethylation.subtype))
    df.all$glioma.DNAmethylation.subtype[idx] <- df.all$glioma.gcimp.model[idx]
    return(list(final.classification = data.frame(Sample = df.all$samples,Final_classification = df.all$glioma.DNAmethylation.subtype),
                model.classifications = df.all,
                model.probabilities = df.prob))
}
