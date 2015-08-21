#' @title Differentially expression analysis (DEA) using limma package.
#' @description Differentially expression analysis (DEA) using limma package.
#' @param Pdatatable
#' @param PathFolder
#' @param FC.cut
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma makeContrasts
#' @importFrom limma contrasts.fit
#' @importFrom limma toptable
#' @examples
#' \dontrun{
#' to add example
#' }
#' @export
#' @return List of list with tables in 2 by 2 comparison
#' of the top-ranked genes from a linear model fitted by DEA's limma
TCGAanalyze_DEA_Affy <- function(Pdatatable, PathFolder, FC.cut = 0.01){

    f <- factor(Pdatatable$Disease)
    groupColors<-names(table(f))

    tmp<-matrix(0,length(groupColors),length(groupColors))
    colnames(tmp) <- groupColors
    rownames(tmp) <- groupColors
    tmp[upper.tri(tmp)] <- 1


    sample_tab<-Pdatatable
    f <- factor(Pdatatable$Disease)
    design <- model.matrix(~0+f)
    colnames(design) <- levels(f)
    fit <- lmFit(AffySet, design) ## fit is an object of class MArrayLM.

    groupColors <- names(table(Pdatatable$Disease))

    CompleteList<-vector("list",sum(tmp))

    k<-1

    for( i in 1: length(groupColors)){
        col1 <- colnames(tmp)[i]
        for( j in 1: length(groupColors)){
            col2 <- rownames(tmp)[j]

            if( i!=j ){

                if(tmp[i,j]!=0){


                    Comparison <- paste(col2,"-",col1,sep="")

                    if(i==4 && j==6){ Comparison <- paste(col1,"-",col2,sep="") }
                    if(i==5 && j==6){ Comparison <- paste(col1,"-",col2,sep="") }

                    print( paste(i, j, Comparison,"to do..." ))

                    cont.matrix <- makeContrasts(I=Comparison,levels=design)

                    fit2 <- contrasts.fit(fit, cont.matrix)
                    fit2 <- eBayes(fit2)



                    sigI <- topTable(fit2,coef=1, adjust.method="BH", sort.by="B", p.value = 0.05, lfc = FC.cut, number = 50000)

                    sigIbis <- sigI[order(abs(as.numeric(sigI$logFC)), decreasing=T),]
                    names(CompleteList)[k]<-gsub("-","_",Comparison)
                    CompleteList[[k]]<-sigIbis
                    k<-k+1
                }
            }
        }

    }

    return(CompleteList)
    }
