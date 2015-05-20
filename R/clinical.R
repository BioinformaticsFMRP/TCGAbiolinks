#' @title SampleTypes
#' @description
#'   SampleTypes
#' @param barcode barcode
#' @param typesample typesample
#' @export
#' @return SampleTypes
SampleTypes <- function(barcode, typesample){
  table.code <- c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','20','40','50','60','61')
  names(table.code) <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC","TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB","CELL","XP","XCL")

  if(is.element(typesample,names(table.code))){
    barcode <- barcode[grep(table.code[typesample], substr(barcode, 14, 15))]
  }else{
    return("Error message: sample type does not exist")
  }
  return(barcode)
}

#' @title MultiSampleTypes
#' @description
#'   MultiSampleTypes
#' @param barcode barcode
#' @param typesample typesample
#' @export
#' @return MultiSampleTypes
MultiSampleTypes <- function(barcode,typesample){
  # Tumor AND Solid Tissue Normal NOT FROM THE SAME PATIENTS
  table.code <- c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','20','40','50','60','61')
  names(table.code) <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC","TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB","CELL","XP","XCL")

  if(sum(is.element(typesample,names(table.code))) == length(typesample)) {

    string <- substr(barcode, 14, 15)
    barcode.all <- NULL
    for(sample.i in typesample){
      barcode.all <- union(barcode.all, barcode[grep(table.code[sample.i], string)])
    }
    return(barcode.all)
  }else{
    return("Error message: one or more sample types do not exist")
  }
}

#' @title MatchedCoupledSampleTypes
#' @description
#'   MatchedCoupledSampleTypes
#' @param barcode barcode
#' @param typesample typesample
#' @export
#' @return MultiSampleTypes
MatchedCoupledSampleTypes <- function(barcode,typesample){
  # Tumor AND Solid Tissue Normal FROM THE SAME PATIENTS
  table.code <- c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','20','40','50','60','61')
  names(table.code) <- c("TP","TR","TB","TRBM","TAP","TM","TAM","THOC","TBM","NB","NT","NBC","NEBV","NBM","CELLC","TRB","CELL","XP","XCL")
  if(length(typesample)!=2){
    return("Error message: exactly two types need to be provided")
  }

  if(sum(is.element(typesample,names(table.code))) == length(typesample)) {

    string <- substr(barcode, 14, 15)
    barcode.1 <- barcode[grep(table.code[typesample[1]], string)]
    barcode.2 <- barcode[grep(table.code[typesample[2]], string)]


    barcode.common <- intersect(substr(barcode.1,1,13), substr(barcode.2,1,13))
    if(length(barcode.common) > 0){
      return(union(barcode.1[grep(barcode.common,barcode.1)],barcode.2[grep(barcode.common,barcode.2)]))
    }else{
      return("Error message: there exist no matched samples")
    }
  } else{
    return("Error message: one or more sample types do not exist")
  }
}

#' @title stage_BRCA
#' @description
#'   stage_BRCA
#' @param barcode barcode
#' @param stage stage
#' @param clinical_patient_data clinical_patient_data
#' @export
#' @return stage_BRCA
stage_BRCA <- function(barcode, stage, clinical_patient_data){
  table.stages <- c("Stage I$|Stage IA$|Stage IB$", "Stage I$", "Stage IA$", "Stage IB$", "Stage II$|Stage IIA$|Stage IIB$", "Stage II$",
                    "Stage IIA$", "Stage IIB$", "Stage III$|Stage IIIA$|Stage IIIB$|Stage IIIC$", "Stage III$", "Stage IIIA$", "Stage IIIB$", "Stage IIIC$", "Stage IV$")
  names(table.stages) <- c("stage_IX", "stage_I", "stage_IA", "stage_IB", "stage_IIX", "stage_IIA", "stage_IIB", "stage_IIIX", "stage_IIIA",
                           "stage_IIIB", "stage_IIIC", "stage_IV")

  if(is.element(stage, names(table.stages))) {
    clinical_patient_data<-as.data.frame(clinical_patient_data)
    stage.i <- clinical_patient_data[grep(table.stages[stage], clinical_patient_data$ajcc_pathologic_tumor_stage), ]
    stage.i <- stage.i[,"bcr_patient_barcode"]
    samples <- substr(barcode, 1, 12)
    barcode <- intersect(samples,stage.i)
  }else{
    return("Error message: stage does not exist")
  }
  return(barcode)
}

#' @title gender_BRCA
#' @description
#'   gender_BRCA
#' @param barcode barcode
#' @param gender gender
#' @param clinical_patient_data clinical_patient_data
#' @export
#' @return stage_BRCA
gender_BRCA <- function(barcode, gender, clinical_patient_data){

  if(is.element(gender,c("MALE", "FEMALE"))){
    clinical_patient_data<-as.data.frame(clinical_patient_data)
    s.gender <- as.data.frame(clinical_patient_data)[grep(paste0("^", gender,"$"), clinical_patient_data$gender), ][,"bcr_patient_barcode"]
    samples<-substr(barcode, 1, 12)
    #find common patients between FEMALE and barcode data
    barcode<-intersect(samples,s.gender)
  }else{
    return("Error message gender doesn't exist")
  }

  return(barcode)
}

#' @title ER_status_BRCA
#' @description
#'   ER_status_BRCA
#' @param barcode barcode
#' @param ER ER
#' @param clinical_patient_data clinical_patient_data
#' @export
#' @return ER_status_BRCA
ER_status_BRCA <- function(barcode,ER, clinical_patient_data){
  ## ER should be "Positive" or "Negative"
  # consider only barcode and ER status
  if(is.element(ER, c("Positive", "Negative"))){
    status <- as.data.frame(clinical_patient_data)[grep(paste0("^",ER,"$"), clinical_patient_data$er_status_by_ihc), ][,"bcr_patient_barcode"]
    samples<-substr(barcode, 1, 12)
    #find common patients between ER status and barcode data
    barcode<-intersect(samples,status)
    return(barcode)
  }else{
    return("Error message: ER status does not exist")
  }
}

#' @title PR_status_BRCA
#' @description
#'   PR_status_BRCA
#' @param barcode barcode
#' @param PR PR
#' @param clinical_patient_data clinical_patient_data
#' @export
#' @return PR_status_BRCA
PR_status_BRCA  <- function(barcode,PR, clinical_patient_data){
  ## PR should be "Positive" or "Negative"

  if(is.element(PR, c("Positive", "Negative"))){
    #for breast cancer
    status <- as.data.frame(clinical_patient_data)[grep(paste0("^", PR, "$"), clinical_patient_data$pr_status_by_ihc), ][,"bcr_patient_barcode"]
    samples<-substr(barcode, 1, 12)
    #find common patients between PR status and barcode data
    barcode<-intersect(samples,status)
  }else{
    return("Error message: PR status does not exist")
  }

  return(barcode)

}

#' @title HER_status_BRCA
#' @description
#'   HER_status_BRCA
#' @param barcode barcode
#' @param HER HER
#' @param clinical_patient_data clinical_patient_data
#' @export
#' @return HER_status_BRCA
HER_status_BRCA  <- function(barcode, HER, clinical_patient_data){
  if(is.element(HER, c("Positive", "Negative"))){
    clinical_patient_data<-as.data.frame(clinical_patient_data)
    #for breast cancer HER+
    status <- as.data.frame(clinical_patient_data)[grep(paste0("^",HER,"$"), clinical_patient_data$her2_status_by_ihc), ][,"bcr_patient_barcode"]
    samples<-substr(barcode, 1, 12)
    #find common patients between HER+ e barcode data
    barcode<-intersect(samples,status)
    return(barcode)
  }else{
    return("Error message: HER status does not exist")
  }


}

#' @title clinical_data_site_cancer
#' @description
#'   clinical_data_site_cancer
#' @param cancer cancer
#' @export
#' @return clinical_data_site_cancer
clinical_data_site_cancer <-function(cancer){
  return(paste0("https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/", cancer,"/bcr/biotab/clin/"))
}

#' @title clinic
#' @description
#'   clinic
#' @param cancer cancer
#' @param clinical_data_type clinical_data_type
#' @export
#' @importFrom RCurl getURL
#' @return clinic
clinic <- function(cancer,clinical_data_type){

  URL <- paste0(clinical_data_site_cancer (cancer), "nationwidechildrens.org_", clinical_data_type, "_", cancer, ".txt")
  writeLines(getURL(URL,ssl.verifypeer = FALSE), file(paste0(clinical_data_type,".txt")))

  return(file(paste0(clinical_data_type,".txt")))
}
