\name{TCGADownload}
\alias{TCGADownload}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TCGADownload
}
\description{

Download data of samples from TCGA specified by tumor,type, species and platform type.

}
\usage{
TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder = "",
            PlatformType, nsample = 0, listSample = 0)
}

\arguments{
  \item{Tumor}{ a character string indicating the cancer type for which to download data. Options include ACC, BLCA, BRCA, CESC, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML, LGG, LIHC, LUAD, LUSC, OV, PAAD, PRAD, READ, SARC, SKCM, STAD, THCA, UCEC, UCS. Look at https://tcga-data.nci.nih.gov/tcga/ for Available Cancer Types. 

}
  \item{Type}{
character containing biological data type to analyze
}
  \item{Species}{
a character string containing the assay platform from which you want to get the data
}
  \item{PlatformAndAssociatedData}{
data frame 615 observations of 12 variables, indicating the different characteristics of the data e.g. tumour, type, species.  
}
  \item{downloadFolder}{
path of the directory to store the acquired data matrix files
}
  \item{PlatformType}{
character containing platform type from which you want to get the data
}
  \item{nsample}{
number of samples to be analyzed, default all samples
}
  \item{listSample}{
  a character that indicates TCGA barcodes of the patients/samples to analyze
}
}
\details{
TCGA download retrieves and stores data of samples belonging to the specified cancer type and measured by the specified assay platform.
}
\value{
TCGADownload returns a matrix with Samples in columns with barcode and probeID in rows.
ProbeID and values changes with platforms.

\item{humanmethylation27}{If platform is humanmethylation27 returns values as Hybridization.RE}
\item{illuminahiseq_rnaseq}{If platform is illuminahiseq_rnaseq returns values as raw_counts}
\item{agilentg4502a_07_3}{If platform is agilentg4502a_07_3 returns values as Hybridization.REF}
\item{illuminahiseq_rnaseqv2}{If platform is illuminahiseq_rnaseqv2 returns values as raw_counts}
\item{humanmethylation450}{If platform is humanmethylation450 returns values as Hybridization.REF}
\item{illuminaga_mirnaseq}{If platform is illuminaga_mirnaseq returns values as read_count}
\item{illuminahiseq_mirnaseq}{If platform is illuminahiseq_mirnaseq returns values as read_count}
\item{genome_wide_snp_6}{If platform is genome_wide_snp_6 returns values as Segment_Mean}
\item{mda_rppa_core}{If platform is mda_rppa_core returns values as Protein Expression}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Antonio Colaprico, Luciano Garofano, Claudia Cava, Gianluca Bontempi, Michele Ceccarelli
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


library(TCGAbiolinks)
## it downloads illuminahiseq_rnaseq data of breast cancer samples and stores in downloadFolder

Tumor <- "BRCA"
downloadFolder <- "~/Data"
setwd(downloadFolder)
Type = "mRNA"
Species <- "RNASeq"   
PlatformType <- "illuminahiseq_rnaseq"
BRCA_rnaseq <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)


## it downloads illuminahiseq_rnaseqv2 data of breast cancer samples and stores in downloadFolder

Type = "mRNA"
Species <- "RNASeqV2" 
PlatformType <- "illuminahiseq_rnaseqv2" 
BRCA_rnaseqv2 <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)

## it downloads agilentg4502a_07_3 data of breast cancer samples and stores in downloadFolder

Type = "mRNA"
Species <- "Exp-Gene" 
PlatformType <- "agilentg4502a_07_3"
BRCA_agilent <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)


## it downloads humanmethylation27 data of breast cancer samples and stores in downloadFolder

Type = "Methylation"
Species <- "Methyl"
PlatformType <- "humanmethylation27" 
BRCA_methylation27 <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)


## it downloads humanmethylation450 data of breast cancer samples and stores in downloadFolder

Type = "Methylation"
Species <- "Methyl"
PlatformType <- "humanmethylation450" 
BRCA_methylation450 <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)


## it downloads illuminaga_mirnaseq data of breast cancer samples and stores in downloadFolder

Type = "miRNA"
Species <- "miRNASeq"
PlatformType <- "illuminaga_mirnaseq"
BRCA_ga_mirnaseq <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)


## it downloads illuminahiseq_mirnaseq data of breast cancer samples and stores in downloadFolder

Type = "miRNA"
Species <- "miRNASeq"
PlatformType <- "illuminahiseq_mirnaseq"
BRCA_hiseq_mirnaseq <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)


## it downloads genome_wide_snp_6 data of breast cancer samples and stores in downloadFolder

Type = "SNP"
Species <- "CNV (SNP Array)"
PlatformType <- "genome_wide_snp_6"  
BRCA_genome_wide_snp_6 <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)


## it downloads illuminahiseq_dnaseqc data of breast cancer samples and stores in downloadFolder

Species <- "CNV (Low Pass DNASeq)"
PlatformType <- "illuminahiseq_dnaseqc"  
BRCA_hiseq_dnaseqc <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)


## it downloads mda_rppa_core data of breast cancer samples and stores in downloadFolder

Type = "Protein"
Species <- "Exp-Protein"
PlatformType <- "mda_rppa_core"
BRCA_protein <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)


## it downloads illuminaga_dnaseq data of breast cancer samples and stores in downloadFolder

Type = "Exome"
Species <- "Somatic Mutation"
PlatformType <- "illuminaga_dnaseq"
BRCA_mutation <- TCGADownload(Tumor, Type, Species, PlatformAndAssociatedData, downloadFolder,PlatformType)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
