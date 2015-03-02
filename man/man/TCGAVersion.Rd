\name{TCGAVersion}
\alias{TCGAVersion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TCGAVersion
}
\description{
TCGAVersion downloads automatically the latest version of data samples from TCGA specified by tumor and platform type.
}
\usage{
TCGAVersion(Tumor, PlatformType, PlatformAndAssociatedData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tumor}{
a character string indicating the cancer type for which to download data. Options include ACC, BLCA, BRCA, CESC, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML, LGG, LIHC, LUAD, LUSC, OV, PAAD, PRAD, READ, SARC, SKCM, STAD, THCA, UCEC, UCS. Look at https://tcga-data.nci.nih.gov/tcga/ for Available Cancer Types. 
}
  \item{PlatformType}{
character containing platform type from which you want to get the data
}
  \item{PlatformAndAssociatedData}{
data frame 615 observations of 12 variables, indicating the different characteristics of the data e.g. tumour, type, species.  
}
}
\details{
}
\value{
TCGAVersion returns a matrix with Samples in columns with barcode and probeID in rows.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
