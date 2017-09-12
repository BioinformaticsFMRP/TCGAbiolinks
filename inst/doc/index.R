## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, eval=FALSE----------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("TCGAbiolinks")

## ----message=FALSE, warning=FALSE, eval=FALSE----------------------------
#  devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')

## ----message=FALSE, warning=FALSE, include=TRUE--------------------------
library(TCGAbiolinks)
library(dplyr)
library(DT)

## ------------------------------------------------------------------------
version
packageVersion("TCGAbiolinks")

