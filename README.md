# TCGAbiolinks

### How do I get set up? ###

* Installing dependencies
```R
install.packages(c("downloader","RCurl","httr","devtools","stringr",
                    "exactRankTests","XML","GGally","parallel","ggplot2",
                    "survival","biomaRt","rjson","ggbiplot","dnet","igraph",
                    "rvest"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges","edgeR","EDASeq","biomaRt","supraHex"))
devtools::install_github("vqv/ggbiplot")
```

* Creating roxygen documentation
```r
devtools::document()
```
* Building the package
```r
devtools::build()
```
 
* Summary of set up
```r
install.packages(path_to_package, repos = NULL, type="source")
```

### Main structure of the repository ###
| Folder  | Description |
| ------------- | ------------- |
| R	  | main R files
| man	| Manual files (can be created by roxygen with devtools::document())
| demo	| Example of how to run the code
| inst	| Files that should be acessed by the installed package
| DESCRIPTION	| Package description
| NAMESPACE	| Package namespace (can be created by roxygen with devtools::document())
| README.md | Project readme
