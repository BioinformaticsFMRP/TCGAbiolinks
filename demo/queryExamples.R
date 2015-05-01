# query examples
query <- TCGAQuery(tumor="gbm",added.since="01/01/2013",added.up.to = "06/01/2013")
query <- TCGAQuery(tumor = "gbm", platform = "bio")
query <- TCGAQuery(tumor = "gbm", platform = "HumanMethylation450",level="3")
query <- TCGAQuery(listSample = "TCGA-61-1743-01A-01D-*")
query <- TCGAQuery(listSample = "TCGA-61-1743-01A-01D-0649-04")
