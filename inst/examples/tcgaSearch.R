query <- tcgaQuery(tumor = "gbm",
                    added.since = "01/01/2013",
                    added.up.to = "06/01/2013")

query <- tcgaQuery(tumor = c("gbm","lgg"),
                    platform = c("HumanMethylation450","HumanMethylation27"))

query <- tcgaQuery(tumor = "gbm", platform = "bio")

query <- tcgaQuery(tumor = "gbm", platform = "HumanMethylation450",
                    level = "3")

query <- tcgaQuery(samples = "TCGA-61-1743-01A-01D-*")

query <- tcgaQuery(samples = "TCGA-61-1743-01A-01D-0649-04")

query <- tcgaQuery(samples = "TCGA-61-1743-01A-01D-0649-04", level = 3)

query <- tcgaQuery(samples = "TCGA-61-1743-01A-01D-0649-04",
                    tumor = "OV", platform = "CGH-1x1M_G4447A")
