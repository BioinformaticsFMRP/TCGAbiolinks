query <- tcgaSearch(tumor = "gbm",
                    added.since = "01/01/2013",
                    added.up.to = "06/01/2013")

query <- tcgaSearch(tumor = "gbm", platform = "bio")

query <- tcgaSearch(tumor = "gbm", platform = "HumanMethylation450",
                    level = "3")

query <- tcgaSearch(listSample = "TCGA-61-1743-01A-01D-*")

query <- tcgaSearch(listSample = "TCGA-61-1743-01A-01D-0649-04")

query <- tcgaSearch(listSample = "TCGA-61-1743-01A-01D-0649-04", level = "3")

query <- tcgaSearch(listSample = "TCGA-61-1743-01A-01D-0649-04",
                    tumor = "OV", platform = "CGH-1x1M_G4447A")
