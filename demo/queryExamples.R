# query examples

# user      system   elapsed
# 2.188     0.421    19.934
query <- TCGAQuery(tumor="gbm",added.since="01/01/2013",added.up.to = "06/01/2013")

# user      system   elapsed
# 0.085     0.028     1.063
query <- TCGAQuery(tumor = "gbm", platform = "bio")

# user      system   elapsed
# 0.063     0.024     1.495
query <- TCGAQuery(tumor = "gbm", platform = "HumanMethylation450",level="3")

#  user      system   elapsed
#  1.763     0.704    28.846
query <- TCGAQuery(listSample = "TCGA-61-1743-01A-01D-*")

#  user      system   elapsed
#  0.244     0.106     4.261
query <- TCGAQuery(listSample = "TCGA-61-1743-01A-01D-0649-04")

# user      system   elapsed
# 0.239     0.090     3.987
query <- TCGAQuery(listSample = "TCGA-61-1743-01A-01D-0649-04", level="3")

# user      system   elapsed
# 0.260     0.131     7.577
query <- TCGAQuery(listSample = "TCGA-61-1743-01A-01D-0649-04", tumor = "OV", platform = "CGH-1x1M_G4447A")
