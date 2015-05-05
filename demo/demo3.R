samples <- c( "TCGA-02-0010*","TCGA-02-0014*", "TCGA-02-0058*", "TCGA-02-0114*", "TCGA-02-2483*", "TCGA-06-0128*",
              "TCGA-06-0129*", "TCGA-06-1805*", "TCGA-06-2570*", "TCGA-06-5417*", "TCGA-06-6389*", "TCGA-06-A7TL*",
              "TCGA-12-0827*", "TCGA-14-1821*", "TCGA-16-0850*", "TCGA-19-2629*", "TCGA-26-5133*", "TCGA-27-2521*",
              "TCGA-CS-6665*", "TCGA-DB-5277*", "TCGA-DB-A4XB*", "TCGA-DH-5142*", "TCGA-DH-A7UT*", "TCGA-DU-6396*",
              "TCGA-DU-6408*", "TCGA-DU-7010*", "TCGA-DU-A5TP*", "TCGA-DU-A7TG*", "TCGA-E1-5304*", "TCGA-E1-A7YE*",
              "TCGA-E1-A7YI*", "TCGA-E1-A7YK*", "TCGA-E1-A7YV*", "TCGA-E1-A7YY*", "TCGA-FG-A6J3*", "TCGA-FG-A87N*",
              "TCGA-HT-7477*", "TCGA-HT-7601*", "TCGA-HT-7686*", "TCGA-HT-7689*", "TCGA-HT-7880*", "TCGA-HT-7902*",
              "TCGA-HT-8018*", "TCGA-HT-8106*", "TCGA-HT-8111*", "TCGA-HT-8563*", "TCGA-HT-A5R7*", "TCGA-HT-A614*",
              "TCGA-HT-A618*", "TCGA-HT-A61A*", "TCGA-HW-8319*", "TCGA-IK-7675*", "TCGA-P5-A5EU*", "TCGA-QH-A65S*",
              "TCGA-RY-A83Z*", "TCGA-S9-A6TV*", "TCGA-S9-A6WI*", "TCGA-S9-A7IS*", "TCGA-S9-A7R7*", "TCGA-TM-A84F*",
              "TCGA-TM-A84I*", "TCGA-TQ-A7RM*", "TCGA-TQ-A7RR*", "TCGA-02-0011*", "TCGA-02-0024*", "TCGA-02-0027*",
              "TCGA-02-0034*", "TCGA-02-0047*", "TCGA-02-0060*", "TCGA-02-0069*", "TCGA-02-0107*", "TCGA-02-0113*",
              "TCGA-06-0121*", "TCGA-06-0139*", "TCGA-06-0141*", "TCGA-06-0142*", "TCGA-06-0650*", "TCGA-06-0875*",
              "TCGA-06-0881*", "TCGA-06-0882*", "TCGA-06-1086*", "TCGA-06-1801*", "TCGA-06-1806*", "TCGA-06-2566*",
              "TCGA-06-2569*", "TCGA-06-5410*", "TCGA-06-5416*", "TCGA-06-5858*", "TCGA-06-6391*", "TCGA-06-6698*",
              "TCGA-06-A5U0*", "TCGA-12-0821*", "TCGA-12-0822*", "TCGA-12-0826*", "TCGA-12-1088*", "TCGA-14-0736*",
              "TCGA-14-0740*", "TCGA-14-0871*", "TCGA-14-1455*", "TCGA-14-3477*", "TCGA-15-1447*", "TCGA-19-1385*",
              "TCGA-19-1388*", "TCGA-19-1389*", "TCGA-19-5956*", "TCGA-26-1440*", "TCGA-28-2510*", "TCGA-28-2512*",
              "TCGA-28-5218*", "TCGA-32-1973*", "TCGA-32-1980*", "TCGA-32-2498*", "TCGA-76-4934*", "TCGA-CS-6669*",
              "TCGA-DB-A75P*", "TCGA-DH-5140*", "TCGA-DU-6392*", "TCGA-DU-6404*", "TCGA-DU-A7TB*", "TCGA-FG-5963*",
              "TCGA-FG-8181*", "TCGA-FG-8189*", "TCGA-HT-7469*", "TCGA-HT-7680*", "TCGA-HT-7691*", "TCGA-HT-7854*",
              "TCGA-HT-7857*", "TCGA-HT-8015*", "TCGA-HT-8019*", "TCGA-HT-8107*", "TCGA-HT-8558*", "TCGA-HT-8564*",
              "TCGA-P5-A5EY*", "TCGA-P5-A5F6*", "TCGA-QH-A6CS*", "TCGA-S9-A6UA*", "TCGA-S9-A89V*", "TCGA-TM-A84C*",
              "TCGA-TM-A84J*", "TCGA-VM-A8C9*")

query <- TCGAQuery(tumor = "LGG", platform = "HumanMethylation450",level="3",listSample = samples)
query <- TCGAQuery(tumor = "LGG", platform = "HumanMethylation27",level="3",listSample = samples)
query <- TCGAQuery(tumor = "GBM", platform = "HumanMethylation450",level="3",listSample = samples)
query <- TCGAQuery(tumor = "GBM", platform = "HumanMethylation27",level="3",listSample = samples)

TCGADownload(query,path="dataDemo3")

