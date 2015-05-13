.createDirectory <-
function(base){
  i="";
  while(file.exists(paste(base, i, sep=""))){
    if(i==""){
      i=1;
    }else{
      i=i+1;
    }
  }
  toDir = paste(base, i, sep="")
  dir.create(toDir)
  
  toDir
}
.DownloadData_fromURL <-
function(url, sep = ",", header = TRUE){
  require(httr)
  #handle_find(url)
  #handle_reset(url)  
  request <- GET(url)
  stop_for_status(request)
  handle <- textConnection(content(request, as = 'text'))
  on.exit(close(handle))
  read.table(handle, sep = sep, header = header)
}
.DownloaDmageTAB <-
function(Description,TumorDataList, keySpecies,startK, stopK, typeProtein = F ){
  Description2 <- paste(Description, keySpecies, sep = "")
  Description_i_ord <- paste(Description2, "?C=M;O=D", sep = "")
  x <- .DownloadURL(Description_i_ord)
  if(length(x)!=10){
    siteNewLevel <- .FindGrepSite(x,Key="mage-tab",Description2)
    siteNewLevelSdrf <- DownloadSdrf(siteNewLevel)
    tmp2 <-  siteNewLevelSdrf$Comment..TCGA.Barcode.
    
    if(typeProtein==T){
      siteNewLevelDesign <- DownloadTypeFile(siteNewLevel,"design")
      tmp2 <- siteNewLevelDesign$Sample.description
    }
    
    tmp2 <- tmp2[grep("TCGA",tmp2)]
    NumberSample <- length(unique(substr(tmp2, startK, stopK)))
    msgOUT <-  paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " ,  " .n samples ", NumberSample, sep="")
    print(msgOUT)
    SampleTmp <- unique(substr(tmp2, startK, stopK))
    idx<- which(names(TumorDataList) == unique(tmp4$Platform))
    TumorDataList[[idx]] <- SampleTmp
  }
  return(TumorDataList)
}
.DownloaDmageTAB_sdrf <-
function(Description,keySpecies,KeyGrep1 = "mage-tab", KeyGrep2 = "sdrf"){
  Description2 <- paste(Description, keySpecies, sep = "")
  Description_i_ord <- paste(Description2, "?C=M;O=D", sep = "")
  x <- .DownloadURL(Description_i_ord)
  if(length(x)!=10){
    siteNewLevel <- .FindGrepSite(x,Key=KeyGrep1,Description2)
    x <- .DownloadURL(siteNewLevel)
    x2 <- x[grep(KeyGrep2,x)]
    
    x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
    x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
    site3 <- paste(siteNewLevel, x2,sep="" )
    #site4 <- paste(keySpecies,unlist(strsplit(site3,keySpecies))[2],sep="")
    site4 <- unlist(strsplit(site3,keySpecies))[2]
    
    print(site4)
    return(site4)
  } else{return("")}
  
}
.DownloadManifest <-
function(siteNewLevel){
  site3 <- paste(siteNewLevel, "MANIFEST.txt",sep="")
  x <- .DownloadURL(site3)
  writeLines(x, "x2.txt" )
  tmp2 <- read.table("x2.txt", quote="\"", stringsAsFactors = F)[,2]
  tmp2 <- tmp2[ nchar(tmp2) > 20 ]
  return(tmp2)
}
.DownloadSdrf <-
function(siteNewLevel){
  x <- .DownloadURL(siteNewLevel)
  x2 <- x[grep("sdrf",x)]
  x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
  x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
  site3 <- paste(siteNewLevel, x2,sep="" )
  
  x <- .DownloadURL(site3)
  writeLines(x, "x2.txt" )
  tmp2 <- as.data.frame(read.delim("x2.txt",stringsAsFactors = F))
  return(tmp2)
}
.DownloadTypeFile <-
function(siteNewLevel,keyDown){
  x <- .DownloadURL(siteNewLevel)
  x2 <- x[grep(keyDown,x)]
  x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
  x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
  site3 <- paste(siteNewLevel, x2,sep="" )
  
  x <- .DownloadURL(site3)
  writeLines(x, "x2.txt" )
  tmp2 <- as.data.frame(read.delim("x2.txt",stringsAsFactors = F))
  return(tmp2)
}

.DownloadURL <-
  function(Site){
    # setInternet2(use = TRUE)
    Site <- URLencode(Site)
    x=  getURL(Site, ssl.verifypeer = FALSE)
    
    x <- unlist(strsplit(x,"\n"))
    
    return(x)
    
  }    
    
.FindGrepSite <-
function(x,Key,Description){
  x2 <- x[grep(Key, x)]
  
  if( Key != "sdrf"){ x2 <- x2 [- grep("tar.gz", x2)][1] }
  
  
  x2  <- as.matrix(sapply(strsplit(x2, ">"), function(y) y[2]))
  x2  <- as.matrix(sapply(strsplit(x2, "<"), function(y) y[1]))
  site2 <- paste(Description, x2,sep="" )
  return(site2)
}
.Random.seed <-
c(403L, 16L, 1230691973L, 505388895L, 2058239974L, -1291062672L, 
-170568301L, 946504241L, 375120096L, -1398924930L, -360014839L, 
1508174907L, -1079954934L, -556385812L, -1765564817L, 423432613L, 
1202162012L, 1779479858L, 612830253L, 921715431L, -1126586354L, 
-1703544344L, -6022165L, 680526985L, -1744224232L, -1225171130L, 
-1736933983L, -915942349L, 1855256994L, 1905433748L, 1619516183L, 
-1612671731L, 706421796L, 1757246474L, -1384060555L, -1456230193L, 
-561208362L, -659924096L, -1363955773L, 384117089L, -701759536L, 
-2054326482L, -2046825767L, -1807898389L, -1126123174L, 974062364L, 
-1042360801L, 872975989L, 485296908L, 60253922L, 2124778429L, 
-973689801L, 2092515038L, -240215528L, -2072205957L, 1450794265L, 
-1889754712L, -237931178L, -561928751L, -265142909L, 437044242L, 
1949527972L, 780691303L, 829469885L, -282270668L, -182319206L, 
1109414885L, 1656835967L, 1034440454L, 1839958160L, 475783475L, 
-957573423L, -1747360640L, 255903070L, -38290775L, -756054693L, 
-293783702L, -194576884L, -1064470961L, 1692375813L, 1430698940L, 
-934662510L, 584035725L, -2126357049L, -2056577362L, 1009039112L, 
1295983883L, -787966999L, -589994824L, 1438454630L, -1355098111L, 
-1697642861L, -1802176190L, 1455199988L, -2068774153L, 837620461L, 
-1343017596L, 1190921130L, -273155819L, 1567142831L, 1546308534L, 
2064208352L, -161828573L, 1752737345L, -1437316816L, 898120526L, 
-358080071L, 1202256843L, -2040057734L, 1485902204L, 1316692415L, 
-1440194923L, -74799956L, 2073169538L, -1809988003L, 793569623L, 
-29232066L, 1262278776L, -115899557L, -1215758663L, 289395464L, 
-1963054666L, 428848881L, 1791975715L, -8195214L, 649476548L, 
-1477115513L, -84821155L, 1574434132L, -563562630L, 2026673861L, 
-593504225L, -647282778L, 2101985584L, -79995949L, 553091569L, 
-828897504L, -1629531714L, -1863884343L, -133048965L, 1747731274L, 
-1383562324L, 16935599L, -1349789595L, -989177828L, -1361987982L, 
139954413L, 298830119L, -1753172658L, 2081527208L, -1711132245L, 
1541753545L, 1070545624L, -719042682L, 615275105L, 187335667L, 
-1573202334L, 1799418068L, -1465484841L, 977234893L, -1542971548L, 
181930826L, -463217995L, -364991985L, -678532330L, -1948924480L, 
-460253821L, -577004639L, -1759501424L, 340672622L, 793806361L, 
1683646123L, -797238630L, 1525319260L, -683355297L, -323891915L, 
386694348L, -245234142L, 1228863613L, -2036022921L, -2090242402L, 
324590936L, 1302269115L, -772410791L, -457385496L, -1037837290L, 
-1430870127L, -1292720061L, 2099711826L, 1761817572L, 2055391783L, 
-603569539L, -625067148L, 1936545626L, -7024219L, -1814639169L, 
-1778897338L, 311000272L, -339503373L, 1046679569L, 1724388928L, 
166935070L, 670938345L, 910212635L, 1752399402L, -2038881972L, 
-58291441L, -600531387L, 832022524L, 1568898642L, 1127818189L, 
1405881735L, -972371858L, 1813286984L, 1061560907L, 1502216105L, 
1694813944L, 588466726L, -969874623L, 1935280339L, -858510206L, 
-376534348L, -730227145L, -1467431635L, -1089622460L, -889080982L, 
-1629970475L, 771098991L, -433470346L, -749804396L, 1705622658L, 
-684752960L, -1213738692L, 1915943568L, 832418850L, -1315236856L, 
-1701648244L, 1360733212L, 455884434L, 1720042624L, -141946860L, 
959192360L, 1902550586L, -644666160L, 1200972908L, -324545868L, 
1336942610L, 1700509664L, -212813876L, 56085984L, 2110694946L, 
1286865192L, -552453876L, 1787586108L, -265432846L, -791096592L, 
2029429572L, -438582696L, -99087814L, -730341648L, 94377148L, 
195950100L, -1007986494L, 1006856288L, -1870531204L, 1180890320L, 
1050844194L, 364954280L, -509868148L, 1772499388L, -732195790L, 
438744896L, 755714228L, 690387784L, 945330874L, -1004824816L, 
1845764332L, -1742277548L, -1907908846L, -367895648L, -603824180L, 
38681440L, 1925361858L, 892896808L, -1999247956L, 1937235772L, 
-1706244238L, -161850576L, 1089118756L, -1730553800L, 250332378L, 
850487696L, -935393540L, -1315262572L, 1771801410L, -1824275968L, 
626891324L, -67568048L, -1849140574L, 2107464520L, 1240105484L, 
670256348L, -658884014L, 387180672L, -528172204L, 1703623144L, 
1051437818L, -874842288L, -1759888340L, -854084236L, -461027694L, 
1682610784L, 949068428L, -1733657760L, -1622981918L, 411246760L, 
435486924L, 609421820L, -1847849038L, 866563568L, 113747396L, 
986242904L, 876530554L, 1890809264L, 1675491772L, 1497442388L, 
640617410L, -1877261536L, 1711153084L, 1177709136L, 1529011682L, 
1073983848L, 1997229388L, 1441073404L, 1334551282L, 1986057984L, 
418580084L, -1211348984L, -1153713734L, 878614480L, 1064974828L, 
2099744660L, 130216658L, 1811910240L, 169272652L, 1899861280L, 
-699223422L, -169973272L, -1180710932L, -865634756L, 263343730L, 
696075120L, 1912752420L, 300184120L, -216243814L, -1501296176L, 
-332333636L, -1233008492L, -2110514046L, 222750144L, 626172732L, 
1908929680L, 524501794L, 602344968L, -1936216820L, 1544992156L, 
-548741614L, -540963968L, -2123939052L, -1317768920L, -68096966L, 
1085418448L, -665451924L, -743310796L, 719190674L, 593534304L, 
401189196L, 1322178528L, 248911266L, -1818141272L, -2032507252L, 
-1231872580L, 24440434L, 2107399152L, 78739012L, 526440536L, 
956962874L, 27838576L, 221509692L, 2075105556L, -1247001918L, 
1774522080L, -1281652612L, 1677422288L, -1184613982L, 1943629864L, 
1390975756L, -1636068932L, 2144514866L, 1551334080L, 958252084L, 
1401624264L, -1035479878L, -35031408L, -1825009940L, -193574316L, 
-1929665902L, 1346060576L, 1465108940L, -1243004192L, -707331902L, 
528495656L, -1410631764L, 563995324L, -490116110L, -397210064L, 
798836644L, -1367331400L, -755964070L, 616214672L, 744858492L, 
606802196L, -1109399230L, 1975576832L, 737504700L, 1415792976L, 
1281624354L, -85368248L, 263858188L, -1368037796L, 1752697298L, 
-1108169216L, -1566720684L, 647498984L, -857470726L, -839876912L, 
286247724L, -1118623500L, 1515037074L, 27960672L, 1381651468L, 
-822080288L, 704426594L, -1398676312L, 1702806604L, 1532510076L, 
553335090L, 1285735024L, -626412732L, 1664859096L, -1728730246L, 
1197987888L, -1915027524L, 1787382356L, 770469954L, 430319264L, 
-656520824L, -195033295L, -1934075997L, 64241572L, 1072257618L, 
-933326425L, -1620204047L, 387281206L, 876145108L, -1631420571L, 
919826287L, 1184068008L, 912664326L, -1778597325L, -1474944875L, 
-1410170638L, -385312448L, -1768446295L, 1689100059L, 292980972L, 
1389298938L, -965784273L, 2075048249L, 1976945486L, 1870377372L, 
472886029L, -2117264009L, -126477952L, -1051377698L, 636905643L, 
-2052768851L, -62078950L, 534107864L, 72054593L, 1081649139L, 
1103882772L, 912815010L, -278172489L, 1183712897L, 44473990L, 
965142212L, 1642498037L, 1945300575L, -726464328L, 1939434902L, 
-681787485L, -2056383643L, -1874423166L, -482451472L, 1750329753L, 
-1402746997L, 1327092348L, 179855530L, -1603222817L, 1443680361L, 
2119567998L, 906728044L, 1115852477L, -1236167801L, -1296427152L, 
-1417028082L, -1166490821L, 395247645L, 1993611466L, -1048263704L, 
-1329065775L, -960312125L, -393309244L, -294045582L, 517835079L, 
1800516881L, -932425578L, 548730356L, 810338181L, 1753642703L, 
463782856L, 278688358L, -643257197L, -175774731L, 1422288338L, 
-127838688L, 1084143625L, 1442397499L, -1305562996L, -1734853350L, 
2048611215L, 519646425L, -1372388498L, 1249088252L, 2049941613L, 
1091431447L, -10375712L, 2122496446L, -1333629045L, 1926918157L, 
-202082246L, 690927736L, -51178463L, -580278573L, -1135290892L, 
987707010L, -1229876457L, 184451553L, -1229181146L, 2011265572L, 
-1128268587L, 1871096959L, -324083432L, -1917677514L, 1124026051L, 
1192438533L, -1562472798L, 1016233104L, -775648199L, 913818859L, 
-815163940L, 1371872266L, 1081932543L, -525749431L, 19024094L, 
-29903220L, -208027043L, 1832763495L, -631054960L, 793669038L, 
461322331L, -1847939011L, 352364202L, -1412021432L, 640518001L, 
2086383843L, -1700206620L, 1262403346L, -340598041L, 1686930609L, 
1685939958L, 1909829908L, 380626469L, -1194099921L, 1769827560L, 
-329921082L, 1701985395L, 838817493L, 1714810034L, -1318051712L, 
-1143075735L, 719001819L, -665940948L, 1681650746L, -1423482129L, 
-652698247L, 428462606L, -856972324L, -356053939L, 193666999L, 
-1629912768L, 1090250782L, -2130862357L, 1501777517L, -1023568166L, 
258852504L, 1386769665L, 1518949043L, 2075794900L, 1504378082L, 
1234965879L, 1517980865L, -1561169850L, -208339068L, -2106102609L
)
.TCGAQuery2 <-
function(Tumor,siteTCGA){
  
  TumorData <- matrix(0, 1, 8)
  colnames(TumorData) <- c("Total", "Exome", "SNP", "Methylation", "mRNA", "miRNA", "Clinical","Protein")
  rownames(TumorData) <- Tumor
  LevelsPlatforms <- c("Level_1", "Level_2", "Level_3")
  
  #for ( i in 2: ncol(TumorData)){
  
  
  i<-5
  Type <-  colnames(TumorData)[i]
  tmp <- PlatformAndAssociatedData[PlatformAndAssociatedData$Type%in% Type,]
  
  tmp2a <- tmp[ grep(tolower(Tumor),tmp$Tumor),]
  
  if( nrow(tmp2a)!=0){
    
    Species <- unique(tmp2a$Species)
    
    
    
    for( j in 1:length(Species)){
      SpecieCurr <- Species[j]
      tmp3 <- tmp2a[tmp2a$Species %in% SpecieCurr, ]
      Centers <- unique(tmp3$Center)
      
      for( k in 1:length(Centers)){
        CenterCurr <- Centers[k]
        tmp3b <- tmp3[tmp3$Center %in% CenterCurr,]
        
        Platforms <- unique(tmp3b$Platform)
        
        for( q in 1:length(Platforms)){
          
          PlatformCurr <- Platforms[q]
          tmp4 <- tmp3b[tmp3b$Platform %in% PlatformCurr, ]
          
          key1<- paste(unique(tmp4$CenterType), unique(tmp4$Center), unique(tmp4$Platform), sep="/")
          Description <- paste(siteTCGA, tolower(Tumor), "/",key1, sep="")
          
          if( length(grep("agilentg4502a_07",unique(tmp4$Platform))) > 0   || unique(tmp4$Platform) ==  "ht_hg-u133a" || unique(tmp4$Platform) ==  "hg-u133_plus_2" || unique(tmp4$Platform) ==  "illuminaga_mrna_dge" ){
            Description <- paste(Description, "/transcriptome/", sep = "")
            Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
            x <- DownloadURL(Description_i_ord)
            
            
            
            for( w in 1:length(LevelsPlatforms)){
              
              siteNewLevel <- .FindGrepSite(x,Key=LevelsPlatforms[w],Description)
              tmp2 <- DownloadManifest(siteNewLevel)
              
              
              if(length(grep("agilentg4502a_07",unique(tmp4$Platform))) > 0){
                NumberSample <- length(unique(substr(tmp2, 1, 23)))
              }
              
              if(unique(tmp4$Platform) ==  "ht_hg-u133a"){
                NumberSample <- length(unique(substr(tmp2, 44, 58)))
              }
              
              if(unique(tmp4$Platform) ==  "hg-u133_plus_2" || unique(tmp4$Platform) ==  "illuminaga_mrna_dge"){
                NumberSample <- length(unique(substr(tmp2, 1, 16)))
              }
              
              message <- paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " ,LevelsPlatforms[w] ,  " .n samples ", NumberSample, sep="")
              print(message)
              
            } #end for LevelsPlatform
          } #end platform Exp-Gene
          
          
          if(unique(tmp4$Platform) ==  "huex-1_0-st-v2"){
            Description <- paste(Description, "/exon/", sep = "")
            Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
            x <- DownloadURL(Description_i_ord)
            siteNewLevel <- .FindGrepSite(x,Key="mage-tab",Description)
            siteNewLevelSdrf <- DownloadSdrf(siteNewLevel)
            tmp2 <-  siteNewLevelSdrf$Comment..TCGA.Barcode.
            tmp2 <- tmp2[grep("TCGA",tmp2)]
            NumberSample <- length(unique(substr(tmp2, 1, 16)))
            message <- paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " , " .n samples ", NumberSample, sep="")
            print(message)
          }
          
          
          
          if(unique(tmp4$Platform) ==  "illuminahiseq_rnaseq" || unique(tmp4$Platform) ==  "illuminaga_rnaseq"){
            Description <- paste(Description, "/rnaseq/", sep = "")
            Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
            x <- DownloadURL(Description_i_ord)
            siteNewLevel <- .FindGrepSite(x,Key=LevelsPlatforms[3],Description)
            tmp2 <- DownloadManifest(siteNewLevel)
            NumberSample <- length(unique(substr(tmp2, 13, 29)))
            message <- paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " ,LevelsPlatforms[3] ,  " .n samples ", NumberSample, sep="")
            print(message)
          }
          
          
          
          if(unique(tmp4$Platform) ==  "illuminahiseq_rnaseqv2" || unique(tmp4$Platform) ==  "illuminaga_rnaseqv2" ){
            Description <- paste(Description, "/rnaseqv2/", sep = "")
            Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
            x <- DownloadURL(Description_i_ord)
            siteNewLevel <- .FindGrepSite(x,Key=LevelsPlatforms[3],Description)
            tmp2 <- DownloadManifest(siteNewLevel)
            
            NumberSample <- length(unique(substr(tmp2, 13, 29)))
            message <- paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " ,LevelsPlatforms[3] ,  " .n samples ", NumberSample, sep="")
            print(message)
          }
          
          
          if(unique(tmp4$Platform) ==  "illuminahiseq_totalrnaseqv2"){
            Description <- paste(Description, "/totalrnaseqv2/", sep = "")
            Description_i_ord <- paste(Description, "?C=M;O=D", sep = "")
            x <- DownloadURL(Description_i_ord)
            siteNewLevel <- .FindGrepSite(x,Key=LevelsPlatforms[3],Description)
            tmp2 <- DownloadManifest(siteNewLevel)
            
            NumberSample <- length(unique(substr(tmp2, 9, 44)))
            message <- paste(Type, " ", SpecieCurr, " ",   CenterCurr, " ", unique(tmp4$Platform) , " " ,LevelsPlatforms[3] ,  " .n samples ", NumberSample, sep="")
            print(message)
          }
          
          
          
        } # end platform
      } #end for Centers
    } # end for species
    
    
    
    #}
  }
  
}
