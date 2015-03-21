# @title .onAttach
# @description  Load required data into gloval enviroment
# @keywords internal
.onAttach <- function (libname, pkgname){
  load(file = system.file("extdata/PlatformMat.rda",
                          package="TCGAbiolinks"),
       .GlobalEnv)
}

createDir <- function(base){
  i="";
  while(file.exists(paste(base, i, sep=""))){
    if(i==""){
      i=1;
    }else{
      i=i+1;
    }
  }
  toDir = paste0(base, i)
  dir.create(toDir, showWarnings = F, recursive = T, mode = "0777")
  toDir
}
DownloadHTML <- function(url){
  bo2 = T
  count <- 0
  handle_find(url)
  while(bo2){
    request = try(GET(url, timeout(100)), silent = T)
    if( class(request) == "try-error"){
      Sys.sleep(1)
      bo2 = T
      count = count + 1
      handle_find(url)
      if(count%%10==0) print(paste("Reconnection attempt #",count,sep=""))
    }else if(count>=200){
      stop("Connetion limit exceded. Check your internet connection and your proxy settings.
           If you are downloading very big files (proteins for example) you should add the proper variable.
           Take a look to the documentation. If the problem persists please contact the mantainers.")
    }else{
      bo2 = F
    }
  }
  u<-read.table(textConnection(content(request, as = 'text')), sep = ",", header = T)

  return(deparse(u))
}
GrepSite <- function(x,Key){
  x <- x[grep(Key, x)]
  x = sapply(strsplit(x, ">"), function(y) y[2])
  x = sapply(strsplit(x, "<"), function(y) y[1])
  x <- x[grep("/", x)]
  if(length(grep("lost",x))!=0) x <- x[-grep("lost", x)]
  return(x)
}
