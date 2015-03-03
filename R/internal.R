# @title .onAttach
# @description  Load required data into gloval enviroment
# @keywords internal
.onAttach <- function (libname, pkgname){
  load(file = system.file("extdata/PlatformAndAssociatedData_WebSites3.RData", 
                          package="TCGAbiolinks"), 
       .GlobalEnv) 
}
