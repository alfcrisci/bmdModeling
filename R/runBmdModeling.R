
#' Run the BmdModeling Application
#' @return no return value
#' @import shiny
#' @importFrom devtools install_github
#' @export
runBmd <- function(){
  
  # List all packages used in the shiny app but not in R functions
  requiredPackages <- c("rmarkdown", "V8", "shinyjs", "png", "grid")
  newPackages <- requiredPackages[!(requiredPackages %in% installed.packages()[,"Package"])]
  
  if(length(newPackages) > 0) {
    
    install.packages(newPackages)
    
  } 
  
  if(!requireNamespace("shinysky")){
    
    devtools::install_github("AnalytixWare/ShinySky")
    
  }
  
  
  tmpDir <- tempdir()
  
  setwd(tmpDir)
  
  # Copy server.R and ui.R (not folder www)
  uiDir <- system.file("ui", package = "bmdModeling")
  uiFiles <- list.files(path = uiDir, full.names = TRUE)
  uiFiles <- uiFiles[!grepl("www", uiFiles)]
  
  sapply(uiFiles, function(x){
        file.copy(from = x, to = file.path(tmpDir, basename(x)),
            overwrite = TRUE)}
  )
  
  # Make www directory and copy its files
  if (!dir.exists(file.path(tmpDir, "www"))) {
    
    dir.create(path = file.path(tmpDir, "www"))
    
  }
  
  wwwFiles <- list.files(path = file.path(uiDir, "www"), full.names = TRUE)
  
  sapply(wwwFiles, function(x){
        file.copy(from = x, to = file.path(tmpDir, "www", basename(x)),
            overwrite = TRUE)}
  )
  
  
  runApp(appDir = tmpDir)
  
}
