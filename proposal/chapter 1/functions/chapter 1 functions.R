#functions to use for example R-manuscripts project

## Load host-dependent directory environment
    winos <- ifelse(grepl("windows", Sys.info()['sysname'], ignore.case=T), 1, 0)
    # specify absolute file path to file_dir_params.R
    if(winos==1) source("C:/mypath/myproject/Chapter 1/functions/file_dir_params.R")      #Change to your path
    if(winos==0) source("~/mypath/myproject/Chapter 1/functions/file_dir_params.R")       #Change to your path 
    rm(winos)

####################

## Call any functions outside of this script
  source(paste0(prjfuns, "theme_Publication.R"))
  source(paste0(prjfuns, "plot.fpt.r"))
  source(paste0(prjfuns, "multiplot.r"))
  #x <- readRDS(paste0(src, "somefile.Rds"))                                      #example code for reference with .rds files
  #source(chartr("/","\\", "/mypath/myproject/functions/multiplot.r"))            #example code if windows is causing issues with 
                                                                                  #  file paths, not solved by file_dir_params 

## When compiling multiple .R or .Rmd files, it is best to clear your workspace after each compilation
## Call at the beginning of each file
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices",
                      "package:utils","package:datasets","package:methods",
                      "package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}


## Include any functions used in your project below

foobar <- function(foo) {
  print("Hello World!")
  print(paste("foo is", foo)
}

foo_sq <-function(foo) {
  foo*foo
  }

