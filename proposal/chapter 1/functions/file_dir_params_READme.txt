If using 'file_dir_params.R' to specify file paths for collaboration or using across multiple OS platforms:
  1. change the file paths within file_dir_params.R to specify your project path
  2. Add the following code to the beginning of either:
      a. A single project functions Rscript that you call in all documents in your project (e.g. chapter 1 functions.R)
      b. ALL files within your project that include file imports or exports (i.e. do not stand alone)
  
    ## start code
    winos <- ifelse(grepl("windows", Sys.info()['sysname'], ignore.case=T), 1, 0)
    # specify absolute file path to file_dir_params.R
    if(winos==1) source("C:/mypath/myproject/Chapter 1/functions/file_dir_params.R")
    if(winos==0) source(""~/mypath/myproject/Chapter 1/functions/file_dir_params.R")
    rm(winos)
    ## end code
   
