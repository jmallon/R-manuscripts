#--begin file_dir_params.R script--#
### Credit goes to Derek Darves
### Step by step explanations availble at http://derekyves.github.io/2016/05/10/codeshare.html

#----------------------------#
# Make a new environment:
fdirs <- new.env()

# Function to standardize host OS name
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else {
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
fdirs$computeros <- get_os()

## Declare the root project and data directory:
show  <- 1
#build <- 1                                                   #use if you are creating the directories for the first time

if(grepl("windows", fdirs$computeros)==F){
  fdirs$prjdir <- "~/Documents/GitHub/R-manuscripts/"
}else{                                                        #using windows machine
  fdirs$prjdir <- "C:/USERNAME/Documents/GitHub/R-manuscripts/"
}

# Add some child objects to the fdirs environment: (Example, a results and a functions subfolder
fdirs$prjdta   <- paste0(fdirs$prjdir, "data/")
if(show==1 & interactive()) cat("\nprjdta =", fdirs$prjdta) 
#if(build==1) system(paste0("mkdir -p ", fdirs$prjrslts)) 

fdirs$prjfuns   <- paste0(fdirs$prjdir, "functions/")
if(show==1 & interactive()) cat("\nprjfuns =", fdirs$prjfuns)
#if(build==1) system(paste0("mkdir -p ", fdirs$prjfuns)) 

fdirs$prjrslts   <- paste0(fdirs$prjdta, "results/")        #not used here but could be used to save output
if(show==1 & interactive()) cat("\nprjrslts =", fdirs$prjrslts) 
#if(build==1) system(paste0("mkdir -p ", fdirs$prjrslts)) 


# Attach the new environment (and safely reload if already attached):
while("fdirs" %in% search())
  detach("fdirs")
attach(fdirs)
