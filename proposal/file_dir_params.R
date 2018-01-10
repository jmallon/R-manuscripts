#--begin file_dir_params.R script--#
#from  http://derekyves.github.io/2016/05/10/codeshare.html

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
#build <- 1

if(grepl("windows", fdirs$computeros)==F){
  fdirs$prjdir <- "~/Box Sync/R/Chapter 1/Scripts/"
  fdirs$prjdta <- "./Data/"
}else{
  fdirs$prjdir <- "C:/Users/Julie/Box Sync/R/Chapter 1/Scripts/"
  fdirs$prjdta <- "./Data/" #update for windows 
}

# Add some child objects to the fdirs environment:
fdirs$prjrslts   <- paste0(fdirs$prjdta, "Results/")
if(show==1 & interactive()) cat("\nprjrslts =", fdirs$prjrslts) #?
#if(build==1) system(paste0("mkdir -p ", fdirs$prjrslts)) #?

fdirs$prjfuns   <- paste0(fdirs$prjdir, "Functions/")
if(show==1 & interactive()) cat("\nprjfuns =", fdirs$prjfuns) #?
#if(build==1) system(paste0("mkdir -p ", fdirs$prjfuns)) #?


# Attach the new environment (and safely reload if already attached):
while("fdirs" %in% search())
  detach("fdirs")
attach(fdirs)

#x <- readRDS(paste0(src, "somefile.Rds"))