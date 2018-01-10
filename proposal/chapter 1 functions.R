#functions to use for chapter 1 vulture stopover analysis

## Load host-dependent directory environment
winos <- ifelse(grepl("windows", Sys.info()['sysname'], ignore.case=T), 1, 0)
if(winos==1) source("C:/Users/Julie/Box Sync/R/Chapter 1/Scripts/Functions/file_dir_params.R")
if(winos==0) source("~/Box Sync/R/Chapter 1/Scripts/Functions/file_dir_params.R")
rm(winos)
####################


#x <- readRDS(paste0(src, "somefile.Rds"))
source(paste0(prjfuns, "theme_Publication.R"))
source(paste0(prjfuns, "plot.fpt.r"))
source(paste0(prjfuns, "multiplot.r"))

#source(chartr("/","\\", "/Users/Julie/Documents/R/functions/multiplot.r"))
#source("/Users/Julie/R/functions/multiplot.r")

detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}




regularize_gappy_track <- function(x,y,datetime) {
  old_df <- data.frame(x,y,datetime)
  #time_idf if the difference between current and next time
  old_df[,'time_dif'] <- c(as.numeric(diff(old_df[,'datetime'])),NA)
  new_df <- data.frame(x=double(),y=double()) #empty dataframe to fill
  for (i in 1:(nrow(old_df))) {
    if(is.na(old_df[i,'time_dif'])) #include the last observation
      new_df <- rbind(new_df,old_df[i,c('x','y')])
    else if (old_df[i,"time_dif"] == 0) #remove duplicate observations
      next
    else if (old_df[i,'time_dif'] > 1 )
    {
      xx <- old_df[i:(i+1),'x']
      yy <- old_df[i:(i+1),'y']
      timedif <- old_df[i,'time_dif'] #i.e. 2
      tt = c(1,timedif+1)
      df <- data.frame(tt,xx,yy)
      
      mod <- lm(xx ~ tt,df)
      xs <- predict(mod,data.frame(tt=1:(timedif) ) )
      
      mod <- lm(yy ~ tt,df)
      ys <- predict(mod,data.frame(tt=1:(timedif) ) )
      
      new_df <- rbind(new_df,data.frame(x=xs,y=ys))
    }
    else #time_dif==1
      new_df <- rbind(new_df,old_df[i,c('x','y')])
  }
  new_df$OBJECTID <- 1:nrow(new_df)
  new_df
}




ScanTrack <- function(time,x,y=NULL, col=NULL, cex=NULL, ...)
{
  if(is.null(y)) if(is.complex(x)){y <- Im(x); x <- Re(x)} else if(ncol(x) == 2){y <- x[,2]; x <- x[,1]}
  
  if(is.null(col)) col=rgb(0,0,0,.5) 
  if(is.null(cex)) cex = 0.5
  
  layout(rbind(c(1,2), c(1,3)))
  plot(x,y,asp=1, type="o", pch=19, col=col, cex=cex, ...)
  plot(time,x, type="o", pch=19, col=col, xaxt="n", xlab="", cex=cex, ...)
  plot(time,y, type="o", pch=19, col=col, cex=cex, ...)
}

ScanTrack4<- function(time,x,y=NULL, col=NULL, cex=NULL, ...)
{
  if(is.null(y)) if(is.complex(x)){y <- Im(x); x <- Re(x)} else if(ncol(x) == 2){y <- x[,2]; x <- x[,1]}
  
  if(is.null(col)) col=rgb(0,0,0,.5) 
  if(is.null(cex)) cex = 0.5
  
  layout(rbind(c(1,2), c(1,3),c(1,4)))
  plot(x,y,asp=1, type="o", pch=19, col=col, cex=cex, ...)
  plot(time,x, type="o", pch=19, col=col, xaxt="n", xlab="", cex=cex, ...)
  plot(time,y, type="o", pch=19, col=col, cex=cex, ...)
  }


# last observation moved forward
# replaces all NA values with last non-NA values
na.lomf <- function(x) {
  
  na.lomf.0 <- function(x) {
    non.na.idx <- which(!is.na(x))
    if (is.na(x[1L])) {
      non.na.idx <- c(1L, non.na.idx)
    }
    rep.int(x[non.na.idx], diff(c(non.na.idx, length(x) + 1L)))
  }
  
  dim.len <- length(dim(x))
  
  if (dim.len == 0L) {
    na.lomf.0(x)
  } else {
    apply(x, dim.len, na.lomf.0)
  }
}

DaysfromStart <- function(timestamp) {
  (year(timestamp) - min(year(timestamp)) + yday(timestamp)/365)*365}

is.even <- function(x) x %% 2 == 0

TimeDiff<-function(timestamp){
  difftime(timestamp[-1], timestamp[-length(timestamp)], unit = "hour") %>% as.numeric  }

GetDisplacements<-function(timestamp,X,Y,n){
  dt <-TimeDiff(timestamp)
  dx <-diff(X)/dt
  dy <-diff(Y)/dt
  idxy <- which(abs(dx) > n|abs(dy) > n) #problems
  idxy
}


getDiffTable <- function(e, plotme = FALSE){
  Z <- with(e, x + 1i*y)
  Z1 <- Z[-length(Z)]
  Z2 <- Z[-1]
  dZ <- diff(Z)
  S <- Mod(dZ)
  dT <- with(e, difftime(DateTime[-1], DateTime[-length(DateTime)], unit = "hour")) %>% as.numeric
  Speed <- S/dT/1000
  local.mid <- with(e, DateTime[-nrow(e)] + 
                      difftime(DateTime[-1], DateTime[-length(DateTime)])/2)
  Hour <- hour(local.mid)
  year<-year(local.mid)
  day<-yday(local.mid)
  day.year<-interaction(day,year)
  e2 <- data.frame(ID = e$Tag[1], Z1, Z2, dZ, S, dT, Speed, 
                   local.mid,Hour,day.year)
  return(e2)
} 

long2UTM <- function(long) { (floor((long + 180)/6) %% 60) + 1}

#Half year 'season' function
halfyear<-function(df) {
  ifelse(df$day<200,paste(as.character(df$year),"A"), paste(as.character(df$year),"B")) #break down tracks into spring and fall by year
}

#FocalTrack function
FocalTrack<- function(df){
  df$halfyr<-halfyear(df)
  ddply(df,.(halfyr,day), function(df) c(x=mean(df$X)/1000, y = mean(df$Y)/1000))
}

SelectFPT<-function(fpt.obj,x,plot=TRUE){ 
  #x is threshold
  #only works for an equal number of crossup and cross down
  f1 <- fpt.obj$fpt[-nrow(fpt.obj)]
  f2 <- fpt.obj$fpt[-1]
  crossup <- which(f1 < x & f2 > x)
  crossdown <- which(f1 > x & f2 < x)
    if(plot==TRUE){
      plot(fpt.obj, type="l")
      abline(v = fpt.obj$time[crossup], col=2)
      abline(v = fpt.obj$time[crossdown], col=3)
    }
  t.start <- fpt.obj$time[crossup]
  t.end <- fpt.obj$time[crossdown]
  z<-sort(c(t.start,t.end))
  a<-diff(z)
  units(a)<-"hours"
  z.remove<-z[which(a<5)]
  z<-z[!z %in% z.remove]
  
  #if it is uneven, remove the time from the smallest interval
  if(!is.even(length(z))){
    z<-z[-which(diff(z)==min(diff(z)))[1]]
  }
  
    stops<-as.data.frame(matrix(z,ncol=2,byrow=T))
    stops<-data.frame("t.start"=as.POSIXct(stops$V1,origin=origin),
                      "t.end"=as.POSIXct(stops$V2,origin=origin))
    stops$test<-difftime(stops$t.start,stops$t.end,units="hours")
    stops
  }

LavielleEdits<-function(data,mig.TF,seg.length){
  #Use lavielle to edit migration or nonmigration segments 
  segs<-NULL
  df2<-data[data$migration==mig.TF,]
  df2$hr<-hour(df2$study.local.timestamp)
  for(i in unique(df2$cluster)){
    df3<-subset(df2,cluster==i&hr<20&hr>6)
    a=nrow(df3)
    #k=round(a/seg.length/2,digits=0)
    k=round(a/seg.length,digits=0)-2
      df.traj <- with(df3,as.ltraj(data.frame(X,Y), as.POSIXct(study.local.timestamp), id=tag.local.identifier))
      
      #k number of segments, of a length = seg.length
      ltest<-lavielle(df.traj,seg.length,k)
      lt <- findpath(ltest,  k,plotit=FALSE)
      #x<-lapply(lt, function (x) lapply(x,mean,na.rm=T))
      #x<-as.data.frame(matrix(unlist(x), length(lt), 10, byrow=T, dimnames=list(names(x), names(x[[1]]))),colnames=dimnames(x)[[2]])
      x<-data.frame("R2n"=unlist(lapply(lt,function(x) tail(x$R2n,n=1))),"n"=summary(lt)[3])
      x$mR2n<-x$R2n/x$n

      if(mig.TF==1){
       keep<-which(x$mR2n<1e+8)
      }else{
        keep<-which(x$mR2n>1e+8)
      }
      
      keep<-which(x$mR2n<1e+4)
      keep<-keep[which(keep %in% c(1,k))]
      segs<-rbind(segs,summary(lt)[keep,5:6])
    }
  segs
}


LavEditDF<-function(LavFindPath.output,data){
  #edit the original data frame with results from LavielleEdits
  start<-findInterval(data$study.local.timestamp,LavFindPath.output$date.begin)
  end<-findInterval(data$study.local.timestamp,(LavFindPath.output$date.end))
  phase<-start-end
  phase
}

#R2n= the squared distance from the first point of the trajectory

LavielleStopover<-function(data,seg.length=10,k=10){
  #Use lavielle to edit migration or nonmigration segments 
  df3<-data
    df.traj <- with(df3,as.ltraj(data.frame(X,Y), as.POSIXct(study.local.timestamp), id=tag.local.identifier))
    #if(nrow(df3)<seg.length*k){
    #  k=round(nrow(df3)/seg.length)-1
    #  if(k<3){
    #    seg.length=5
    #    k=10
    #    if(nrow(df3)<seg.length*k){
    #      k=round(nrow(df3)/seg.length)-1
    #    }
    #  }
    #}
    .lavielle<-lavielle(df.traj,seg.length,k)
    k2<-chooseseg(.lavielle,output="opt",draw=FALSE)
    if(k2==1){
      .lavielle<-lavielle(df.traj,5,k+1)
      k2<-chooseseg(.lavielle,output="opt",draw=FALSE)
    }
    if(k2>1){
      lt <- findpath(.lavielle,  k2,plotit=FALSE)
      
      x<-data.frame("R2n"=unlist(lapply(lt,function(x) tail(x$R2n,n=1))),"n"=summary(lt)[3])
      x$mR2n<-x$R2n/x$n
      
      keep<-which(x$mR2n<1e+9&x$R2n<9e+10)
      segs<-summary(lt)[keep,5:6]
      segs
    }else{
      segs<-NA
      segs
    }
}

#ChooseK<-function(.lavielle){
#  x<-chooseseg(.lavielle,draw=FALSE)
#  k<-which(abs(diff(x$Jk))==min(abs(diff(x$Jk))))
#}
  



segDiffs <- function(e){
  Z <- with(e, X + 1i*Y)
  Z1 <- Z[length(Z)]
  Z2 <- Z[1]
  dZ <- Z1-Z2
  Dist <- Mod(dZ)
  dT <- with(e, difftime(study.local.timestamp[length(study.local.timestamp)], study.local.timestamp[1], unit = "hour")) %>% as.numeric
  Speed <- Dist/dT
  e2 <- data.frame(Dist, dT, Speed)
  return(e2)
}

suncalc<-function(d,time,Lat=48.1442,Long=-122.7551){
  ## d is the day of year
  ## Lat is latitude in decimal degrees
  ## Long is longitude in decimal degrees (negative == West)
  
  ##This method is copied from:
  ##Teets, D.A. 2003. Predicting sunrise and sunset times.
  ##  The College Mathematics Journal 34(4):317-321.
  
  ## At the default location the estimates of sunrise and sunset are within
  ## seven minutes of the correct times (http://aa.usno.navy.mil/data/docs/RS_OneYear.php)
  ## with a mean of 2.4 minutes error.
  
  ## Function to convert degrees to radians
  rad<-function(x)pi*x/180
  
  ##Radius of the earth (km)
  R=6378
  
  ##Radians between the xy-plane and the ecliptic plane
  epsilon=rad(23.45)
  
  ##Convert observer's latitude to radians
  L=rad(Lat)
  
  ## Calculate offset of sunrise based on longitude (min)
  ## If Long is negative, then the mod represents degrees West of
  ## a standard time meridian, so timing of sunrise and sunset should
  ## be made later.
  timezone = -4*(abs(Long)%%15)*sign(Long)
  
  ## The earth's mean distance from the sun (km)
  r = 149598000
  
  theta = 2*pi/365.25*(d-80)
  
  z.s = r*sin(theta)*sin(epsilon)
  r.p = sqrt(r^2-z.s^2)
  
  t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))
  
  ##a kludge adjustment for the radius of the sun
  that = t0+5 
  
  ## Adjust "noon" for the fact that the earth's orbit is not circular:
  n = 720-10*sin(4*pi*(d-80)/365.25)+8*sin(2*pi*d/365.25)
  
  ## now sunrise and sunset are:
  sunrise = (n-that+timezone)/60
  sunset = (n+that+timezone)/60
  
  require(lubridate)
  daylight = ifelse(hour(time)>=sunrise,ifelse(hour(time)<=sunset,1,0),0)
  
  return(daylight)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
