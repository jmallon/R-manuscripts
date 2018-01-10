## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
#rm(list=ls())

## Load host-dependent directory environment
winos <- ifelse(grepl("windows", Sys.info()['sysname'], ignore.case=T), 1, 0)
if(winos==1) source("C:/Users/Julie/Box Sync/R/Chapter 1/Scripts/Functions/file_dir_params.R")
if(winos==0) source("~/Box Sync/R/Chapter 1/Scripts/Functions/file_dir_params.R")
rm(winos)

## ----libraries,include=FALSE---------------------------------------------
library(here)
source(paste0(prjfuns, "chapter 1 functions.r"))
#load(here(paste0(prjdta, "dist_season_annotate.robj")))

library(dplyr)
#library(MASS)
#library(scales)
#library(zoo)
#library(raster)
library(ggplot2)
#library(gridExtra)

## ----Raw Data,eval=FALSE,echo=FALSE--------------------------------------
## load(here(paste0(prjdta,"stopovers_with_all_weather2.robj")))
## 
##   #add annotations about where along the route it was
## sprng <- subset(stops_full_weather2,Season=="Spring") %>% as.data.frame
##   mins <- aggregate(study.local.timestamp~ID_migration, sprng ,min)[,]
##   maxs <- aggregate(study.local.timestamp~ID_migration, sprng ,max)[,]
##   mins<-merge(sprng,mins) %>% subset(select=c(ID_migration,X,Y))
##   names(mins)<-c("ID_migration","winter.X","winter.Y")
##   maxs<-merge(sprng,maxs) %>% subset(select=c(ID_migration,X,Y))
##   names(maxs)<-c("ID_migration","summer.X","summer.Y")
##   sprng = plyr::join(sprng,mins,by='ID_migration')
##   sprng = plyr::join(sprng,maxs,by='ID_migration')
## 
## fall <- subset(stops_full_weather2,Season=="Fall") %>% as.data.frame
##   mins <- aggregate(study.local.timestamp~ID_migration, fall ,min)
##   maxs <- aggregate(study.local.timestamp~ID_migration, fall ,max)
##   mins<-merge(fall,mins) %>% subset(select=c(ID_migration,X,Y))
##   names(mins)<-c("ID_migration","summer.X","summer.Y")
##   maxs<-merge(fall,maxs) %>% subset(select=c(ID_migration,X,Y))
##   names(maxs)<-c("ID_migration","winter.X","winter.Y")
##   fall = plyr::join(fall,mins,by='ID_migration')
##   fall = plyr::join(fall,maxs,by='ID_migration')
## 
## stop.weather2.dists<-rbind(sprng,fall)
## 
## dist.winter <- function(x){
##   my_vector <- vector()
##   for( i in 1:length(x$X)){
##     d <- pointDistance(cbind(x[i,'X'], x[i,'Y']),
##                        cbind(x[i,'winter.X'], x[i,'winter.Y']), lonlat=FALSE)
##     my_vector <- c(my_vector, d)
##   }
##   my_vector
## }
## dist.summer <- function(x){
##   my_vector <- vector()
##   for( i in 1:length(x$X)){
##     d <- pointDistance(cbind(x[i,'X'], x[i,'Y']),
##                        cbind(x[i,'summer.X'], x[i,'summer.Y']), lonlat=FALSE)
##     my_vector <- c(my_vector, d)
##   }
##   my_vector
## }
## 
## stop.weather2.dists$winter.dist<-dist.winter(stop.weather2.dists)
## stop.weather2.dists$summer.dist<-dist.summer(stop.weather2.dists)
## stop.weather2.dists<-stop.weather2.dists %>% group_by(ID_migration) %>%
##                         mutate(summer.dist.p=summer.dist/max(summer.dist),
##                                winter.dist.p=winter.dist/max(winter.dist))
## stop.weather2.dists<-as.data.frame(stop.weather2.dists)
## stop.weather2.dists$dist.category<-with(stop.weather2.dists,ifelse(summer.dist.p<0.25,"summer",ifelse(summer.dist.p<0.75,"middle","winter")))
## save(stop.weather2.dists,file=paste0(prjdta, "raw_weather_distances_df.robj"))
## 
## raw_weather_data<-stop.weather2.dists[,c(37:38,42:45,47:48,50:51,54:61,67:70, 23:25,1:8,15,52, 77:79)]
## raw_weather_data<-subset(raw_weather_data,roost==0)
## raw_weather_data<-subset(raw_weather_data,!is.na(NDVI_16d)&!is.na(thermal))
## annotations2<-raw_weather_data[,c(20:38)]
## raw_weather_data2<-raw_weather_data[,1:19]
## raw_weather_data2<- data.frame(lapply(raw_weather_data2, function(x) scale(x, center = F)))
## annotated_raw_data<-data.frame(raw_weather_data2,annotations2)
## save(annotated_raw_data,file=paste0(prjdta, "dist_season_raw_weather_annotate.robj"))

## ----Stopover Freq Speed1,eval=TRUE,echo=FALSE,warning=FALSE,message=FALSE----
load(here(paste0(prjdta, "dist_season_raw_weather_annotate.robj")))
load(here(paste0(prjdta, "raw_weather_distances_df.robj")))
stop.weather2.dists <- stop.weather2.dists[!duplicated(stop.weather2.dists),]
stop.weather2.dists<-subset(stop.weather2.dists,ID_migration!=257)
pop.summary<-stop.weather2.dists %>% group_by(pop,ID_migration) %>% filter(roost==0)  %>% mutate(dist=max(winter.dist),p.stop=mean(stopover),duration=sum(dT,na.rm=T),speed=dist[]/duration[]/1000) %>% group_by(pop,ID_migration,stop.cat) %>% mutate(cat.dur=sum(dT,na.rm=T)/duration) %>% summarize(prop.tort=cat.dur[1],dist=dist[1],p.stop=p.stop[1],duration=duration[1],speed=speed[1],p.tort.stops=prop.tort[]/p.stop[]) %>% group_by(pop,ID_migration) %>% filter(stop.cat=="tortuous") #missing migrations without any tort. keep Na instead?

#pop.sum2<-pop.summary %>% group_by(pop) %>% summarize(dist=mean(dist),stopover=mean(stopover),duration=mean(duration),speed=dist[]/duration[]/1000)

## ----Stopover Freq Speed2,eval=TRUE,message=FALSE,warning=FALSE,echo=FALSE----
colrs <- c("aura"="white", "meridionalis"="black", "ruficollis"="red")
p.speed<-ggplot(pop.summary,aes(y=p.stop,x=speed))+geom_point(aes(fill=pop),color="black",shape=21)+
      stat_smooth(method = "glm", method.args = list(family = "binomial"), size = 1,se=FALSE)+
theme_bw()+theme(legend.position = c(0.85,0.85),
                   legend.background = element_rect(color="gray", size=0.5, 
                                                    linetype="solid"),
                    # remove the vertical grid lines
                    panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                    # explicitly set the horizontal lines (or they will disappear too)
                     panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank())+
  scale_fill_manual(values = colrs, name="Race",
                         breaks=c("aura", "meridionalis", "ruficollis"),
                         labels=c("Aura", "Meridionalis", "Ruficollis"))+
  xlab("Mean Hourly Speed (km/hr)")+ylab("Mean Proportion of Time at Stopovers")
p.speed

#ggplot(pop.summary,aes(y=speed,x=p.tort.stops))+geom_point(aes(color=pop))+
#      stat_smooth(method = "glm", method.args = list(family = "binomial"), size = 1,se=FALSE)+
#theme_bw() goes past 1

speed.binom<-glm(p.stop~speed,data=pop.summary,family=binomial)
#summary(speed.binom)
speed.coef<-summary(speed.binom)$coef[2,] %>% as.numeric
1 - summary(speed.binom)$deviance/summary(speed.binom)$null.deviance #R^2
library(ResourceSelection)
hoslem.test(pop.summary$p.stop, fitted(speed.binom)) #fits really well!

## ----Raw Models,eval=TRUE,echo=FALSE,message=FALSE,warning=FALSE---------
library(MASS)
#with RAW data
stop.weather2.dists$non.tort.stops<-with(stop.weather2.dists,ifelse(stopover==1,  ifelse(stop.cat=="tortuous",0,1),0))

fit.global <- lda(non.tort.stops ~ accum_precip_m+ thermal+tailwind+tot_atm_water+
              temp+ orographic+ downward_LW+ downward_SW,
                   data=stop.weather2.dists, na.action = "na.omit")

#center, then divide within group standard deviation
colz<-c(38,44, 50:51,54,56,61,67)
c.data<-stop.weather2.dists %>% mutate_at(colz, function(x) c(scale(x, center = T,scale=F))) 
grp.stdev<-c.data %>% group_by(non.tort.stops) %>% summarise_at(colz,function(x) sd(x, na.rm=T))
c.data0<-subset(c.data,non.tort.stops==0)
c.data1<-subset(c.data,non.tort.stops==1)
c.data0<- c.data0 %>%   mutate(temp=temp/grp.stdev$temp[1],orographic=orographic/grp.stdev$orographic[1],
         thermal=thermal/grp.stdev$thermal[1],
         tot_atm_water=tot_atm_water/grp.stdev$tot_atm_water[1],
         downward_LW=downward_LW/grp.stdev$downward_LW[1],      
         accum_precip_m=accum_precip_m/grp.stdev$accum_precip_m[1],
         tailwind=tailwind/grp.stdev$tailwind[1],
         downward_SW=downward_SW/grp.stdev$downward_SW[1])
c.data1<- c.data1 %>%   mutate(temp=temp/grp.stdev$temp[2],orographic=orographic/grp.stdev$orographic[2],
         thermal=thermal/grp.stdev$thermal[2],
         tot_atm_water=tot_atm_water/grp.stdev$tot_atm_water[2],
         downward_LW=downward_LW/grp.stdev$downward_LW[2],      
         accum_precip_m=accum_precip_m/grp.stdev$accum_precip_m[2],
         tailwind=tailwind/grp.stdev$tailwind[2],
         downward_SW=downward_SW/grp.stdev$downward_SW[2])
scaled.data<-rbind(c.data0,c.data1)
  #experience
  exper<-scaled.data %>% group_by(individual.local.identifier,ID_migration) %>% summarize(y=length(unique(year))) %>% mutate(experience=as.factor(cumsum(y)))
  scaled.data<-merge(scaled.data,exper[,c(1,2,4)],by=c("individual.local.identifier","ID_migration"))

#refit with scaled data
fit.global2 <- lda(non.tort.stops ~ accum_precip_m+ thermal+tailwind+tot_atm_water+
              temp+ orographic+ downward_LW+ downward_SW,
                   data=scaled.data)
          ax1.loadings<- fit.global2$scaling[,1]
          
          
  plda <- predict(object = fit.global2,
                newdata = scaled.data)
global.df2<-data.frame(LD1=plda$x[,1],PP=plda$posterior[,2],
                pop=scaled.data$pop,
                Season=scaled.data$Season,
                stop.id=scaled.data$stop.id,
                       stopover=scaled.data$non.tort.stops,
                       dist=scaled.data$dist.category,
                stop.cat=scaled.data$stop.cat,experience=scaled.data$experience)

## ----Raw Plots Pop,eval=TRUE,echo=FALSE----------------------------------
ld.glm <- glm(stopover ~ LD1*pop, data=global.df2, family=binomial)

XLD <- seq(from=min(global.df2$LD1,na.rm=T), to=max(global.df2$LD1,na.rm=T), by=.1)
  loadings.axis<- ax1.loadings
  loadings.axis[which(ax1.loadings>0)]<-0.9
  loadings.axis[which(ax1.loadings<0)]<-0.1
  
#jpeg('../Plots/tol_race.jpg')

plot(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",LD1=XLD),type="response"),
type="o",pch=21,bg="white",ylim=c(0,1),ylab="Probability of Stopover",xlab="Multivariate Weather")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",LD1=XLD),type="response"),
type="o",pch=21,bg="black")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",LD1=XLD),type="response"),
type="o",pch=21,bg="red")

mtext("precipitation",side=1,col="red",outer=T,line=-1,at=as.numeric( loadings.axis[1]))
mtext("thermal",side=1,col="red",outer=T,line=-1,at=as.numeric( loadings.axis[2]))
mtext("tailwind",side=1,col="red",outer=T,line=-2,at=as.numeric( loadings.axis[3]))
mtext("total atm. water",side=1,col="red",outer=T,line=-2,at=as.numeric( loadings.axis[4]))
mtext("temperature",side=1,col="red",outer=T,line=-3,at=as.numeric( loadings.axis[5]))
mtext("orographic",side=1,col="red",outer=T,line=-4,at=as.numeric( loadings.axis[6]))
mtext("downward LW",side=1,col="red",outer=T,line=-3,at=as.numeric( loadings.axis[7]))
mtext("downward SW",side=1,col="red",outer=T,line=-5,at=as.numeric( loadings.axis[8]))


legend(x="topleft",pch=21,col=c("black","black","black"),
       pt.bg=c("white","black","red"),
legend=c("aura","meridionalis","ruficollis"))

p.tol.race <- recordPlot()
dev.off()

## ----Raw Plots Season,echo=FALSE-----------------------------------------
par(mfrow=c(2,2))
ld.glm <- glm(stopover ~ LD1*Season, data=global.df2, family=binomial)

plot(XLD,
predict(ld.glm,newdata=data.frame(Season="Spring",LD1=XLD),type="response"),
type="o",pch=21,bg="white",ylab="Probability of Stopover",
xlab="Multivariate Weather",ylim=c(0,1),main="species")
points(XLD,
predict(ld.glm,newdata=data.frame(Season="Fall",LD1=XLD),type="response"),
type="o",pch=21,bg="black")
legend(x="topleft",pch=21,col=c("black","black"),pt.bg=c("white","black"),
legend=c("Spring","Fall"))



###Population and Season interaction

ld.glm <- glm(stopover ~ LD1*pop*Season, data=global.df2, family=binomial)


#aura
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1),
main="aura")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1))


#ruficollis
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",ylab="",ylim=c(0,1),
main="ruficollis")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1))


#meridionalis
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1),
main="meridionalis")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1))

p.tol.season <- recordPlot()
dev.off()

## ------------------------------------------------------------------------
ld.glm <- glm(stopover ~ LD1*pop*experience, data=global.df2, family=binomial)

#aura
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",LD1=XLD,experience="1"),
        type="response"),type="o",pch=21,bg="white",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1),
main="aura")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",experience="2",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",experience="3",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1))

#ruficollis
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",LD1=XLD,experience="1"),
        type="response"),type="o",pch=21,bg="white",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1),
main="ruficollis")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",experience="2",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",experience="3",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1))

#meridionalis
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",LD1=XLD,experience="1"),
        type="response"),type="o",pch=21,bg="white",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1),
main="meridionalis")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",experience="2",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",experience="3",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",
ylab="Probability of Stopover",xlab="Multivariate Weather",ylim=c(0,1))

#need to 1) make experience years instead of migrations and 2) properly subset so only birds that fit all the years are in the model


## ----Raw Plots Dist,echo=FALSE-------------------------------------------
par(mfrow=c(2,2))
ld.glm <- glm(stopover ~ LD1*dist, data=global.df2, family=binomial)

plot(XLD,
predict(ld.glm,newdata=data.frame(dist="summer",LD1=XLD),type="response"),
type="o",pch=21,bg="white",ylab="",ylim=c(0,1),main="species")
points(XLD,
predict(ld.glm,newdata=data.frame(dist="middle",LD1=XLD),type="response"),
type="o",pch=21,bg="black")
points(XLD,
predict(ld.glm,newdata=data.frame(dist="winter",LD1=XLD),type="response"),
type="o",pch=21,bg="red")
legend(x="topleft",pch=21,col=c("black","black","black"),pt.bg=c("white","black","red"),
legend=c("summer","middle","winter"))


###Population and Distance interaction
ld.glm <- glm(stopover ~ LD1*pop*dist, data=global.df2, family=binomial)

plot(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",dist="summer",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",ylab="",ylim=c(0,1),
main="meridionalis")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",dist="middle",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",ylab="",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",dist="winter",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",ylab="",ylim=c(0,1))

plot(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",dist="summer",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",ylab="",ylim=c(0,1),
main="aura")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",dist="middle",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",ylab="",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",dist="winter",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",ylab="",ylim=c(0,1))

plot(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",dist="summer",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",ylab="",ylim=c(0,1),
main="ruficollis")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",dist="middle",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",ylab="",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",dist="winter",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",ylab="",ylim=c(0,1))

p.tol.dist <- recordPlot()
dev.off()


## ----Raw Full Interactions meridionalis, eval=TRUE,echo=FALSE------------
par(mfrow=c(1,2))
ld.glm <- glm(stopover ~ LD1*pop*dist*Season, data=global.df2, family=binomial)

#spring
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",dist="summer",
                                  Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",ylab="",ylim=c(0,1),
main="meridionalis - spring")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",dist="middle",
                                  Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",ylab="",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",dist="winter",
                                  Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",ylab="",ylim=c(0,1))
legend(x="topleft",pch=21,col=c("black","black","black"),pt.bg=c("white","black","red"),
legend=c("summer","middle","winter"))
#FALL
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",dist="summer",
                                  Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",ylab="",ylim=c(0,1),
main="meridionalis - Fall")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",dist="middle",
                                  Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",ylab="",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="meridionalis",dist="winter",
                                  Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",ylab="",ylim=c(0,1))
legend(x="topleft",pch=21,col=c("black","black","black"),pt.bg=c("white","black","red"),
legend=c("summer","middle","winter"))

p.tol.m <- recordPlot()
dev.off()

## ----aura,eval=TRUE,echo=FALSE-------------------------------------------
par(mfrow=c(1,2))
#spring
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",dist="summer",
                                  Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",ylab="",ylim=c(0,1),
main="aura - spring")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",dist="middle",
                                  Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",ylab="",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",dist="winter",
                                  Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",ylab="",ylim=c(0,1))
legend(x="topleft",pch=21,col=c("black","black","black"),pt.bg=c("white","black","red"),
legend=c("summer","middle","winter"))

#fall
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",dist="summer",
                                  Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",ylab="",ylim=c(0,1),
main="aura - Fall")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",dist="middle",
                                  Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",ylab="",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="aura",dist="winter",
                                  Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",ylab="",ylim=c(0,1))
legend(x="topleft",pch=21,col=c("black","black","black"),pt.bg=c("white","black","red"),
legend=c("summer","middle","winter"))

p.tol.a <- recordPlot()
dev.off()

## ----ruficollis,eval=TRUE,echo=FALSE-------------------------------------
par(mfrow=c(1,2))
#spring
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",dist="summer",
                                  Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",ylab="",ylim=c(0,1),
main="ruficollis - spring")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",dist="middle",
                                  Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",ylab="",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",dist="winter",
                                  Season="Spring",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",ylab="",ylim=c(0,1))
legend(x="topleft",pch=21,col=c("black","black","black"),pt.bg=c("white","black","red"),
legend=c("summer","middle","winter"))

#fall
plot(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",dist="summer",
                                  Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="white",ylab="",ylim=c(0,1),
main="ruficollis - Fall")
points(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",dist="middle",
                                  Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="black",ylab="",ylim=c(0,1))
points(XLD,
predict(ld.glm,newdata=data.frame(pop="ruficollis",dist="winter",
                                  Season="Fall",LD1=XLD),
        type="response"),type="o",pch=21,bg="red",ylab="",ylim=c(0,1))
legend(x="topleft",pch=21,col=c("black","black","black"),pt.bg=c("white","black","red"),
legend=c("summer","middle","winter"))

p.tol.r <- recordPlot()
dev.off()


## ----output,eval=FALSE,include=FALSE-------------------------------------
## library(knitr)
## knit('Scripts/risk aversion.Rmd', output='Scripts/risk aversion.R', tangle=TRUE)

