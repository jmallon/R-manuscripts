## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
#detachAllPackages()
#rm(list=ls())

## Load host-dependent directory environment
winos <- ifelse(grepl("windows", Sys.info()['sysname'], ignore.case=T), 1, 0)
if(winos==1) source("C:/Users/Julie/Box Sync/R/Chapter 1/Scripts/Functions/file_dir_params.R")
if(winos==0) source("~/Box Sync/R/Chapter 1/Scripts/Functions/file_dir_params.R")
rm(winos)

## ----libraries,include=FALSE---------------------------------------------
library(here)
source(paste0(prjfuns, "chapter 1 functions.r"))

## ----create change.weather2 df,eval=FALSE,include=FALSE------------------
## #library(lubridate)
## change.weather<-subset(stops_full_weather2,!is.na(tot_atm_water))[,c(1:38,40,42:48,50:61,63,67:70)]
## #change.weather$tot_atm_water<-as.numeric(change.weather$tot_atm_water)
## #change.weather<-change.weather %>%  mutate(hour=hour(local.time))
## #detach(package:lubridate)
## change.weather<-change.weather[order(change.weather$ID_migration,change.weather$local.time),]
## columns<-c(36:47,50:59) # minus 1 than I want because grouping removes a column as an option??
## change.weather2<-change.weather %>% group_by(ID_migration) %>% mutate_at(columns,function(x) c(NA,diff(x))) %>% as.data.frame() #find the difference between each observation
## columns<-c(37:48,51:60) # minus 1 than I want because grouping removes a column as an option??
## 
## change.weather2<-data.frame(change.weather2[,c(1:35,49,61:63)],sapply(change.weather2[,columns],function(.x) .x/change.weather2$dT)) #make this a rate, over time
## 
## library(here)
## save(change.weather2,file=paste0(here(prjdta,"change_in_weather.robj")))

## ----where peak change occurs,eval=FALSE,include=FALSE-------------------
## change.weather2<-change.weather2[order(change.weather2$ID_migration,change.weather2$local.time),]
## change.weather2<-change.weather2 %>% group_by(stop.id) %>%
##   mutate(ts.stopover=as.numeric(local.time-min(local.time))/3600)
## change.weather.s<-change.weather2 %>% group_by(stop.id) %>%
##   mutate(ts.stopover=ts.stopover/max(ts.stopover,na.rm=T))
## #create a start column (1 and 0), and a diffstart column where all cells are difference from the next start
## #same for end column
## 
## cw2<-subset(change.weather.s,stopover==1)

## ----create dflag df,eval=FALSE,include=FALSE----------------------------
## #summarize each stopover to use as reference
## lag.id.df<-change.weather.s %>% filter(stopover>0) %>% group_by(ID_migration,stop.id)%>% summarise(n=n(),start=local.time[1],stop=local.time[length(local.time)])
## lag.id.df$lag.id<-lag.id.df$stop.id
## 
## #function to find points between start and stop (with a 6 hr buffer)
## intervalDF<-function(times,data){
##     start<-findInterval(data2$local.time,(times$start-21600))
##     end<-findInterval(data2$local.time,(times$stop+21600))
##     phase<-start-end
##     phase
##   }
## 
## data1<-change.weather.s
## dflag<-NULL
## ids<-unique(lag.id.df$ID_migration)
## for(j in ids){
##   data2<-subset(data1,ID_migration==j)
##   lag.sub<-subset(lag.id.df,ID_migration==j)
##   for(i in 1:nrow(lag.sub)){
##       data3<-data2
##       data3$lagseg<-intervalDF(lag.sub[i,],data2)
##       data4<-subset(data3,lagseg==1)
##       data4$lag.id<-lag.sub$stop.id[i] #unique id
##       dflag<-rbind(data4,dflag) #if lagsegs overlap, then they will appear twice, with unique lag ids
##   }
## }
## 
## #now, find time difference from each point to the start of the lag id and the end of the lag id
## dflag<-merge(dflag,lag.id.df[,c(4:6)],by="lag.id") #create start and end times cols
## dflag$diff.start<-with(dflag,as.numeric(local.time-start)/3600)
## dflag$diff.end<-with(dflag,as.numeric(local.time-stop)/3600)
## 
## library(here)
## save(dflag,file=here(paste0(prjdta, "timeseries.robj"))) #raw weather data
## 
## 

## ----echo=FALSE,warning=FALSE,message=FALSE------------------------------
library(dplyr)
load(here(paste0(prjdta,"change_in_weather.robj")))
load(here(paste0(prjdta, "timeseries.robj")))

n.tort<-dflag[which(dflag$stop.cat=="tortuous"),'stop.id'] %>% unique
n.tort.df<-dflag[which(!(dflag$lag.id %in% n.tort)),] #remove tortuous stopovers to get that weather signal more clearly

library(dplyr)
  library(data.table)
vars=c("surface_roughness", "temp" , "wind_speed" , "high_veg", "DEM" , "orographic",  "NDVI_16d",  "rough_length", "dew_temp" ,  "surface_pressure" , "tailwind" , "thermal"  , "tot_atm_water",
       "downward_LW" , "downward_SW" ,  "sensible_HF","latent_HF", "boundary_height", "accum_precip_m" )

#v.start<-n.tort.df[,c(37:38,40:62,67)] %>% filter(complete.cases(.),diff.start<7) %>%  group_by(diff.start)  %>% summarize_at(vars,mean)

v.start<-n.tort.df[,c(37:38,40:62,67:68)] %>% filter(complete.cases(.),diff.start<7,diff.start>-8) %>%  group_by(diff.start)  %>% summarize_at(vars,mean)

library(tidyr)
library(ggplot2)

v.start %>%
  gather(-diff.start, key = "var", value = "value")%>% 
  ggplot(aes(y = value, x = diff.start)) +
    geom_smooth(se=FALSE) +
    facet_wrap(~var, scales = "free") +
    geom_vline(aes(xintercept=0),col="red")+
    theme_bw()+ggtitle("Start of Stopover")



#Start of stopovers associated with:
#Landscape: increase in elevation `r mds[5,3]` and high vegetation `r mds[4,3]` andndvi `r mds[8,3]` hours from the start of the stopover

#Weather: maximum dew temp `r mds[10,3]`, orographic `r mds[7,3]`, total atmospheric water `r mds[14,3]`, and precipitation `r mds[6,3]` hours from the start of a stopover

#Minimum  surface pressure `r mds[11,4]` , tailwinds `r mds[12,4]` ,temperature `r mds[2,4]`, thermal  intensity`r mds[13,4]`, and wind speed `r mds[3,4]` hours from the start of a stopover


## ----echo=FALSE,warning=FALSE,message=FALSE------------------------------


#v.end<-dflag[,c(37:38,40:62,68)] %>% filter(complete.cases(.),diff.end>-7) %>% group_by(diff.end) %>% summarize_at(vars,mean) 


v.end<-n.tort.df[,c(37:38,40:62,67:68)] %>% filter(complete.cases(.),diff.end>-7) %>% group_by(diff.end) %>% summarize_at(vars,mean) 

#plot all together
v.end %>%
  gather(-diff.end, key = "var", value = "value")%>% 
  ggplot(aes(y = value, x = diff.end)) +
    geom_smooth(se=FALSE) +
    facet_wrap(~var, scales = "free") +
  geom_vline(aes(xintercept=0),col="red")+
    theme_bw()+ggtitle("End of Stopover")

#install.packages("tidyr",repos="https://cran.rstudio.com/")



#End of stopovers associated with:
#no maxima or minima of landscape variables

#Weather: maximum orographic `r eds[5,3]`, dew temp `r eds[2,3]`,  temperature `r eds[10,3]`, wind speed `r eds[3,3]`, tailwind `r eds[12,3]`,  surface pressure `r eds[8,3]` and thermal intensity `r eds[12,3]` hours from the end of a stopover

#Minimum  total atmospheric water `r eds[14,4]`, and precipitation `r eds[6,4]` hours from the end of the stpover

## ----difference between start and end variables, echo=FALSE,message=FALSE----
d1<-v.start %>%
  gather(-diff.start, key = "var", value = "value") %>% rename(diff.hr=diff.start,start.val=value)
d2<-v.end %>%
  gather(-diff.end, key = "var", value = "value") %>% rename(diff.hr=diff.end,end.val=value)
d12<-merge(d1,d2,by=c("diff.hr","var"))
d12$diff.val<-d12$start.val-d12$end.val

d12 %>%
  ggplot(aes(y = diff.val, x = diff.hr)) +
    geom_smooth(se=FALSE) +
    facet_wrap(~var, scales = "free") +
  geom_vline(aes(xintercept=0),col="red")+
    theme_bw()+ggtitle("Difference between start and end of stopovers")

## ----start of stopover model, echo=FALSE,message=FALSE-------------------
library(lme4)

lag.start<-subset(n.tort.df,diff.start<7)
lag.s<-lag.start[,39:62]
lag.s<- data.frame(lapply(lag.s, function(x) scale(x, center = T)))
lag.start<-cbind(as.data.frame(lag.start[,c(1:38,63:68)]),lag.s)


lag.before<-subset(lag.start,diff.start<1&!is.na(thermal)&!is.na(NDVI_16d))
lag.after<-subset(lag.start,diff.start>-1&!is.na(thermal))
                  
b1<-  lmer(diff.start~#tot_atm_water+
            # surface_pressure+
             thermal+#temp+ #wind_speed+
             tailwind+   #NDVI_16d+
             dew_temp+high_veg+#surface_roughness+
             downward_LW+#downward_SW+ sensible_HF+#latent_HF+
             boundary_height+accum_precip_m+
             (1|individual.local.identifier),data=lag.before)
#drop1(b1)

summary(b1)$coef
#AIC(b1,b2,b3,b4) verified best random structure



## ----eval=T,echo=FALSE---------------------------------------------------

lag.end<-subset(n.tort.df,diff.end>-7)
lag.e<-lag.end[,39:62]
lag.e<- data.frame(lapply(lag.e, function(x) scale(x, center = T)))
lag.end<-cbind(as.data.frame(lag.end[,c(1:38,63:68)]),lag.e)

lag.e.before<-subset(lag.end,diff.end<1&!is.na(thermal))

b2<-  lmer(diff.end~tot_atm_water+
            surface_pressure+
             temp+ wind_speed+
             dew_temp+high_veg+
             sensible_HF+latent_HF+
             accum_precip_m++(1|individual.local.identifier),data=lag.e.before)
#drop1(b2)

summary(b2)$coef


## ----eval=FALSE,echo=FALSE-----------------------------------------------
##   la1<-subset(lag.e.after,pop=="aura")
##   lb1<-subset(lag.e.before,pop=="aura")
##   la2<-subset(lag.after,pop=="aura")
##   lb2<-subset(lag.before,pop=="aura")
## 
##   #lam1<-lmer(diff.end~surface_pressure+precip_fraction+thermal+temp+wind_speed+ tailwind+dew_temp+high_veg+
##              #  (1|individual.local.identifier),data=la1,na.action="na.omit")
##   #drop1(lam1)
## 
##   am2<-lmer(diff.end~tot_atm_water+surface_pressure+precip_fraction+thermal+temp+wind_speed+
##                orographic+rough_length+
##                (1|individual.local.identifier),data=lb1,na.action="na.omit")
##  # drop1(am2)
## 
##     #lam3<-lmer(diff.start~tot_atm_water+surface_pressure+precip_fraction+thermal+temp+ +dew_temp+
##               # (1|individual.local.identifier),data=la2,na.action="na.omit")
##     #drop1(lam3)
## 
##     am1<-lmer(diff.start~surface_pressure+precip_fraction+thermal+temp+wind_speed+ dew_temp+high_veg+
##                (1|individual.local.identifier),data=lb2,na.action="na.omit")
##     drop1(lam4)
## 
##  #meridionalis
##   lm2<-subset(lag.e.before,pop=="meridionalis")
##   lm1<-subset(lag.before,pop=="meridionalis")
## 
## 
##       mm1<-lmer(diff.start~surface_pressure+precip_fraction+thermal+temp+wind_speed+
##                tailwind+NDVI_16d+high_veg+
##                (1|individual.local.identifier),data=lm1,na.action="na.omit")
##       drop1(mm1)
##     mm2<-lmer(diff.end~tot_atm_water+surface_pressure+precip_fraction+temp+wind_speed+
##                tailwind+dew_temp+rough_length+
##                (1|individual.local.identifier),data=lm2,na.action="na.omit")
##     drop1(mm2)
## 
## 
##      #ruficollis
##   lr2<-subset(lag.e.before,pop=="ruficollis")
##   lr1<-subset(lag.before,pop=="ruficollis")
##             rm1<-lmer(diff.start~tot_atm_water+surface_pressure+thermal+temp+
##                +NDVI_16d+dew_temp+high_veg+
##                (1|individual.local.identifier),data=lr1,na.action="na.omit")
##             drop1(rm1)
##             rm2<-lmer(diff.end~precip_fraction+temp+wind_speed+
##                +dew_temp+
##                (1|individual.local.identifier),data=lr2,na.action="na.omit")
##             drop1(rm2)
##    # sjt.lmer(rm1, rm2)
## 
##     #all before start
##    # sjt.lmer(b1,am1,mm1,rm1, show.header = TRUE, string.dv="Before start of stopovers", depvar.labels  =c("Species", "Aura", "Meridionalis", "Ruficollis"))
## 
##         #all before end
##   #  sjt.lmer(b2,am2,mm2,rm2, show.header=TRUE,string.dv="Before end of stopovers", depvar.labels  =c("Species", "Aura", "Meridionalis", "Ruficollis"))

## ---- eval=FALSE,echo=FALSE----------------------------------------------
## lag.id.df<-stop.weather %>% filter(stopover>0) %>% group_by(ID_migration,stop.id)%>% summarise(n=n(),start=local.time[1],stop=local.time[length(local.time)])
## lag.id.df$lag.id<-lag.id.df$stop.id
## 
## #function to find points between start and stop (with a 6 hr buffer)
## intervalDF<-function(times,data){
##     start<-findInterval(data2$local.time,(times$start-28800))
##     end<-findInterval(data2$local.time,(times$stop+28800))
##     phase<-start-end
##     phase
##   }
## 
## data1<-stop.weather
## sw.diff.h<-NULL
## ids<-unique(lag.id.df$ID_migration)
## for(j in ids){
##   data2<-subset(data1,ID_migration==j)
##   lag.sub<-subset(lag.id.df,ID_migration==j)
##   for(i in 1:nrow(lag.sub)){
##       data3<-data2
##       data3$lagseg<-intervalDF(lag.sub[i,],data2)
##       data4<-subset(data3,lagseg==1)
##       data4$lag.id<-lag.sub$stop.id[i] #unique id
##       sw.diff.h<-rbind(data4,sw.diff.h) #if lagsegs overlap, then they will appear twice, with unique lag ids
##   }
## }
## 
## #now, find time difference from each point to the start of the lag id and the end of the lag id
## sw.diff.h<-plyr::join(sw.diff.h,lag.id.df[,c(4:6)],by="lag.id") #create start and end times cols
## sw.diff.h$diff.start<-with(sw.diff.h,as.numeric(local.time-start)/3600)
## sw.diff.h$diff.end<-with(sw.diff.h,as.numeric(local.time-stop)/3600)
## 
## 
## end.obs<-sw.diff.h[,37:58] %>% filter(diff.end %in% c(-4,0,4)) %>%
## mutate(h = case_when(diff.end == -4 ~ "before",
##                      diff.end == 0 ~"end",
##                   diff.end == 4~"after",
##                    TRUE~as.character(NA)))
## 
## start.obs<-sw.diff.h[,37:58] %>% filter(diff.start %in% c(-4,0,4)) %>%
## mutate(h = case_when(diff.start == -4 ~ "before",
##                      diff.start == 0 ~"start",
##                   diff.start == 4~"after",
##                    TRUE~as.character(NA)))
## 
## 
## #model difference between before and start
## so.s<-start.obs[,1:15]
## so.s<- data.frame(lapply(so.s, function(x) scale(x, center = T)))
## so.s<-cbind(as.data.frame(start.obs[,c(16:23)]),so.s)
## 
## bs<-so.s %>% filter(complete.cases(.),h=="start"|h=="before") %>%  mutate(h1=case_when(h == "before"~0,
##                      h == "start"~1,TRUE~2))
## as<-so.s %>% filter(complete.cases(.),h=="start"|h=="after") %>%  mutate(h1=case_when(h == "after"~1,
##                      h == "start"~0,TRUE~2))
## #global model
## #glm(h1 ~ surface_roughness + temp + DEM +  precip_fraction+NDVI_16d+  dew_temp+ orographic +
##  #                 surface_pressure+  thermal+tot_atm_water,
##   #                 data=subset(as,pop=="meridionalis"),  na.action="na.omit")

## ----eval=FALSE,echo=FALSE-----------------------------------------------
## pa1 <- glm(h1 ~   temp + DEM +  precip_fraction+  dew_temp +
##                   +  thermal,
##                    data=subset(bs,pop=="aura"),  na.action="na.omit")
## drop1(pa1)[which.min(drop1(pa1)$AIC),]
## summary(pa1)
## 
## paa1 <- glm(h1 ~ temp  +  precip_fraction +  thermal,
##                    data=subset(as,pop=="aura"),  na.action="na.omit")
## drop1(paa1)[which.min(drop1(paa1)$AIC),]
## summary(paa1)

## ----eval=FALSE,echo=FALSE-----------------------------------------------
## pm1 <- glm(h1 ~  temp +  precip_fraction+  dew_temp +  thermal,
##                    data=subset(bs,pop=="meridionalis"),  na.action="na.omit")
## drop1(pm1)[which.min(drop1(pm1)$AIC),]
## summary(pm1)
## 
## pma1 <- glm(h1 ~         temp +  precip_fraction +  thermal+tot_atm_water,
##                    data=subset(as,pop=="meridionalis"),  na.action="na.omit")
## drop1(pma1)[which.min(drop1(pma1)$AIC),]
## summary(pma1)

## ----eval=FALSE,echo=FALSE-----------------------------------------------
## pr1 <- glm(h1 ~    temp +  dew_temp +
##                   surface_pressure+  thermal,
##                    data=subset(bs,pop=="ruficollis"),  na.action="na.omit")
## drop1(pr1)[which.min(drop1(pr1)$AIC),]
## 
## pra1 <- glm(h1 ~  thermal+tot_atm_water,
##                    data=subset(as,pop=="ruficollis"),  na.action="na.omit")
## drop1(pra1)[which.min(drop1(pra1)$AIC),]
## summary(pra1)

## ----output,eval=FALSE,include=FALSE-------------------------------------
## library(knitr)
## here()
## knit(here('Scripts/start and end of stopovers.Rmd'), output=here('Scripts/start_end_stopovers.R'), tangle=TRUE)

