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

## ----categorizing stopovers, echo=FALSE, include=FALSE, warning=FALSE,message=FALSE,eval=FALSE----
## load(here(paste0(prjdta,"stopovers_with_all_weather.robj")))
## library(adehabitatLT)
## library(adehabitatHR)
## library(dplyr)
## 
## #group 1
## df_stopovers<-subset(stops_full_weather,stopover==1)
## df_stopovers5<-df_stopovers %>% group_by(stop.id) %>% mutate(n=n()) %>% filter(n>4)
## group1<-df_stopovers[which(!(df_stopovers$stop.id %in% df_stopovers5$stop.id)),'stop.id'] %>% unique() %>% as.numeric()
## #calculate area from minimum convex polygon of the stopover
## centroids <- with(df_stopovers5,
##                   SpatialPointsDataFrame(coords = cbind("x" = X, "y" = Y),
##                                     data = data.frame("stop.id"=stop.id,"ID_migration"=ID_migration)))
## areas<-mcp(centroids[,1],percent=100, unout="km2") %>% as.data.frame()
## hist(areas$area,breaks=20)
## nn<-0.05
## g1<-areas %>% filter(area<=nn)
## g1<-as.numeric(as.character(g1$id))
## group1<-c(group1,g1)
## 
## groups<-areas %>% filter(area>nn)
## groups<-as.numeric(as.character(groups$id))
## df_groups<-df_stopovers5[which(df_stopovers5$stop.id %in% groups),]
## 
## #calcualte distance from adehabitat
## dt<-diff(df_groups$local.time) %>% as.numeric()
## x<-which(dt==0)
## df_groups2<-df_groups[-x,]
## 
##   df.traj <- as.ltraj(data.frame(df_groups2$X,df_groups2$Y), as.POSIXct(df_groups2$local.time), id=df_groups2$stop.id)
##   df_g2_dists<-ld(df.traj)
## 
##   #manually annotated into categories
## not<-c(34,60, 76,102, 138, 116, 120, 184, 190,200,222,234, 244,252,260, 294, 444, 464, 536,612)
## overnight<-c(12,26, 74, 84, 92, 106, 164,166,178,192,202,206,222, 236,250,262,266, 270,
##              284, 298, 310, 318,332,334, 336,358,416, 426,  434, 436,444, 454,464, 482, 558,
##              564,574,600, 604,612, 656, 660,662, 680, 694, 710 )
## 
## linear1<-c(24,32, 44, 50, 70, 82,104,132,136,150, 154,160,170,210,214,220,224,226,234, 238,244,
##            302, 338,344, 364, 366,   412, 440, 452, 458, 476, 488, 490, 498, 500, 504, 508, 524, 526, 532,
##            552, 562,566, 592, 596, 606, 616, 620,640,648, 650,  664, 674,696,  706, 714, 716)
## linear2<-c(46, 48, 118, 410,478, 568,570,  614, 634, 644,658, 688, 690, 692,  708, 720 )
## tort<-c(6, 8,10,18, 20, 22,28,30, 40, 52, 54, 58, 64, 88, 90,94,96,100,130,134,140,
##         142,144,  152,156, 158, 176, 182, 196, 208,218,  248, 258, 272, 274, 276, 286,304,
##         312, 314, 322, 324, 328, 342,346,352, 354, 360,374, 380, 382, 386, 388, 394, 396, 398,
##         402, 404, 424, 432,  438, 442, 462, 474, 486, 492, 502, 506, 520, 522, 528,530, 548, 554,
##         556, 560, 572,578,  588, 588, 590, 608, 610, 636,  666,668,  676, 678, 682, 684, 686,
##         704, 718, 722)
## group1<-c(group1,534, 470,168)
## yy<-c(not,overnight,linear1,linear2,tort,group1)
## 
## 
## df_stopovers[which(!(df_stopovers$stop.id %in% yy)),'stop.id'] %>% unique()
## 
## stops_full_weather2<-stops_full_weather %>% mutate(stop.cat=case_when(
##                                            .$stop.id %in% overnight ~ "overnight",
##                                            .$stop.id %in% linear1 ~ "linear1",
##                                            .$stop.id %in% linear2 ~ "linear2",
##                                            .$stop.id %in% tort ~"tortuous",
##                                            .$stop.id %in% group1 ~"sedentary",
##                                            TRUE ~ "NA"))
## stops_full_weather2<-merge(stops_full_weather2,areas,by.x="stop.id",by.y="id",all.x=T,all.y=F)
## stops_full_weather2[which(stops_full_weather2$stop.id %in% not),'stopover']<-0
## save(stops_full_weather2,file=here(paste0(prjdta,"stopovers_with_all_weather2.robj")))

## ----summarize stopovers by category, echo=FALSE, message=FALSE, include=FALSE----
load(here(paste0(prjdta,"stopovers_with_all_weather2.robj")))
library(dplyr)
nn<-0.05
k1<-stops_full_weather2 %>% filter(stopover==1) %>% summarize(n=length(unique(stop.id)))
k2<-stops_full_weather2 %>% filter(stopover==1,stop.cat=="tortuous") %>% group_by(stop.id) %>% summarize(area=area[1])
stops_full_weather2$tort.size<-ifelse(stops_full_weather2$area>hist(k2$area,breaks=100)$breaks[2],1,0)
tort.lim<-hist(k2$area,breaks=100)$breaks[2]
stops_full_weather2<-stops_full_weather2 %>% filter(stopover==1) %>% 
  mutate(stop.cat2=case_when(.$stop.cat=="tortuous"&.$tort.size==1~"large tortuous",
                            .$stop.cat=="tortuous"&.$tort.size==0~"small tortuous",
                            TRUE~.$stop.cat))
k3<-stops_full_weather2%>% group_by(stop.cat2) %>% summarize(n=length(unique(stop.id)))
k4<-stops_full_weather2 %>% filter(stopover==1) %>% 
  summarize(n.bird=length(unique(individual.local.identifier)),n.migs=length(unique(ID_migration))) #for MS

## ----echo=FALSE,include=FALSE,message=FALSE------------------------------
library(lme4)
df_model<-subset(stops_full_weather2,roost==0)

#sedentary
df_model$response<-with(df_model,ifelse(stop.cat2=="sedentary",1,
                                          ifelse(migration==1,0,NA)))
df_model1<-df_model[which(df_model$ID_migration %in% unique(subset(df_model,response==1)$ID_migration)),]
n01<-df_model1 %>% filter(response==1) %>% summarize(n()) %>% as.numeric()
set.seed(11)
df_model1<-df_model1 %>% group_by(response) %>% sample_n(n01)
df_model2<-df_model1[,c(37:48,50:51,54:67)]
df_model2<- data.frame(lapply(df_model2, function(x) scale(x, center = T)))
df_model3<-cbind(as.data.frame(df_model1[,c(1:35,52, 68:72)]),df_model2)

m1<-  glmer(response~#surface_roughness+
              tot_atm_water+#accum_precip_m+
              thermal+    #temp+ 
              tailwind+    
              orographic+surface_pressure+downward_LW+#downward_SW+ boundary_height+  
              (1|ID_migration),data=df_model3,family="binomial")
#AIC(m1)
print(summary(m1),correlation=T)

#overnight
df_model$response<-with(df_model,ifelse(stop.cat2=="overnight",1,
                                          ifelse(migration==1,0,NA)))
df_model1<-df_model[which(df_model$ID_migration %in% unique(subset(df_model,response==1)$ID_migration)),]
n01<-df_model1 %>% filter(response==1) %>% summarize(n()) %>% as.numeric()
set.seed(11)
df_model1<-df_model1 %>% group_by(response) %>% sample_n(n01)
df_model2<-df_model1[,c(37:48,50:51,54:67)]
df_model2<- data.frame(lapply(df_model2, function(x) scale(x, center = T)))
df_model3<-cbind(as.data.frame(df_model1[,c(1:35,52, 68:72)]),df_model2)

m2<-  glmer(response~#surface_roughness+tot_atm_water+accum_precip_m+thermal+
              temp + #surface_pressure+#orographic+tailwind+
              downward_SW+ #downward_LW+
              #boundary_height+  
              (1|ID_migration),data=df_model3,family="binomial")
#AIC(m2)
print(summary(m2),correlation=T)

#lienar1
df_model$response<-with(df_model,ifelse(stop.cat2=="linear1",1,
                                          ifelse(migration==1,0,NA)))
df_model1<-df_model[which(df_model$ID_migration %in% unique(subset(df_model,response==1)$ID_migration)),]
n01<-df_model1 %>% filter(response==1) %>% summarize(n()) %>% as.numeric()
set.seed(11)
df_model1<-df_model1 %>% group_by(response) %>% sample_n(n01)
df_model2<-df_model1[,c(37:48,50:51,54:67)]
df_model2<- data.frame(lapply(df_model2, function(x) scale(x, center = T)))
df_model3<-cbind(as.data.frame(df_model1[,c(1:35,52, 68:72)]),df_model2)

m3<-  glmer(response~#surface_roughness+
              tot_atm_water+
              accum_precip_m+#thermal+
              #temp+  tailwind 
              +orographic+
              surface_pressure+
             # downward_SW+ 
              boundary_height+  #downward_LW+
              (1|ID_migration),data=df_model3,family="binomial")
#AIC(m3)
print(summary(m3),correlation=T)

#lienar2
df_model$response<-with(df_model,ifelse(stop.cat2=="linear2",1,
                                          ifelse(migration==1,0,NA)))
df_model1<-df_model[which(df_model$ID_migration %in% unique(subset(df_model,response==1)$ID_migration)),]
n01<-df_model1 %>% filter(response==1) %>% summarize(n()) %>% as.numeric()
set.seed(11)
df_model1<-df_model1 %>% group_by(response) %>% sample_n(n01)
df_model2<-df_model1[,c(37:48,50:51,54:67)]
df_model2<- data.frame(lapply(df_model2, function(x) scale(x, center = T)))
df_model3<-cbind(as.data.frame(df_model1[,c(1:35,52, 68:72)]),df_model2)

m4<-  glmer(response~surface_roughness+tot_atm_water+
              #accum_precip_m+#thermal+temp+ tailwind +#orographic+
              surface_pressure+
              downward_SW+ #boundary_height+  downward_LW+
              (1|ID_migration),data=df_model3,family="binomial")
#AIC(m4)
print(summary(m4),correlation=T)

#small tortuous
df_model$response<-with(df_model,ifelse(stop.cat2=="small tortuous",1,
                                          ifelse(migration==1,0,NA)))
df_model1<-df_model[which(df_model$ID_migration %in% unique(subset(df_model,response==1)$ID_migration)),]
n01<-df_model1 %>% filter(response==1) %>% summarize(n()) %>% as.numeric()
set.seed(11)
df_model1<-df_model1 %>% group_by(response) %>% sample_n(n01)
df_model2<-df_model1[,c(37:48,50:51,54:67)]
df_model2<- data.frame(lapply(df_model2, function(x) scale(x, center = T)))
df_model3<-cbind(as.data.frame(df_model1[,c(1:35,52, 68:72)]),df_model2)

m5<-  glmer(response~#surface_roughness+
              #tot_atm_water+
              accum_precip_m+#thermal+
              temp+ 
             # tailwind 
              +orographic+
              #surface_pressure+
              downward_SW+ 
              #boundary_height+ #downward_LW+
              (1|ID_migration),data=df_model3,family="binomial")
#AIC(m5)
print(summary(m5),correlation=T)

#large tortuous
df_model$response<-with(df_model,ifelse(stop.cat2=="large tortuous",1,
                                          ifelse(migration==1,0,NA)))
df_model1<-df_model[which(df_model$ID_migration %in% unique(subset(df_model,response==1)$ID_migration)),]
n01<-df_model1 %>%  filter(response==0) %>% summarize(n())  %>% as.numeric() 
#this is reversed because I have more stopover points than migration points
set.seed(11)
df_model1<-df_model1 %>% group_by(response) %>% sample_n(n01)
df_model2<-df_model1[,c(37:48,50:51,54:67)]
df_model2<- data.frame(lapply(df_model2, function(x) scale(x, center = T)))
df_model3<-cbind(as.data.frame(df_model1[,c(1:35,52, 68:72)]),df_model2)

m6<-  glmer(response~surface_roughness+
              tot_atm_water+
              accum_precip_m+thermal+
              #temp+  tailwind +
              orographic+#surface_pressure+
             downward_SW+ 
              boundary_height+ # downward_LW+
              (1|ID_migration),data=df_model3,family="binomial")
#AIC(m6)
print(summary(m6),correlation=T)


## ----echo=FALSE,eval=FALSE-----------------------------------------------
##  library(sjPlot)
## 
##     sjt.glmer(m1,m2,m3,m4,m5,m6,string.dv="Response", depvar.labels  =c("Sedentary", "Overnight","Single Linear", "Mulitple Linear", "Small Tortuous", "Large Tortuous"))

## ----echo=FALSE----------------------------------------------------------

df_model$response<-with(df_model,ifelse(stop.cat2=="sedentary",1,
                                        ifelse(stop.cat2=="large tortuous",0,NA)))
df_m_stp<-subset(df_model,stopover==1)
df_m_stp1<-df_m_stp[,c(37:48,50:51,54:67)]
df_m_stp1<- data.frame(lapply(df_m_stp1, function(x) scale(x, center = T)))
df_m_stp2<-cbind(as.data.frame(df_m_stp[,c(1:35,52, 68:72)]),df_m_stp1)

train.data<-subset(df_m_stp2,!is.na(response))
test.data<-subset(df_m_stp2,is.na(response))

n02<-train.data %>%  group_by(response) %>% summarize(n()) 
set.seed(11)
train.data1<-train.data %>% group_by(response) %>% sample_n(as.numeric(n02[2,2]))


m.train<-  glmer(response~surface_roughness+
              tot_atm_water+ accum_precip_m+thermal+
             # tailwind + #orographic+
                surface_pressure+
              #downward_SW+  boundary_height+ downward_LW+
              (1|ID_migration),data=train.data1,family="binomial")
#AIC(m.train)
print(summary(m.train),correlation=T)


  pred <- predict(object = m.train,
                newdata = train.data, allow.new.level=T)
  train.results<-data.frame(fit=pred,
                stop.cat=train.data$stop.cat2,
                stop.id=train.data$stop.id)
  train.results$fit01<-ifelse(train.results$fit<.5,0,1)

  pred <- predict(object = m.train,
                newdata = test.data, allow.new.level=T)
test.results<-data.frame(fit=pred,
                stop.cat=test.data$stop.cat2,
                stop.id=test.data$stop.id)
test.results$fit01<-ifelse(test.results$fit<.5,0,1)

 r1<-table(test.results$fit01,test.results$stop.cat)
 r2<-table(test.results$stop.cat)
library(scales)

## ----eval=FALSE,include=FALSE--------------------------------------------
## library(knitr)
## knit('Scripts/weather vs feeding.Rmd', output='Scripts/weather vs feeding.R', tangle=TRUE)

