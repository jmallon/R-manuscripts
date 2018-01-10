## ------------------------------------------------------------------------
library(here)
here.<-here()
#library(rgdal)
library(ggplot2)
#library(stringr)
#library(PBSmapping)
#library(maps)
#library(data.table)
#library(maptools)
#library(grid)
#library(sp)

## ----echo=FALSE,warning=FALSE,message=FALSE------------------------------
winos <- ifelse(grepl("windows", Sys.info()['sysname'], ignore.case=T), 1, 0)

if(winos==1) load("C:/Users/Julie/Box Sync/R/Chapter 1/Data/stopovers_with_all_weather2.robj")
if(winos==0) load("~/Box Sync/R/Chapter 1/Data/stopovers_with_all_weather2.robj")
      
      

#Plot both types
library(ggmap)
library(rworldmap)
library(rworldxtra)
library(scales)

latlimits <- c(-50, 55) 
longlimits <- c(-120, -40) 

newmap <- getMap(resolution = "high")
map.americas <- subset(newmap,REGION=='North America'|REGION=='South America and the Caribbean')
americas.points <- fortify(map.americas)
americas.points$region <- americas.points$id

americas.df <- americas.points[,c("long","lat","group", "region")]

## ----eval=FALSE----------------------------------------------------------
## 
## library(grDevices)
## telem_map<-ggplot() +
##   geom_polygon(data = americas.df, aes(x = long, y = lat,group=group),fill="grey90") +
##   scale_y_continuous(breaks = (-6:6) * 10,
##                      labels = c("60{\degree}S", "50{\degree}S", "expression(40~ degree~S)","30ºS",
##                                 "20ºS", "10ºS", "0ºN", "10ºN",
##                                 "20ºN","30ºN", "40ºN", "50ºN", "60ºN")) +
##   scale_x_continuous(breaks = (-12:-4) * 10,
##                      labels = c("120ºW","110ºW", "100ºW", "90ºW",
##                                 "80ºW", "70ºW", "60ºW","50ºW", "40ºW"))+
##   labs(y="",x="")+
##   coord_map("azequalarea", orientation = c(0, -90, 0))+
## geom_point(data=stops_full_weather2, aes(x=location.long, y=location.lat,color=pop), alpha = .5,size=1.4,pch=16)+  #scale_color_gradient2(high = "#d7191c",mid="#fdae61",low = "#2c7bb6",midpoint = 15)+
##  # scale_colour_gradientn(colours = c( c2),guide_legend(title="Fat (g) used \n per day"))+
## 
##   theme(panel.background = element_rect(fill = "white", colour = "grey"),
##         panel.grid.major = element_line(colour = "grey90"),
##         panel.grid.minor = element_blank(),
##         axis.ticks = element_blank(),
##         axis.text.x = element_text (size = 10, vjust = 0),
##         axis.text.y = element_text (size = 10, hjust = 1.3),
##         legend.position=c(.2, .3)) + coord_cartesian(xlim = longlimits, ylim = latlimits) #+ggtitle("Daily Energy Use")

## ----eval=FALSE,include=FALSE--------------------------------------------
##  library(knitr)
##  knit(here('Scripts/acopian telem map.Rmd'), output=here('Scripts/acopian telem map.R'), tangle=TRUE)

