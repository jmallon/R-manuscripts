#sourced example

#in-line references
n=rnorm(1)

#table 2
library(ggplot2)
m1<- glm(price ~ carat + depth + clarity, 
         data = diamonds)
m2<- glm(price ~ carat + depth + clarity, 
         data = subset(diamonds, cut =="Ideal"))
m3<- glm(price ~ carat + depth + clarity, 
         data = subset(diamonds, cut =="Very Good"))
m4<- glm(price ~ carat + depth + clarity, 
         data = subset(diamonds, cut =="Fair"))


#Figure 1 
library(ggmap)
library(rworldmap)
library(rworldxtra)
library(scales)

latlimits <- c(-50, 10) 
longlimits <- c(-80, -40) 

newmap <- getMap(resolution = "high")
map.americas <- subset(newmap,REGION=='South America and the Caribbean')
americas.points <- fortify(map.americas)
americas.points$region <- americas.points$id

americas.df <- americas.points[,c("long","lat","group", "region")]

set.seed(12)
points.df<-data.frame(lat = c(sample(-20:10,5),sample(-30:-40,3)), #FIX
                      lon = c(sample(-70: -40, 5),sample(-60:-70,3)),
                      group = c("A", "B"))


#Figure 2
p <- ggplot(mpg, aes(class, hwy))
p.mpg<- p + geom_boxplot()

