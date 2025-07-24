library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)
library(sf)
library(ggmap)
library(maps)
library(patchwork)

#https://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
setwd("~/Desktop/urbanurchins")


#Mussel watch sites####
M<-read.csv("musselwatch/musselwatch.csv", header=TRUE, sep=",")

##Historical Exposure Sites ####
#subset only Twin Points, Treasure Island, Cabrillo and White Pt
HS<-read.csv("ch2_sites.csv", header=TRUE, sep=",")



# set wd
# packages
library(sf)
library(dplyr)
library(ggrepel)
library(ggspatial)
library(leaflet)
library(ggplot2)
library(maps) 
library(tools)
library(rnaturalearth)
library(forcats)
library(lubridate)
library(tidyr)
##################

# read in our site data
#######################
###ALL SITES#####
sites<-read.csv("metadata/urbanurchins_metadata_thinned.csv", header=TRUE, sep=",")

sites19<-sites[!duplicated(sites$Site),]

LAonly<-sites19[ which(sites19$region=='LA'),]
SDonly<-sites19[ which(sites19$region=='SD'),]
Viconly<-sites19[ which(sites19$region=='Vic'),]
#######################

# create maps
#############
# prelim visualization of sites of interest in google maps
m <- leaflet()
m <- addTiles(m)
m <- addCircleMarkers(m, long=long, lat=lat, radius =2, opacity = 1, 
                      label = sites$name, labelOptions = labelOptions(noHide = T))
m 

# set object
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)
world_crop <- st_crop(world, c(xmin =-117.2 , xmax = 110, ymin = -60, ymax = 60))
class(world)
class(world_crop)

# add state borders for the sake of clarity
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
head(states)
states <- cbind(states, st_coordinates(st_centroid(states)))
states$ID <- toTitleCase(states$ID)
head(states)
##just map##
map<-
  ggplot(data = world) +
  theme_bw() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = states, fill = 'antiquewhite1') +
  coord_sf(xlim = c(-127, -112), ylim = c(32, 49), expand = FALSE) + ##FIT TO REGION OF INTEREST
  theme(plot.title = element_text(size = 24), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "aliceblue"), legend.position = 'right')
ggsave("map_only.png", map, width=30, height=30, units = "cm")


#three paired cities ####
three<-read.csv("3sites.csv", header=TRUE, sep=",")
three<-three %>% 
    unite(dev_region, c(region, dev), sep="_", remove=FALSE)
#sites3<-
ggplot(data = world) +
    theme_bw() + 
    geom_sf(data = world_crop, fill = 'antiquewhite1') +
    geom_sf(data = states, fill = 'antiquewhite1') +
    geom_point(data = three, aes(x=long, y=lat, color=dev, size=5)) +
    #geom_text_repel(data = three, aes(x = long, y = lat, label = dev_region),
                    #size = 10, nudge_x = c(-1.25,1.25),fontface = "bold")+
    coord_sf(xlim = c(-127, -112), ylim = c(32, 49), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+ 
    theme(plot.title = element_text(size = 24), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
          panel.background = element_rect(fill = "aliceblue"), legend.position = 'right')

ggsave("map_3_notext.png", sites3, width=30, height=30, units = "cm")

##everything
all<-ggplot(data = world) +
  theme_bw() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = states, fill = 'antiquewhite1') +
  geom_point(data = sites19, aes(x=Longitude, y=Latitude, color=Dev, size=2)) +
  geom_text_repel(data = sites19, aes(x = Longitude, y = Latitude, label = Site_ID),
                  size = 10, nudge_x = c(-1.25,1.25),fontface = "bold")+
  coord_sf(xlim = c(-127, -112), ylim = c(32, 49), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+  theme(plot.title = element_text(size = 24), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "aliceblue"), legend.position = 'right')

ggsave("map_all.png", all, width=30, height=30, units = "cm")

##zoom on Vic
V<-
ggplot(data = world) +
  theme_bw() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = world, fill = 'antiquewhite1') +
  geom_point(data = Viconly, aes(x=Longitude, y=Latitude, color=Dev, size=5)) +
  coord_sf(xlim = c(-125, -122), ylim = c(47.5, 49), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+   theme(plot.title = element_text(size = 24), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
                                                                                       panel.background = element_rect(fill = "aliceblue"), ,legend.position="none")

##zoom on LA
L<-
  ggplot(data = world) +
  theme_bw() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = states, fill = 'antiquewhite1') +
  geom_point(data = LAonly, aes(x=Longitude, y=Latitude, color=Dev, size=5)) +
  coord_sf(xlim = c(-118.5, -117.5), ylim = c(33.5, 34), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+  theme(plot.title = element_text(size = 24), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"), ,legend.position="none")

## zoom on SD                                                                                                                                                                           
S<-
ggplot(data = world) +
  theme_bw() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = states, fill = 'antiquewhite1') +
  geom_point(data = SDonly, aes(x=Longitude, y=Latitude, color=Dev, size=5)) +
  coord_sf(xlim = c(-118, -117), ylim = c(32.6, 33.1), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+   theme(plot.title = element_text(size = 24), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
                                                                                      panel.background = element_rect(fill = "aliceblue"),legend.position="none")
mini_maps<-(V/L/S)

ggsave("maps_mini.png", mini_maps, width=20, height=40, units = "cm")

