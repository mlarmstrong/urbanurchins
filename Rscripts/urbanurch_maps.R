library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)
library(sf)
library(ggmap)
library(maps)
library(patchwork)
library(ggrepel)
library(ggspatial)
library(leaflet)
library(tools)
library(rnaturalearth)
library(forcats)
library(lubridate)
library(tidyr)
#################

#https://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
setwd("~/Desktop/urbanurchins")

#figure theme to keep things consistent
theme_bigbox <- function(base_size = 11, base_family = '') {
  theme_classic() %+replace% 
    theme(text = element_text(size = 20), 
          panel.border = element_rect(color = "white", fill = NA, size = 1)) 
}


#Mussel watch sites####
M<-read.csv("preliminary_work/musselwatch/musselwatch.csv", header=TRUE, sep=",")

#just sites relevant to my dataset
MWS_data<-read.csv("MWS_waste_storm.csv", header=TRUE, sep=",")

#MW station rank
station_rankMWS<-read.csv("musselwatch_stationrank.csv",header=TRUE, sep=",")

##all treatment plant data
treatment_plants<-read.csv("treatmentplants_latlong.csv", header=TRUE, sep=",")

vic_plant<-treatment_plants %>% 
  filter(Region=="Vic")

##only MWS near my sites
scb_mws<-read.csv("Musselwatch_scbstations.csv", header=TRUE, sep=",")
#make variables continuous for plotting
scb_mws$AP <- as.numeric(scb_mws$AP)

##Historical Exposure Sites ####
#subset only Twin Points, Treasure Island, Cabrillo and White Pt
HS<-read.csv("ch2_sites.csv", header=TRUE, sep=",")
#

# read in our site data
#######################
###ALL SITES#####
sites<-read.csv("metadata/urbanurchins_metadata_thinned.csv", header=TRUE, sep=",")

sites19<-sites[!duplicated(sites$Site),]

LAonly<-sites19[ which(sites19$region=='LA'),]
SDonly<-sites19[ which(sites19$region=='SD'),]
Viconly<-sites19[ which(sites19$region=='Vic'),]

SCB<-sites19[ which(sites19$region=='LA' | sites19$region=='SD'),]
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
#map<-
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
  theme_bigbox() + 
    geom_sf(data = world_crop, fill = 'antiquewhite1') +
    geom_sf(data = states, fill = 'antiquewhite1') +
    geom_point(data = three, aes(x=long, y=lat, color=dev, size=5)) +
    #geom_text_repel(data = three, aes(x = long, y = lat, label = dev_region),
                    #size = 10, nudge_x = c(-1.25,1.25),fontface = "bold")+
    coord_sf(xlim = c(-127, -112), ylim = c(32, 49), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+ 
  labs(x = "Longitude", y = "Latitude") +
    theme(plot.title = element_text(size = 12), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
          panel.background = element_rect(fill = "aliceblue"), legend.position = 'none')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("map_3_notext.png", sites3, width=30, height=30, units = "cm")

##zoom on Vic
#V<-
ggplot(data = world) +
  theme_bigbox() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = world, fill = 'antiquewhite1') +
  geom_point(data = Viconly, aes(x=Longitude, y=Latitude, color=Dev, size=5)) +
  coord_sf(xlim = c(-125, -122), ylim = c(47.5, 49), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+   theme(plot.title = element_text(size = 12), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
                                                                                       panel.background = element_rect(fill = "aliceblue"), ,legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
##zoom on LA
#L<-
  ggplot(data = world) +
  theme_bigbox() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = world, fill = 'antiquewhite1') +
  geom_point(data = LAonly, aes(x=Longitude, y=Latitude, color=Dev, size=5)) +
  coord_sf(xlim = c(-118.5, -117.5), ylim = c(33.5, 34), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+  theme(plot.title = element_text(size = 12), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"), ,legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## zoom on SD                                                                                                                                                                           
#S<-
ggplot(data = world) +
  theme_bigbox() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = world, fill = 'antiquewhite1') +
  geom_point(data = SDonly, aes(x=Longitude, y=Latitude, color=Dev, size=5)) +
  coord_sf(xlim = c(-118, -117), ylim = c(32.6, 33.1), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+   theme(plot.title = element_text(size = 12), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
                                                                                      panel.background = element_rect(fill = "aliceblue"),legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mini_maps<-(V/L/S)
ggsave("maps_mini.png", mini_maps, width=20, height=40, units = "cm")

##fig 1
library(cowplot)
# make a blank plot that can receive a label for structure plot later
blank_plot <- ggplot() + theme_void() + theme_bigbox() +coord_fixed(ratio = 2) 

fig1 <- plot_grid(
  sites3,
  mini_maps,
  blank_plot,
  pca +theme(legend.position="none"),
  ncol = 4,
  rel_widths = c(1, 1, 0.5, 1),
  align = "lr",
  axis = "t",            # top and bottom axes aligned
  labels = "AUTO",
  label_size = 30,
  label_fontface = "plain"
)
ggsave("mapsfig1.png", fig1, width=60, height=30, units = "cm")


##LA and SD data
ggplot(data = world) +
  theme_bigbox() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = world, fill = 'antiquewhite1') +
  geom_point(data = SCB, aes(x=Longitude, y=Latitude, color=Dev, size=5)) +
  coord_sf(xlim = c(-118.5, -117), ylim = c(32.6, 34), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+   theme(plot.title = element_text(size = 12), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
                                                                                     panel.background = element_rect(fill = "aliceblue"),legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))                                                                                     

##add my sites + X's for Mussel Watch Sites
suppmap<-ggplot(data = world) +
  theme_bigbox() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = world, fill = 'antiquewhite1') +
  geom_point(data = SCB, aes(x=Longitude, y=Latitude, size=3,color= Dev)) +
  coord_sf(xlim = c(-118.5, -117), ylim = c(32.6, 34), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" , "urban" = "#665d4b"))+   
  geom_point(shape=4, data = MWS_data, aes(x=MWS.Long, y=MWS.Lat, size=3))+
  theme(plot.title = element_text(size = 12), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "aliceblue"),legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))    

ggsave("suppmap_MWS.png", suppmap, width=20, height=20, units = "cm")

###old maps ####
##add musselwatch points & shapes for sewage vs stormwater vs both vs none?
ggplot(data = world) +
  theme_bigbox() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = world, fill = 'antiquewhite1') +
  geom_point(data = MWS_data, aes(x=MWS.Long, y=MWS.Lat, size=3,color=Urban, shape=MWS.exposure)) +
  coord_sf(xlim = c(-118.5, -117), ylim = c(32.6, 34), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_shape_manual(values=c("wastewater"=17, "stormwater only"=15))+ #triangle for wastewater
  scale_color_manual(values = c("nonurban"="#99d1ec" , "urban" = "#665d4b"))+   
  geom_point(shape=4, data = MWS_data, aes(x=wastewater.plant.Longitude, y=wastewater.plant.Latitude, size=3))+
  theme(plot.title = element_text(size = 12), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
                                                                                     panel.background = element_rect(fill = "aliceblue"),legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))    


#vic treatment plant map
#V<-
ggplot(data = world) +
  theme_bigbox() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = world, fill = 'antiquewhite1') +
  geom_point(data = Viconly, aes(x=Longitude, y=Latitude, color=Dev, size=3)) +
  coord_sf(xlim = c(-125, -122), ylim = c(47.5, 49), expand = FALSE) + ##FIT TO REGION OF INTEREST
  scale_color_manual(values = c("nonurban"="#99d1ec" ,"urban" = "#665d4b"))+   
  geom_point(shape=4, data = vic_plant, aes(x=wastewater.plant.Longitude, y=wastewater.plant.Latitude, size=3))+
  theme(plot.title = element_text(size = 12), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
                                                                                     panel.background = element_rect(fill = "aliceblue"), ,legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###ranking of mussel watch sites

ggplot(data = world) +
  theme_bigbox() + 
  geom_sf(data = world_crop, fill = 'antiquewhite1') +
  geom_sf(data = world, fill = 'antiquewhite1') +
  geom_point(shape=20, data = scb_mws, aes(x=long, y=lat, size=3, color= AP)) +
  coord_sf(xlim = c(-118.5, -117), ylim = c(32.6, 34), expand = FALSE) + ##FIT TO REGION OF INTEREST
  geom_point( data = SCB, aes(x=Longitude, y=Latitude, size=3, shape=Dev))+
  scale_shape_manual(values = c("nonurban"= 3, "urban" = 4))+   
  theme(plot.title = element_text(size = 12), panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5),
        panel.background = element_rect(fill = "aliceblue"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))    
