####################################################
####### North Carolina Map - Stability Trial ####### 
####################################################

####### Author: Vinicius Garnica
####### Date: Mar, 2024

### Load packages ------------------------------------------------------------------------------------ 
pacman::p_load(ggplot2,   # must have pacman installed!
               rgdal,
               sf,
               ggspatial,
               ggsn,
               cowplot,
               ggrepel,
               elevatr,
               soilDB,
               tidyverse)


### Set the path for workflow
rm(list = ls())
theme_set(theme_bw())

### Load data ------------------------------------------------------------------------------------ 

nc = st_read(dsn = "C:/Users/Garnica/Documents/GitHub/SNB_stability/data/map/nc_eco_l3.shp") # download file
nc = st_transform(nc, crs = 4326) # transform to EPSG 4326 defines a full coordinate reference system

### Download remaining data ------------------------------------------------------------------------------------ 
nc_eco = st_as_sf(nc) # convert to sf object
nc_state = map_data("state", "north carolina") # get state boundaries
nc_county = tigris::counties(state = "North Carolina", cb = TRUE) %>% # get county boundaries
  st_as_sf() %>%
  st_transform(st_crs(nc))

### Make color selections ------------------------------------------------------------------------------------ 
order = c("Blue Ridge", "Piedmont","Southeastern Plains","Middle Atlantic Coastal Plain") # Specify the order of the regions
region_colors = c("#33a02c", "#F0C35E","#E9D2B4","#98BCD9") # select colors for ecoregions

### Add data points for experiments ------------------------------------------------------------------------------------ 
sites=data.frame(Sites=c("PY22","MR22","SB22","RW22", "CL22", "KS22",
                       "PY23","UN23","SB23","LB23", "OX23", "KS23",
                       "BE24","AL24","SB24","LB24", "RO24", "KS24"),
               Year = c(2022,2022,2022, 2022,2022,2022,
                        2023,2023,2023,2023,2023,2023,
                        2024,2024,2024,2024,2024,2024),
               Lat=c(35.849167,34.949337,35.696489,34.511219,35.671203,35.378138,
                     35.849213,35.089407,35.699578,34.713276,36.35593184,35.375757,
                     35.559115,34.859825,35.697913,34.754557,36.322885,35.377746),
               Lon=c(-76.655278,-80.424468,-80.623680,-79.3061264,-78.511441,-77.559817,
                     -76.659167,-80.468470,-80.620807,-78.9671584,-78.54890126219952,-77.560750,
                     -76.568492,-80.512979,-80.624366,-79.0343698, -78.8776242,-77.557221))

sites= sites %>% 
  mutate(Region= case_when(Sites=="AL24" ~'○ Piedmont',
                           str_detect(Sites, 'KS') ~'■ Southeastern Plains',
                           Sites=="BE24" ~'× Middle Atlantic Coastal Plain',
                           Sites=="RO24" ~'○ Piedmont',
                           Sites=="OX23" ~'○ Piedmont',
                           Sites=="CL22" ~'○ Piedmont',
                           str_detect(Sites, 'SB') ~ '○ Piedmont',
                           Sites=="RW22" ~'■ Southeastern Plains',
                           str_detect(Sites, 'UN') ~ '○ Piedmont',
                           str_detect(Sites, 'MR') ~ '○ Piedmont',
                           str_detect(Sites, 'LB') ~ '■ Southeastern Plains',
                           str_detect(Sites, 'PY') ~ '× Middle Atlantic Coastal Plain'))


# Make a map for legend ------------------------------------------------------------------------------------ 
p1=ggplot() + 
  geom_sf(data=nc_eco,aes(fill = US_L3NAME),linewidth=0.1,color = "black",alpha=0.4)+
  geom_sf(data=nc_county,linewidth=0.05,color = "grey50",alpha=0.4)+
  scale_fill_manual(values = setNames(region_colors, order))+
  scale_shape_manual(values=c(1, 15))+
  labs(fill = "Ecoregion") +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_blank(),
        plot.title.position = "plot",
        legend.position = c(0.15, 0.22), 
        legend.justification = "center",
        legend.box.just = "center",
        legend.margin = margin(0,0,0,0),
        legend.text = element_text(size = 8))

# Extract the colour legend - leg1
legend1 = get_legend(p1)
legend1_grid = cowplot::plot_grid(legend1, align = "v", nrow = 1)

### Make a plot with no legend ------------------------------------------------------------------------------------ 
p3=ggplot() + 
  geom_sf(data=nc_eco,aes(fill = US_L3NAME),linewidth=0.1,color = "black",alpha=0.4)+
  geom_sf(data=nc_county,linewidth=0.05,color = "grey50",alpha=0.4)+
  geom_point(data =  sites, aes(x = Lon, y = Lat,shape=as.factor(Region)),size = 2,position = position_dodge2(w = 0.4)) +
  scale_fill_manual(values = setNames(region_colors, order))+
  geom_label_repel(data =  sites, aes(x = Lon, y = Lat,label = Sites),box.padding = 0.7,position = position_dodge2(width = 0.4),segment.color = 'grey50',)+
  scale_shape_manual(values=c(1, 4,15))+
  xlab("Longitude") + ylab("Latitude")  +
  annotation_scale(location = "br", width_hint = 0.3) + 
  north(nc,location = "topleft",scale = 0.15)+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_blank(),
        plot.title.position = "plot",
        legend.position =  "none")

### Final plot ------------------------------------------------------------------------------------ 
eco=p3 + annotation_custom(ggplotGrob(legend1_grid))
ggsave("results/plots/locations.tiff",  width = 10, height = 5, units = "in", dpi = 300)


### Get elevation data ------------------------------------------------------------------------------------ 
elevation = get_elev_point(locations = data.frame(x=sites$Lon,y=sites$Lat), prj = "EPSG:4326", src = "epqs")
data = data.frame(sites,elevation$elevation)
names(data)[4] = "longitude"
names(data)[3] = "latitude"
names(data)[1] = "site ID"


### Get soil data ------------------------------------------------------------------------------------ 
soil = fetchSoilGrids(data,loc.names = c("site ID", "latitude", "longitude"))
soil_dt = left_join(data,soil@horizons %>% dplyr::filter(label == "5-15")%>% dplyr::select(id,matches("mean")),by=c("site ID" = "id")) # obtain soil characteristics for 5-15 cm

write.csv(soil_dt,"soil_dt.csv")



