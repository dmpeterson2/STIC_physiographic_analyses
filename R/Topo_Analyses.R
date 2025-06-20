#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Topo_Analyses
#Coder: Delaney Peterson (dmpeterson2@crimson.ua.edu)
#Date Finalized: 6/15/2025
#Purpose: Bring in all necessary data, and perform initial analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup workspace --------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#clear workspace 
remove(list=ls())

#load packages
library(FedData)
library(sf)
library(stars)
library(tidyverse)
library(terra)
library(whitebox)
library(MultiscaleDTM)
library(mapview)

#turn off spherical geometry
sf_use_s2(FALSE)

#Create temp dir
temp_dir <- "C:\\R\\temp_dir\\"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gather data ------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 1 - bring in the site locations file 
sites <- read_csv("C:\\R\\data_dir\\data\\site_locations.csv") %>% 
  #drop any erroneous columns
  filter(!is.na(siteId))

#convert the csv to a shapefile using sf
sites_sf <- sites %>% st_as_sf(coords = c("long", "lat"), agr = "siteId",
                               crs = '+proj=longlat +datum=WGS84 +no_defs') %>% 
  #reproject
  st_transform(sites, crs = '+proj=utm +zone=16 +datum=NAD83 +units=m +no_defs')


## 2 - read in the watershed boundaries delineated (found in the ENVI files)
tal_shed <- st_read("C:\\R\\data_dir\\data\\watershed_TAL.shp") 
whr_shed <- st_read("C:\\R\\data_dir\\data\\watershed_WHR.shp") 
prf_shed <- st_read("C:\\R\\data_dir\\data\\watershed_PRF.shp") 


## 3 - read in the DEMs from USGS 
tal_DEM <- raster("C:\\R\\data_dir\\data\\croppedDEM_TAL.tif")
prf_DEM <- raster("C:\\R\\data_dir\\data\\croppedDEM_PRF.tif")
whr_DEM <- raster("C:\\R\\data_dir\\data\\croppedDEM_WHR.tif")


## 4 - use whitebox to gather topographic variables 

## TALLADEGA -----

writeRaster(tal_DEM, paste0(temp_dir, "dem_tal.tif"), overwrite = T)

wbt_fast_almost_gaussian_filter(
  input = "dem_tal.tif",
  output = "dem_smooth_tal.tif",
  sigma = 1.8,
  wd = temp_dir)

wbt_fill_single_cell_pits(
  dem = "dem_smooth_tal.tif",
  output = "dem_tal_fill.tif",
  wd = temp_dir)

wbt_breach_depressions(
  dem = "dem_tal_fill.tif",
  output = "dem_tal_breach.tif",
  wd = temp_dir)

wbt_d8_pointer(
  dem = "dem_tal_breach.tif",
  output = "flowdir_tal.tif",
  wd = temp_dir)

wbt_d8_flow_accumulation(
  input = "dem_tal_breach.tif",
  #you have to use the DEM not the flow direction 
  output = "flowaccum_tal.tif",
  wd = temp_dir
)

#stream network analysis
 
wbt_extract_streams(
  flow_accum = "flowaccum_tal.tif",
  output = "streams_hires_tal.tif",
  threshold = 10000,
  wd = temp_dir
)

wbt_raster_streams_to_vector(
  d8_pntr = "flowdir_tal.tif",
  streams = "streams_hires_tal.tif",
  output = "streams_hires_tal_fix.shp", 
  wd = temp_dir
)

# now calculate the variables of interest

wbt_d8_flow_accumulation(
  input = "dem_tal_breach.tif",
  output = "tal_uaa.tif",
  out_type="catchment area", 
  wd = temp_dir
)

wbt_length_of_upstream_channels(
  d8_pntr = "flowdir_tal.tif",
  streams = "tal_streams_hires_fix.tif",
  output = "tal_upstream_length.tif",
  wd = temp_dir)

wbt_aspect(
  dem = "dem_tal_breach.tif",
  output = "tal_aspect.tif",
  wd = temp_dir
)

wbt_total_curvature(
  dem = "dem_tal_fill.tif",
  output = "tal_curvature.tif",
  wd = temp_dir
)

tal_tpi <- TPI(tal_DEM, w = 5, shape = "circle")

# now bring the rasters into the R environment

tal_aspect <- raster(paste0(temp_dir, "tal_aspect.tif"))
tal_curvature <- raster(paste0(temp_dir, "tal_curvature.tif"))
tal_upstream <- raster(paste0(temp_dir, "tal_upstream_length.tif"))

#now extract the values for each sensor
tal_sensors <- sites %>% 
  filter(watershed == "TAL") %>% 
  dplyr::select(siteId, long, lat)

tal_tpi_st <- st_as_stars(tal_tpi)
#extract the data!
tal_sensor_tpi <- st_extract(tal_tpi_st, subset(sites_sf, watershed == "TAL")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  dplyr::select(siteId, tpi) %>% 
  unique()

tal_aspect_st <- st_as_stars(tal_aspect)
tal_sensor_aspect <- st_extract(tal_aspect_st, subset(sites_sf, watershed == "TAL")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  mutate(aspect = tal_aspect) %>% 
  dplyr::select(siteId, aspect) %>% 
  unique()

tal_curvature_st <- st_as_stars(tal_curvature)
tal_sensor_curv <- st_extract(tal_curvature_st, subset(sites_sf, watershed == "TAL")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  mutate(total_curvature = tal_curvature) %>% 
  dplyr::select(siteId, total_curvature) %>% 
  unique()

tal_upstream_st <- st_as_stars(tal_upstream)
tal_sensor_upstream <- st_extract(tal_upstream_st, subset(sites_sf, watershed == "TAL")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  mutate(Upstream_Length_m = tal_upstream_length) %>% 
  dplyr::select(siteId, Upstream_Length_m) %>% 
  unique()

#merge them all into the topo file you will use later (paired with the tal_envi from HydroShare)
tal_topo <- left_join(tal_envi, tal_sensor_upstream) %>% 
  dplyr::select(siteId, drainage_area_m, Upstream_Length_m) %>% 
  left_join(., tal_sensor_aspect) %>% 
  left_join(., tal_sensor_curv) %>% 
  left_join(., tal_sensor_tpi) %>% 
  mutate(drainage_density = Upstream_Length_m/drainage_area_m,
         drainage_density_km = (Upstream_Length_m/1000)/(drainage_area_m/1000000))


wbt_slope(
  dem = "dem_tal_fill.tif",
  output = "slope_tal.tif",
  units = "percent",
  wd = temp_dir
)

# do some quick and dirty analyses of the patterns within our watershed

tal_slope <- raster(paste0(temp_dir, "tal_slope.tif"))
tal_slope <- crop(tal_slope, tal_shed)
tal_slope <- mask(tal_slope, tal_shed)
tal_meanslope <- cellStats(tal_slope, stat = "mean")
tal_minslope <- cellStats(tal_slope, stat = "min")
tal_maxslope <- cellStats(tal_slope, stat = "max")



# now repeat for the other two watersheds

## PAINT ROCK -----

writeRaster(prf_DEM, paste0(temp_dir, "dem_prf.tif"), overwrite = T)

wbt_fast_almost_gaussian_filter(
  input = "dem_prf.tif",
  output = "dem_smooth_prf.tif",
  sigma = 1.8,
  wd = temp_dir)

wbt_fill_single_cell_pits(
  dem = "dem_smooth_prf.tif",
  output = "dem_prf_fill.tif",
  wd = temp_dir)

wbt_breach_depressions(
  dem = "dem_prf_fill.tif",
  output = "dem_prf_breach.tif",
  wd = temp_dir)

wbt_d8_pointer(
  dem = "dem_prf_breach.tif",
  output = "flowdir_prf.tif",
  wd = temp_dir)

wbt_d8_flow_accumulation(
  input = "dem_prf_breach.tif",
  #you have to use the DEM not the flow direction for some reason
  output = "flowaccum_prf.tif",
  wd = temp_dir
)

wbt_extract_streams(
  flow_accum = "flowaccum_prf.tif",
  output = "streams_hires_prf.tif",
  threshold = 60000,
  wd = temp_dir
)

wbt_raster_streams_to_vector(
  d8_pntr = "flowdir_prf.tif",
  streams = "streams_hires_prf.tif",
  output = "streams_hires_prf.shp", 
  wd = temp_dir
)

# variables of interest

wbt_length_of_upstream_channels(
  d8_pntr = "flowdir_prf.tif",
  streams = "streams_hires_prf.tif",
  output = "prf_upstream_length.tif",
  wd = temp_dir)

wbt_aspect(
  dem = "dem_prf_breach.tif",
  output = "prf_aspect.tif",
  wd = temp_dir
)

wbt_total_curvature(
  dem = "dem_prf_fill.tif",
  output = "prf_curvature.tif",
  wd = temp_dir
)

prf_tpi <- TPI(prf_DEM, w = 5, shape = "circle")

prf_aspect <- raster(paste0(temp_dir, "prf_aspect.tif"))
prf_curvature <- raster(paste0(temp_dir, "prf_curvature.tif"))
prf_upstream <- raster(paste0(temp_dir, "prf_upstream_length.tif"))


#now extract the values
prf_sensors <- sites %>% 
  filter(watershed == "PRF") %>% 
  dplyr::select(siteId, long, lat)

prf_tpi_st <- st_as_stars(prf_tpi)
#extract the data!
prf_sensor_tpi <- st_extract(prf_tpi_st, subset(sites_sf, watershed == "PRF")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  dplyr::select(siteId, tpi) %>% 
  unique()

prf_aspect_st <- st_as_stars(prf_aspect)
prf_sensor_aspect <- st_extract(prf_aspect_st, subset(sites_sf, watershed == "PRF")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  mutate(aspect = prf_aspect) %>% 
  dplyr::select(siteId, aspect) %>% 
  unique()

prf_curvature_st <- st_as_stars(prf_curvature)
prf_sensor_curv <- st_extract(prf_curvature_st, subset(sites_sf, watershed == "PRF")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  mutate(total_curvature = prf_curvature) %>% 
  dplyr::select(siteId, total_curvature) %>% 
  unique()

prf_upstream_st <- st_as_stars(prf_upstream)
prf_sensor_upstream <- st_extract(prf_upstream_st, subset(sites_sf, watershed == "PRF")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  mutate(Upstream_Length_m = prf_upstream_length) %>% 
  dplyr::select(siteId, Upstream_Length_m) %>% 
  unique()

prf_topo <- left_join(prf_envi, prf_sensor_upstream) %>% 
  dplyr::select(siteId, drainage_area_m, Upstream_Length_m) %>% 
  left_join(., prf_sensor_aspect) %>% 
  left_join(., prf_sensor_curv) %>% 
  left_join(., prf_sensor_tpi) %>% 
  mutate(drainage_density = Upstream_Length_m/drainage_area_m,
         drainage_density_km = (Upstream_Length_m/1000)/(drainage_area_m/1000000))


wbt_slope(
  dem = "dem_prf_fill.tif",
  output = "slope_prf.tif",
  units = "percent",
  wd = temp_dir)

prf_slope <- raster(paste0(temp_dir, "prf_slope.tif"))
prf_slope <- crop(prf_slope, prf_shed)
prf_slope <- mask(prf_slope, prf_shed)
prf_meanslope <- cellStats(prf_slope, stat = "mean")
prf_minslope <- cellStats(prf_slope, stat = "min")
prf_maxslope <- cellStats(prf_slope, stat = "max")



## SHAMBLEY CREEK -----

writeRaster(whr_DEM, paste0(temp_dir, "dem_whr.tif"), overwrite = T)

wbt_fast_almost_gaussian_filter(
  input = "dem_whr.tif",
  output = "dem_smooth_whr.tif",
  sigma = 1.8,
  wd = temp_dir)

wbt_fill_single_cell_pits(
  dem = "dem_smooth_whr.tif",
  output = "dem_whr_fill.tif",
  wd = temp_dir)

wbt_breach_depressions(
  dem = "dem_whr_fill.tif",
  output = "dem_whr_breach.tif",
  wd = temp_dir)

wbt_d8_pointer(
  dem = "dem_whr_breach.tif",
  output = "flowdir_whr.tif",
  wd = temp_dir)

wbt_d8_flow_accumulation(
  input = "dem_whr_breach.tif",
  #you have to use the DEM not the flow direction for some reason
  output = "flowaccum_whr.tif",
  wd = temp_dir
)

wbt_extract_streams(
  flow_accum = "flowaccum_whr.tif",
  output = "streams_hires_whr.tif",
  threshold = 12000,
  wd = temp_dir
)

wbt_raster_streams_to_vector(
  d8_pntr = "flowdir_whr.tif",
  streams = "streams_hires_whr.tif",
  output = "streams_hires_whr.shp", 
  wd = temp_dir
)

# variables of interest

wbt_length_of_upstream_channels(
  d8_pntr = "flowdir_whr.tif",
  streams = "streams_hires_whr.tif",
  output = "whr_upstream_length.tif",
  wd = temp_dir)

wbt_aspect(
  dem = "dem_whr_breach.tif",
  output = "whr_aspect.tif",
  wd = temp_dir
)

wbt_total_curvature(
  dem = "dem_whr_fill.tif",
  output = "whr_curvature.tif",
  wd = temp_dir
)

whr_tpi <- TPI(whr_DEM, w = 5, shape = "circle")

whr_aspect <- raster(paste0(temp_dir, "whr_aspect.tif"))
whr_curvature <- raster(paste0(temp_dir, "whr_curvature.tif"))
whr_upstream <- raster(paste0(temp_dir, "whr_upstream_length.tif"))


#now extract the values
whr_sensors <- sites %>% 
  filter(watershed == "WHR") %>% 
  dplyr::select(siteId, long, lat)

whr_tpi_st <- st_as_stars(whr_tpi)
#extract the data!
whr_sensor_tpi <- st_extract(whr_tpi_st, subset(sites_sf, watershed == "WHR")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  dplyr::select(siteId, tpi) %>% 
  unique()

whr_aspect_st <- st_as_stars(whr_aspect)
whr_sensor_aspect <- st_extract(whr_aspect_st, subset(sites_sf, watershed == "WHR")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  mutate(aspect = whr_aspect) %>% 
  dplyr::select(siteId, aspect) %>% 
  unique()

whr_curvature_st <- st_as_stars(whr_curvature)
whr_sensor_curv <- st_extract(whr_curvature_st, subset(sites_sf, watershed == "WHR")) %>% 
  st_join(sites_sf) %>% 
  as_tibble() %>% 
  mutate(total_curvature = whr_curvature) %>% 
  dplyr::select(siteId, total_curvature) %>% 
  unique()

whr_upstream_st <- st_as_stars(whr_upstream)
#make a buffer for these 
whr_buf <- st_buffer(subset(sites_sf, watershed == "WHR"), dist = 4)

buffer_fun <- function(n) {
  buf <- whr_buf[n,]
  
  value <- terra::extract(whr_upstream, buf, fun = mean, na.rm = T, df = T)
  
  value$ID <- buf$siteId
  
  value 
}

whr_sensor_upstream <- buffer_fun(n = 1:nrow(whr_buf)) %>% 
  mutate(Upstream_Length_m = whr_upstream_length,
         siteId = ID) %>% 
  dplyr::select(siteId, Upstream_Length_m) %>% 
  unique()


whr_topo <- left_join(whr_envi, whr_sensor_upstream) %>% 
  dplyr::select(siteId, drainage_area_m, Upstream_Length_m) %>% 
  left_join(., whr_sensor_aspect) %>% 
  left_join(., whr_sensor_curv) %>% 
  left_join(., whr_sensor_tpi) %>% 
  mutate(drainage_density = Upstream_Length_m/drainage_area_m,
         drainage_density_km = (Upstream_Length_m/1000)/(drainage_area_m/1000000))

wbt_slope(
  dem = "dem_whr_fill.tif",
  output = "slope_whr.tif",
  units = "percent",
  wd = temp_dir
)
whr_slope <- raster(paste0(temp_dir, "whr_slope.tif"))
whr_slope <- crop(whr_slope, whr_shed)
whr_slope <- mask(whr_slope, whr_shed)
whr_meanslope <- cellStats(whr_slope, stat = "mean")
whr_minslope <- cellStats(whr_slope, stat = "min")
whr_maxslope <- cellStats(whr_slope, stat = "max")

