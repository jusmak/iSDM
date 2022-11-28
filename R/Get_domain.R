Get_domain = function(wd, species) {
  
  #get whole study area
  r = raster('Data/Study_area/Study_area_laea.tif')
  
  #Survey locations, get coordinates of the survey sites
  #transform the coordinates of the records to plane coordinate system
  new_table = read.csv(paste0('Data/PO_observations_non_thinned/', species, ".csv"))
  colnames(new_table) = c('lat', 'lon')
  new_table = new_table[,c(2,1)]
  
  #remove duplicates
  new_table=new_table[!duplicated(coordinates(new_table)),]
  
  coordinates(new_table) = c("lon", "lat")
  proj4string(new_table) = CRS("+init=epsg:4326") # WGS 84
  CRS.new = crs(r)
  points_temp = spTransform(new_table, CRS.new)
  
  #Load the expert range map
  range_temp = readOGR(paste0('Data/Range_shp/', species, '.shp'))
  
  #change projection
  range_temp = sp::spTransform(range_temp, CRSobj = CRS.new)

  if (species == 'Chaetocercus_bombus' | species == 'Amazilia_tzacatl' | species == 'Florisuga_mellivora') {
    range_temp_rec = as(extent(range_temp), "SpatialPolygons")
    proj4string(range_temp_rec) = CRS.new
    domain_tmp = try(st_as_sf(points_temp) %>% st_union(st_as_sf(range_temp_rec)) %>% st_buffer(dist = 750000) %>% st_union())
  } else {
    domain_tmp = try(st_as_sf(points_temp) %>% st_union(st_as_sf(range_temp)) %>% st_buffer(dist = 750000) %>% st_union())
  }
  
  #transform domain to a raster, use study area raster as a spatial reference
  domain_tmp_sp = as(domain_tmp,  'Spatial')
  domain = rasterize(domain_tmp_sp, r)
  
  #mask with study area
  domain = raster::mask(domain, r)
  
  writeRaster(domain, filename = paste0('Data/Domains_new/', species, '.tif'), overwrite = T)
  
}
  