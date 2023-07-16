# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Data extraction for the CPR analysis


# ----------------------------- LOAD PACKAGES ---------------------------

library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(rgbif)
library(sp)
library(raster)
library(ncdf4)
library(rnaturalearth)
library(ggspatial)
library(vegan)
library(sdmpredictors)


# -----------------------------------------------------------------------


# ---------------------- EXTRACT DIVERSITY DATA -------------------------
# extract diversity data via the rgbif package: internet connection required
# this extraction takes some time

# get a list of the included species
# only include Copepoda that are identified up to species level
spec_list <- read.csv('CPR_SpeciesList.csv',header=TRUE,row.names=NULL) %>% 
				dplyr::filter(Class=='Copepoda' & !is.na(Species))

# add the GBIF species key to the list
spec_list$speciesKey <- 0

for(i in 1:length(spec_list$speciesKey)){
	spec_list$speciesKey[i] <- rgbif::name_backbone(name=paste(spec_list$Genus[i],spec_list$Species[i])) %>% pull(usageKey) 
}

# initialise the count container
# each observation is reported with coordinates and time
dat_count_all <- matrix(nrow=0,ncol=4) %>% as.data.frame()
names(dat_count_all) <- c('lon','lat','key','year') 

for(i in 1:length(spec_list$speciesKey)){

	print(paste(i,'of',length(spec_list$speciesKey))) 

	dat_count_sub <- rgbif::occ_data(datasetKey='6d56415d-b007-4273-9c74-bcd6b2467434',speciesKey=spec_list$speciesKey[i],limit=100)$data 
	if(!is.null(dat_count_sub)){
		dat_count_sub <- dat_count_sub %>% dplyr::select(c(decimalLongitude,decimalLatitude,speciesKey,year))
		dat_count_all <- rbind(dat_count_all,dat_count_sub)
	} else{
		print('0')
	}

}

# problem with Metridia lucens: 
# species is in dataset, but is not recognised by rgbif
# add manually to the dataset
# downloaded from: https://doi.org/10.17031/1629
metr <- read.csv('CPR_Metridia.csv',header=TRUE,sep='\t') %>% dplyr::select(c(decimalLongitude,decimalLatitude,speciesKey,year))
dat_count_all <- rbind(dat_count_all,metr)

# fix names
names(dat_count_all) <- c('lon','lat','key','year')
# write out combined Copepoda data
write.csv(dat_count_all,'CPR_combined_div.csv')

# -----------------------------------------------------------------------


# ------------- CREATE SUBAREA NORTH ATLANTIC ---------------------------
# dimensions of the grid are arbitrary to get a uniform distribution of samples

# get count data
# only keep counts with geographic location
dat_all_counts <- read.csv('CPR_combined_div.csv',header=TRUE,row.names=1) %>% dplyr::filter(!is.na(lon) & !is.na(lat))

# define the Lambert Conformal Conic projection
crs_atl <- CRS('+proj=lcc +lon_0=-20 +lat_1=40 +lat_2=60')

# transform to spatial points
points <- dat_all_counts %>% sf::st_as_sf(.,crs=4326,coords=c('lon','lat')) %>% sf::st_transform(crs=crs_atl)

# define the arbitrary grid extent
ex <- sf::st_sf(geom=sf::st_sfc(sf::st_point(c(-2200000,5436156)),sf::st_point(c(-2200000,7500000)),
						sf::st_point(c(1850000,5436156)),sf::st_point(c(1850000,7500000)))) %>% 
		sf::st_bbox() %>% raster::extent()

# create the simulation grid raster
r <- raster::raster(ncol=45,nrow=14,crs=crs_atl,ext=ex)
# set dummy values for plotting
values(r) <- 0

# create land layer for plotting
shp_land <- rnaturalearth::ne_countries(scale='medium') %>% sf::st_as_sf() %>% sf::st_geometry() %>%
				sf::st_crop(., sf::st_sf(geom=sf::st_sfc(sf::st_point(c(-65,33)),sf::st_point(c(-65,70)),
									sf::st_point(c(20,33)),sf::st_point(c(20,70)))) %>% sf::st_bbox()) %>%
				sf::st_transform(crs=crs_atl)
# crop land layer to raster size
shp_atl_land <- shp_land %>% sf::st_crop(r) 

# create dataframe with centre coordinates of the grid cells
r_df_out <- as.data.frame(r, xy = TRUE, na.rm = TRUE) %>% dplyr::select(-c(layer))

# plot the simulation grid
p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=r_df_out,aes(x=x,y=y),fill=NA,color='black') + theme_bw() +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank())

# uncomment to get the png file
#p %>% ggsave('CPR_grid.png',.,device='png',width=15,height=10,units='cm')

# polygonise the raster grid
r_poly <- r %>% raster::rasterToPolygons() %>% as(.,'SpatialPolygons') %>% sf::st_as_sf(.)

# create containers for the surface area, circumference and the land cover
area_vec <- rep(0,nrow(r_df_out))
peri_vec <- rep(0,nrow(r_df_out))
land_vec <- rep(0,nrow(r_df_out))

# for each grid cell
for(i in 1:length(r_poly[[1]])){
	# get intersection with land
	a <- (sf::st_intersection(r_poly[[1]][i],shp_atl_land) %>% sf::st_area()) / (r_poly[[1]][i] %>% sf::st_area())
	if(length(a) == 0){
		# if no land, NA
		a <- NA
	}

	# surface area in square kilometres
	area_vec[i] <- (r_poly[[1]][i] %>% sf::st_area() )/1e6
	# transform to line in order to calculate circumference (in kilometres)
	peri_vec[i] <- (r_poly[[1]][i] %>% st_cast('MULTILINESTRING') %>% st_length()) / 1e3
	land_vec[i] <- a
}

r_df_out <- r_df_out %>% dplyr::mutate(area=area_vec,peri=peri_vec,land_perc=land_vec)

# plot land cover ratios
p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=r_df_out,aes(x=x,y=y,fill=land_perc),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='land percentage') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_land_perc.png',.,device='png',width=18,height=12,units='cm')

# write out base grid information
write.csv(r_df_out,'CPR_grid_base.csv')

# -----------------------------------------------------------------------


# ------------------- EXTRACT ENVIRONMENTAL VARIABLES -------------------
# note: this analysis takes some time due to the size of the raster layers
# ignore warnings
# problems of (RAM) memory could arise due to the size of rasters

# extract temperature layer from Bio-ORACLE; internet connection required
temp <- sdmpredictors::load_layers('BO2_tempmean_ss') %>% 
			raster::crop(.,sf::st_sf(geom=sf::st_sfc(sf::st_point(c(-65,33)),sf::st_point(c(-65,70)),
												sf::st_point(c(20,33)),sf::st_point(c(20,70)))) %>% sf::st_bbox())
# write out temperature raster
temp %>% raster::writeRaster('CPR_temp.tif',overwrite=TRUE)

# extract phytoplankton layer from Bio-ORACLE; internet connection required
phyto <- sdmpredictors::load_layers('BO2_carbonphytomean_ss') %>% 
			raster::crop(.,sf::st_sf(geom=sf::st_sfc(sf::st_point(c(-65,33)),sf::st_point(c(-65,70)),
												sf::st_point(c(20,33)),sf::st_point(c(20,70)))) %>% sf::st_bbox())
# write out phytoplankton raster
phyto %>% raster::writeRaster('CPR_phyto.tif',overwrite=TRUE)

# unzip environment layers
# downloadable from: https://drive.google.com/drive/folders/1HJ-QKm95RALnsqqGJZJo9WAxAD1Eoao4?usp=drive_link
unzip('CPR_env.zip')

# get the current velocity layers
# downloaded from: https://doi.org/10.48670/moi-00016
velo_x <- raster::stack('CPR_velo_x.nc') %>% mean() 
velo_y <- raster::stack('CPR_velo_y.nc') %>% mean()
# calculate kinetic energy
velo <- sqrt(velo_x*velo_x + velo_y*velo_y)

# transform environmental rasters to points
# definition crs_atl -> line 87
temp_p <- temp %>% raster::projectRaster(crs=crs_atl) %>% raster::rasterToPoints() %>% as.data.frame() %>% 
					dplyr::rename(temp=BO2_tempmean_ss) %>% sf::st_as_sf(coords=c('x','y'),crs=crs_atl)
phyto_p <- phyto %>% raster::projectRaster(crs=crs_atl) %>% raster::rasterToPoints() %>% as.data.frame() %>% 
					dplyr::rename(phyto=BO2_carbonphytomean_ss) %>% sf::st_as_sf(coords=c('x','y'),crs=crs_atl)
velo_p <- velo %>% raster::projectRaster(crs=crs_atl) %>% raster::rasterToPoints() %>% as.data.frame() %>% 
					dplyr::rename(velo=layer) %>% sf::st_as_sf(coords=c('x','y'),crs=crs_atl)

# initialise containers for the environmental grid values
# definition for r_poly -> line 124
temp_vector <- rep(NULL,length(r_poly[[1]]))
phyto_vector <- rep(NULL,length(r_poly[[1]]))
velo_vector <- rep(NULL,length(r_poly[[1]]))

# for each grid cell
for(i in 1:length(r_poly[[1]])){
	# intersection between environmental layer and simulation grid
	p_sub_temp <- sf::st_intersection(temp_p,r_poly[[1]][i])
	if(length(p_sub_temp[[1]])==0){
		t <- NA
	} else{
		# temperature of simulation cell = mean of all points/cells of environmental variable
		# ! APPROXIMATION
		t <- mean(p_sub_temp$temp)
	}
	temp_vector[i] <- t

	p_sub_phyto <- sf::st_intersection(phyto_p,r_poly[[1]][i])
	if(length(p_sub_phyto[[1]])==0){
		ph <- NA
	} else{
		# phytoplankton concentration of simulation cell = mean of all points/cells of environmental variable
		# ! APPROXIMATION
		ph <- mean(p_sub_phyto$phyto)
	}
	phyto_vector[i] <- ph

	p_sub_velo <- sf::st_intersection(velo_p,r_poly[[1]][i])
	if(length(p_sub_velo[[1]])==0){
		v <- NA
	} else{
		# current velocity of simulation cell = mean of all points/cells of environmental variable
		# ! APPROXIMATION
		v <- mean(p_sub_velo$velo)
	}
	velo_vector[i] <- v
}

# read in base grid information
r_df_out <- read.csv('CPR_grid_base.csv',header=TRUE,row.names=1) %>%
				dplyr::select(c(x,y,land_perc))

r_df_out <- r_df_out %>% dplyr::mutate(temp=temp_vector,phyto=phyto_vector,velo=velo_vector)

# write out environmental grid information
write.csv(r_df_out,'CPR_grid_env.csv')

# plot the environmental variables
# definition shp_atl_land -> line 108
p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=r_df_out,aes(x=x,y=y,fill=temp),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='temperature (°C)') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_temp.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=r_df_out,aes(x=x,y=y,fill=phyto),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name=expression(phytoplankton~concentration~(μmol~m^'-3'))) +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_phyto.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=r_df_out,aes(x=x,y=y,fill=velo),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='current velocity (m/s)') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_velo.png',.,device='png',width=18,height=12,units='cm')

# -----------------------------------------------------------------------


# ---------------------- MAP DIVERSITY TO GRID --------------------------

# make plot of all the samples in the grid
# collect combined count data
dat_all_counts <- read.csv('CPR_combined_div.csv',header=TRUE,row.names=1) %>% dplyr::filter(!is.na(lon) & !is.na(lat))
# get unique sampling locations
# definition crs_atl -> line 87
# definition r -> line 98 
dat_unique_counts <- dat_all_counts %>% dplyr::select(c(lon,lat)) %>% unique() %>% 
					sf::st_as_sf(.,crs=4326,coords=c('lon','lat')) %>% sf::st_transform(crs=crs_atl) %>% sf::st_crop(r) %>% sf::st_geometry()

# definition shp_atl_land -> line 108
p <- ggplot() + geom_sf(data=shp_atl_land) +
		geom_sf(data=dat_unique_counts,size=0.25) + theme_bw() + 
		theme(panel.background = element_rect(fill = "#d0dfff"),
				panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',colour = "grey")) +
		coord_sf(expand=FALSE) + 
		ggspatial::annotation_north_arrow(location = "tl",which_north = "true", pad_x = unit(0.5, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.1,'cm'), pad_x=unit(0.5,'cm'))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_samples.png',.,device='png',width=18,height=12,units='cm')

# make diversity grid
dat_samples_div <- dat_all_counts %>% dplyr::group_by(lon,lat,key) %>% dplyr::summarise(count=length(key)) %>%
					tidyr::pivot_wider(values_from=count,names_from=key)
dat_samples_div[is.na(dat_samples_div)] <- 0

# transform to spatial points and crop to grid
points <- dat_samples_div %>% sf::st_as_sf(.,crs=4326,coords=c('lon','lat')) %>% sf::st_transform(crs=crs_atl) %>% 
			sf::st_crop(r)
# remove species not present in grid
xy <- points %>% sf::st_coordinates() %>% as.data.frame() %>% dplyr::rename(x=X,y=Y)
dat_samples_div_sub <- points %>% sf::st_drop_geometry()
dat_samples_div_sub <- dat_samples_div_sub[,colSums(dat_samples_div_sub)!=0]
points_sub <- cbind(xy,dat_samples_div_sub) %>% sf::st_as_sf(coords=c('x','y'),crs=crs_atl)

# write out diversity matrix for samples
write.csv(cbind(xy,dat_samples_div_sub),'CPR_samples_div.csv')

# initialise output container for diversity and composition for simulation grid
# definition r_poly -> line 124
df_spec_out <- matrix(nrow=0,ncol=length(names(dat_samples_div_sub))) %>% as.data.frame()

# for each grid cell
for(i in 1:length(r_poly[[1]])){
	print(i)
	# get total number of observations for each species in grid cell
	inter <- sf::st_intersection(points_sub,r_poly[[1]][i]) %>% sf::st_drop_geometry() %>% 
				colSums() %>% as.data.frame() %>% t()
	rownames(inter) <- i
	# add to container
	df_spec_out <- rbind(df_spec_out,inter)
}

# combine with coordinates of grid cells and calculate normalised diversity
r_df_out <- read.csv('CPR_grid_base.csv',header=TRUE,row.names=1) %>% dplyr::select(c(x,y))
r_df_out <- cbind(r_df_out,df_spec_out)
r_df_out <- r_df_out %>% dplyr::mutate(diversity=r_df_out[,-c(1,2)] %>% vegan::decostand(method='pa') %>% rowSums())
r_df_out$diversity <- r_df_out$diversity / max(r_df_out$diversity,na.rm=TRUE)
# NA = no counts in the grid cell
r_df_out$diversity[r_df_out$diversity==0] <- NA
r_df_out <- r_df_out[,c('x','y','diversity',names(r_df_out)[3:(ncol(r_df_out)-1)])]

# write out grid diversity
write.csv(r_df_out,'CPR_grid_div.csv')

# plot the grid diversity
p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=r_df_out,aes(x=x,y=y,fill=diversity),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='diversity') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_div.png',.,device='png',width=18,height=12,units='cm')

# -----------------------------------------------------------------------


# ------------------ MAKE FINAL SIMULATION GRID -------------------------

# read environment grid and add diversity
dat_env_div <- read.csv('CPR_grid_env.csv',header=TRUE,row.names=1) %>% 
				dplyr::select(c(x,y,land_perc,temp,phyto,velo)) %>%
				dplyr::mutate(diversity=read.csv('CPR_grid_div.csv',header=TRUE,row.names=1) %>% dplyr::pull(diversity))

# to eliminate edge effects, only consider cells with less than 80% land cover
# if necessary, also discard cells with no count data
# all 'non-active' cells = NA
dat_env_div$land_perc[is.na(dat_env_div$land_perc)] <- 0
dat_env_div$land_perc[dat_env_div$land_perc > 0.8 | is.na(dat_env_div$diversity)] <- NA
dat_env_div$temp[is.na(dat_env_div$land_perc)] <- NA
dat_env_div$phyto[is.na(dat_env_div$land_perc)] <- NA
dat_env_div$velo[is.na(dat_env_div$land_perc)] <- NA
dat_env_div$diversity[is.na(dat_env_div$land_perc)] <- NA

# add whether the cell is active or not; land < 0.8 = active = 1; land > 0.8 = not active = 0 
dat_env_div_sim <- dat_env_div %>% dplyr::mutate(active=1)
dat_env_div_sim$active[is.na(dat_env_div_sim$land_perc)] <- 0

# define the geometry of the simulation
# each cell will remember whether it has a specific active neighbour in its neighbourhood
# this effectively creates the boundary conditions

# list of latitudes and longitudes
lat_list <- unique(dat_env_div_sim$y)
lon_list <- unique(dat_env_div_sim$x)

# . x .
# . o .
# . . .
dat_env_div_sim$upper <- 1 
# . . .
# . o .
# . x .
dat_env_div_sim$lower <- 1
# . . .
# . o x
# . . .
dat_env_div_sim$right <- 1
# . . .
# x o .
# . . .
dat_env_div_sim$left <- 1
# x . .
# . o .
# . . .
dat_env_div_sim$upper_left <- 1
# . . x
# . o .
# . . .
dat_env_div_sim$upper_right <- 1
# . . .
# . o .
# x . .
dat_env_div_sim$lower_left <- 1
# . . .
# . o .
# . . x
dat_env_div_sim$lower_right <- 1

# cells at top edge have no upper neighbours
dat_env_div_sim$upper[dat_env_div_sim$y == max(lat_list)] <- 0
# cells at bottom edge have no lower neighbours
dat_env_div_sim$lower[dat_env_div_sim$y == min(lat_list)] <- 0
# cells at right edge have no right neighbours
dat_env_div_sim$right[dat_env_div_sim$x == max(lon_list)] <- 0
# cells at left edge have no left neighbours
dat_env_div_sim$left[dat_env_div_sim$x == min(lon_list)] <- 0
# and combinations for corners
dat_env_div_sim$upper_left[dat_env_div_sim$y == max(lat_list) | dat_env_div_sim$x == min(lon_list)] <- 0
dat_env_div_sim$upper_right[dat_env_div_sim$y == max(lat_list) | dat_env_div_sim$x == max(lon_list)] <- 0
dat_env_div_sim$lower_left[dat_env_div_sim$y == min(lat_list) | dat_env_div_sim$x == min(lon_list)] <- 0
dat_env_div_sim$lower_right[dat_env_div_sim$y == min(lat_list) | dat_env_div_sim$x == max(lon_list)] <- 0

# check whether upper and lower neighbours are active
for(i in 1:length(lat_list)){
		index <- which(dat_env_div_sim$y == lat_list[i])
		if(i > 1){
			index_up <- which(dat_env_div_sim$y == lat_list[i-1])
			dat_env_div_sim$upper[index] <- dat_env_div_sim$active[index_up]
		}

		if(i < length(lat_list)){
			index_low <- which(dat_env_div_sim$y == lat_list[i+1])
			dat_env_div_sim$lower[index] <- dat_env_div_sim$active[index_low]
		}
}

# check whether left and right neighbours are active
for(i in 1:length(lon_list)){
	index <- which(dat_env_div_sim$x == lon_list[i])
	if(i > 1){
		index_left <- which(dat_env_div_sim$x == lon_list[i-1])
		dat_env_div_sim$left[index] <- dat_env_div_sim$active[index_left]
	}

	if(i < length(lon_list)){
		index_right <- which(dat_env_div_sim$x == lon_list[i+1])
		dat_env_div_sim$right[index] <- dat_env_div_sim$active[index_right]
	}
}

# check whether corner neighbours are active
for(i in 1:length(lat_list)){
	for(j in 1:length(lon_list)){
		index <- which(dat_env_div_sim$y == lat_list[i] & dat_env_div_sim$x == lon_list[j])
		if(i > 1 && j > 1){
			index_upper_left <- which(dat_env_div_sim$y == lat_list[i-1] & dat_env_div_sim$x == lon_list[j-1])
			dat_env_div_sim$upper_left[index] <- dat_env_div_sim$active[index_upper_left]
		}

		if(i > 1 && j < length(lon_list)){
			index_upper_right <- which(dat_env_div_sim$y == lat_list[i-1] & dat_env_div_sim$x == lon_list[j+1])
			dat_env_div_sim$upper_right[index] <- dat_env_div_sim$active[index_upper_right]
		}

		if(i < length(lat_list) && j > 1){
			index_lower_left <- which(dat_env_div_sim$y == lat_list[i+1] & dat_env_div_sim$x == lon_list[j-1])
			dat_env_div_sim$lower_left[index] <- dat_env_div_sim$active[index_lower_left]
		}

		if(i < length(lat_list) && j < length(lon_list)){
			index_lower_right <- which(dat_env_div_sim$y == lat_list[i+1] & dat_env_div_sim$x == lon_list[j+1])
			dat_env_div_sim$lower_right[index] <- dat_env_div_sim$active[index_lower_right]
		}
	}
}

# if the cell itself is not active, also indicate it has no active neighbours
dat_env_div_sim$upper[dat_env_div_sim$active == 0] <- 0
dat_env_div_sim$lower[dat_env_div_sim$active == 0] <- 0
dat_env_div_sim$right[dat_env_div_sim$active == 0] <- 0
dat_env_div_sim$left[dat_env_div_sim$active == 0] <- 0
dat_env_div_sim$upper_left[dat_env_div_sim$active == 0] <- 0
dat_env_div_sim$upper_right[dat_env_div_sim$active == 0] <- 0
dat_env_div_sim$lower_left[dat_env_div_sim$active == 0] <- 0
dat_env_div_sim$lower_right[dat_env_div_sim$active == 0] <- 0

# plot these geometries
# definition shp_atl_land -> line 108
p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(active)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='active') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_active.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(upper)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='upper') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_upper.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(lower)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='lower') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_lower.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(left)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='left') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_left.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(right)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='right') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_right.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=shp_atl_land) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(lower_left)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='lower left') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('CPR_grid_lower_left.png',.,device='png',width=18,height=12,units='cm')

# transform degrees Celsius to Kelvin
dat_env_div_sim$temp <- dat_env_div_sim$temp + 273.15
# normalise the other environmental variables: min = 1
dat_env_div_sim$phyto <- dat_env_div_sim$phyto / min(dat_env_div_sim$phyto,na.rm=TRUE)
dat_env_div_sim$velo <- dat_env_div_sim$velo / min(dat_env_div_sim$velo,na.rm=TRUE)

# write out final simulation grid
write.csv(dat_env_div_sim,'CPR_grid_env_div_sim.csv')