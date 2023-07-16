# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Data extraction for the NABB analysis


# ----------------------------- LOAD PACKAGES ---------------------------
library(tidyverse)
library(rnaturalearth)
library(raster)
library(sf)
sf_use_s2(FALSE)
library(sp)
library(ggspatial)
library(vegan)
# -----------------------------------------------------------------------


# ---------------------------- UNZIP FILES ------------------------------

# unzip the main file: 2022Release_Nor.zip / Environment 
# downloadable from: https://drive.google.com/drive/folders/1HJ-QKm95RALnsqqGJZJo9WAxAD1Eoao4?usp=drive_link
unzip('NABB_data.zip')

# unzip the occurrence data 
# downloaded from: https://doi.org/10.5066/P97WAZE5
unzip('2022Release_Nor.zip',files='States.zip')
# unzip the occurrence data for each state
unzip('States.zip')

# unzip the overview of sampling routes
unzip('2022Release_Nor.zip',files='Routes.zip')
unzip('Routes.zip')

# unzip the species list
unzip('2022Release_Nor.zip',files='SpeciesList.txt')

# read in the sampling route information and create unique ID
# unique ID = state + route number
routes <- read.csv('routes.csv',header=T) %>% 
			mutate(routeID=paste0(StateNum,'_',Route))

# -----------------------------------------------------------------------


# ------------------------- PLOT ALL SAMPLES ----------------------------

# Canada and USA shapefile
shp <- rnaturalearth::ne_countries(country=c('canada','united states of america'),scale='medium')
# crop to data extent
shp_extent <- raster::crop(shp,raster::extent(min(routes$Longitude),max(routes$Longitude),min(routes$Latitude),max(routes$Latitude)))

# use WGS84 for plotting sampling locations in Canada and USA 
# (union leaves ugly lines)
shp_sf <- sf::st_union(sf::st_as_sf(shp_extent)) %>% sf::st_transform(.,crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
# get sample locations as points
smpl_points <- sp::SpatialPoints(coords=routes %>% dplyr::select(c(Longitude,Latitude)),
									proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) %>%
				sf::st_as_sf()
# include Mexico for the plot
shp_mex <- rnaturalearth::ne_countries(country='mexico',scale='medium') %>% sf::st_as_sf() %>%
			sf::st_transform(.,crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# make the plot
p <- ggplot() + geom_sf(data=shp_sf) + geom_sf(data=shp_mex) +
		geom_sf(data=smpl_points,size=0.25) + theme_bw() + 
		theme(panel.background = element_rect(fill = "#d0dfff"),
				panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',colour = "grey")) +
		coord_sf(expand=FALSE) + scale_x_continuous(limits = c(-179, -50)) + scale_y_continuous(limits=c(20,85),breaks=seq(30,80,10)) +
		ggspatial::annotation_north_arrow(location = "tl",which_north = "true", pad_x = unit(0.5, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,	height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.1,'cm'), pad_x=unit(0.5,'cm'))

# uncomment to get the png file
#p %>% ggsave('NABB_smpl_points.png',.,device='png',width=15,height=10,units='cm')

# ----------------------------------------------------------------------- 


# --------------------CREATE SUBAREA CONTIGUOUS USA ---------------------

# get contiguous USA shapefile
USA <- rnaturalearth::ne_countries(country='united states of america',scale='medium',returnclass='sf') %>% sf::st_geometry() %>%
		sf::st_crop(xmin=-130,xmax=-50,ymin=20,ymax=55)

# transform to Albers projection (equal area)
USA_t <- USA %>% sf::st_transform(crs=5070)

# create raster for simulation and data extraction
r <- raster::raster(ncol=45,nrow=14,ext=raster::extent(sf::st_bbox(USA_t)),crs=crs(USA_t))
# set dummy value for plotting
values(r) = 0

# transform to dataframe for ggplot
r_df_plot <- as.data.frame(r, xy = TRUE, na.rm = TRUE)

# plot the raster grid
p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=r_df_plot,aes(x=x,y=y),fill=NA,color='black') + theme_bw() +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank())

# uncomment to get the png file
#p %>% ggsave('NABB_grid.png',.,device='png',width=15,height=10,units='cm')

# alternative to extract coordinates
xy <- r %>% raster::rasterToPoints() %>% as.data.frame() %>% dplyr::select(c(x,y))

# polygonise the raster for intersection -> list of polygons 
r_poly <- r %>% raster::rasterToPolygons() %>% as(.,'SpatialPolygons') %>% sf::st_as_sf(.)

# initialise output dataframe with surface area, perimeter and land cover for each raster cell
r_df_out <- matrix(ncol=3,nrow=length(r_poly[[1]])) %>% as.data.frame()
names(r_df_out) <- c('area','peri','land_perc')

# r_poly[[1]][1] = cell 1 of raster
for(i in 1:length(r_poly[[1]])){
	# get intersection with land
	# error: great lakes are counted as 'land'
	a <- (sf::st_intersection(r_poly[[1]][i],USA_t) %>% sf::st_area()) / (r_poly[[1]][i] %>% sf::st_area())
	if(length(a) == 0){
		# if no land, NA
		a <- NA
	}

	# transform m to km
	r_df_out$area[i] <- (r_poly[[1]][i] %>% sf::st_area() )/1e6
	# transform to line in order to calculate perimeter
	# transform m to km
	r_df_out$peri[i] <- (r_poly[[1]][i] %>% sf::st_cast('MULTILINESTRING') %>% sf::st_length()) / 1e3
	# store land percentage
	r_df_out$land_perc[i] <- a
}

# combine with coordinates
r_df_out <- cbind(xy,r_df_out)
# write out the base simulation grid
write.csv(r_df_out,'NABB_grid_base.csv')

# plot land cover ratios
p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=r_df_out,aes(x=x,y=y,fill=land_perc),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='land percentage') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_land_perc.png',.,device='png',width=18,height=12,units='cm')

# ----------------------------------------------------------------------- 


# ---------------------- EXTRACT ENVIRONMENTAL VARIABLES ----------------
# note: this analysis takes some time due to the size of the raster layers
# ignore warnings
# problems of (RAM) memory could arise due to the size of rasters

# get environmental rasters and reproject to Albers projection

# temperature; downloaded from: https://prism.oregonstate.edu/normals/
r_temp <- raster::raster('Environment/PRISM_tmean_30yr_normal_4kmM4_annual_asc.asc') %>% raster::projectRaster(.,crs=5070)
# precipitation; downloaded from: https://prism.oregonstate.edu/normals/
r_prec <- raster::raster('Environment/PRISM_ppt_30yr_normal_4kmM4_annual_asc.asc') %>% raster::projectRaster(.,crs=5070)
# bathymetry; downloaded from: https://www.gebco.net/data_and_products/gridded_bathymetry_data/
r_bath <- raster::raster('Environment/GEBCO_bathymetry.tif') %>% raster::projectRaster(.,crs=5070) %>% raster::crop(.,r_prec)
# forest biomass; downloaded from: https://data.fs.usda.gov/geodata/rastergateway/biomass/
r_fore <- raster::raster('Environment/conus_forest_biomass_mg_per_ha.img')

# transform to points
# this is an approximation. However, the resolution of the environmental data >> the one of the simulation grid. 
# the error will be small
point_temp <- r_temp %>% raster::rasterToPoints() %>% as.data.frame() %>% 
				rename(temp=PRISM_tmean_30yr_normal_4kmM4_annual_asc) %>% sf::st_as_sf(coords=c('x','y'),crs=5070)
point_prec <- r_prec %>% raster::rasterToPoints() %>% as.data.frame() %>% 
				rename(prec=PRISM_ppt_30yr_normal_4kmM4_annual_asc) %>% sf::st_as_sf(coords=c('x','y'),crs=5070)

# change all bathymetric values below zero to NA
r_bath_NA <- r_bath %>% raster::reclassify(.,c(-Inf,0,NA))
# create a gradient layer out of the bathymetry layer
# from: https://gis.stackexchange.com/questions/260013/vector-gradient-of-a-raster-image-in-r
r_grad <- sqrt( (r_bath_NA %>% raster::focal(.,matrix(c(-1/2,0,1/2),nrow=3)) )^2 +  (r_bath_NA %>% raster::focal(.,matrix(c(-1/2,0,1/2),ncol=3)) )^2)
# crop forest layer to contiguous USA
# definition USA_t -> line 91
r_fore_crop <- raster::mask(r_fore,as(USA_t,'Spatial'))

# initialise the environmental containers for the simulation grid
temp_vector <- rep(NULL,length(r_poly[[1]]))
prec_vector <- rep(NULL,length(r_poly[[1]]))
grad_vector <- rep(NULL,length(r_poly[[1]]))
fore_vector <- rep(NULL,length(r_poly[[1]]))

# for each raster cell
# definition r_poly -> line 115
for(i in 1:length(r_poly[[1]])){
	print(paste0('cell ',i,' of ',length(r_poly[[1]])))
	# TEMPERATURE
	# intersection between environmental layer and simulation grid
	p_sub_temp <- sf::st_intersection(point_temp,r_poly[[1]][i])
	# if intersection with cell is empty, set value to NA
	if(length(p_sub_temp[[1]])==0){
		t <- NA
	} else{
		# value of simulation cell = mean of all points/cells of environmental variable
		# ! APPROXIMATION
		t <- mean(p_sub_temp$temp)
	}
	temp_vector[i] <- t

	# PRECIPITATION
	p_sub_prec <- sf::st_intersection(point_prec,r_poly[[1]][i])
	if(length(p_sub_prec[[1]])==0){
		pc <- NA
	} else{
		pc <- mean(p_sub_prec$prec)
	}
	prec_vector[i] <- pc
	
	# GRADIENT
	# first crop gradient raster to cell size because resolution too high
	r_grad_sub <- raster::crop(r_grad,sf::st_bbox(r_poly[[1]][i]))
	# transform cropped raster to points
	point_grad_sub <- r_grad_sub %>% raster::rasterToPoints() %>% as.data.frame() %>% 
					rename(grad=layer) %>% sf::st_as_sf(coords=c('x','y'),crs=5070)
	# if temperature or precipitation are NA, gradient is also NA
	if(is.na(pc) | is.na(t)){
		gr <- NA
	} else {
		gr <- mean(point_grad_sub$grad,na.rm=TRUE)
	}
	grad_vector[i] <- gr

	# FOREST BIOMASS
	# first crop forest raster to cell size because resolution too high
	r_fore_sub <- raster::crop(r_fore_crop,sf::st_bbox(r_poly[[1]][i]))
	# transform cropped raster to points
	point_fore_sub <- r_fore_sub %>% raster::rasterToPoints() %>% as.data.frame() %>% rename(fore=Layer_1) 
	# if temperature or precipitation are NA, forest is also NA
	if(is.na(pc) | is.na(t)){
		fo <- NA
	} else {
		fo <- mean(point_fore_sub$fore,na.rm=TRUE)
	}
	fore_vector[i] <- fo

}

# read in the earlier created grid and add environmental variables
r_df_out <- read.csv('NABB_grid_base.csv',header=TRUE,row.names=1) %>% 
		mutate(temp=temp_vector,prec=prec_vector,grad=grad_vector,fore=fore_vector)
# adjust for edge effects (for example, the Great Lakes)
r_df_out$land_perc[is.na(r_df_out$temp)] <- NA
r_df_out$land_perc[is.na(r_df_out$prec)] <- NA
r_df_out$land_perc[is.na(r_df_out$grad)] <- NA
r_df_out$land_perc[is.na(r_df_out$fore)] <- NA
r_df_out$temp[is.na(r_df_out$land_perc)] <- NA
r_df_out$prec[is.na(r_df_out$land_perc)] <- NA
r_df_out$grad[is.na(r_df_out$land_perc)] <- NA
r_df_out$fore[is.na(r_df_out$land_perc)] <- NA

# plot the environmental variables on the simulation grid
# definition USA_t -> line 91
# TEMPERATURE
p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=r_df_out,aes(x=x,y=y,fill=temp),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='temperature (°C)') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_temp.png',.,device='png',width=18,height=12,units='cm')

# PRECIPITATION
p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=r_df_out,aes(x=x,y=y,fill=prec),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='precipitation (mm)') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_prec.png',.,device='png',width=18,height=12,units='cm')

# GRADIENT
p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=r_df_out,aes(x=x,y=y,fill=grad),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='vertical gradient (m)') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_grad.png',.,device='png',width=20,height=12,units='cm')

# FOREST BIOMASS
p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=r_df_out,aes(x=x,y=y,fill=fore),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='forest biomass (Mg/ha)') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_fore.png',.,device='png',width=20,height=12,units='cm')

# write out updated grid
write.csv(r_df_out %>% dplyr::select(c(x,y,land_perc,temp,prec,grad,fore)),'NABB_grid_env.csv')

# ----------------------------------------------------------------------- 


# ----------------------- EXTRACT SPECIES LIST --------------------------
# clean up the species list (a total mess)

# read files with species
spec_list <- read.delim('SpeciesList.txt',skip=14,header=FALSE)

# initialise output species list: AOU = unique species identifier
df_spec <- matrix(ncol=9,nrow=nrow(spec_list)) %>% as.data.frame()
names(df_spec) <- c('Seq','AOU','English','French','Spanish','Order','Family','Genus','Species')

# for each row in list = for each species
for(j in 1:nrow(spec_list)){
	# extract row; 
	# since R version 4.3.0 -> convert explicitly to UTF-8 characters
	t <- iconv(spec_list[j,], from = "ISO-8859-1", to = "UTF-8")
	# split based on spaces
	splt <- strsplit(t,split=' ')[[1]]
	# counters for bookkeeping
	count <- 0
	count2 <- 0
	for(i in 1:length(splt)){
		c <- splt[i]
		# add new split: '!' on locations based on messy file
		if(c != '' && count!=0 && (count2 < 1 || count2 > 2)){
			splt[i-1] <- '!'
			count <- 0
			count2 <- count2 + 1
		} else if(c != '' && count==0 && (count2 > 0 && count2 < 3)){
			splt[i] <- paste0('!',splt[i])
			count2 <- count2 + 1
		} else if(c == ''){
			count <- count + 1
		}
	}
	# create string with new seperators
	str_new <- paste(splt,collapse='')
	# resplit following new seperators
	splt <- strsplit(str_new,split='!')[[1]][-1]
	# manual adjustments for inconsistent format
	if(j == 28){
		splt <- c(splt[1:5],'Anseriformes',splt[6:8])
	}
	if(j == 128){
		splt <- c(splt[1:4],'not',splt[5:8])
	}
	if(j == 477){
		splt <- c(splt[1:4],'not',splt[5:8])
	}
	if(j == 514){
		splt <- c(splt[1:4],'not',splt[5:8])
	}
	if(j == 691){
		splt <- c(splt[1:4],'not',splt[5:8])
	}

	df_spec[j,] <- splt
}

df_spec <- df_spec %>% dplyr::select(c(Seq,AOU,Order,Family,Genus,Species))
# write out new species list in csv format
write.csv(df_spec,'NABB_SpeciesList.csv')

# ----------------------------------------------------------------------- 


# ------------------ COMBINE DIVERSITY DATA -----------------------------

# list of all zipped state files
files <- list.files(path="States", pattern="*.zip", full.names=TRUE, recursive=FALSE)

# unzip all these files in 'States' directory
for(file in files){
	unzip(file,exdir='States')
}

# list of all state csv files
files <- list.files(path="States", pattern='*.csv', full.names=TRUE, recursive=FALSE)

# initialise container for species counts
dat_all_counts <- matrix(ncol=6,nrow=0) %>% as.data.frame()
names(dat_all_counts) <- c('AOU','CountryNum','StateNum','SpeciesTotal','Longitude','Latitude')
# container for 2021 subset
dat_all_counts_2021 <- matrix(ncol=6,nrow=0) %>% as.data.frame()
names(dat_all_counts_2021) <- c('AOU','CountryNum','StateNum','SpeciesTotal','Longitude','Latitude')

# for each state
for(file in files){
	# read all occurrence data
	dat <- read.csv(file,header=T) %>%
			# create unique route ID
			dplyr::mutate(routeID=paste0(StateNum,'_',Route)) %>%
			dplyr::select(c(AOU,SpeciesTotal,routeID,CountryNum,StateNum)) %>%
			# group by species and route and calculate total number of observations
			# CountryNum and StateNum needed to subset later
			dplyr::group_by(AOU,routeID) %>% 
			dplyr::summarise(SpeciesTotal=sum(SpeciesTotal),CountryNum=median(CountryNum),StateNum=median(StateNum))

	# same calculations for 2021 subset
	dat_2021 <- read.csv(file,header=T) %>% dplyr::filter(Year %in% c(2021)) %>%
			dplyr::mutate(routeID=paste0(StateNum,'_',Route)) %>%
			dplyr::select(c(AOU,SpeciesTotal,routeID,CountryNum,StateNum)) %>%
			dplyr::group_by(AOU,routeID) %>% 
			dplyr::summarise(SpeciesTotal=sum(SpeciesTotal),CountryNum=median(CountryNum),StateNum=median(StateNum))

	# link routeID to longitude and latitude in routes dataframe;
	# definition routes -> line 44
	dat_lonlat <- dat %>% dplyr::inner_join(.,routes %>% dplyr::select(c(Longitude,Latitude,routeID)))
	dat_lonlat_2021 <- dat_2021 %>% dplyr::inner_join(.,routes %>% dplyr::select(c(Longitude,Latitude,routeID)))
	# add to final container 
	dat_all_counts <- rbind(dat_all_counts,dat_lonlat)
	dat_all_counts_2021 <- rbind(dat_all_counts_2021,dat_lonlat_2021)
}

# ungroup
dat_all_counts <- dat_all_counts %>% dplyr::ungroup() 
dat_all_counts_2021 <- dat_all_counts_2021 %>% dplyr::ungroup() 
# only USA without Alaska
dat_all_counts_USA <- dat_all_counts %>% dplyr::filter(CountryNum==840 & StateNum != 3) 
dat_all_counts_2021_USA <- dat_all_counts_2021 %>% dplyr::filter(CountryNum==840 & StateNum != 3) 
# calculate how many observations are discarded
attrition <- 1 - nrow(dat_all_counts_USA)/nrow(dat_all_counts)
attrition_2021 <- 1 - nrow(dat_all_counts_2021_USA)/nrow(dat_all_counts_2021)
# write out the combined data
write.csv(dat_all_counts,'NABB_combined_div.csv')
write.csv(dat_all_counts_2021,'NABB_combined_div_2021.csv')

# ----------------------------------------------------------------------- 


# ----------------- MAKE ROUTES DATASET ---------------------------------

# read files with combined data
# subset for contiguous USA
dat_all_counts <- read.csv('NABB_combined_div.csv',header=TRUE,row.names=1) %>% dplyr::filter(CountryNum==840 & StateNum != 3) 
# same calculations for 2021 subset
dat_all_counts_2021 <- read.csv('NABB_combined_div_2021.csv',header=TRUE,row.names=1) %>% dplyr::filter(CountryNum==840 & StateNum != 3) 

# create species counts (column) for each route (row)
dat_routes <- dat_all_counts %>% tidyr::pivot_wider(values_from=SpeciesTotal,names_from=AOU) %>% 
					dplyr::select(-c(routeID,CountryNum,StateNum))
dat_routes_2021 <- dat_all_counts_2021 %>% tidyr::pivot_wider(values_from=SpeciesTotal,names_from=AOU) %>% 
					dplyr::select(-c(routeID,CountryNum,StateNum))
# transform NA to 0
dat_routes[is.na(dat_routes)] <- 0
dat_routes_2021[is.na(dat_routes_2021)] <- 0

# write out the routes data
write.csv(dat_routes,'NABB_routes_div.csv')
write.csv(dat_routes_2021,'NABB_routes_div_2021.csv')

# plot/map diversity for all samples
# calculate normalised diversity
dat_routes_div <- dat_routes %>% dplyr::mutate(diversity=dat_routes[,-c(1,2)] %>% vegan::decostand(method='pa') %>% rowSums()) %>% 
					dplyr::select(c(Longitude,Latitude,diversity)) %>%
					dplyr::mutate(diversity=diversity/max(diversity))

# transform to crs 5070
dat_routes_plot <- cbind(dat_routes %>% sf::st_as_sf(.,coords=c('Longitude','Latitude'),crs=4326) %>% sf::st_transform(crs=5070) %>%
							sf::st_coordinates() %>% as.data.frame() %>% dplyr::rename(x=X,y=Y) , 
						dat_routes_div %>% dplyr::select(c(diversity)))

# definition USA_t -> line 91
p <- ggplot() + geom_sf(data=USA_t) + geom_point(data=dat_routes_plot,aes(x=x,y=y,col=diversity),size=0.75) + theme_bw() +
	scale_colour_viridis_c(option='magma') +
	ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
	ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
	theme(axis.title=element_blank())

# uncomment to get the png file
#p %>% ggsave('NABB_routes_div.png',.,device='png',width=18,height=10,units='cm')

# subset for Passeriformes
spec_list <- read.csv('NABB_SpeciesList.csv',header=TRUE)
AOU_pass <- spec_list %>% dplyr::filter(Order=='Passeriformes') %>% pull(AOU) 
dat_routes_PASS <- dat_routes[,names(dat_routes) %in% c('Longitude','Latitude',AOU_pass)]
# delete routes without counts (if any)
dat_routes_PASS <- dat_routes_PASS[rowSums(dat_routes_PASS[,-c(1,2)]) != 0,]

# calculate Passeriformes diversity
dat_routes_PASS_div <- dat_routes_PASS %>% dplyr::mutate(diversity=dat_routes_PASS[,-c(1,2)] %>% vegan::decostand(method='pa') %>% rowSums()) %>% 
					dplyr::select(c(Longitude,Latitude,diversity)) %>%
					dplyr::mutate(diversity=diversity/max(diversity))

# transform to crs 5070
dat_routes_PASS_plot <- cbind(dat_routes_PASS %>% sf::st_as_sf(.,coords=c('Longitude','Latitude'),crs=4326) %>% sf::st_transform(crs=5070) %>%
							sf::st_coordinates() %>% as.data.frame() %>% dplyr::rename(x=X,y=Y) , 
						dat_routes_PASS_div %>% dplyr::select(c(diversity)))

# definition USA_t -> line 91
p <- ggplot() + geom_sf(data=USA_t) + geom_point(data=dat_routes_PASS_plot,aes(x=x,y=y,col=diversity),size=0.75) + theme_bw() +
	scale_colour_viridis_c(option='magma') +
	ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
	ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
	theme(axis.title=element_blank())

# uncomment to get the png file
#p %>% ggsave('NABB_routes_PASS_div.png',.,device='png',width=18,height=10,units='cm')

# ----------------------------------------------------------------------- 


# ------------------- MAP DIVERSITY TO GRID -----------------------------

# read the routes data
dat_routes <- read.csv('NABB_routes_div.csv', header=TRUE, row.names=1)  
# remove 'X' from column names
names(dat_routes)[-c(1,2)] <- sub('X','',names(dat_routes)[-c(1,2)])

# create spatial points for each route and transform coordinates
dat_points <- dat_routes %>% sf::st_as_sf(.,coords=c('Longitude','Latitude'),crs=4326) %>%
					sf::st_transform(crs=5070)

# initialise output container for diversity and composition for simulation grid
# definition r_poly -> line 115
dat_grid_div <- matrix(nrow=0,ncol=ncol(dat_routes)) %>% as.data.frame()
# subset for Passeriformes and do the same calculations
# definition AOU_pass -> line 492
dat_grid_PASS_div <- matrix(nrow=0,ncol=length(AOU_pass)+2) %>% as.data.frame()

# for each grid cell
for(i in 1:length(r_poly[[1]])){
	# get total number of observations for each species in grid cell
	inter <- sf::st_intersection(dat_points,r_poly[[1]][i]) %>% sf::st_drop_geometry() %>%
				colSums()
	inter <- inter %>% as.data.frame() %>% t()
	# get cell coordinates
	# definition xy -> line 112
	inter <- cbind(xy[i,],inter)
	rownames(inter) <- i
	dat_grid_div <- rbind(dat_grid_div,inter)

	inter_PASS <- inter[,names(inter) %in% c('x','y',AOU_pass)]
	
	dat_grid_PASS_div <- rbind(dat_grid_PASS_div,inter_PASS)
}

# calculate the diversity for each grid cell
dat_grid_div <- dat_grid_div %>% dplyr::mutate(diversity=dat_grid_div[,-c(1,2)] %>% vegan::decostand(method='pa') %>% rowSums())
# rearrange columns
dat_grid_div <- dat_grid_div[,c('x','y','diversity',names(dat_grid_div)[3:(length(names(dat_grid_div))-1)])]
dat_grid_div$diversity[dat_grid_div$diversity == 0] <- NA
# normalise
dat_grid_div$diversity <- dat_grid_div$diversity / max(dat_grid_div$diversity,na.rm=TRUE)

dat_grid_PASS_div <- dat_grid_PASS_div %>% dplyr::mutate(diversity=dat_grid_PASS_div[,-c(1,2)] %>% vegan::decostand(method='pa') %>% rowSums())
dat_grid_PASS_div <- dat_grid_PASS_div[,c('x','y','diversity',names(dat_grid_PASS_div)[3:(length(names(dat_grid_PASS_div))-1)])]
dat_grid_PASS_div$diversity[dat_grid_PASS_div$diversity == 0] <- NA
dat_grid_PASS_div$diversity <- dat_grid_PASS_div$diversity / max(dat_grid_PASS_div$diversity,na.rm=TRUE)

# write out grid diversity dataset
write.csv(dat_grid_div,'NABB_grid_div.csv')
write.csv(dat_grid_PASS_div,'NABB_grid_PASS_div.csv')

# plot the empircal diversity
# definition USA_t -> line 91
p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=dat_grid_div,aes(x=x,y=y,fill=diversity),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='diversity') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_div.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=dat_grid_PASS_div,aes(x=x,y=y,fill=diversity),color='black') + theme_bw() +
		scale_fill_viridis_c(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='diversity') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_PASS_div.png',.,device='png',width=18,height=12,units='cm')

# ----------------------------------------------------------------------- 


# ------------ MAKE FINAL SIMULATION GRID -------------------------------

# read environmental grid and add diversity
dat_env_div <- read.csv('NABB_grid_env.csv',header=TRUE,row.names=1) %>% 
				dplyr::select(c(x,y,land_perc,temp,prec,grad,fore)) %>%
				dplyr::mutate(diversity=read.csv('NABB_grid_div.csv',header=TRUE,row.names=1) %>% dplyr::pull(diversity),
								diversity_PASS=read.csv('NABB_grid_PASS_div.csv',header=TRUE,row.names=1) %>% dplyr::pull(diversity))

# to eliminate edge effects, only consider cells with more than 20% land cover
# if necessary, also discard cells with no count data
# all 'non-active' cells = NA
dat_env_div$land_perc[dat_env_div$land_perc < 0.2 | is.na(dat_env_div$diversity)] <- NA
dat_env_div$temp[is.na(dat_env_div$land_perc)] <- NA
dat_env_div$prec[is.na(dat_env_div$land_perc)] <- NA
dat_env_div$grad[is.na(dat_env_div$land_perc)] <- NA
dat_env_div$fore[is.na(dat_env_div$land_perc)] <- NA
dat_env_div$diversity[is.na(dat_env_div$land_perc)] <- NA
dat_env_div$diversity_PASS[is.na(dat_env_div$land_perc)] <- NA

# add whether the cell is active or not; land > 0.2 = active = 1; land < 0.2 = not active = 0 
dat_env_div_sim <- dat_env_div %>% mutate(active=1)
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
# definition USA_t -> line 91
p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(active)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='active') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_active.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(upper)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='upper') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_upper.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(lower)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='lower') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_lower.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(left)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='left') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_left.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(right)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='right') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_right.png',.,device='png',width=18,height=12,units='cm')

p <- ggplot() + geom_sf(data=USA_t) + geom_tile(data=dat_env_div_sim,aes(x=x,y=y,fill=factor(lower_left)),color='black') + theme_bw() +
		scale_fill_viridis_d(option='magma',na.value=rgb(1,1,1,0),alpha=0.65,name='lower left') +
		ggspatial::annotation_north_arrow(location = "tr",which_north = "true", pad_x = unit(2.4, "cm"), pad_y = unit(0.2, "cm"),
			style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm")) +
		ggspatial::annotation_scale(location = "bl", width_hint=0.3, text_cex=0.6, pad_y=unit(0.08,'cm'), pad_x=unit(0.5,'cm')) +
		theme(axis.title=element_blank(), legend.title=element_text(size=10))

# uncomment to get the png file
#p %>% ggsave('NABB_grid_lower_left.png',.,device='png',width=18,height=12,units='cm')

# transform degrees Celcius to Kelvin
dat_env_div_sim$temp <- dat_env_div_sim$temp + 273.15
# normalise other environmental variables: lowest value = 1
dat_env_div_sim$prec <- dat_env_div_sim$prec / min(dat_env_div_sim$prec,na.rm=TRUE)
dat_env_div_sim$grad <- dat_env_div_sim$grad / min(dat_env_div_sim$grad,na.rm=TRUE)
dat_env_div_sim$fore <- dat_env_div_sim$fore / min(dat_env_div_sim$fore,na.rm=TRUE)

# write final simulation grid
write.csv(dat_env_div_sim,'NABB_grid_env_div_sim.csv')

# ----------------------------------------------------------------------- 