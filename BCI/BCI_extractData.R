# Master Thesis IMBRSea
# The Physics of Biodiversity: 
# exploring the dynamics behind spatial biodiversity patterns
#
# contact: kobe.simoens@imbrsea.eu
# date: 01/08/2023
#
# Data extraction for the BCI analysis


# ----------------------------- LOAD PACKAGES ---------------------------
library(tidyverse)
library(ggspatial)
# -----------------------------------------------------------------------


# ----------------------- UNZIP DATA FILES ------------------------------
# unzip the downloaded data files
# downloaded from: https://doi.org/10.15146/5xcp-0d46

unzip('BCI_data.zip')
unzip('bci.tree.zip',files='bci.tree8.rdata')

# only use the final census
# change the final number (1-8) for another census
load('bci.tree8.rdata')
# only retain trees that were still alive during the census
dat_combined <- bci.tree8 %>% dplyr::filter(status=='A') %>% dplyr::select(c(treeID,sp,gx,gy))

# plot the study area and the trees
p <- dat_combined %>% ggplot() + geom_point(aes(x=gx,y=gy),cex=0.001) + theme_bw() +
							xlab('x (m)') + ylab('y (m)') + xlim(0,1050) +
							ggspatial::annotation_north_arrow(location = "tr", pad_x = unit(0.05, "cm"), pad_y = unit(0.2, "cm"),
													style = north_arrow_fancy_orienteering,height = unit(1.5, "cm"), width = unit(1.5, "cm"))

# uncomment to get the png file							
p %>% ggsave('CPR_samples.png',.,device='png',width=20,height=10,units='cm')

# -----------------------------------------------------------------------


# -------------------- CREATE SAMPLES -----------------------------------
# artificially create non-overlapping subsamples of certain size

# function to create the samples
# input: size of a sample in metres: kx = horizontal; ky = vertical
# please use integer fractions of the total sizes
# output: csv file with diversity matrix

extractSamples <- function(kx,ky){

	# get list of species
	# definition dat_combined -> line 29
	spec_list <- unique(dat_combined$sp)

	# transform continuous coordinate to discrete coordinate
	# discrete coordinate represents the order of the discrete sample
	# for example: a tree with continuous horizontal coordinate 496 m 
	# is located in horizontal cell 19 (starting from 0)
	# create unique cell identifier
	dat_combined <- dat_combined %>% dplyr::mutate(xc=floor(gx/kx),yc=floor(gy/ky)) %>% 
										dplyr::mutate(cell=paste0(xc,'_',yc))

	# calculate total number of counts per species per sample
	df_sum <- dat_combined %>% dplyr::group_by(cell,sp) %>% 
								dplyr::summarise(count=length(sp)) %>% dplyr::ungroup()


	# create diversity matrix
	df_new <- df_sum %>% pivot_wider(values_from=count,names_from=sp)  
	df_new[is.na(df_new)] <- 0

	# calculate the centre of the sample coordinates
	df_new$x <- NA
	df_new$y <- NA

	for(i in 1:nrow(df_new)){
		df_new$x[i] <- (str_split(df_new$cell[i],pattern='_')[[1]][1] %>% strtoi())*kx + kx/2
		df_new$y[i] <- (str_split(df_new$cell[i],pattern='_')[[1]][2] %>% strtoi())*ky + ky/2
	}

	df_new <- df_new[,c('x','y',names(df_new[2:(ncol(df_new)-2)]))]

	# write out the diversity matrix for the artificial samples
	write.csv(df_new,paste0('BCI_grid_',kx,'_',ky,'.csv'))

}

# horizontal size of sample:
kx <- 25 # metres
# vertical size of sample:
ky <- 25 # metres
# create the samples:
extractSamples(kx,ky)

# -----------------------------------------------------------------------