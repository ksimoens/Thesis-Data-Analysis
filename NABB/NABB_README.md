# Thesis-Data-Analysis
## Master Thesis IMBRSea

The Physics of Biodiversity:
exploring the dynamics behind spatial biodiversity patterns  
<br/>
contact: kobe.simoens@imbrsea.eu  
date: 01/08/2023  

---

## Data extraction for the NABB analysis

---

### NABB_extractData.R

R-script to perform the data extraction and manipulation.  
Maps are created as well.  
Constructed in R version 4.3 ("Already Tomorrow").  
Read comments in the code to avoid problems with memory usage.

---

### NABB_data.zip

zip-file with all the downloaded data on environment and diversity.  
The file is unzipped automatically in **NABB_extractData.R**.  
Therefore, **NABB_extractData.R** has to be placed in the same directory as the zip-file.  
The entire data manipulation only requires **NABB_extractData.R** and the zip-file.  
<br/>
The zip-file contains:  

- **2022Release_Nor.zip**: the NABB diversity data  
The data is downloaded from https://doi.org/10.5066/P97WAZE5 .  
Details on the data format can be found under this link.

- **Environment**: directory with environmental data:  

	- *PRISM_tmean_30yr_normal_4kmM4_annual_asc.asc* + related files: contiguous USA temperature data  
	The data is downloaded from https://prism.oregonstate.edu/normals/ .  
	Details on the data format can be found under this link.

	- *PRISM_ppt_30yr_normal_4kmM4_annual_asc.asc* + related files: contiguous USA precipitation data  
	The data is downloaded from https://prism.oregonstate.edu/normals/ .  
	Details on the data format can be found under this link.

	- *GEBCO_bathymetry.tif*: contiguous USA bathymetry data  
	The data is downloaded from https://www.gebco.net/data_and_products/gridded_bathymetry_data/ .
	Details on the data format can be found under this link.  

	- *conus_forest_biomass_mg_per_ha.img*: contiguous USA forest biomass data  
	The data is downloaded from https://data.fs.usda.gov/geodata/rastergateway/biomass/ .  
	Details on the data format can be found under this link. 

---

### NABB_combined_div.csv

Count data of all species and all routes combined in one csv file.  
The csv file is generated in **NABB_extractData.R**.  
The columns are:

- **AOU**: 'American Ornithologists' Union' species identifier  
Details on the taxonomy can be found under https://doi.org/10.5066/P97WAZE5 .

- **routeID**: unique route identifier  
Created by combining the state number with the state route number.  

- **SpeciesTotal**: the total number of individuals of a particular species found in a particular route

- **CountryNum**: country identifier  
Details can be found under https://doi.org/10.5066/P97WAZE5 .

- **StateNum**: state identifier  
Details can be found under https://doi.org/10.5066/P97WAZE5 .

- **Longitude**: longitudinal coordinate of the route starting location  
Coordinate system: WGS84

- **Latitude**: latitudinal coordinate of the route starting location  
Coordinate system: WGS84

---

### NABB_combined_div_2021.csv

Exactly the same information as in **NABB_combined_div.csv**.  
The only difference is:

- **SpeciesTotal**: the total number of individuals of a particular species found in a particular route in 2021

---

### NABB_routes_div.csv

The diversity matrix for all sampling routes.  
The csv file is generated in **NABB_extractData.R**.  
The columns are:

- **Longitude**: longitudinal coordinate of the route starting location  
Coordinate system: WGS84

- **Latitude**: latitudinal coordinate of the route starting location  
Coordinate system: WGS84

- **#AOU#**:  
All remaining columns represent a unique species.  
Column names are the unique AOU identifiers for the species.  
Values are the total number of individuals of a particular species found in a particular route.

---

### NABB_routes_div_2021.csv

The diversity matrix for all sampling routes in 2021.  
Exactly the same information as in **NABB_routes_div.csv**.  
Values now represent the total number of individuals of a particular species found in a particular route in 2021.

---

### NABB_grid_base.csv

Basic information on the NABB simulation grid.  
The csv file is generated in **NABB_extractData.R**.  
The columns are:

- **x**: longitudinal coordinate of the centre of the grid cell  
Coordinate system: NAD83 Conus Albers

- **y**: latitudinal coordinate of the centre of the grid cell  
Coordinate system: NAD83 Conus Albers

- **area**: surface area of the grid cell  
Unit: square kilometres (km<sup>2</sup>)

- **peri**: circumference of the grid cell  
Unit: kilometres (km)

- **land_perc**: fraction of the grid cell surface area covered by land  
NA-value = no land in the grid cell

---

### NABB_grid_env.csv

Environmental data mapped onto the simulation grid.  
The csv file is generated in **NABB_extractData.R**.  
The columns are:

- **x**: longitudinal coordinate of the centre of the grid cell  
Coordinate system: NAD83 Conus Albers

- **y**: latitudinal coordinate of the centre of the grid cell  
Coordinate system: NAD83 Conus Albers

- **land_perc**: fraction of the grid cell surface area covered by land  
NA-value = no land in the grid cell

- **temp**: mean temperature value in the grid cell  
Unit: degrees Celsius (<sup>o</sup>C)  
NA-value = no land in the grid cell

- **prec**: mean precipitation value in the grid cell  
Unit: millimetres (mm)  
NA-value = no land in the grid cell

- **grad**: mean altitudinal gradient value in the grid cell  
Unit: metres (m)  
NA-value = no land in the grid cell

- **fore**: mean forest biomass value in the grid cell  
Unit: megagrammes per hectare (Mg/ha)  
NA-value = no land in the grid cell

---

### NABB_SpeciesList.csv

Reworking the species list for subsetting the count data.  
The csv file is generated in **NABB_extractData.R**.  
The list contains all species in the NABB dataset.  
The columns are:

- **Seq**: identifier within the NABB dataset

- **AOU**: 'American Ornithologists' Union' species identifier  
Details on the taxonomy can be found under https://doi.org/10.5066/P97WAZE5 .

- **Order**: ordo classification

- **Family**: familia classification

- **Genus**: genus classification

- **Species**: species classification

---

### NABB_grid_div.csv

The diversity matrix for the simulation grid.  
The csv file is generated in **NABB_extractData.R**.  
The columns are:

- **x**: longitudinal coordinate of the centre of the grid cell  
Coordinate system: NAD83 Conus Albers

- **y**: latitudinal coordinate of the centre of the grid cell  
Coordinate system: NAD83 Conus Albers

- **diversity**: normalised total number of species found in the grid cell  
max = 1  
NA-value = no counts in the grid cell 

- **#AOU#**:  
All remaining columns represent a unique species.  
Column names are the unique AOU identifiers for the species.  
Values are the total number of individuals of a particular species found in a particular grid cell. 

---

### NABB_grid_PASS_div.csv

Exactly the same information as in **NABB_grid_div.csv**.  
However, the dataset is restricted to species belonging to the order of Passeriformes.  
The differences are:

- **diversity**: normalised number of Passeriformes species found in the grid cell  
max = 1  
NA-value = no counts in the grid cell

- **#AOU#**:  
Only Passeriformes species are included in the columns.  
Values are the total number of individuals of a particular Passeriformes species found in a particular grid cell. 

---

### NABB_grid_env_div_sim.csv

Final NABB simulation grid that is used in the simulation of the Mechanistic Model.  
The csv file is generated in **NABB_extractData.R**.  
The columns are:

- **x**: longitudinal coordinate of the centre of the grid cell  
Coordinate system: NAD83 Conus Albers

- **y**: latitudinal coordinate of the centre of the grid cell  
Coordinate system: NAD83 Conus Albers

- **land_perc**: fraction of the grid cell surface area covered by land  
NA-value = land fraction smaller than 0.2

- **temp**: normalised mean temperature value in the grid cell  
min = 1   
NA-value = land fraction smaller than 0.2

- **prec**: normalised mean precipitation value in the grid cell  
min = 1   
NA-value = land fraction smaller than 0.2

- **grad**: normalised mean altitudinal gradient value in the grid cell  
min = 1   
NA-value = land fraction smaller than 0.2

- **fore**: normalised mean forest biomass value in the grid cell  
min = 1   
NA-value = land fraction smaller than 0.2

- **diversity**: normalised total number of species found in the grid cell  
max = 1    
NA-value = land fraction smaller than 0.2

- **diversity_PASS**: normalised number of Passeriformes species found in the grid cell  
max = 1  
NA-value = land fraction smaller than 0.2

- **active**: is grid cell active in the simulation?  
0-value: inactive = land fraction smaller than 0.2  
1-value: active = land fraction larger than 0.2  

- **upper**: activity of the upper neighbour  

- **lower**: activity of the lower neighbour

- **right**: activity of the right neighbour

- **left**: activity of the left neighbour

- **upper_left**:activity of the upper-left neighbour

- **upper_right**: activity of the upper-right neighbour

- **lower_left**: activity of the lower-left neighbour

- **lower_right**: activity of the lower-right neighbour

---
