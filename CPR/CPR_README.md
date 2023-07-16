# Thesis-Data-Analysis
## Master Thesis IMBRSea

The Physics of Biodiversity:
exploring the dynamics behind spatial biodiversity patterns  
<br/>
contact: kobe.simoens@imbrsea.eu  
date: 01/08/2023  

---

## Data extraction for the CPR analysis

---

### CPR_extractData.R

R-script to perform the data extraction and manipulation.  
Maps are created as well.  
Constructed in R version 4.3 ("Already Tomorrow").  
Read comments in the code to avoid problems with memory usage.

---

### CPR_env.zip

zip-file with all the downloaded data on environment.   
Downloadable from: https://drive.google.com/drive/folders/1HJ-QKm95RALnsqqGJZJo9WAxAD1Eoao4?usp=drive_link   
The file is unzipped automatically in **CPR_extractData.R**.  
Therefore, **CPR_extractData.R** has to be placed in the same directory as the zip-file. 
<br/>
The zip-file contains:  

- *CPR_temp.tif*: North Atlantic Sea Surface Temperature data  
The data is downloaded from Bio-ORACLE using the R package sdmpredictors.  
Details on the data format can be found here: https://doi.org/10.1111/j.1466-8238.2011.00656.x .

- *CPR_phyto.tif*: North Atlantic phytoplankton concentration data  
The data is downloaded from Bio-ORACLE using the R package sdmpredictors.  
Details on the data format can be found here: https://doi.org/10.1111/j.1466-8238.2011.00656.x .

- *CPR_velo_x.nc*: North Atlantic zonal current velocity data  
The data is downloaded from Copernicus: https://doi.org/10.48670/moi-00016  .
Details on the data format can be found under this link.

- *CPR_velo_x.nc*: North Atlantic meridional current velocity data  
The data is downloaded from Copernicus: https://doi.org/10.48670/moi-00016  .
Details on the data format can be found under this link.

---

### CPR_SpeciesList.csv

List of planktonic species collected in the CPR sampling campaigns.  
The list was compiled manually for the purpose of this thesis.  
The columns are:  

- **Species**: species classification

- **Genus**: genus classification

- **Family**: familia classification 

- **Order**: ordo classification

- **Class**: classis classification

- **Fylum**: fylum classification

- **count**: number of entries of the species in the CPR data set

---

### CPR_Metridia.csv

Count data of *Metridia lucens* species.  
Due to an unkown error in the rgbif data extraction, the data for this species is extracted directly from:  
https://doi.org/10.17031/1629  
The columns are defined in the documentation on the rgbif package.  
The file is downloaded in a tab-separated format.  
The extracted columns are:

- **decimalLongitude**: the longitudinal sample coordinate  
Coordinate system: WGS84

- **decimalLatitude**: the latitudinal sample coordinate  
Coordinate system: WGS84

- **speciesKey**: unique GBIF species identifier

- **year**: year in which the sample was collected 

---

### CPR_combined_div.csv

Count data of all species and all samples combined in one csv file.  
The csv file is generated in **CPR_extractData.R**.  
The columns are:

- **lon**: the longitudinal sample coordinate   
Coordinate system: WGS84

- **lat**: the latitudinal sample coordinate   
Coordinate system: WGS84

- **key**: unique GBIF species identifier

- **year**: year in which the sample was collected 

---

### CPR_grid_base.csv

Basic information on the CPR simulation grid.  
The csv file is generated in **CPR_extractData.R**.  
The columns are:

- **x**: longitudinal coordinate of the centre of the grid cell  
Coordinate system: Lambert Conformal Conic as defined in **CPR_extractData.R**.

- **y**: latitudinal coordinate of the centre of the grid cell  
Coordinate system: Lambert Conformal Conic as defined in **CPR_extractData.R**.

- **area**: surface area of the grid cell  
Unit: square kilometres (km<sup>2</sup>)

- **peri**: circumference of the grid cell  
Unit: kilometres (km)

- **land_perc**: fraction of the grid cell surface area covered by land  
NA-value = no land in the grid cell

---

### CPR_grid_env.csv

Environmental data mapped onto the simulation grid.  
The csv file is generated in **CPR_extractData.R**.  
The columns are:

- **x**: longitudinal coordinate of the centre of the grid cell  
Coordinate system: Lambert Conformal Conic as defined in **CPR_extractData.R**.

- **y**: latitudinal coordinate of the centre of the grid cell  
Coordinate system: Lambert Conformal Conic as defined in **CPR_extractData.R**.

- **land_perc**: fraction of the grid cell surface area covered by land  
NA-value = no land in the grid cell

- **temp**: mean sea surface temperature value in the grid cell  
Unit: degrees Celsius (<sup>o</sup>C)  
NA-value = no data = land cell

- **phyto**: mean phytoplankton concentration value in the grid cell  
Unit: micromoles per cubic metre (μmol/m<sup>3</sup>)  
NA-value = no data = land cell

- **velo**: mean horizontal current kinetic energy value in the grid cell  
Unit: metres per second (m/s)  
NA-value = no data = land cell

---

### CPR_grid_div.csv

The diversity matrix for the simulation grid.  
The csv file is generated in **CPR_extractData.R**.  
The columns are:

- **x**: longitudinal coordinate of the centre of the grid cell  
Coordinate system: Lambert Conformal Conic as defined in **CPR_extractData.R**.

- **y**: latitudinal coordinate of the centre of the grid cell  
Coordinate system: Lambert Conformal Conic as defined in **CPR_extractData.R**.

- **diversity**: normalised total number of species found in the grid cell  
max = 1  
NA-value = no counts in the grid cell 

- **#key#**:  
All remaining columns represent a unique species.  
Column names are the unique GBIF identifiers for the species.  
Values are the total number of entries of a particular species found in a particular grid cell.  
This number is only indicative and cannot be used as a real abundance. 

---

### CPR_grid_env_div_sim.csv

Final CPR simulation grid that is used in the simulation of the Mechanistic Model.  
The csv file is generated in **CPR_extractData.R**.  
The columns are:

- **x**: longitudinal coordinate of the centre of the grid cell  
Coordinate system: Lambert Conformal Conic as defined in **CPR_extractData.R**.

- **y**: latitudinal coordinate of the centre of the grid cell  
Coordinate system: Lambert Conformal Conic as defined in **CPR_extractData.R**.

- **land_perc**: fraction of the grid cell surface area covered by land  
NA-value = land fraction higher than 0.8 OR no species counts

- **temp**: mean temperature value in the grid cell  
Units: Kelvin (K)  
NA-value = land fraction higher than 0.8 OR no species counts

- **phyto**: normalised mean phytoplankton concentration value in the grid cell  
min = 1   
NA-value = land fraction higher than 0.8 OR no species counts

- **velo**: normalised mean horizontal current kinetic energy value in the grid cell  
min = 1   
NA-value = land fraction higher than 0.8 OR no species counts

- **diversity**: normalised total number of species found in the grid cell  
max = 1    
NA-value = land fraction higher than 0.8 OR no species counts

- **active**: is grid cell active in the simulation?  
0-value: inactive = land fraction larger than 0.8  
1-value: active = land fraction smaller than 0.8  

- **upper**: activity of the upper neighbour  

- **lower**: activity of the lower neighbour

- **right**: activity of the right neighbour

- **left**: activity of the left neighbour

- **upper_left**:activity of the upper-left neighbour

- **upper_right**: activity of the upper-right neighbour

- **lower_left**: activity of the lower-left neighbour

- **lower_right**: activity of the lower-right neighbour

---
