# Thesis-Data-Analysis
## Master Thesis IMBRSea

The Physics of Biodiversity:
exploring the dynamics behind spatial biodiversity patterns  
<br/>
contact: kobe.simoens@imbrsea.eu  
date: 01/08/2023  

---

## Data extraction for the BCI analysis

---

### BCI_extractData.R

R-script to perform the data extraction and manipulation.  
Maps are created as well.  
Constructed in R version 4.3 ("Already Tomorrow").

---

### BCI_data.zip

zip-file with all the downloaded data on diversity.  
The data is downloaded from https://doi.org/10.15146/5xcp-0d46 .     
The file is unzipped automatically in **BCI_extractData.R**.  
Therefore, **BCI_extractData.R** has to be placed in the same directory as the zip-file.  
<br/>
The zip-file contains:  

- **bci.tree.zip**: BCI tree census data of all 8 censuses  
Details on the data format can be found here: https://doi.org/10.15146/5xcp-0d46 .  
	- *bci.treeX.rdata* (X = 1-8): loads the data frame into R  
	Extracted columns are:  
		- **treeID**: unique individual tree identifier  

		- **sp**: unique species identifier to which the individual belongs  

		- **gx**: horizontal coordinate  
		Unit: metres (m)  
		Reference: southwest corner of the plot  

		- **gy**: vertical coordinate  
		Unit: metres (m)  
		Reference: southwest corner of the plot

	- *bci.spptable.rdata*: loads information on species  
	Not used in the data extraction; for reference only

---

### BCI_grid_25_25.csv

Diversity matrix for artificial samples.  
Horizontal size of a sample = 25 metres.  
Vertical size of a sample = 25 metres.  
The csv file is generated in **BCI_extractData.R**.  
The columns are:

- **x**: horizontal coordinate of sample centres  
Unit: metres (m)  
Reference: southwest corner of the plot  

- **y**: vertical coordinate of sample centres  
Unit: metres (m)  
Reference: southwest corner of the plot  

- **#sp#**:  
All remaining columns represent a unique species.  
Column names are the unique species identifiers for the species.    
Values are the total number of tree individuals of a particular species found in a particular sample.
---

### BCI_grid_10_10.csv

Exactly the same as for **BCI_grid_25_25.csv** but:  
Horizontal size of a sample = 10 metres.  
Vertical size of a sample = 10 metres. 

---