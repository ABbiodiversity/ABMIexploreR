<!-- badges: start -->

![In Development](https://img.shields.io/badge/Status-In%20Development-yellow) 
![Languages](https://img.shields.io/badge/Languages-R-blue)

<!-- badges: end -->

# ABMIexploreR

> Custom Predictions using ABMI species distribution models in Alberta

## Overview

## Installation

``` r
if (!require("remotes"))
    install.packages("remotes")
remotes::install_github("ABbiodiversity/ABMIexploreR")
```

The bioclimatic rasters are downloaded alongside the package, but need to be loaded into memory before users can make predictions.

``` r
library(ABMIexploreR)

abmi_load_bioclimatic()

```

## Lookup tables

### Species information

The lookup table providing information on taxonomic group, model region (e.g., vegetation and or soil), and model fit can be extracted from the abmi_species() function.

``` r
species.lookup <- abmi_species()
str(species.lookup)
# 'data.frame':	909 obs. of  21 variables:
#  $ SpeciesID           : chr  "WETO" "CATO" "BCFR" "WOFR" ...
#  $ ScientificName      : chr  "Bufo boreas" "Bufo hemiophrys" "Pseudacris maculata" "Rana sylvatica" ...
#  $ TSNID               : chr  "773513" "773521" "99005527" "775117" ...
#  $ CommonName          : chr  "Western Toad" "Canadian Toad" "Boreal Chorus Frog" "Wood Frog" ...
#  $ Rank                : chr  "Species" "Species" "Species" "Species" ...
#  $ Taxon               : chr  "Amphibians" "Amphibians" "Amphibians" "Amphibians" ...
#  $ Nonnative           : chr  NA NA NA NA ...
#  $ Occurrences         : int  117 53 636 421 163 1906 289 53 2123 100 ...
#  $ nSites              : int  117 53 636 421 102 956 194 46 948 48 ...
#  $ UseAvailabilityNorth: logi  FALSE FALSE FALSE FALSE TRUE TRUE ...
#  $ UseAvailabilitySouth: logi  FALSE FALSE FALSE FALSE TRUE TRUE ...
#  $ ModelNorth          : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#  $ ModelSouth          : logi  FALSE FALSE TRUE TRUE TRUE TRUE ...
#  $ LinkHabitat         : chr  "logit" "logit" "logit" "logit" ...
#  $ LinkSpclim          : chr  "logit" "logit" "logit" "logit" ...
#  $ AUCSpclimNorth      : num  0.862 0.917 0.764 0.722 0.808 ...
#  $ AUCVegHFNorth       : num  0.742 0.728 0.731 0.64 0.755 ...
#  $ AUCVegHFSpclimNorth : num  0.819 0.796 0.753 0.658 0.846 ...
#  $ AUCSpclimSouth      : num  NA NA NA NA 0.808 ...
#  $ AUCVegHFSouth       : num  NA NA NA NA 0.661 ...
#  $ AUCVegHFSpclimSouth : num  NA NA NA NA 0.751 ...
```

Here is the number of species by taxonomic group available in the package:

``` r
data.frame(table(species.lookup$Taxon))
#           Var1 Freq
#     Amphibians    4
#          Birds  124
#     Bryophytes  118
#        Lichens  145
#          Mites  108
# VascularPlants  434
```

### Bioclimatic information
The bioclimatic variables used in our hierarchical models can be viewed as both a lookup table, or visualized to see the spatial distribution of these variables across Alberta.

``` r
# Bioclimatic data
str(abmi_climate())
# 'data.frame':	20 obs. of  2 variables:
#  $ Column    : chr  "wS" "wN" "Easting" "Northing" ...
#  $ Definition: chr  "Weighting proportion for the soil models when both soil and vegetation models are available." "Weighting proportion for the vegetation models when both vegetation and soil models are available." "Easting coordinate (m)" "Northing coordinate (m)" ...

plot_climate(variable = "MAT")

```

### Landcover information
We can also view the lookup tables that define the landcover variables used in the two types of models available in this package (vegetation and soil). 

``` r
# Vegetation data
str(abmi_landcover(landcover = "Vegetation"))
# 'data.frame':	91 obs. of  5 variables:
#  $ Type    : chr  "Vegetation" "Vegetation" "Vegetation" "Vegetation" ...
#  $ Variable: chr  "Climate" "WhiteSpruceR" "WhiteSpruce1" "WhiteSpruce2" ...
#  $ Name    : chr  "Climate" "WhiteSpruceR" "WhiteSpruce1" "WhiteSpruce2" ...
#  $ Label   : chr  NA "White Spruce 0-9" "White Spruce 10-19" "White Spruce 20-39" ...
#  $ Color   : chr  NA "#9A9723" "#9A9723" "#9A9723" ...

# Soil data
str(abmi_landcover(landcover = "Soil"))
# 'data.frame':	24 obs. of  5 variables:
#  $ Type    : chr  "Soil" "Soil" "Soil" "Soil" ...
#  $ Variable: chr  "Climate" "paspen" "Loamy" "SandyLoam" ...
#  $ Name    : chr  "Climate" "pAspen" "Loamy" "SandyLoam" ...
#  $ Label   : chr  NA "pAspen" "Loamy" "Sandy / Loamy" ...
#  $ Color   : chr  NA "#9A9723" "#DAD157" "#663301" ...

```
## Species coefficients
The coefficients for both the vegetation and soil based models can be extracted into a data frame for single species or all species. Built in plotting functions can be called to visualize the coefficients. There are 100 bootstrapped iterations of these coefficients. However, these functions only produce the median bootstrap run that best represents the species.

``` r

# Extraction of single species coeffcients for vegetation models
single.species <- coefficient_extraction(species = "Amblystegium.serpens", model = "Vegetation")
str(single.species)
# Named num [1:90] 0.0774 0.1193 0.1885 0.2209 0.2188 ...
#  - attr(*, "names")= chr [1:99] "WhiteSpruceR" "WhiteSpruce1" "WhiteSpruce2" "WhiteSpruce3" ...

# Extraction of single species coeffcients for soil models
single.species <- coefficient_extraction(species = "Amblystegium.serpens", model = "Soil")
str(single.species)
# Named num [1:22] 0.0478 0.0414 0.0252 0.0458 0.0285 ...
# - attr(*, "names")= chr [1:23] "Loamy" "SandyLoam" "RapidDrain" "ClaySub" ...

# Extraction of all species coeffcients          
all.species <- coefficient_extraction(species = "All", model = "Vegetation")
str(all.species)
# num [1:667, 1:99] 0.02262 0.01237 0.00528 0.01527 0.07737 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : chr [1:667] "BCFR" "CATO" "WETO" "WOFR" ...
#  ..$ : chr [1:99] "WhiteSpruceR" "WhiteSpruce1" "WhiteSpruce2" "WhiteSpruce3" 

# Visualization of the vegetation models
plot_coef(species = "Amblystegium.serpens", 
          model = "Vegetation")

# Visualization of the soil models
plot_coef(species = "Amblystegium.serpens", 
          model = "Soil")


```

## Generating Species Predictions

You need to define the species ID and (`abmi_species()$SpeciesID`) the bootstrap ID (`i`). The bootstrap ID can be between 1 and 100. Users can also define a bootstrap value of 0, which selects the representative median bootstrap iteration for that species

``` r
## define species and bootstrap id
spp <- "Amblystegium.serpens"
i <- 0
```

### Example data

We provided an example data set that shows you how to organize the landcover data:

``` r
## example data to see what is needed and how it is formatted
load(system.file("extdata", "simulated-data.rda", package="ABMIexploreR"))

## Spatial Coordinates
## Formatted in EPSG: 3400
str(simulated.data$simulated.climate)
# 'data.frame':	400 obs. of  2 variables:
#  $ X: num  624879 625878 626877 627876 628875 ...
#  $ Y: num  6091950 6091950 6091950 6091950 6091950 ...

## veg+HF composition data matrix
colnames(simulated.data$simulated.vegetation)
#  [1] "DeciduousR"     "Deciduous1"     "Deciduous2"     "Deciduous3"     "Deciduous4"  
#  [6] "Deciduous5"     "Deciduous6"     "Deciduous7"     "Deciduous8"     "MixedwoodR"   
# [11] "Mixedwood1"     "Mixedwood2"     "Mixedwood3"     "Mixedwood4"     "Mixedwood5"   
# [16] "Mixedwood6"     "Mixedwood7"     "Mixedwood8"     "PineR"          "Pine1"        
# [21] "Pine2"          "Pine3"          "Pine4"          "Pine5"          "Pine6"        
# [26] "Pine7"          "Pine8"          "WhiteSpruceR"   "WhiteSpruce1"   "WhiteSpruce2" 
# [31] "WhiteSpruce3"   "WhiteSpruce4"   "WhiteSpruce5"   "WhiteSpruce6"   "WhiteSpruce7" 
# [36] "WhiteSpruce8"   "TreedBogR"      "TreedBog1"      "TreedBog2"      "TreedBog3"    
# [41] "TreedBog4"      "TreedBog5"      "TreedBog6"      "TreedBog7"      "TreedBog8"    
# [46] "TreedFenR"      "TreedFen1"      "TreedFen2"      "TreedFen3"      "TreedFen4"    
# [51] "TreedFen5"      "TreedFen6"      "TreedFen7"      "TreedFen8"      "TreedSwamp"   
# [56] "GrassHerb"      "Shrub"          "GraminoidFen"   "Marsh"          "ShrubbyBog"   
# [61] "ShrubbyFen"     "ShrubbySwamp"   "Bare"           "SnowIce"        "Water"        
# [66] "Urban"          "Rural"          "Industrial"     "Mine"           "Wellsites"    
# [71] "EnSoftLin"      "EnSeismic"      "HardLin"        "TrSoftLin"      "Crop"         
# [76] "RoughP"         "TameP"          "HWater"         "CCDeciduousR"   "CCDeciduous1" 
# [81] "CCDeciduous2"   "CCDeciduous3"   "CCDeciduous4"   "CCMixedwoodR"   "CCMixedwood1" 
# [86] "CCMixedwood2"   "CCMixedwood3"   "CCMixedwood4"   "CCPineR"        "CCPine1"      
# [91] "CCPine2"        "CCPine3"        "CCPine4"        "CCWhiteSpruceR" "CCWhiteSpruce1"
# [96] "CCWhiteSpruce2" "CCWhiteSpruce3" "CCWhiteSpruce4"

## soil+HF composition data matrix
colnames(simulated.data$simulated.soil)
#  [1] "ClaySub"     "Other"       "RapidDrain"  "Loamy"       "SandyLoam"  
#  [6] "ThinBreak"   "Blowout"     "SoilUnknown" "Water"       "Urban"      
# [11] "Rural"       "Industrial"  "Mine"        "Wellsites"   "EnSoftLin"  
# [16] "EnSeismic"   "HardLin"     "TrSoftLin"   "Crop"        "RoughP"     
# [21] "TameP"       "HWater"      "HFor" 
```

### Spatial data extraction

``` r
## Define the spatial coordinates and 
spatial.locations <- as.matrix(simulated.data$simulated.climate)
spatial.locations <- vect(x = spatial.locations[, c("X", "Y")],
                          type = "points")
crs(spatial.locations) <- "EPSG:3400"

# Extract the appropriate bioclimatic data
# Users can define unique cell names with the cell.id option
climate.input <- extract_climate(spatial.grid = spatial.locations, 
                          cell.id = 1:400,
                          reproject = FALSE)
                          
```

### Predictions

The prediction functions require composition data, i.e. giving the areas or proportions of
different landcover types (columns) in a spatial unit (rows). The corresponding relative abundance values will be returned in a matrix format:

``` r
## Define the compositional data
veg.data <- simulated.data$simulated.vegetation
soil.data <- simulated.data$simulated.soil

## Make single species prediction
model.output <- species_predict(species = spp, 
                                   veg = veg.data, 
                                   soil = soil.data, 
                                   climate = climate.input, 
                                   modified = FALSE, 
                                   boot = i)

```

### Prediction blending

For applicable species, models can be generated using both vegetation and soil based information.
These two models can be blended together in areas where the predictions have spatial overlap.

``` r
## Blend the predictions
blend.pred <- blend_predict(climate = climate.input, 
                            veg = model.output$Vegetation, 
                            soil = model.output$Soil)

```
### Prediction visualization

For applicable species, models can be generated using both vegetation and soil based information.
These two models can be blended together in areas where the predictions have spatial overlap.

``` r
## Visualize Blend the predictions
species.pred <- simulated.data$simulated.climate
species.pred$Species <- blend.pred
species.pred <- rast(x = species.pred,
                          type = "xyz")
crs(species.pred) <- "EPSG:3400"

plot_species(spat.raster = species.pred,
             variable = "Species")

```

### Prediction uncertainty

We can also generate uncertainty estimates for each species:

``` r
bootstrap.matrix <- NULL

for (i in 1:20) {
model.output <- species_predict(species = spp, 
                                   veg = veg.data, 
                                   soil = NULL, 
                                   climate = climate.input, 
                                   modified = FALSE, 
                                   boot = i)
    bootstrap.matrix <- cbind(bootstrap.matrix, model.output$Vegetation)
}

t(apply(bootstrap.matrix[1:5,], 1, quantile, c(0.5, 0.05, 0.95)))
##         50%        5%       95%
## 1 0.5139610 0.4541572 0.5672434
## 2 0.2078696 0.1587900 0.2419546
## 3 0.4180285 0.3991523 0.4348132
## 4 0.3762806 0.3444155 0.4120810
## 5 0.2520864 0.2271224 0.2793064
```

## Coefficient Modification

We recommend that users do not modify the coefficients that are available in this package. However, there are circumstances where users may wish to fix specific habitat coefficients to 0. 
If the coefficients are modified, users will need to include the `TRUE` flag for the `modified` option in the `species_predict` function. Coefficient adjustments are applied to all bootstrap iterations and across all species.

``` r

# Create the adjusted coefficients
coefficient_adjustment(model = "Vegetation", coef = "Crop", value = 0)
model.output <- species_predict(species = spp, 
                                   veg = veg.data, 
                                   soil = soil.data, 
                                   climate = climate.input, 
                                   modified = TRUE, 
                                   boot = 0)
                                   
blend.pred <- blend_predict(climate = climate.input, 
                            veg = model.output$Vegetation, 
                            soil = model.output$Soil)

# Visualize the new predictions                              
species.pred <- simulated.data$simulated.climate
species.pred$Species <- blend.pred
species.pred <- rast(x = species.pred,
                          type = "xyz")
crs(species.pred) <- "EPSG:3400"

plot_species(spat.raster = species.pred,
             variable = "Species")
             
```

## Parallel Processing

The `ABMIexploreR` package can be implemented in parallel processing. In order to do this, the `abmi_load_bioclimatic` and `extract_climate` need to be initialized for each core so climate information can be passed to the `species_predict` function. 

``` r

# Load parallel R package
library(parallel)

# Define the species list (lets use all of the bird models)
species.list <- abmi_species()
species.list <- species.list[species.list$Taxon == "Birds", "SpeciesID"]
names(species.list) <- species.list

# Define the number of bootstraps
boot.iter <- 1:10

# Initialize cores
n.clusters <- 10
core.input <- makeCluster(n.clusters)
clusterExport(core.input, c())
clusterEvalQ(core.input, {
  
  library(ABMIexploreR)
  library(Matrix)
  library(terra)
  
  # Load the landcover data
  load(system.file("extdata", "simulated-data.rda", package="ABMIexploreR"))
  veg.data <- simulated.data$simulated.vegetation
  soil.data <- simulated.data$simulated.soil
    
  # Load bioclimatic data
  abmi_load_bioclimatic()
  
  # Generate the spatial points
  spatial.locations <- as.matrix(simulated.data$simulated.climate)
  spatial.locations <- vect(x = spatial.locations[, c("X", "Y")],
                            type = "points")
  crs(spatial.locations) <- "EPSG:3400"


  # Generate the matching spatial grid 
  climate.grid <- extract_climate(spatial.grid = spatial.locations, 
                            cell.id = rownames(veg.data), 
                            reproject = FALSE)
  
})

# Generate median bootstrap prediction across multiple species
species.predictions <- parLapply(cl = core.input, 
                                 X = species.list, 
                                 function(spp) species_predict(species = spp, 
                                                               veg = veg.data, 
                                                               soil = soil.data, 
                                                               climate = climate.grid, 
                                                               modified = FALSE, 
                                                               boot = 0))
                                                               
str(species.predictions)
## List of 124
##  $ AlderFlycatcher            :List of 2
##   ..$ Vegetation: Named num [1:400] 0.0742 0.0056 0.0627 0.0997 0.0523 ...
##   .. ..- attr(*, "names")= chr [1:400] "1" "2" "3" "4" ...
##   ..$ Soil      : Named num [1:304] 0.0477 0.0413 0.0199 0.0369 0.0525 ...
##   .. ..- attr(*, "names")= chr [1:304] "1" "2" "21" "22" ...                                                             
# Or the 100 bootstrap predictions across a single species
species.predictions <- parLapply(cl = core.input, 
                                 X = boot.iter, 
                                 function(boot) species_predict(species = "Amblystegium.serpens", 
                                                               veg = veg.data, 
                                                               soil = soil.data, 
                                                               climate = climate.grid, 
                                                               modified = FALSE, 
                                                               boot = boot))
str(species.predictions)
## List of 10
##  $ :List of 2
##   ..$ Vegetation: Named num [1:400] 0.517 0.198 0.417 0.376 0.252 ...
##   .. ..- attr(*, "names")= chr [1:400] "1" "2" "3" "4" ...
##   ..$ Soil      : Named num [1:304] 0.581 0.523 0.28 0.469 0.547 ...
##   .. ..- attr(*, "names")= chr [1:304] "1" "2" "21" "22" ...                                                              
stopCluster(core.input)
            
```