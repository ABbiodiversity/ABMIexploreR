# Extract climate

Function for extracting the relevant bioclimatic variables

## Usage

``` r
extract_climate(spatial.grid, cell.id = NULL, reproject = FALSE)
```

## Arguments

- spatial.grid:

  Terra spatial vector of points that the user wants to extract
  bioclimatic information for

- cell.id:

  Vector of cell ID used to link the bioclimatic and landcover data.

- reproject:

  Logical; If TRUE, will call terra::project to reproject to the
  appropriate CRS (NAD83 / Alberta 10-TM (Forest) (EPSG:3400))
