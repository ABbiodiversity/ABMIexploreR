# Mammal predictions

Function for calculating the species level predictions for mammals

## Usage

``` r
.mammal_predicts(
  species,
  taxon = NULL,
  veg = NULL,
  soil = NULL,
  climate = NULL,
  boot = 0,
  modified = FALSE
)
```

## Arguments

- species:

  Unique Species ID.

- taxon:

  Taxon associated with the unique Species ID.

- veg:

  Matrix of standardized vegetation information (row sum equals 1). If
  NULL and model is available, no prediction is generated.

- soil:

  Matrix of standardized soil information (row sum equals 1). If NULL
  and model is available, no prediction is generated.

- climate:

  Matrix of the bioclimatic variables that is matched to the vegetation
  and/or soil matrices.

- boot:

  Defines bootstrap iteration to generate the predictions from.
