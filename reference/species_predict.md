# Species predictions

Wrapper function for calculating the species level predictions

## Usage

``` r
species_predict(
  species,
  veg = NULL,
  soil = NULL,
  climate = NULL,
  modified = FALSE,
  boot = 0
)
```

## Arguments

- species:

  Unique Species ID.

- veg:

  Matrix of standardized vegetation information (row sum equals 1). If
  NULL and model is available, no prediction is generated.

- soil:

  Matrix of standardized soil information (row sum equals 1). If NULL
  and model is available, no prediction is generated.

- climate:

  Matrix of the bioclimatic variables that is matched to the vegetation
  and/or soil matrices.

- modified:

  Logical; Default behaviour (FALSE) means the original species
  coefficients are used. TRUE indicators the function should use the
  coefficients created through the coefficient_adjustment function.

- boot:

  Defines bootstrap iteration to generate the predictions from. If 0,
  uses the median bootstrap run unique for each species.
