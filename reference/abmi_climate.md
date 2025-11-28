# Provide users with bioclimatic lookup table

This function provides the lookup table of the bioclimatic covariates
used in the climate model.

## Usage

``` r
abmi_climate()
```

## Value

A dataframe with all of the coefficients used in the climate model. Each
taxonomic group uses a unique set of these variables.

## Examples

``` r
if (FALSE) { # \dontrun{
climate.lookup <- abmi_climate()
} # }
```
