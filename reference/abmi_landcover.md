# Internal lookup tables

This function provides the lookup table of the vegetation, soil, and
footprint definitions

## Usage

``` r
abmi_landcover(landcover)
```

## Arguments

- landcover:

  Define 'Vegetation' or 'Soil' landcover classes.

## Value

A dataframe with all of the coefficients used in the joint climate and
landcover model.

## Provide users with access to internal lookup tables.

NA

## Examples

``` r
if (FALSE) { # \dontrun{
landcover.lookup <- abmi_landcover(landcover = "Vegetation")
} # }
```
