# Truncate prediction

Function for truncating the predictions from the species_predicts and
blend_predict function. Users can use truncate based on the current
landscape predictions (current) or truncate two time periods using the
same threshold.

## Usage

``` r
truncate_predict(species, current = NULL, reference = NULL)
```

## Arguments

- species:

  Unique Species ID.

- current:

  Blended or single model type predictions (e.g., vegetation or soil)
  based on current landscape conditions

- reference:

  Blended or single model type predictions (e.g., vegetation or soil)
  based on reference landscape conditions
