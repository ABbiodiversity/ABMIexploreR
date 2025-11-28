# Coefficient Adjustment

This function will adjust the stored coefficients to user defined
values. In general, we advice not adjusting these coefficients. However,
there are specific situations where adjusting stored coefficients to 0
(e.g., Crop, Urban Industrial, Hard Linear) may be appropriate.

## Usage

``` r
coefficient_adjustment(model, coef = NULL, value = NULL)
```

## Arguments

- model:

  Vector defined as "Vegetation", "Soil", or "All" so the appropriate
  model is adjusted.

- coef:

  Vector of characters that the user wishes to modify.

- value:

  User defined value to be assigned to the specified coefficients.
  Adjustment is made across all bootstraps.
