# Fill Unknown Ages

There are often parts of forested landcover datasets where the age of
the forest is unknown. When this occurs, these areas must be assigned an
age. This function assigns ages to these areas based on the average
forest age distributions of Natural Subregions in Alberta. These
distributions have been shifted such that the minimum age is 80 years,
assuming that the areas of unknown age have not been harvested or
burned, and therefore are likely to be at least 80 years old.

## Usage

``` r
fill_unknown_ages(land.cover, age.distribution = "age.old.nsr.maltman")
```

## Arguments

- land.cover:

  The dataframe with unknown ages that needs to be updated.

- age.distribution:

  The data file that contains average age distributions by natural
  subregion.

## Value

A data frame with missing ages filled in.

## Details

Note that the dataframe should have all the necessary columns, in order.
If the columns are not found, an error will be returned. If they are out
of order,the function will try to put them back in the correct order.
The correct order for each veg type follows this pattern: PineR, Pine1,
..., Pine8.
