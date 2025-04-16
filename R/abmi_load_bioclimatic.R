#' Bioclimatic data
#'
#' Load the bioclimatic raster into memory

#' @param dir Directory path.
#'
#' @return
#'
#' `abmi_load_bioclimatic Loads the bioclimatic data into the current session.
#'
#' @name bioclimatic
NULL


#' @export
#' @rdname bioclimatic
abmi_load_bioclimatic <- function(dir=NULL) {
    
    fn <- system.file("extdata", "bioclimatic-data.tif", package="ABMIexploreR")
    ABMIexploreR:::.msg("Loading bioclimatic data")
    bioclimatic.raster <- try(terra::rast(fn))
    assign(x = ".bioclimatic.raster", 
           value = bioclimatic.raster, 
           envir = .GlobalEnv)
    
}
