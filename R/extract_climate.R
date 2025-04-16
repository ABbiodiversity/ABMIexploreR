#' Extract climate
#'
#' @description Function for extracting the relevant bioclimatic variables
#'
#' @param spatial.grid Terra spatial vector of points that the user wants to extract bioclimatic information for 
#' @param cell.id Vector of cell ID used to link the bioclimatic and landcover data.
#' @param reproject Logical; If TRUE, will call terra::project to reproject to the appropriate CRS (NAD83 / Alberta 10-TM (Forest) (EPSG:3400))
#'
#' @import terra
#'
#' @export
#'

extract_climate <- function(spatial.grid, 
                            cell.id = NULL,
                            reproject = FALSE) {
    
    ##########################
    # Check for matching CRS #
    ##########################
    
    if(reproject == TRUE) {
        
        spatial.grid <- terra::project(spatial.grid, 
                                       crs(.bioclimatic.raster))
        
    }
    
    ###############################
    # Extract climate information #
    ###############################
    
    bioclimatic.data <- terra::extract(.bioclimatic.raster, 
                                      spatial.grid, 
                                      ID = FALSE) 
    
    if(nrow(bioclimatic.data) != length(cell.id)) {
        
        stop("Length of cell.id does not match the number of rows in spatial.grid.")
        
    } 
    
    # Assign row names
    rownames(bioclimatic.data) <- cell.id
    
    return(bioclimatic.data)
    
}
