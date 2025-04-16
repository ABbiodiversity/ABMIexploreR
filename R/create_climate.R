#' Create bioclimatic data
#'
#' @description Function for generating the taxon specific bioclimatic data.
#'
#' @param climate.grid Matrix containing the stored bioclimatic variables. Rownames are used as the cell identifier.
#' @param taxon Taxonomic group associated with the focal species of interest. 
#'
#' @export
#'

.create_climate <- function(climate.grid, taxon) {
    
    ##########################
    # Check for valid inputs #
    ##########################
    
    if(!(taxon %in% c("Amphibians", "Birds", "Bryophytes",
                      "Lichens", "Mammals", "Mites", "VascularPlants"))) {
        
        stop("Need to have a valid taxonomic group (Amphibians, Birds, Bryophytes, 
         Lichens, Mammals, Mites, VascularPlants).")
        
    }
    
    # Determine which pixels fall within the vegetation and soil model regions
    UseN <- rownames(climate.grid)[climate.grid$wN > 0]
    UseS <- rownames(climate.grid)[climate.grid$wS > 0]
    
    # Create unique climate grids 
    if (taxon == "Mammals") {
        
        # Need to write
        
    } 
    
    if (taxon != "Mammals") {
        
        climate.summary <- with(climate.grid, cbind(
            Intercept = 1,
            Easting = Easting,
            Northing = Northing,
            Easting2 = Easting2,
            Northing2 = Northing2,
            EastingNorthing = EastingNorthing,
            pAspen = pAspen,
            MAT = MAT,
            FFP = FFP,
            MWMT = MWMT,
            EMT = EMT,
            MAP = MAP,
            PET = PET,
            CMD = CMD,
            TD = TD,
            MAPPET = MAPPET,
            CMDMAT = CMDMAT,
            MAT2 = MAT2,
            MWMT2 = MWMT2,
            bio9 = bio9,
            bio15 = bio15))
        
        rownames(climate.summary) <- rownames(climate.grid)
        climate.summary <- list(vegetation = climate.summary[UseN, ],
                                soil = climate.summary[UseS, ])
        
    } 
    
    gc()
    
    return(climate.summary)
    
}
