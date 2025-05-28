#' Truncate prediction
#'
#' @description Function for truncating the predictions from the species_predicts and blend_predict function. Users can use truncate based on the current landscape predictions (current) or truncate two time periods using the same threshold. 
#'
#' @param species Unique Species ID.
#' @param current Blended or single model type predictions (e.g., vegetation or soil) based on current landscape conditions
#' @param reference Blended or single model type predictions (e.g., vegetation or soil) based on reference landscape conditions
#'
#' @export
#'

truncate_predict <- function(species, current = NULL, reference = NULL) {
    
    # Call the species lookup table function
    species.lookup <- abmi_species()
    
    # Isolate the relevant prediction threshold 
    threshold <- species.lookup[species.lookup$SpeciesID %in% species, "Threshold"]
    
    if(!is.null(current) & !is.null(reference)) {
        
        # Calculate the value based on the threshold considering both current and reference landscapes
        prediction.threshold <- quantile(c(current, reference), threshold)
        
        current <- ifelse(current >= prediction.threshold,
                          prediction.threshold,
                          current)
        
        reference <- ifelse(reference >= prediction.threshold,
                            prediction.threshold,
                            reference)
        
        
    } else {
        
        # Calculate the value based on the threshold considering only the current landscape
        prediction.threshold <- quantile(current, threshold)
        current <- ifelse(current >= prediction.threshold,
                          prediction.threshold,
                          current)
        
    }
    
    # Return the values
    return(list(current = current,
                reference = reference))
}
