#' Mammal predictions
#'
#' @description Function for calculating the species level predictions for mammals
#'
#' @param species Unique Species ID. 
#' @param taxon Taxon associated with the unique Species ID.
#' @param veg Matrix of standardized vegetation information (row sum equals 1). If NULL and model is available, no prediction is generated.
#' @param soil Matrix of standardized soil information (row sum equals 1). If NULL and model is available, no prediction is generated.
#' @param climate Matrix of the bioclimatic variables that is matched to the vegetation and/or soil matrices.
#' @param boot Defines bootstrap iteration to generate the predictions from.
#'
#' @import Matrix
#'
#' @export
#'

.mammal_predicts <- function(species, taxon = NULL, veg = NULL, 
                                  soil = NULL, climate = NULL, boot = 0,
                                  modified = FALSE) {
    
    # If modified coefficients are defined by the user, use them
    if(modified) {
        
        COEFS <- .species.coefs.mod
        
    } else {
        
        COEFS <- species.coefs
        
    }
    
    # # If boot equals zero, pull the median bootstrap value
    # # Currently we only have a single model. Therefore, we fix boot to 1 regardless of the user defined option
    # if (boot == 0) {
    #     
    #     # Call the species lookup table function
    #     species.lookup <- abmi_species()
    #     
    #     # Isolate the relevant bootstrap iteration
    #     boot <- species.lookup[species.lookup$SpeciesID %in% species, "Bootstrap"]
    #     
    # } 
    
    boot <- 1
    
    # Define vegetation and soil coefficients
    if(species %in% rownames(COEFS[[taxon]]$Climate)) {
        
        species.climate <- COEFS[[taxon]]$Climate[species,,boot]
        
    } else {
        
        species.climate <- NULL
        
    }
    
    if(species %in% rownames(COEFS[[taxon]]$Vegetation$PA)) {
        
        species.veg <- COEFS[[taxon]]$Vegetation$PA[species,,boot]
        species.veg.agp <- COEFS[[taxon]]$Vegetation$AGP[species,,boot]
        species.veg.ta <- COEFS[[taxon]]$Vegetation$TA[species,,boot]
        
    } else {
        
        species.veg <- NULL
        species.veg.agp <- NULL
        species.veg.ta <- NULL
        
    }
    
    if(species %in% rownames(COEFS[[taxon]]$Soil$PA)) {
        
        species.soil <- COEFS[[taxon]]$Soil$PA[species,,boot]
        species.soil.agp <- COEFS[[taxon]]$Soil$AGP[species,,boot]
        species.soil.ta <- COEFS[[taxon]]$Soil$TA[species,,boot]
        
    } else {
        
        species.soil <- NULL
        species.soil.agp <- NULL
        species.soil.ta <- NULL
        
    }
    
    # Define link and inverse link functions
    # These are separate for the total abundance and climate
    inv_link_pa <- inv_link_function[[taxon]]$Binomial
    link_pa <- link_function[[taxon]]$Binomial
    
    inv_link_abundance <- inv_link_function[[taxon]]$Gamma
    link_abundance <- link_function[[taxon]]$Gamma
    
    # Define blank object to be returned to user
    veg.pred <- NULL
    soil.pred <- NULL
    
    #######################################
    # Create the global bioclimatic model #
    #######################################
    
    # Standardize the climate data
    climate <- climate$Climate
    climate.global <- as.matrix(climate[, names(species.climate)])
    climate.coef <- species.climate[colnames(climate.global)]
    
    # Predict space/climate component
    climate.global <- matrix(inv_link_pa(drop(climate.global %*% climate.coef)), ncol = 1,
                           dimnames = list(rownames(climate.global), "Climate"))
    
    # Truncate climate prediction
    climate.global <- ifelse(climate.global >= quantile(climate.global, 0.99),
                           quantile(climate.global, 0.99),
                           climate.global)
    
    #################################################
    # Perform the separate veg and soil predictions #
    #################################################
    
    # Identify if a vegetation prediction can be made
    if(!is.null(veg) & !is.null(species.veg)) {
        
        # Subset the climate data to the region of interest
        climate.pred <- climate.global[rownames(veg), ]
        
        # Use this to predict the joint climate contribution
        # Presence Absence
        climate.pred <- (climate.pred * species.veg["Climate"])
        
        # Using these prediction, create a matrix and get the climate adjusted veg coefficients
        climate.matrix <- matrix(climate.pred, nrow = nrow(veg), ncol = ncol(veg))
        
        # Standardize the vegetation coefficients
        # PA
        veg.coef <- species.veg[colnames(veg)]
        
        # AGP is derived from Total Adundance / Presence Absence
        # This is because the Total Abundance and Presence Absence coefficients
        # have been standardized using a Observed/Predicted.
        # The stored AGP coefficients have not been adjusted but are required to
        # isolate the pAspen coefficient in the south model.
        veg.coef.agp <- inv_link_abundance(species.veg.ta[names(veg.coef)]) / inv_link_pa(veg.coef)
        
        # Multiple AGP by PA to calculate the new climate modified joint coefficients
        veg.coef <- t(inv_link_pa(t(climate.matrix) + veg.coef) * veg.coef.agp)
        
        # Prediction
        veg.pred <- rowSums(veg * veg.coef)
        
    }
    
    # Identify if soil prediction can be made
    if(!is.null(soil) & !is.null(species.soil)) {
        
        # Calculate the pAspen component
        paspen.global <- matrix(climate[, "pAspen"], ncol = 1,
                                dimnames = list(rownames(climate), "pAspen"))
        
        # Subset the climate and pAspen data to the region of interest
        climate.pred <- climate.global[rownames(soil), ]
        paspen.global <- paspen.global[rownames(soil), ]
        
        # Use these to predict the joint climate contribution
        climate.pred <- (climate.pred * species.soil["Climate"])
        
        # PA pAspen
        paspen.pa.pred <- (paspen.global * species.soil["pAspen"])
        
        # AGP pAspen
        paspen.agp.pred <- inv_link_abundance(paspen.global * species.soil.agp["pAspen"])
        
        # Create the joint PA climate and paspen
        climate.matrix <- matrix(climate.pred + paspen.pa.pred, nrow = nrow(soil), ncol = ncol(soil))
        
        # Create the AGP paspen matrix
        paspen.agp.matrix <- matrix(paspen.agp.pred, nrow = nrow(soil), ncol = ncol(soil))
        
        # Standardize the soil coefficients
        # PA
        soil.coef <- species.soil[colnames(soil)]
        
        # AGP is derived from Total Adundance / Presence Absence
        # This is because the Total Abundance and Presence Absence coefficients
        # have been standardized using a Observed/Predicted.
        # The stored AGP coefficients have not been adjusted but are required to
        # isolate the pAspen coefficient in the south model.
        soil.coef.agp <- link_abundance(inv_link_abundance(species.soil.ta[names(soil.coef)]) / inv_link_pa(soil.coef))
        
        # Multiple AGP by PA to calculate the new climate modified joint coefficients
        # We add the agp.paspen value to modified agp coefficients otherwise
        # if agp paspen == 0 we get 0
        soil.coef <- t(inv_link_pa(t(climate.matrix) + soil.coef) * 
                           inv_link_abundance(t(paspen.agp.matrix) * soil.coef.agp))
        
        # Prediction
        soil.pred <- rowSums(soil * soil.coef)
        
        # Garbage collect
        gc()
        
    }
    
    # Return the predictions to user
    return(list(Vegetation = veg.pred, 
                Soil = soil.pred))
    
}
