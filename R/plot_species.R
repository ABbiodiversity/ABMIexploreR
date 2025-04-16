#' Plot species
#'
#' @description Function for visualizing the predictions from species_predicts
#'
#' @param spat.raster SpatRaster object of the predictions.
#' @param variable Attribute within the SpatRaster object for plotting.
#'
#' @import ggplot2
#' @import MetBrewer
#' @import terra
#' @import tidyterra
#'
#' @export
#'

plot_species <- function(spat.raster = NULL,
                         variable = NULL) {
    
    ####################################
    # Check for matching variable name #
    ####################################
    
    if(variable %in% names(spat.raster)) {
        
        plot.species <- ggplot() +
            geom_spatraster(data = spat.raster, aes_string(fill = variable),
                            maxcell = length(values(spat.raster))) +
            scale_fill_gradientn(name = paste0(variable), 
                                 colors = rev(met.brewer(name = "Hiroshige", n = 100, type = "continuous")), 
                                 guide = "colourbar",
                                 na.value = NA) +
            theme_light() +
            theme(axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size=18),
                  axis.text.y = element_text(size=18),
                  panel.grid.major.y = element_blank(),
                  legend.text = element_text(size=14),
                  legend.title = element_text(size=16),
                  legend.key.size = unit(1, "cm"),
                  axis.line = element_line(colour = "black"),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))
        
        return(plot.species)
        
    } else {
        
        return(paste0("Provide variable (", variable, ") not found in the spat.raster object."))
        
    }
    
}
