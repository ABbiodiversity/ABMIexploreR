#' Plot species coefficient
#'
#' @description Function for visualizing the species coefficients
#'
#' @param species Unique Species ID. 
#' @param model Vector defined as "vegetation" or "soil" that provides the relevant coefficients.
#'
#' @import ggplot2
#' @import MetBrewer
#'
#' @export
#'

plot_coef <- function(species,
                      model = NULL) {
    
    # Call the coefficient extraction figure to isolate the species and model of interest
    model.coef <- coefficient_extraction(species = species,
                                         model = model)
    
    # Vegetation Plot
    if(model == "Vegetation") {
        
        # Isolate the vegetation lookup table
        vegetation.template <- landcover.names$Vegetation
        vegetation.template <- vegetation.template[!is.na(vegetation.template$Label), ]
        
        # Isolate the main coefficients
        main.coef <- data.frame(Name = names(model.coef),
                                Coef = model.coef,
                                Label = factor(names(model.coef),
                                               levels = names(model.coef)))
        
        main.coef$Color <- vegetation.template$Color[match(main.coef$Name, vegetation.template$Name)]
        
        # Define maximum value
        max.value <- max(c(main.coef$Coef))
        
        # Define the harvest area coefficients and label
        cc.coef <- main.coef[grep("CC", main.coef$Name), ]
        cc.coef$Name <- gsub("CC", "", cc.coef$Name)
        cc.coef$Label <- factor(cc.coef$Label, levels = cc.coef$Label)  
        main.coef <- main.coef[!is.na(main.coef$Color), ]
        main.coef <- droplevels(main.coef)
        
        cc.coef <- data.frame(Class = c(rep("WhiteSpruce", 5), rep(NA, 4),
                                        rep("Pine", 5), rep(NA, 4),
                                        rep("Deciduous", 5), rep(NA, 4),
                                        rep("Mixedwood", 5), rep(NA, 4),
                                        rep(NA, 29)),
                              Coef = c(cc.coef$Coef[1:5], rep(NA, 4),
                                      cc.coef$Coef[6:10], rep(NA, 4),
                                      cc.coef$Coef[11:15], rep(NA, 4),
                                      cc.coef$Coef[16:20], rep(NA, 4),
                                      rep(NA, 29)))
        
        # Create the plot
        coef.plot <- print(ggplot(data = main.coef, aes(x = Label, y = Coef, fill = Color)) +
            geom_bar(stat = "identity", fill = main.coef$Color) +
            scale_x_discrete(labels = main.coef$Label) +
            geom_point(aes(x = main.coef$Name, y = cc.coef$Coef), show.legend = FALSE) +
            geom_line(aes(x = main.coef$Name, y = cc.coef$Coef, 
                          group = cc.coef$Class, linetype = "dotted"), 
                      size = 1, show.legend = FALSE) +
            scale_linetype_manual(values=c("dotted")) +
            guides(scale = "none") + 
            labs(x = "Coefficient", y = "Relative Abundance") +
            ggtitle(species) +
            theme_light() +
            coord_cartesian(ylim = c(0, max.value), clip = "off") +
            theme(axis.title = element_text(size=16),
                  axis.text.x = element_text(size=16, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size=16),
                  title = element_text(size=12),
                  legend.text = element_text(size=16),
                  legend.title = element_blank(),
                  axis.line = element_line(colour = "black"),
                  panel.border = element_rect(colour = "black", fill=NA, size=1)))
        
    }
    
    # Soil Plot
    if(model == "Soil") {
        
        # Isolate the soil lookup table
        soil.template <- landcover.names$Soil
        soil.template <- soil.template[!is.na(soil.template$Label), ]
        
        # Standardize the coefficients
        main.coef <- data.frame(Name = names(model.coef),
                                Coef = model.coef,
                                Label = factor(names(model.coef),
                                               levels = names(model.coef)))
        
        main.coef$Color <- soil.template$Color[match(main.coef$Name, soil.template$Name)]
        main.coef <- main.coef[!is.na(main.coef$Color), ]
        
        max.value <- max(c(main.coef$Coef))
        
        # Create the plot
        coef.plot <- print(ggplot(data = main.coef, aes(x = Label, y = Coef, fill = Color)) +
            geom_bar(stat = "identity", fill = main.coef$Color) +
            scale_x_discrete(labels = main.coef$Label) +
            guides(scale = "none") + 
            labs(x = "Coefficient", y = "Relative Abundance") +
            ggtitle(species) +
            theme_light() +
            coord_cartesian(ylim = c(0, max.value), clip = "off") +
            theme(axis.title = element_text(size=16),
                  axis.text.x = element_text(size=16, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size=16),
                  title = element_text(size=12),
                  legend.text = element_text(size=16),
                  legend.title = element_blank(),
                  axis.line = element_line(colour = "black"),
                  panel.border = element_rect(colour = "black", fill=NA, size=1)))
        
    }
    
    return(coef.plot)
    
}
