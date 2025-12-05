#' Fill Unknown Ages
#'
#' There are often parts of forested landcover datasets where the age of the
#' forest is unknown. When this occurs, these areas must be assigned an age.
#' This function assigns ages to these areas based on the average forest age
#' distributions of Natural Subregions in Alberta. These distributions have been
#' shifted such that the minimum age is 80 years, assuming that the areas of
#' unknown age have not been harvested or burned, and therefore are likely to be
#' at least 80 years old.
#' 
#' Note that the dataframe should have all the necessary columns, in order. If
#' the columns are not found, an error will be returned. If they are out of 
#' order,the function will try to put them back in the correct order. The
#' correct order for each veg type follows this pattern: 
#' PineR, Pine1, ..., Pine8.
#'
#' @param land.cover The dataframe with unknown ages that needs to be updated.
#' @param age.distribution The data file that contains average age distributions by natural subregion.
#'
#' @return A data frame with missing ages filled in.
#'
#' @import data.table
#'
#' @export
#'
fill_unknown_ages <-  function(land.cover, age.distribution="age.old.nsr.maltman") {
    veg.types <- c("Pine", "WhiteSpruce", "Mixedwood", "Deciduous", "TreedBog")
    
    # Load age distributions
    load('data/ages-by-nsr.rda')
    all.distributions <- c('age.all.abmi', "age.all.maltman", "age.nsr.abmi", "age.nsr.maltman", "age.old.all.abmi", "age.old.all.maltman", "age.old.nsr.abmi", "age.old.nsr.maltman", "ages.list", "nsr.age.proportions")
    
    age.dist <- eval(parse(text=age.distribution))
    
    nsr.age.distributions <- age.dist$reference 
    
    # Type checking
    
    if (!("data.frame" %in% class(land.cover) || 
          "data.table" %in% class(land.cover))) {
        stop("land.cover must be a data.frame or data.table object.")
    }
    if (!any(class(land.cover) == "data.table") ) {
        land.cover <- data.table::data.table(land.cover)
    }
    
    # Check for expected columns
    veg.col.names <- ABMIexploreR::landcover.names$Vegetation$Name
    all.expected.cols <- unique(c("NSRNAME", veg.col.names))
    all.expected.cols <- all.expected.cols[!(all.expected.cols %in% c("Climate", "pAspen", "SoilUnknown", "HWater"))]
    
    missing.cols <- setdiff(all.expected.cols, names(land.cover))
    
    if (length(missing.cols) > 0) {
        stop(paste0("Expected columns are missing: ", 
                    paste(missing.cols, collapse = ", "),
                    ". Please check the dataframe and ensure all columns are present."))
    }
    
    
    # Ensure columns are in expected order
    extras <- setdiff(names(land.cover), all.expected.cols)
    data.table::setcolorder(land.cover, c(all.expected.cols, extras))
    
    # Clean the age distributions to ensure the names of the columns match
    message("Updating age distributions...")
    
    # Replace "Decid" with "Deciduous" in column names
    if (any("Decid" == names(nsr.age.distributions))){
        names(nsr.age.distributions) <- gsub("Decid", "Deciduous", 
                                             names(nsr.age.distributions))
    }
    if (any("DecidR" == names(nsr.age.distributions$Deciduous))){
        names(nsr.age.distributions$Deciduous) <- 
            gsub("Decid", "Deciduous", 
                 names(nsr.age.distributions$Deciduous))
    }
    
    # Replace "Spruce" with "WhiteSpruce"
    if (any("Spruce" == names(nsr.age.distributions))){
        names(nsr.age.distributions) <- gsub("Spruce", "WhiteSpruce", 
                                             names(nsr.age.distributions))
    }
    
    if (any("SpruceR" == names(nsr.age.distributions$WhiteSpruce))){
        names(nsr.age.distributions$WhiteSpruce) <- 
            gsub("Spruce", "WhiteSpruce", 
                 names(nsr.age.distributions$WhiteSpruce))
    }
    
    message("Reclassifying areas of unknown age...")
    
    for (vt in veg.types) {
        if (any(paste0(vt, "U") == colnames(land.cover))) {
            message(paste0("Updating ages for ", vt, "..."))
            
            # Get column indices
            veg.type.first.col <- match(paste0(vt, "R"), 
                                        colnames(land.cover))
            veg.type.last.col <- match(paste0(vt, "8"), 
                                       colnames(land.cover))
            nsr.name.col <- match("NSRNAME", colnames(land.cover))
            veg.type.unknown.age.col = match(paste0(vt, "U"), 
                                             colnames(land.cover))
            
            # Check if there are any unknown areas
            if (sum(land.cover[[veg.type.unknown.age.col]]) > 0) {
                # Get a subset of the data table
                cols <- c(nsr.name.col, veg.type.unknown.age.col)
                unknown.veg <- land.cover[get(paste0(vt, "U")) > 0, 
                                          ..cols]
                names(unknown.veg) <- c("NSRNAME", "U")
                
                # Grab the assumed distribution
                raw.nsr.distribution <- data.table::as.data.table(
                    nsr.age.distributions[[vt]], keep.rownames="NSRNAME"
                )
                #print(head(raw.nsr.distribution))
                
                # Join the distributions to the cells by NSRNAME
                nsr.distribution.by.cell <- 
                    data.table::merge.data.table(unknown.veg[,1], 
                                                 raw.nsr.distribution, 
                                                 by="NSRNAME", all.x=TRUE)
                
                # Separate the text and numbers
                nsr.names <- nsr.distribution.by.cell[, 1, with = FALSE]
                distributions <- nsr.distribution.by.cell[, -c(1, 2, 12), 
                                                          with = FALSE]
                
                # Get estimated areas
                distributions.by.cell <- distributions * unknown.veg$U
                
                # Update the dataframe with the original + new area
                column.names <- names(distributions)
                rows.to.update <- land.cover[[paste0(vt, "U")]] > 0
                for (col.name in column.names) {
                    land.cover[rows.to.update, (col.name) := 
                                   get(col.name) + 
                                   distributions.by.cell[[col.name]]]
                }
                
                # Set original unknown areas to zero
                land.cover[, (paste0(vt, "U")) := 0]
                
                
            } else {
                print(paste0(vt, " does not need updating."))
            }
        }
    }
    return(as.data.frame(land.cover))
}
