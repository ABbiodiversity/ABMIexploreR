
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
#' @param land_cover_dataframe The dataframe with unknown ages that needs to be
#' updated.
#' @param nsr_age_distributions The data file that contains average age
#' distributions by natural subregion.
#'
#' @return A data frame with missing ages filled in.
#'
#' @import data.table
#'
#' @export
#'
fill_unknown_ages <-  function(land_cover_dataframe, age_distribution="age.old.nsr.maltman") {
    veg_types <- c("Pine", "WhiteSpruce", "Mixedwood", "Deciduous", "TreedBog")
    
    # Load age distributions
    load('data/ages-by-nsr.rda')
    all_distributions <- c('age.all.abmi', "age.all.maltman", "age.nsr.abmi", "age.nsr.maltman", "age.old.all.abmi", "age.old.all.maltman", "age.old.nsr.abmi", "age.old.nsr.maltman", "ages.list", "nsr.age.proportions")
    to_remove <- all_distributions[!(all_distributions == age_distribution)]
    rm(list=to_remove, envir=.GlobalEnv)
    
    age_dist <- eval(parse(text=age_distribution))
    
    nsr_age_distributions <- age_dist$reference 
    
    # Type checking
    
    if (!("data.frame" %in% class(land_cover_dataframe) || 
          "data.table" %in% class(land_cover_dataframe))) {
        stop("land_cover_dataframe must be a data.frame or data.table object.")
    }
    
    # Check for expected columns
    missing_cols <- c()
    for (vt in veg_types){
        message(paste0("Checking for ", vt, " columns."))
        expected_cols <- paste0(vt, c("R", 1:8))
        missing <- setdiff(expected_cols, names(land_cover_dataframe))
        missing_cols <- c(missing_cols, missing)
    }
    
    if (length(missing_cols) > 0) {
        stop(paste0("Expected columns are missing: ", 
                    paste(missing_cols, collapse = ", "),
                    ". Please check the dataframe and ensure all columns are present."))
    }
    
    if (!(any("NSRNAME" == names(land_cover_dataframe)))){
        stop("NSRNAME not found in dataframe, please add natural subregion data for each observation.")
    }
    
    # Ensure columns are in expected order
    for (vt in veg_types) {
        # Expected columns for this veg type
        expected_cols <- paste0(vt, c("R", 1:8))
        # Find current positions
        current_positions <- match(expected_cols, names(land_cover_dataframe))
        # If not in order or not together, reorder
        if (!all(c(diff(current_positions) == 1, !is.na(current_positions)))) {
            message(paste0("Re-ordering ", vt, " columns."))
            # Move columns to correct positions
            other_cols <- setdiff(names(land_cover_dataframe), expected_cols)
            land_cover_dataframe <- land_cover_dataframe[, c(other_cols, 
                                                             expected_cols)]
        }
    }
    
    # Clean the age distributions to ensure the names of the columns match
    message("Updating age distributions...")
    nsr_age_distributions <- nsr_age_distributions$reference
    
    # Replace "Decid" with "Deciduous" in column names
    if (any("Decid" == names(nsr_age_distributions))){
        names(nsr_age_distributions) <- gsub("Decid", "Deciduous", 
                                             names(nsr_age_distributions))
    }
    if (any("DecidR" == names(nsr_age_distributions$Deciduous))){
        names(nsr_age_distributions$Deciduous) <- 
            gsub("Decid", "Deciduous", 
                 names(nsr_age_distributions$Deciduous))
    }
    
    # Replace "Spruce" with "WhiteSpruce"
    if (any("Spruce" == names(nsr_age_distributions))){
        names(nsr_age_distributions) <- gsub("Spruce", "WhiteSpruce", 
                                             names(nsr_age_distributions))
    }
    
    if (any("SpruceR" == names(nsr_age_distributions$WhiteSpruce))){
        names(nsr_age_distributions$WhiteSpruce) <- 
            gsub("Spruce", "WhiteSpruce", 
                 names(nsr_age_distributions$WhiteSpruce))
    }
    
    message("Reclassifying areas of unknown age...")
    if (!any(class(land_cover_dataframe) == "data.table") ) {
        land_cover_dataframe <- data.table::data.table(land_cover_dataframe)
    }
    for (vt in veg_types) {
        if (any(paste0(vt, "U") == colnames(land_cover_dataframe))) {
            message(paste0("Updating ages for ", vt, "..."))
            
            # Get column indices
            veg_type_first_col <- match(paste0(vt, "R"), 
                                        colnames(land_cover_dataframe))
            veg_type_last_col <- match(paste0(vt, "8"), 
                                       colnames(land_cover_dataframe))
            nsr_name_col <- match("NSRNAME", colnames(land_cover_dataframe))
            veg_type_unknown_age_col = match(paste0(vt, "U"), 
                                             colnames(land_cover_dataframe))
            
            # Check if there are any unknown areas
            if (sum(land_cover_dataframe[[veg_type_unknown_age_col]]) > 0) {
                # Get a subset of the data table
                cols <- c(nsr_name_col, veg_type_unknown_age_col)
                unknown_veg <- land_cover_dataframe[get(paste0(vt, "U")) > 0, 
                                                    ..cols]
                names(unknown_veg) <- c("NSRNAME", "U")
                
                # Grab the assumed distribution
                raw_nsr_distribution <- data.table::as.data.table(
                    nsr_age_distributions[[vt]], keep.rownames="NSRNAME"
                    )
                print(head(raw_nsr_distribution))
                
                # Join the distributions to the cells by NSRNAME
                nsr_distribution_by_cell <- 
                    data.table::merge.data.table(unknown_veg[,1], 
                                                 raw_nsr_distribution, 
                                                 by="NSRNAME", all.x=TRUE)
                
                # Separate the text and numbers
                nsr_names <- nsr_distribution_by_cell[, 1, with = FALSE]
                distributions <- nsr_distribution_by_cell[, -c(1, 2, 12), 
                                                          with = FALSE]
                
                # Get estimated areas
                distributions_by_cell <- distributions * unknown_veg$U
                
                # Update the dataframe with the original + new area
                column_names <- names(distributions)
                rows_to_update <- land_cover_dataframe[[paste0(vt, "U")]] > 0
                for (col_name in column_names) {
                    land_cover_dataframe[rows_to_update, (col_name) := 
                                             get(col_name) + 
                                             distributions_by_cell[[col_name]]]
                }
                
                # Set original unknown areas to zero
                land_cover_dataframe[, (paste0(vt, "U")) := 0]
                
                
            } else {
                print(paste0(vt, " does not need updating."))
            }
        }
    }
    return(as.data.frame(land_cover_dataframe))
}
