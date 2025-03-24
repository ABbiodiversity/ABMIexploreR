#' Internal helper functions
#'
#' Provide SDM coefficients
#'
#' @description Helper functions to list species climate, vegetation, soil, and footprint coefficients
#'
#' @param x Input
#'

#' @name helper
NULL

#' @export
#' @rdname helper
abmi_species <- function() {
    if (!.loaded())
        stop("Use abmi_load_coefs() to load coefs")
    out <- NULL
    for (taxon in names(.abmi1$COEFS)) {
        tmp <- .abmi1$COEFS[[taxon]]$species
        tmp <- tmp[tmp$ModelNorth | tmp$ModelSouth,]
        out <- rbind(out, tmp)
    }
    out
}

#' Provide landcover definitions
#'
#' @description Helper functions to list species vegetation, soil, and footprint definitions
#' @param landcover Users can define Vegetation or Soil landcover classes
#' 
#' @export
#' @rdname helper
abmi_landcover <- function(landcover) {
    
    return(.landcover.names[[landcover]])
    
}

#' Provide climate definitions
#'
#' @description Helper functions to list bioclimatic definitions
#' 
#' @export
#' @rdname helper
abmi_climate <- function() {
    
    return(.bioclimatic.lookup)
    
}

#' Handle package options
#'
#' @description Helper functions to attach or detech package options to the current space
#' 
#' @export
#' @rdname helper
#' 
.onAttach <- function(libname, pkgname){
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields=c("Version", "Date"))
    packageStartupMessage(paste(pkgname, ver[1], "\t", ver[2]))
    if (is.null(getOption("ABMIexploreR"))) {
        options("ABMIexploreR" = list(
            url="https://github.com/beallen/ABMIexploreR-data/v2025",
            dir=".",
            verbose=1))
    }
    invisible(NULL)
}

.onUnload <- function(libpath){
    options("ABMIexploreR" = NULL)
    invisible(NULL)
}