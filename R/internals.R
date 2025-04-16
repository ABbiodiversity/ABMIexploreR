## internal functions -----------------

## hidden environment to store coefs in
.abmi1 <- new.env(parent=emptyenv())
## store object for subset of the grid and species
#.abmi1=ABMIexploreR:::.abmi1

## path to COEFS
.path <- function() {
    getOption("ABMIexploreR")$url
}

## check if coefs are loaded
.loaded <- function() {
    length(names(.abmi1)) > 0
}

## check if user wants messages
.verbose <- function() {
    x <- getOption("ABMIexploreR")$verbose
    !is.null(x) && x > 0
}

## print messages if verbose
.msg <- function(x, level=1) {
    if (.verbose()) {
        if (level >= getOption("ABMIexploreR")$verbose) {
            pf <- switch(level,
                "1"="[INFO] ",
                "2"="[NOTE] ",
                "3"="[WARN] ",
                "4"="")
            msg <- paste0(pf, x, "\n")
            if (level >= 4) {
                stop(msg, call.=FALSE)
            } else {
                cat(msg)
            }
            utils::flush.console()
        }
    }
    invisible()
}


