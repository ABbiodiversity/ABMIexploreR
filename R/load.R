#' Coefficients
#'
#' Download, load, unload coefficients.

#' @param dir Directory path.
#' @param ... Arguments parameters passed to [utils::download.file()].
#'
#' @return
#'
#' `abmi_dowload_coefs` dowloads a file as a side effect.
#'
#' `abmi_load_coefs` loads the coefficients into the current session.
#'
#' `abmi_unload_coefs` unloads the coefficients from the current session.
#'
#' @seealso [utils::download.file()].
#'
#' @name coefs
NULL

#' @export
#' @rdname coefs
abmi_dowload_coefs <- function(dir=NULL, ...) {
    if (is.null(dir))
        dir <- getOption("ABMIexploreR")$dir
    fn <- file.path(getOption("ABMIexploreR")$url, "COEFS.RData")
    fo <- file.path(dir, "COEFS.RData")
    if (file.exists(fo))
        .msg("File already exists", 2)
    .msg("Downloading coefs")
    utils::download.file(
        url=fn,
        destfile=fo, ...)
}

#' @export
#' @rdname coefs
abmi_load_coefs <- function(dir=NULL) {
    if (.loaded())
        invisible(NULL)
    if (is.null(dir))
        dir <- getOption("ABMIexploreR")$dir
    fn <- file.path(dir, "COEFS.RData")
    .msg("Loading coefs")
    out <- try(load(fn, envir=.abmi1))
    if (inherits(out, "try-error"))
        stop("Use abmi_download_coefs() to download coefs")
    invisible(out)
}

#' @export
#' @rdname coefs
abmi_unload_coefs <- function() {
    if (.loaded())
        .msg("Unloading coefs")
    rm(list=ls(envir=.abmi1), envir=.abmi1)
}
