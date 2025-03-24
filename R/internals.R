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

## transform climate variables
.make_clim <- function(x, taxon) {
    if (taxon == "birds") {
        z <- with(x, cbind(
            pWater_KM=pWater,
            pWater2_KM=pWater^2,
            xPET=(PET - 0) / 800,
            xMAT=(MAT - 0) / 6,
            xAHM=(AHM - 0) / 50,
            xFFP=(FFP - 0) / 130,
            xMAP=(MAP - 0) / 2300,
            xMWMT=(MWMT - 0) / 20,
            xMCMT=(MCMT - 0) / 25,
            xY=(POINT_Y - 54.1) / 2.43,
            xX=(POINT_X - (-113.3)) / 2.12))
        z <- cbind(z,
            xY2=z[,"xY"]^2,
            xX2=z[,"xX"]^2,
            `xFFP:xMAP`=z[,"xFFP"]*z[,"xMAP"],
            `xMAP:xPET`=z[,"xMAP"]*z[,"xPET"],
            `xAHM:xMAT`=z[,"xAHM"]*z[,"xMAT"],
            `xX:xY`=z[,"xX"]*z[,"xY"])
    } else {
        LAT <- pmin(x$POINT_Y, 56.5)
        z <- with(x, cbind(
            Intercept=1,
            Lat=LAT,
            Long=POINT_X,
            AHM=AHM,
            PET=PET,
            FFP=FFP,
            MAP=MAP,
            MAT=MAT,
            MCMT=MCMT,
            MWMT=MWMT,
            Lat2=LAT^2,
            Long2=POINT_X^2,
            LatLong=POINT_X*LAT,
            MAPPET=MAP*x$PET,
            MATAHM=MAT*x$AHM,
            MAPFFP=MAP*x$FFP,
            MAT2=MAT^2,
            MWMT2=MWMT^2))
        if (taxon == "mammals") {
            z <- cbind(z,
                Lat3=z[,"Lat"]^3,
                Lat2Long2=z[,"Lat"]^2 * z[,"Long"]^2,
                LongMAT=z[,"Long"] * z[,"MAT"])
            z[,"MAT2"] <- z[,"MAT"]*(z[,"MAT"]+10)
            z <- cbind(z, x[,c("PeaceRiver", "NSR1CentralMixedwood",
                    "NSR1DryMixedwood", "NSR1Foothills",
                    "NSR1Mountain", "NSR1North", "NSR1Parkland", "NSR1Shield")])
        }
    }
    rownames(z) <- rownames(x)
    as.matrix(z)
}

.get_taxon <- function(spp) {
    for (i in names(.abmi1$COEFS)) {
        SPP <- rownames(.abmi1$COEFS[[i]]$species)
        if (tolower(spp) %in% tolower(SPP))
            return(i)
    }
    .msg(sprintf("Coefs not available for species %s", spp), 4)
}
