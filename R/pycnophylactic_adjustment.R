#' Adjust a smoothed raster to preserve mass of unsmoothed raster
#'
#' The function adjusts a smoothed raster based on the difference in means 
#' between the original raster and the smoothed raster). Pycnophylactic 
#' interpolation (PP) is applied to the difference, resulting in a smooth 
#' surface where of the errors, but where the mean error within the zone 
#' equals the original mean error. The smoothed input raster is then adjusted 
#' with the new surface. In practice, smoothing this way is a combination of 
#' the two surfaces: the original smoothed raster, and the smoothed error 
#' obtained from PP. The function sets all negative values in the 
#' smoothed raster as 0 The function also keeps cells with value 0 as 0, and 
#' therefore any negative value are set to zero and kept as such.
#' 
#' @param r1 unsmoothed \code{RasterLayer}, from which we get zonal mass
#' @param r2 smoothed \code{RasterLayer}, which' zonal mass we adjust
#' @param zones a \code{RasterLayer} with zones in which adjustment takes place
#' @param adjust_threshold How much the values of r2 are allowed to change in 
#'   the adjustment. A value of 0.5 will allow a maximum adjustment of 50% in 
#'   r2 in order to restore balance with r1. Within range [0,1].
#' @param n number of iterations for pycnophylactic adjustment
#' @param smoothing_matrix neighbourhood for pycnophylactic interpolation. 
#'   Experimental, default 3x3 neighbourhood suggested.
#' @param intensive Are values r1 and r2 densities (\code{TRUE}) or counts 
#'  (\code{FALSE}).
#' @param return_error If \code{TRUE}, returns the pycnophylactic adjustment 
#' layer in addition to the adjusted. 
#' @param verbose Print progress indication?
#' 
#' @return Volume corrected r2 as single \code{RasterLayer}, or , if 
#'   \code{return_error == TRUE}, a two-layer raster with corrected r2 and 
#'   the layer used for the correction.
#'   
#' @export
pycnophylactic_adjustment <- function(r1, r2, 
                                      zones, 
                                      adjust_threshold = 0.5,
                                      n=5, 
                                      smoothing_matrix = matrix(1,3,3),
                                      intensive = TRUE,
                                      return_error = FALSE,
                                      verbose = FALSE) {
    
    threshold <- NULL
    . <- NULL
    orig_err_mean <- NULL
    zone <- NULL
    err  <- NULL
    pixel <- NULL
    smooth_mean <- NULL
    mean_sm <- NULL
    adjust <- NULL
    high_err <- NULL
    
    
    if(verbose) message("Preparing..")
    
    if(intensive) {
        r_zonal_mean <- raster::zonal(r1*raster::area(r1), zones)
        rsm_zonal_mean <- raster::zonal(r2*raster::area(r2), zones)
    } else {
        r_zonal_mean <- raster::zonal(r1, zones)
        rsm_zonal_mean <- raster::zonal(r2, zones)
    }
    
    zonal_diff <- cbind(r_zonal_mean, sm_mean = rsm_zonal_mean[,2])
    zonal_diff <- cbind(zonal_diff, diff = zonal_diff[,3]-zonal_diff[,2])
    
    rzd <- raster::reclassify(zones, zonal_diff[,c(1,4)])
    
    ## TEST WHETHER POSSIBLE
    diff_test <- zonal_diff[,4] / zonal_diff[,3]
    diff_test[is.infinite(diff_test) | is.nan(diff_test)] <- 0 
    zonal_thresholds <- data.table(zone = zonal_diff[,1], 
                                   #threshold = abs(diff_test),
                                   threshold = adjust_threshold,
                                   key="zone")
    
    if(is.null(adjust_threshold)) adjust_threshold <- 0.5
    
    test <- any(abs(diff_test) > 1)
    if(test) stop("adjusting for difference not possible with these zones.")
    diff_test <- any(abs(diff_test) > adjust_threshold)
    
    
    
    # pycnophylactic interpolation
    if(verbose) message("Pycnophylactic interpolation of errors..")
    rzd_i <- as.matrix(rzd)
    ncol <- ncol(rzd)
    nrow <- nrow(rzd)
    # boundary condition
    rleft <- rzd_i[,1]
    rright <- rzd_i[,ncol]
    rtop <- c(NA, rzd_i[1,], NA)
    rbottom <- c(NA, rzd_i[nrow,], NA)
    
    if(intensive) {
        areas <- raster::values(raster::area(rzd))
    } else {
        areas <- 1
    }
    
    ## construct a matrix which is used to fix values outside zones,
    ## and to keep r2 <= 0 at 0
    smooth_r2 <- raster::values(r2) * areas
    zs <- raster::values(zones)
    z <- !is.na(zs) & is.na(smooth_r2)
    zs[z] <- NA
    nochange <- rep(NA, length(smooth_r2))
    below_zero <- smooth_r2 <= 0 
    inds <- which(!zs & !is.na(smooth_r2))
    if(length(inds) != 0) nochange[inds] <- smooth_r2[inds]
    nochange[below_zero] <- 0
    negative <- smooth_r2 < 0 
    nc_inds <- !is.na(nochange)
    
    
    adj <- data.table::data.table(pixel = 1:raster::ncell(zones),
                                  zone = zs,
                                  areas = areas,
                                  r1 = raster::values(r1)*areas,
                                  smooth_r2 = smooth_r2,
                                  negative = negative,
                                  orig_err_mean = raster::values(rzd),
                                  key = c("pixel", "zone", "smooth_r2"))
    adj <- zonal_thresholds[adj, on="zone"]
    adj[, max_err := abs(smooth_r2 * threshold*!negative)]
    
    ## TEST WHETHER POSSIBLE errors
    tst <- adj[, .(max_err = round(sum(max_err, na.rm=TRUE),2),
                   err = round(sum(orig_err_mean, na.rm=TRUE),2)), 
               by=zone]
    tst[, ':='(ratio = round(err/max_err,2),
            test = round(abs(tst$err),5) > round(tst$max_err,5))]
    test <- any(tst$test[-1])
    #test <- FALSE
    if(test) stop("cannot reallocate error with current threshold. Please increase.")
    if(!test && diff_test) warning(paste0("threshold smaller than difference ",
                                          "between r1 and r2. adjustment may ",
                                          "be unstable."))
    
    # ITERATE n TIMES
    pb <- utils::txtProgressBar(min = 0, max = n, style=3) 
    for(i in 1:n) {
        
        # add boundary condition
        rzd_i[nc_inds] <- nochange[nc_inds]
        rzd_i <- cbind(rleft, rzd_i, rright)
        rzd_i <- rbind(rtop, rzd_i, rbottom)
        
        # smooth
        rzd_i <- raster::raster(rzd_i)
        rzd_i <- raster::focal(rzd_i, 
                               w=smoothing_matrix, 
                               fun=mean, 
                               pad=TRUE, 
                               na.rm=TRUE)
        rzd_i <- as.matrix(rzd_i)[c(-1,-(nrow+1)), c(-1,-(ncol+1))]
        rzd_i[nc_inds] <- nochange[nc_inds]
        
        adj <- adj[order(pixel)]
        adj[,smooth_mean := as.vector(t(rzd_i))]
        adj[, mean_sm := mean(smooth_mean, na.rm=TRUE), by = zone]
        adj[, adjust := ifelse(orig_err_mean == 0,
                               1,
                               (mean_sm)/(orig_err_mean))]
        adj[, new_err := smooth_mean/adjust]
        adj[, high_err := abs(new_err) > abs(max_err)]
        
        iterr <- adj$new_err
        if(any(adj$high_err)) {
            
            uniq_zone <- unique(adj$zone)
            output <- list()
            adj <- adj[, fixed_err := NA_real_]
            adj <- adj[order(zone,smooth_r2)]
            for(z in uniq_zone) {
                if(is.na(z)) {
                    next
                } 
                max_err <- adj[zone == z, max_err]
                new_err <- adj[zone == z, new_err]
                r2val <- adj[zone == z, smooth_r2]
                fixed_err <- rep(NA, length(new_err))
                shares <- max_err / sum(abs(max_err))
                share_left <- sum(shares)
                shares_visited <- 0
                total_err <- 0
                
                for(j in seq_along(new_err)) {
                    share_of_err <- shares[j] / (1-shares_visited)
                    sherr <- total_err*share_of_err
                    test <- abs(new_err[j]+sherr) > abs(max_err[j])
                    if(is.na(test)) next
                    pos <- new_err[j] >= 0
                    if(test) {
                        
                        if(pos) {
                            dist_err = new_err[j] - max_err[j];
                            fix_err = max_err[j];
                        }  else {
                            dist_err = new_err[j] + max_err[j];
                            fix_err = -max_err[j];
                        }
                        
                        
                        shares_visited <- shares_visited + shares[j]
                        total_err <- total_err + dist_err
                        
                        fixed_err[j] <- fix_err
                        
                    } else {
                        fixed_err[j] <- new_err[j] + sherr
                        total_err <- total_err - sherr
                        shares_visited <- shares_visited + shares[j]
                    }
                }
                errvec <- fixed_err
                adj[zone == z, fixed_err := errvec]
            }
            adj <- adj[order(pixel)]
            iterr <- adj$fixed_err
        }
        rzd_i <- matrix(iterr, nrow=nrow, ncol=ncol,
                        byrow=TRUE)
        if(verbose) utils::setTxtProgressBar(pb, i)
    }
    
    if(verbose) close(pb)
    
    message("Preparing output..")
    raster::values(rzd) <- rzd_i 
    if(intensive) rzd <- rzd/raster::area(rzd)
    sm_out <- r2 - rzd
    if(return_error) {
        out <- raster::stack(sm_out, r2, rzd)
        names(out) <- c("Adjusted", "Original", "Error")
        return(out)
    } else {
        return(sm_out)
    }
} 


