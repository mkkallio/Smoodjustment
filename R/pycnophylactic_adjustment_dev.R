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
pycnophylactic_adjustment2 <- function(r1, r2, 
                                      zones, 
                                      #adjust_threshold = c(-0.5,0.75),
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
        areas <- raster::values(raster::area(zones))
    } else {
        areas <- 1
    }
    
    # remove zones-pixels with nodata in r2
    orig_r1 <- raster::values(r1) * areas
    smooth_r2 <- raster::values(r2) * areas
    below_zero <- smooth_r2 <= 0
    smooth_r2[below_zero] <- 0
    zs <- raster::values(zones)
    z <- !is.na(zs) & is.na(smooth_r2)
    zs[z] <- NA
    raster::values(zones) <- zs
    
    
    ### prepare
    adj <- data.table::data.table(pixel = 1:raster::ncell(zones),
                                  zone = zs,
                                  areas = areas,
                                  orig_r1 = orig_r1,
                                  smooth_r2 = smooth_r2,
                                  key = c("pixel", "zone", "smooth_r2"))
    
    
    ###########################
    # COMPUTE ZONAL DIFFERENCES BETWEEN r1 AND r2
    if(intensive) {# convert to volume?
        r_zonal_mean <- raster::zonal(r1*raster::area(r1), zones)
        rsm_zonal_mean <- raster::zonal(r2*raster::area(r2), zones)
    } else {
        r_zonal_mean <- raster::zonal(r1, zones)
        rsm_zonal_mean <- raster::zonal(r2, zones)
    }
    zonal_diff <- data.table::data.table(zone = r_zonal_mean[,1],
                                         r1_mean = r_zonal_mean[,2],
                                         sm_mean = rsm_zonal_mean[,2])
    
    zonal_diff[, orig_diff := sm_mean-r1_mean]
    zonal_diff[, ratio := orig_diff/sm_mean]
    zonal_diff[is.na(ratio) | is.nan(ratio) | is.infinite(ratio), ratio := 1]
    #threshold <-  abs(zonal_diff$ratio)
    zonal_diff[, threshold := abs(ratio)*1.25 ]
    
    ## TEST WHETHER REALLOCATION IS POSSIBLE
    test <- any(abs(zonal_diff$ratio) > 1)
    if(test) stop("adjusting for difference not possible with these zones.")
    #diff_test <- any(abs(zonal_diff$ratio) > adjust_threshold)
    
    
    
    
    ##################
    # PREPARE FOR PYCNOPHYLACTIC INTERPOLATION
    adj <- adj[zonal_diff, on=.(zone), ':='(orig_diff = i.orig_diff,
                                         threshold = i.threshold,
                                         r_ratio = i.ratio)]
    adj[, error := ifelse(is.na(orig_r1),
                          smooth_r2,
                          smooth_r2-orig_r1)]
    adj[, error_ratio := abs(error)/sum(abs(error), na.rm=TRUE), by=zone]
    
    ## rather than divide the error equally, do some initial distribution
    ## based on the error
    adj[, rzd := error_ratio*sum(orig_diff, na.rm=TRUE), by=zone]

    ### 
    rzd <- zones
    values(rzd) <- adj$rzd
    rzd_i <- as.matrix(rzd)
    ncol <- ncol(rzd)
    nrow <- nrow(rzd)
    
    # set up boundary condition surrounding the raster
    rleft <- rzd_i[,1]
    rright <- rzd_i[,ncol]
    rtop <- c(NA, rzd_i[1,], NA)
    rbottom <- c(NA, rzd_i[nrow,], NA)
    
    ## construct a matrix which is used to fix values outside zones,
    ## and to keep r2 <= 0 at 0
    nochange <- rep(NA, length(smooth_r2))
    negative <- smooth_r2 < 0 
    # cells which are outside of zones and have values
    inds <- which(!zs & !is.na(smooth_r2))
    if(length(inds) != 0) nochange[inds] <- smooth_r2[inds]
    # which are zero or negative
    nochange[below_zero] <- 0
    nc_inds <- !is.na(nochange)
    nc_inds <- matrix(nc_inds, nrow, ncol, byrow=TRUE)
    below_zero <- matrix(below_zero, nrow, ncol, byrow=TRUE)
    
    
    #####
    # MAXIMUM ERROR
    adj[, ':='(negative = negative)]
    adj[, max_err := abs(smooth_r2 * threshold*!negative)]
    
    
    #####
    ## TEST MAXIMUM ERROR
    tst <- adj[, .(max_err = abs(round(sum(max_err, na.rm=TRUE),2)),
                   err = round(sum(orig_diff, na.rm=TRUE),2)), 
               by=zone]
    tst[, ':='(max_ratio = round(err/max_err,2),
               max_test = round(abs(tst$err),5) > round(tst$max_err,5))]
    test <- any(tst$max_test[-1])
    test <- FALSE
    # if(test) stop("cannot reallocate error with current threshold. Please increase.")
    if(test) {
        warning(paste0("WARNING: cannot reallocate error with current ", 
                        "threshold. Returning the error table."))
        return(tst)
    }
    # if(!test && diff_test) warning(paste0("threshold smaller than difference ",
    #                                       "between r1 and r2. adjustment may ",
    #                                       "be unstable."))
    
    
    ##################
    # ITERATE n TIMES
    if(verbose) message("Pycnophylactic interpolation of errors..")
    pb <- utils::txtProgressBar(min = 0, max = n, style=3) 
    #n <- 1
    #keep_iter <- TRUE
    #while(keep_iter & n <= nmax ) {
    for(i in 1:n) {
        
        # add boundary condition
        rzd_i[below_zero] <- 0
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
        rzd_i[below_zero] <- 0
        
        adj <- adj[order(pixel)]
        adj[,smooth_mean := as.vector(t(rzd_i))]
        adj[, mean_sm := ifelse(is.na(smooth_mean),
                                NA,
                                mean(smooth_mean, na.rm=TRUE)), by = zone]
        adj[, adjust := ifelse(orig_diff == 0,
                               1,
                               (mean_sm)/(orig_diff))]
        adj[, new_err := smooth_mean/adjust]
        adj[, high_err := abs(new_err) > max_err]
        
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
                value <- adj[zone == z, smooth_r2]
                fixed_err <- new_err
                shares <- max_err / sum(abs(max_err))
                #share_left <- sum(shares)
                shares_visited <- 0
                total_err <- 0
                keep <- TRUE
                order <- rev(seq_along(new_err))
                nrev <- 0
                while(keep) {
                    nrev <- nrev + 1
                    order <- rev(order)
                    for(j in order) {
                        share_of_err <- shares[j] / (1-shares_visited)
                        sherr <- total_err*share_of_err
                        pos <- fixed_err[j] >= 0
                        test <- abs(fixed_err[j]+sherr) > abs(max_err[j])
                        test
                        if(is.na(test)) next
                        if(test) {
                            if(pos) {
                                dist_err = fixed_err[j] - max_err[j];
                                fix_err = max_err[j];
                            }  else {
                                dist_err = fixed_err[j] + max_err[j];
                                fix_err = -max_err[j];
                            }
                            new_val <- value[j] - dist_err
                            
                            shares_visited <- shares_visited + shares[j]
                            total_err <- total_err + dist_err
                            
                            fixed_err[j] <- fix_err
                            
                        } else {
                            fixed_err[j] <- fixed_err[j] + sherr
                            total_err <- total_err - sherr
                            shares_visited <- shares_visited + shares[j]
                        }
                    }
                    if(round(total_err,5) == 0 || nrev == 25) {
                        keep <- FALSE
                        break
                    }
                    shares_visited <- 0
                }
                errvec <- fixed_err
                adj[zone == z, fixed_err := errvec]
                #print(z)
            }
            adj <- adj[order(pixel)]
            iterr <- adj$fixed_err
            test <- adj[, .(sum_fix = sum(fixed_err, na.rm=TRUE),
                            sum_orig = sum(rzd, na.rm=TRUE)), by=zone]
            test[, test := round(sum_fix, 5) == round(sum_orig, 5)]
            n_test <- sum(test$test, na.rm=TRUE)
        }
        #n <- n+1
        if(n_test == nrow(test)) break
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
        out <- raster::stack(sm_out, r2, r1, rzd)
        names(out) <- c("Adjusted", "Smoothed", "Unsmoothed", "Error")
        return(out)
    } else {
        return(sm_out)
    }
} 

# 
# zm_runoff <- raster::zonal(w_rs*area(w_rs), zones1)
# zm_runoff_smooth <- raster::zonal(out[[1]]*area(w_rs), zones1)
# all.equal(zm_runoff, zm_runoff_smooth)
# 
# any(values(out[[1]]) < 0, na.rm=TRUE)
# t <- values(out[[1]])
# table(t < 0)
# #writeRaster(test, "testnew_50_abs_err.tif")
# 
# dt <- as.data.table(zm_runoff)
# dt[, sm_mean := zm_runoff_smooth[,2]]
# dt[, test_equal := round(mean, 10) == round(sm_mean,10)]
# dt[, ratio := sm_mean/mean][]
# 
# 
# 
# system.time({
#     test <- pycnophylactic_adjustment_dev(w_rs, wsm_rs, zones, 
#                                       adjust_threshold = 0.5, n=5, 
#                                       return_error = TRUE, verbose=TRUE)
# })
# 
# zm_runoff <- raster::zonal(w_rs*area(w_rs), zones)
# zm_runoff_smooth <- raster::zonal(test[[1]]*area(w_rs), zones)
# all.equal(zm_runoff, zm_runoff_smooth)
# 
# any(values(test[[1]]) < 0, na.rm=TRUE)
# t <- values(test[[1]])
# table(t < 0)
# #writeRaster(test, "testnew_50_abs_err.tif")
# 
# dt <- as.data.table(zm_runoff)
# dt[, sm_mean := zm_runoff_smooth[,2]]
# dt[, test_equal := round(mean, 10) == round(sm_mean,10)]
# dt[, ratio := sm_mean/mean][]
# #tempdt <- dt
# 
# 
# system.time({
#     test <- pycnophylactic_adjustment(w_rs, wsm_rs, vzr, 
#                                           adjust_threshold = 0.9, n=5, 
#                                           return_error = TRUE, verbose=TRUE)
# })
# 
# zm_runoff <- raster::zonal(w_rs*area(w_rs), zones)
# zm_runoff_smooth <- raster::zonal(test[[1]]*area(w_rs), zones)
# all.equal(zm_runoff, zm_runoff_smooth)
# 
# any(values(test[[1]]) < 0, na.rm=TRUE)
# t <- values(test[[1]])
# table(t < 0)
# 
# dt <- as.data.table(zm_runoff)
# dt[, sm_mean := zm_runoff_smooth[,2]]
# dt[, test_equal := round(mean, 10) == round(sm_mean,10)]
# dt[, ratio := sm_mean/mean][]