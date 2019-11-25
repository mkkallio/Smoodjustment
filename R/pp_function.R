
pycnophylactic_smoothing <- function(r1, r2, 
                                     zones, 
                                     adjust_threshold = NULL,
                                     n=50, 
                                     smoothing_matrix = matrix(1,3,3),
                                     intensive = TRUE,
                                     return_error = FALSE,
                                     verbose = FALSE) {
    
    if(verbose) message("Preparing..")
    
    # retain
    if(intensive) {
        r_zonal_mean <- raster::zonal(r1*area(r1), zones)
        rsm_zonal_mean <- raster::zonal(r2*area(r2), zones)
    } else {
        r_zonal_mean <- raster::zonal(r1, zones)
        rsm_zonal_mean <- raster::zonal(r2, zones)
    }
    
    zonal_diff <- cbind(r_zonal_mean, sm_mean = rsm_zonal_mean[,2])
    zonal_diff <- cbind(zonal_diff, diff = zonal_diff[,3]-zonal_diff[,2])
    
    rzd <- raster::reclassify(zones, zonal_diff[,c(1,4)])
    
    ## TEST WHETHER POSSIBLE
    diff_test <- zonal_diff[,4] / zonal_diff[,2]
    if(is.null(adjust_threshold)) adjust_threshold <- max(abs(diff_test))
    
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
        areas <- values(area(rzd))
    } else {
        areas <- 1
    }
    
    adj_tmpl <- tibble::tibble(pixel = 1:ncell(zones),
                               zone = values(zones), 
                               areas = areas,
                               r1 = values(r1) * areas,
                               smooth_r2 = values(r2) * areas,
                               no_change = smooth_r2 <= 0,
                               negative = smooth_r2 < 0,
                               max_err = abs(smooth_r2 *
                                                 adjust_threshold * 
                                                 !negative),
                               orig_err_mean = values(rzd))
    
    nochange <- matrix(adj_tmpl$no_change, nrow=nrow, ncol=ncol,
                       byrow=TRUE)
    
    
    ## TEST WHETHER POSSIBLE errors
    tst <- adj_tmpl %>% 
        dplyr::select(zone, max_err, orig_err_mean) %>%
        dplyr::group_by(zone) %>%
        dplyr::summarise(max_err = sum(max_err, na.rm=TRUE),
                  err = sum(orig_err_mean, na.rm=TRUE),
                  test = abs(err) > max_err)
    test <- any(tst$test)
    if(test) stop("cannot reallocate error with current threshold. Please increase.")
    if(!test && diff_test) warning(paste0("threshold smaller than difference ",
                                          "between r1 and r2. adjustment may ",
                                          "be unstable."))
    
    
    
    pb <- txtProgressBar(min = 0, max = n, style=3) 
    for(i in 1:n) {
        
        # add boundary condition
        rzd_i[nochange] <- 0
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
        rzd_i[nochange] <- 0
        
        adj <- adj_tmpl %>%
            mutate(smooth_mean = as.vector(t(rzd_i))) %>%
            dplyr::group_by(zone) %>%
            dplyr::mutate(mean_sm = mean(smooth_mean, na.rm=TRUE),
                          adjust = ifelse(orig_err_mean == 0,
                                          1,
                                          (mean_sm)/(orig_err_mean)),
                          new_err = smooth_mean/adjust,
                          high_err = abs(new_err) > abs(max_err))
        
        iterr <- adj$new_err
        if(any(adj$high_err)) {
            
            uniq_zone <- unique(adj$zone)
            output <- list()
            for(z in uniq_zone) {
                if(is.na(z)) {
                    order <- adj %>% 
                        dplyr::filter(is.na(zone)) %>%
                        dplyr::select(zone, pixel, areas, orig_err_mean,
                               smooth_r2, max_err, new_err) %>%
                        dplyr::arrange(smooth_r2) %>% 
                        dplyr::mutate(fixed_err = new_err)
                    output <- append(output, list(order))
                    
                    next
                } 
                order <- adj %>% 
                    dplyr::filter(zone == z) %>%
                    dplyr::select(zone, pixel, areas, orig_err_mean,
                           smooth_r2, max_err, new_err, adjust) %>%
                    dplyr::arrange(smooth_r2)  
                
                max_err <- order$max_err
                new_err <- order$new_err
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
                            # if(total_err >= 0) {
                            #     dist_err = new_err[j] + max_err[j];
                            #     fix_err = -max_err[j];
                            # } else {
                            #     dist_err = new_err[j] - max_err[j];
                            #     fix_err = -max_err[j];
                            # }
                            
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
               
                
                order <- order %>%
                    dplyr::mutate(fixed_err = fixed_err)
                output <- append(output, list(order))
            }
            output <- do.call(rbind, output) 
            output <- arrange(output, pixel)
            iterr <- output$fixed_err
        }
        rzd_i <- matrix(iterr, nrow=nrow, ncol=ncol,
                        byrow=TRUE)
        if(verbose) setTxtProgressBar(pb, i)
    }
    
    if(verbose) close(pb)
    
    message("Preparing output..")
    values(rzd) <- rzd_i 
    if(intensive) rzd <- rzd/area(rzd)
    sm_out <- r2 - rzd
    if(return_error) {
        out <- raster::stack(sm_out, r2, rzd)
        names(out) <- c("Adjusted", "Original", "Error")
        return(out)
    } else {
        return(sm_out)
    }
} 
