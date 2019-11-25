
pycnophylactic_smoothing <- function(r1, r2, 
                                     zones, 
                                     #adjust_threshold,
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
                               nazone_val = ifelse(is.na(zone),
                                                   smooth_r2,
                                                   NA),
                               noch_val = ifelse(no_change,
                                                   0,
                                                   nazone_val),
                               max_err = abs(smooth_r2 *
                                                 adjust_threshold * 
                                                 !negative),
                               orig_err_mean = values(rzd))
    
    nochange <- matrix(adj_tmpl$no_change, nrow=nrow, ncol=ncol,
                       byrow=TRUE)
    nochange_val <- matrix(adj_tmpl$noch_val, nrow=nrow, ncol=ncol,
                           byrow=TRUE)
    
    ## TEST WHETHER POSSIBLE
    tst <- adj_tmpl %>% 
        select(zone, max_err, orig_err_mean) %>%
        group_by(zone) %>%
        summarise(max_err = sum(max_err, na.rm=TRUE),
                  err = sum(orig_err_mean, na.rm=TRUE),
                  test = abs(err) > max_err)
    test <- any(tst$test)
    if(test) stop("cannot reallocate error with current threshold. Please increase.")
    
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
                        filter(is.na(zone)) %>%
                        select(zone, pixel, areas, orig_err_mean,
                               smooth_r2, max_err, new_err) %>%
                        arrange(smooth_r2) %>% 
                        mutate(fixed_err = new_err)
                    output <- append(output, list(order))
                    
                    next
                } 
                order <- adj %>% 
                    filter(zone == z) %>%
                    select(zone, pixel, areas, orig_err_mean,
                           smooth_r2, max_err, new_err, adjust) %>%
                    arrange(smooth_r2)  
                
                max_err <- order$max_err
                new_err <- order$new_err
                fixed_err <- rep(NA, length(new_err))
                shares <- max_err / sum(abs(max_err))
                share_left <- sum(shares)
                shares_visited <- 0
                total_err <- 0
                
                recte <- rep(NA, length(new_err))
                recsh <- rep(NA, length(new_err))
                recsherr <- rep(NA, length(new_err))
                recshh <- rep(NA, length(new_err))
                rectest1 <- rep(NA, length(new_err))
                rectest2 <- rep(NA, length(new_err))
                recdist <- rep(NA, length(new_err))
                rect <-  rep(NA, length(new_err))
                
                for(j in seq_along(new_err)) {
                    share_of_err <- shares[j] / (1-shares_visited)
                    sherr <- total_err*share_of_err
                    test <- abs(new_err[j]+sherr) > abs(max_err[j])
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
                        
                        recdist[j] <- dist_err
                        rectest1[j] <- fix_err
                    } else {
                        fixed_err[j] <- new_err[j] + sherr
                        total_err <- total_err - sherr
                        shares_visited <- shares_visited + shares[j]
                    }
                    recte[j] <- total_err
                    recsh[j] <- shares_visited
                    recsherr[j] <- sherr
                    recshh[j] <- share_of_err
                    
                }
                # stopifnot(round(sum(fixed_err, na.rm=TRUE),8) == 
                #               round(sum(new_err, na.rm=TRUE)), 8)
                
                order <- order %>%
                    mutate(fixed_err = fixed_err)
                
                
                order$te <- recte
                order$sh <- recsh
                order$sherr <- recsherr
                order$shoferr <- recshh
                #order$test1 <- rectest1
                #order$test2 <- rectest2
                order$dist <- recdist
                #order$fix <- rectest1
                
                # tibble(o_new_err = old_order$new_err, 
                #        new_err = order$new_err,
                #        max_err = order$max_err,
                #        o_fix_err =  old_order$fixed_err,
                #        fix_err = order$fixed_err,
                #        o_te = old_order$te,
                #        te = order$te,
                #        o_sherr = old_order$sherr,
                #        sherr = order$sherr,
                #        o_dist = old_order$dist,
                #        dist = order$dist) %>%
                #     head(20)
                
                output <- append(output, list(order))
            }
            output <- do.call(rbind, output) 
            test <- output %>% group_by(zone) %>%
                summarise(fix = sum(fixed_err, na.rm=TRUE),
                          new_err = sum(new_err, na.rm=TRUE),
                          test = round(fix,5) == round(new_err, 5)) #%>%
                #print(n=30)
            print(all(test$test))
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
    if(intensive) {
        rzd <- rzd/area(rzd)
    } 
    sm_out <- r2 - rzd
    if(return_error) {
        out <- raster::stack(sm_out, r2, rzd)
        names(out) <- c("Adjusted", "Original", "Error")
        return(out)
    } else {
        return(sm_out)
    }
} 


system.time({
    test <- pycnophylactic_smoothing(w_rs, wsm_rs, 
                                     zones1, n=5, adjust_threshold = 0.5,
                                     return_error = TRUE, verbose=TRUE)
})

zm_runoff <- raster::zonal(w_rs, zones1)
zm_runoff_smooth <- raster::zonal(test[[1]], zones1)
all.equal(zm_runoff, zm_runoff_smooth)

zm_runoff <- raster::zonal(w_rs*area(w_rs), zones1)
zm_runoff_smooth <- raster::zonal(test[[1]]*area(test[[1]]), zones1)
all.equal(zm_runoff, zm_runoff_smooth)

any(values(test[[1]]) < 0, na.rm=TRUE)
t <- values(test[[1]])
table(t < 0)

# writeRaster(test, "for_example.tif")
# 
# for(i in 1:5) {
#     system.time({
#         test <- pycnophylactic_smoothing(w_rs, wsm_rs, 
#                                          zones1, n=i, adjust_threshold = 0.5,
#                                          intensive = FALSE
#                                          return_error = TRUE, verbose=TRUE)
#     })
#     
#     writeRaster(test, paste0("itertest_",i,".tif"))
# }




system.time({
    test <- pycnophylactic_smoothing(w_rs, wsm_rs, zones1, n=5, 
                                     adjust_threshold = 0.5,
                                     return_error = TRUE, verbose=TRUE)
})

zm_runoff <- raster::zonal(w_rs, zones1)
zm_runoff_smooth <- raster::zonal(test[[1]], zones1)
all.equal(zm_runoff, zm_runoff_smooth)

zm_runoff <- raster::zonal(w_rs*area(w_rs), zones1)
zm_runoff_smooth <- raster::zonal(test[[1]]*area(test[[1]]), zones1)
all.equal(zm_runoff, zm_runoff_smooth)

any(values(test[[1]]) < 0, na.rm=TRUE)
t <- values(test[[1]])
table(t < 0)

writeRaster(test, "wateruse_5_iter_mountainslopes.tif")




system.time({
    test <- pycnophylactic_smoothing(runoff_rs, runoff_sm_rs, zones1, 
                                     n=5, 
                                     adjust_threshold = 0.5,
                                     return_error = TRUE, verbose=TRUE)
})

zm_runoff <- raster::zonal(runoff_rs, zones1)
zm_runoff_smooth <- raster::zonal(test[[1]], zones1)
all.equal(zm_runoff, zm_runoff_smooth)

zm_runoff <- raster::zonal(runoff_rs*area(w_rs), zones1)
zm_runoff_smooth <- raster::zonal(test[[1]]*area(test[[1]]), zones1)
all.equal(zm_runoff, zm_runoff_smooth)

any(values(test[[1]]) < 0, na.rm=TRUE)
t <- values(test[[1]])
table(t < 0)

writeRaster(test, "runoff__5_iter_mountainslopes.tif", overwrite=TRUE)
