# smoothing with no negative runoff
library(raster)
library(gstat)
library(rasterVis)
library(dplyr)
library(Rcpp)
sourceCpp("func.cpp")



###############
# Create some artificial data

side <- 100
points <- data.frame(x = rep(1:side, each = side),
                     y = rep(1:side, times = side))


variogram <- vgm(psill=1, range=50, model='Exp')
spatial_model <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=variogram, nmax=20)

prediction <- predict(spatial_model, newdata=points, nsim=1)

gridded(prediction) = ~x+y

runoff <- raster(prediction) 
runoff <- runoff + abs(min(values(runoff)))
levelplot(runoff)

variogram <- vgm(psill=0.5, range=25, model='Sph')
spatial_model <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=variogram, nmax=20)

prediction <- predict(spatial_model, newdata=points, nsim=1)

gridded(prediction) = ~x+y

wateruse <- raster(prediction)
levelplot(wateruse)

#make zones
zones <- aggregate(runoff, 20)
values(zones) <- 1:ncell(zones)
zones <- disaggregate(zones, 20)
levelplot(zones)




#################
#################
#################
#################
r1 <- runoff
r1_nb <- matrix(1,9,9)


r2 <- raster::focal(r1, w=r1_nb, fun=mean, pad=TRUE, na.rm=TRUE)
# make sure some cells are negative
remove <- min(values(r2)) #+ 0.1
#remove <- 0.25
r1 <- r1 - remove

r2 <- raster::focal(r1, w=r1_nb, fun=mean, pad=TRUE, na.rm=TRUE)


#######
###### DONT SMOOTH ADDITIONALLY THOSE CELLS WHICH AREA ALREADY ZERO
###### ALTERNATIVELY CHECK HOW IS EQUAL DISTRIBUTION OF TRIM
#######
n <- 50
smoothing_matrix <- matrix(1,3,3)
pycnophylactic_smoothing <- function(r1, r2, 
                                     zones, 
                                     retain_zones,
                                     n=50, 
                                     smoothing_matrix = matrix(1,3,3),
                                     return_error = FALSE,
                                     verbose = FALSE) {
    
    if(verbose) message("Preparing..")
    
    # retain
    
    r_zonal_mean <- raster::zonal(r1, zones)
    rsm_zonal_mean <- raster::zonal(r2, zones)
    zonal_diff <- cbind(r_zonal_mean, sm_mean = rsm_zonal_mean[,2])
    zonal_diff <- cbind(zonal_diff, diff = zonal_diff[,3]-zonal_diff[,2])
    #zonal_diff <- tibble::as_tibble(zonal_diff)
    
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
    
    areas <- values(area(rzd))
    # zone <- values(zones)
    # r1val <- values(r1)
    # smooth_r2 <- values(r2)
    # orig_err_mean <- values(rzd)
    #r2_mat <- as.matrix(r2)
    adj_tmpl <- tibble::tibble(pixel = 1:ncell(zones),
                               zone = values(zones), 
                          r1 = values(r1),
                          smooth_r2 = values(r2),
                          #no_change = smooth_r2 <= threshold & smooth_r2 > 0,
                          no_change = smooth_r2 == 0,
                          negative = smooth_r2 < 0,
                          max_err = abs(smooth_r2 * threshold * !negative),
                          areas = areas,
                          orig_err_mean = values(rzd))
    
    nochange <- which(adj_tmpl$no_change)
    #nochange_val <- adj_tmpl$smooth_r2[nochange]
    
    
    pb <- txtProgressBar(min = 0, max = n, style=3) 
    for(i in 1:n) {
        
        # add boundary condition
        rzd_i[nochange] <- NA
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
        
        # adjust
        # smooth_mean <- as.vector(t(rzd_i))
        # mean_sm <- aggregate(smooth_mean, list(zone), FUN=mean)
        # colnames(mean_sm)[1] <- "zone"
        # 
        # mean_sm <- merge(data.frame(zone = zone), mean_sm, all.x=TRUE)
        # adjust = ifelse(orig_err_mean == 0,
        #                 1,
        #                 (mean_sm*areas)/(orig_err_mean*areas))
        # new_err = smooth_mean/adjust

        
         adj <- adj_tmpl %>%
             mutate(smooth_mean = as.vector(t(rzd_i))) %>%
            dplyr::group_by(zone) %>%
            #filter(zone == z) %>%
            dplyr::mutate(vol_sm = smooth_mean*areas, 
                          mean_sm = mean(smooth_mean),
                          # adjust = case_when(no_change ~ 1,
                          #                    !no_change ~ mean_sm/orig_err_mean),
                          adjust = ifelse(orig_err_mean == 0,
                                          1,
                                          (mean_sm*areas)/(orig_err_mean*areas)),
                                          #mean_sm/orig_err_mean),
                          new_err = smooth_mean/adjust,
                          high_err = abs(new_err) > abs(max_err),
                          test = any(high_err))
                          #adjusted_r2 = smooth_r2-new_err,
                          #test = round(mean(r1),5) == round(mean(adjusted_r2),5),
                          #test2 = round(mean(orig_err_mean),5) == round(mean(new_err),5))
        
        # adj %>% select(-areas,-vol_sm)
        test <- adj$high_err
        table(test)
        if(any(adj$high_err)) {
            # visited <- adj$high_err | adj$no_change
            # ni <- 0
            # keep <- TRUE
            
            #testg <- adj %>% filter(zone == 1000011328)
            uniq_zone <- unique(adj$zone)
            output <- list()
            print("working...")
            for(z in uniq_zone) {
                if(is.na(z)) {
                    order <- adj %>% 
                               filter(is.na(zone)) %>%
                               select(zone, pixel, smooth_r2, 
                                      max_err, new_err) %>%
                               arrange(smooth_r2) %>% 
                        mutate(fixed_err = new_err)
                    output <- append(output, list(order))
        
                    next
                } 
                cat(paste0(z, "..."))
                order <- adj %>% 
                    filter(zone == z) %>%
                    select(zone, pixel, smooth_r2, max_err, new_err) %>%
                    arrange(smooth_r2) 
                
                max_err <- order$max_err
                new_err <- order$new_err
                fixed_err <- rep(NA, length(new_err))
                
                total_dist <- 0
                pb2 <- txtProgressBar(min = 0, max = length(new_err), style=3)
                for(j in seq_along(new_err)) {
                    test <- abs(new_err[j]+total_dist) > abs(max_err[j])
                    if(is.na(test)) next
                    if(test) {
                        pos <- new_err[j] >= 0
                        
                        
                        if(pos) {
                            dist_err = new_err[j] + total_dist - max_err[j];
                            fix_err = max_err[j];
                        }  else {
                            dist_err = new_err[j] + total_dist + max_err[j];
                            fix_err = -max_err[j];
                        }
                        
                        npix <- length(new_err)-j
                        dist_err = dist_err/npix;
                        total_dist = total_dist + dist_err;

                        fixed_err[j] <- fix_err
                        
                    } else {
                        fixed_err[j] <- new_err[j]+total_dist
                    }
                    setTxtProgressBar(pb2, j)
                }
                #if(sum(fixed_err) )
                order$fixed_err <- fixed_err

                print(round(sum(fixed_err, na.rm=TRUE),5) == 
                          round(sum(new_err, na.rm=TRUE),5))
                print(!any(abs(fixed_err) > max_err))
                output <- append(output, list(order))
            }
            output <- do.call(rbind, output)
            
            
           
            # total_err <- rep(NA, length(new_err))
            # dist <- rep(NA, length(new_err))
            # 
            # testi <- fix_error(max_err, new_err, length(max_err))
            # summary(testi)
            # sum(testi)
            # table(abs(testi) > max_err)
            # 
            # new_err <- 5:-5
            # order <- data.frame(id = 1:length(new_err), err = abs(new_err))
            # order <- arrange(order, err)
            # new_err <- new_err[order$id]
            # max_err <- rep(2, length(new_err))
            # fixed_err <- rep(NA, length(new_err))
            # total_err <- rep(NA, length(new_err))
            # dist <- rep(NA, length(new_err))
            
            # total_dist <- 0
            # pb2 <- txtProgressBar(min = 0, max = length(new_err), style=3)
            # for(i in 1:length(new_err)) {
            #     test <- abs(new_err[i]+total_dist) > abs(max_err[i])
            #     if(test) {
            #         pos <- new_err[i] >= 0
            #         
            #         
            #         if(pos) {
            #             dist_err = new_err[i] + total_dist - max_err[i];
            #             fix_err = max_err[i];
            #         }  else {
            #             dist_err = new_err[i] + total_dist + max_err[i];
            #             fix_err = -max_err[i];
            #         }
            #         
            #         npix <- length(new_err)-i
            #         dist_err = dist_err/npix;
            #         total_dist = total_dist + dist_err;
            #         
            #         # ## define how much of error must be distributed
            #         # if(pos) {
            #         #     dist_err <- (new_err[i]+total_dist)-max_err[i]
            #         #     fix_err <- max_err[i]
            #         # }  else {
            #         #     dist_err <- (new_err[i]+total_dist)-max_err[i]
            #         #     fix_err <- -max_err[i]
            #         # }
            #         # 
            #         # # distribute to all remaining pixels
            #         # npix <- length(new_err)-i
            #         # dist_err <- ifelse(dist_err == 0,
            #         #                    0,
            #         #             if_else(npix == 0,
            #         #                     dist_err,
            #         #                     dist_err/npix))
            #         # total_dist <- total_dist + dist_err
            #         # total_err[i] <- total_dist
            #         # dist[i] <- dist_err
            #         # new_err[(i+1):nrow(order)] <- 
            #         #     new_err[(i+1):nrow(order)] + dist_err
            #         fixed_err[i] <- fix_err
            #         
            #     } else {
            #         fixed_err[i] <- new_err[i]+total_dist
            #     }
            #     setTxtProgressBar(pb2, i)
            # }
            # 
            # close(pb2)
            # #cbind(new_err, max_err, fixed_err)#, total_err, dist)
            # 
            # order$fixed_err <- fixed_err
            # order$dist <- dist
            # order$total_err <- total_err
            # summary(order)
            
            
            
            
            # # Create sample vector
            # X <- c(1:100); print(X)
            # 
            # # Create sample matrix
            # M <- matrix(c(1:100),nrow=10); print(M)
            # 
            # # Set limits
            # minV <- 15; maxV <- 85;
            # 
            # # Limit vector
            # sapply(X, function(y) min(max(y,minV),maxV))
            # 
            # # Limit matrix
            # apply(M, c(1, 2), function(x) min(max(x,minV),maxV))
            # 
            # 
            # while(keep && ni <= nx) {
            #     ni <- ni + 1
                
                
                
                
            
                # adj <- adj %>%
                #     dplyr::ungroup() %>%
                #     dplyr::mutate(negative = high_err,
                #                   visited = visited) %>% 
                #     dplyr::group_by(zone) %>% 
                #     dplyr::mutate(trim = adjusted_r2 * negative,
                #                   sum_trim = sum(trim, na.rm=TRUE),
                #                   new_error_nn = sum(new_err*!visited, na.rm=TRUE),
                #                   ratio = ifelse(visited,
                #                                  0,
                #                                  new_err/new_error_nn),
                #                                  #new_err/sum(!negative)),
                #                   err = ifelse(negative,
                #                                 trim,
                #                                 -sum_trim*ratio),
                #                   trim_ratio = sum_trim*ratio,
                #                   aep = new_err >= 0,
                #                   new_err = case_when(negative && aep ~ new_err+err,
                #                                       negative && !aep ~ new_err-err,
                #                                       !negative && aep ~ new_err-err,
                #                                       !negative && !aep ~ new_err+err),
                #                   high_err = new_err > smooth_r2,
                #                   adjusted_r2 = smooth_r2-new_err,
                #                   test3 = round(mean(r1),10) == round(mean(adjusted_r2),10))#,
                # 
                # adj_r2 <- adj$high_err
                # table(adj_r2)
                # table(adj$test3)
                # new_test <- any(adj_r2)
                # 
                # if(!new_test) { 
                #     keep <- FALSE
                # } else {
                #     visited <- visited | adj_r2
                #     table(visited)
                # }
                # print(ni)
                
            }
            
    
            # adjusted_r2 = smooth_r2-new_err,
            # test = mean(new_err) - mean(new_err)) 
        }
        
        rzd_i <- matrix(dplyr::pull(adj, new_err), nrow=nrow, ncol=ncol,
                        byrow=TRUE)
        if(verbose) setTxtProgressBar(pb, i)
    } 
    if(verbose) close(pb)
    
    message("Preparing output..")
    values(rzd) <- rzd_i 
    sm_out <- r2 - rzd
    if(return_error) {
        out <- raster::stack(sm_out, r2, rzd)
        names(out) <- c("Adjusted", "Original", "Error")
        return(out)
    } else {
        return(sm_out)
    }
}



test <- pycnophylactic_smoothing(r1,r2,zones, return_error=TRUE)

test <- out
zm_runoff <- raster::zonal(r1, zones)
zm_runoff_smooth <- raster::zonal(test[[1]], zones)
all.equal(zm_runoff, zm_runoff_smooth)





## test with real data

runoff <- raster("data/runoff_original_mm.tif")
runoff_smoothed <- raster("data/runoff_smoothed_mm.tif")
zones1 <- raster("data/zones_basins10000km2_mountainslowlands.tif")
crs(zones1) <- crs(runoff)
zones2 <- raster("data/zones_seabasins_mountainslowlands.tif")

#zones1_rs <- resample(zones1, runoff, method="ngb")
runoff_rs <- resample(runoff, zones1, method="ngb")
runoff_sm_rs <- resample(runoff_smoothed, zones1, method="ngb")

nas <- is.na(values(runoff_rs))

#values(runoff_rs) <- values(runoff_rs) * nas

val <- values(zones1)
val[nas] <- NA
values(zones1) <- val


system.time({
    test <- pycnophylactic_smoothing(runoff_rs, runoff_sm_rs, 
                                     zones1, n=75, 
                                     return_error = TRUE, verbose=TRUE)
})

zm_runoff <- raster::zonal(runoff_rs, zones1)
zm_runoff_smooth <- raster::zonal(test[[1]], zones1)
all.equal(zm_runoff, zm_runoff_smooth)

any(round(values(test[[1]]),10) < 0, na.rm=TRUE)
t <- values(test[[1]])
table(t < 0)





## test with real data

w <- raster("data/wateruse_original_mm.tif")
wsm <- raster("data/wateruse_smoothed_mm.tif")

zones1 <- raster("data/zones_basins10000km2_mountainslowlands.tif")
crs(zones1) <- crs(w)
zones2 <- raster("data/zones_seabasins_mountainslowlands.tif")
crs(zones2) <- crs(w)
# zones_sf <- sf::read_sf("data/mountainslopes_vector.gpkg")
# zones <- fasterize::fasterize(zones_sf, zones1, field="id")



#zones1_rs <- resample(zones1, runoff, method="ngb")
w_rs <- resample(w, zones1, method="ngb")
wsm_rs <- resample(wsm, zones1, method="ngb")

nas <- is.na(values(w_rs))
nas2 <- is.na(values(wsm_rs))
nas3 <- is.na(values(zones1))
nas <- nas | nas2 | nas3
#values(runoff_rs) <- values(runoff_rs) * nas

val <- values(zones1)
val[nas] <- NA
values(zones1) <- val

val <- values(w_rs)
val[nas] <- NA
values(w_rs) <- val

val <- values(wsm_rs)
val[nas] <- NA
values(wsm_rs) <- val

table(is.na(values(zones1)) == is.na(values(w_rs)))
table(is.na(values(zones1)) == is.na(values(wsm_rs)))

r1 <- w_rs
r2 <- wsm_rs
zones <- zones1
zones <- zones2
verbose <- TRUE
smoothing_matrix <- matrix(1,3,3)
n = 10
adjust_threshold <- 0.5

system.time({
    test <- pycnophylactic_smoothing(w_rs, wsm_rs, zones1, 
                                     adjust_threshold = 0.35, n=5, 
                                     return_error = TRUE, verbose=TRUE)
})

zm_runoff <- raster::zonal(w_rs, zones1)
zm_runoff_smooth <- raster::zonal(test[[1]], zones1)
all.equal(zm_runoff, zm_runoff_smooth)

any(values(test[[1]]) < 0, na.rm=TRUE)
t <- values(test[[1]])
table(t < 0)

# r1 <- runoff_rs
# r2 <- runoff_sm_rs
# zones <- zones1
# n <- 10
# verbose=TRUE







#######
###### DONT SMOOTH ADDITIONALLY THOSE CELLS WHICH AREA ALREADY ZERO
###### ALTERNATIVELY CHECK HOW IS EQUAL DISTRIBUTION OF TRIM
#######
n <- 50
smoothing_matrix <- matrix(1,3,3)
pycnophylactic_smoothing <- function(r1, r2, 
                                     zones, 
                                     adjust_threshold,
                                     n=50, 
                                     smoothing_matrix = matrix(1,3,3),
                                     return_error = FALSE,
                                     verbose = FALSE) {
    
    if(verbose) message("Preparing..")
    
    # retain
    
    r_zonal_mean <- raster::zonal(r1, zones)
    rsm_zonal_mean <- raster::zonal(r2, zones)
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
    
    areas <- values(area(rzd))
    adj_tmpl <- tibble::tibble(pixel = 1:ncell(zones),
                               zone = values(zones), 
                               r1 = values(r1),
                               smooth_r2 = values(r2),
                               no_change = smooth_r2 == 0,
                               negative = smooth_r2 < 0,
                               max_err = abs(smooth_r2 * 
                                                 adjust_threshold * 
                                                 !negative),
                               areas = areas,
                               orig_err_mean = values(rzd))
    
    nochange <- matrix(adj_tmpl$no_change, nrow=nrow, ncol=ncol,
                       byrow=TRUE)
    
    
    
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
            #filter(zone == z) %>%
            dplyr::mutate(vol_sm = smooth_mean*areas, 
                          mean_sm = mean(smooth_mean, na.rm=TRUE),
                          adjust = ifelse(orig_err_mean == 0,
                                          1,
                                          (mean_sm*areas)/(orig_err_mean*areas)),
                          #mean_sm/orig_err_mean),
                          new_err = smooth_mean/adjust,
                          high_err = abs(new_err) > abs(max_err),
                          test = any(high_err))
        #adjusted_r2 = smooth_r2-new_err,
        #test = round(mean(r1),5) == round(mean(adjusted_r2),5),
        #test2 = round(mean(orig_err_mean),5) == round(mean(new_err),5))
        
        # adj %>% select(-areas,-vol_sm)
        test <- adj$high_err
        iterr <- adj$new_err
        if(any(adj$high_err)) {
            # visited <- adj$high_err | adj$no_change
            # ni <- 0
            # keep <- TRUE
            
            #testg <- adj %>% filter(zone == 1000011328)
            uniq_zone <- unique(adj$zone)
            output <- list()
            #print("working...")
            for(z in uniq_zone) {
                if(is.na(z)) {
                    order <- adj %>% 
                        filter(is.na(zone)) %>%
                        select(zone, pixel, smooth_r2, 
                               max_err, new_err) %>%
                        arrange(smooth_r2) %>% 
                        mutate(fixed_err = new_err)
                    output <- append(output, list(order))
                    
                    next
                } 
               # cat(paste0(z, "..."))
                order <- adj %>% 
                    filter(zone == z) %>%
                    select(zone, pixel, smooth_r2, max_err, new_err) %>%
                    arrange(smooth_r2) 
                
                max_err <- order$max_err
                new_err <- order$new_err
                fixed_err <- rep(NA, length(new_err))
                shares <- max_err / sum(abs(max_err))
                share_left <- sum(shares)
                shares_visited <- 0
                
                # recte <- rep(NA, length(new_err))
                # recsh <- rep(NA, length(new_err))
                # recsherr <- rep(NA, length(new_err))
                # recshh <- rep(NA, length(new_err))
                # rectest1 <- rep(NA, length(new_err))
                # rectest2 <- rep(NA, length(new_err))
                # recdist <- rep(NA, length(new_err))
                # rect <-  rep(NA, length(new_err))
                # 
                total_err <- 0
                #pb2 <- txtProgressBar(min = 0, max = length(new_err), style=3)
                for(j in seq_along(new_err)) {
                    share_of_err <- shares[j] / (1-shares_visited)
                    sherr <- total_err*share_of_err
                    test <- abs(new_err[j]+sherr) > abs(max_err[j])
                    #rect[j] <- test
                    #if(is.na(test)) next
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
                       #share <- shares[j]
                       total_err <- total_err + dist_err
                       
                       # npix <- length(new_err)-j
                       #dist_err = dist_err/npix;
                       #total_err = total_err + dist_err;
                        
                        fixed_err[j] <- fix_err
                        # 
                        # recdist[j] <- dist_err
                        #rectest1[j] <- fix_err
                        #rectest2[j] <-  sum(fixed_err[1:j])
                        
                    } else {
                        fixed_err[j] <- new_err[j] + sherr
                        total_err <- total_err - sherr
                        shares_visited <- shares_visited + shares[j]
                        
                        # recdist[j] <- -sherr
                    }
                    #j <- j+1
                   # recte[j] <- total_err
                   # recsh[j] <- shares_visited
                   # recsherr[j] <- sherr
                   # recshh[j] <- share_of_err
                   
    
                   
    
                }
                
                order$fixed_err <- fixed_err
                # order$te <- recte
                # order$sh <- recsh
                # order$sherr <- recsherr
                # order$shoferr <- recshh
                # #order$test1 <- rectest1
                # #order$test2 <- rectest2
                # order$dist <- recdist
                # order$fix <- rectest1
                # cbind(new_err, fixed_err+recdist, round(new_err,5) == round((fixed_err + recdist),5), rect, rectest1)
                
                
                
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
                                     zones1, n=5, 
                                     adjust_threshold = 0.5,
                                     intensive = TRUE,
                                     return_error = TRUE, verbose=TRUE)
})

zm_runoff <- raster::zonal(w_rs, zones1)
zm_runoff_smooth <- raster::zonal(test[[1]], zones1)
all.equal(zm_runoff, zm_runoff_smooth)

any(values(test[[1]]) < 0, na.rm=TRUE)
t <- values(test[[1]])
table(t < 0)

writeRaster(test, "vikatesti3.tif")
            
system.time({
    test2 <- pycnophylactic_smoothing(w_rs, wsm_rs, 
                                     zones1, n=5, adjust_threshold = 0.5,
                                     return_error = TRUE, verbose=TRUE)
})

writeRaster(test2, "vikatesti6.tif")

zm_runoff <- raster::zonal(w_rs, zones1)
zm_runoff_smooth <- raster::zonal(test2[[1]], zones1)
all.equal(zm_runoff, zm_runoff_smooth)

any(values(test[[1]]) < 0, na.rm=TRUE)
t <- values(test[[1]])
table(t < 0)
