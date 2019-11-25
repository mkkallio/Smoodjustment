# set.seed(20191106)
library(raster)
library(gstat)
library(rasterVis)


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
levelplot(runoff)

variogram <- vgm(psill=0.5, range=25, model='Sph')
spatial_model <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=variogram, nmax=20)

prediction <- predict(spatial_model, newdata=points, nsim=1)

gridded(prediction) = ~x+y

wateruse <- raster(prediction)
levelplot(wateruse)



###############
###############
# Smoothing

# get the difference 
difference <- runoff-wateruse
levelplot(difference)
diff_mean <- focal(difference, w=matrix(1,21,21), fun=mean, pad=TRUE, na.rm=TRUE)
levelplot(diff_mean)


runoff_smooth <- focal(runoff, w=matrix(1,9,9), fun=mean, pad=TRUE, na.rm=TRUE)
levelplot(runoff_smooth)
wateruse_smooth <- focal(wateruse, w=matrix(1,5,5), fun=mean, pad=TRUE, na.rm=TRUE)
levelplot(wateruse_smooth)

difference_smooth <- runoff_smooth-wateruse_smooth
levelplot(difference_smooth)
# diff_smooth_mean <- focal(difference_smooth, w=matrix(1,21,21), 
#                           fun=mean, pad=TRUE, na.rm=TRUE)
# levelplot(diff_smooth_mean)


# diff_diff <- diff_mean - diff_smooth_mean
# summary(values(diff_diff))
# levelplot(diff_diff)

blockdiff <- aggregate(difference_smooth, 20, fun=mean, na.rm=TRUE)-
             aggregate(difference, 20, fun=mean, na.rm=TRUE)

zones <- tibble(zone = 1:ncell(blockdiff), deviation = values(blockdiff))
zonerast <- blockdiff
values(zonerast) <- 1:ncell(blockdiff)
zonerast <- disaggregate(zonerast, 20)
blockdiff <- disaggregate(blockdiff,20)

#pycnophylactic smoothing 
n <- 100
sr <- blockdiff


for(i in 1:n) {
    # smooth
    sr <- focal(sr, w=matrix(1,3,3), fun=mean, pad=TRUE, na.rm=TRUE)
    
    # adjust
    sz <- tibble(zone = values(zonerast), smooth_deviation = values(sr)) %>%
        left_join(zones, by="zone") %>%
        group_by(zone) %>%
        mutate(mean_sm = mean(smooth_deviation),
               adjust = ifelse(deviation == 0, 1, mean_sm/deviation),
               adjusted = smooth_deviation/adjust,
               mean_adjusted = mean(adjusted),
               equal = deviation == mean_adjusted) 
    
    values(sr) <- pull(sz, adjusted)
    
} 
sz %>% summarize(meandev = mean(deviation),
                 meanadj = mean(mean_adjusted))
# check all good


## determine how to adjust smoothed runoff and wateruse

#did smoothing result in higher or lower runoff/wateruse?
rd <- runoff_smooth-runoff
wd <- wateruse_smooth-wateruse


# tbl <- tibble(pixel = 1:ncell(rd), 
#               zone = values(zonerast),
#               rd = values(rd), 
#               wd = values(wd), 
#               adjust_amount = values(sr),
#               rsm = values(runoff_smooth),
#               wsm = values(wateruse_smooth),
#               diff = values(difference)) %>% 
#     mutate(rd_h = rd >=0, 
#            wd_h = wd >=0,  
#            bd_h = rd_h == wd_h,
#            rdwd = abs(rd) >= abs(wd),
#            sr_h = adjust_amount >= 0,
#            # if rd is higher than wd, allocate sr to rd, else wd
#            adj_rd = case_when(bd_h & rdwd ~ FALSE,
#                               bd_h & !rdwd ~ TRUE,
#                               !bd_h & rdwd ~ TRUE,
#                               !bd_h & !rdwd ~FALSE),
#            r_adj = adj_rd*adjust_amount*ifelse(rd_h, -1, 1)*ifelse(bd_h,1,-1),
#            w_adj = (!adj_rd)*adjust_amount*ifelse(wd_h, -1, 1)*ifelse(bd_h,1,-1),
#            diff_sm = rsm-wsm,
#            diff_adj = (rsm + r_adj)-(wsm + w_adj)) %>%
#     select(-rd_h, -wd_h, -bd_h, -rdwd, -sr_h) %>% print(n=30)
#     group_by(zone) %>%
#     summarise(meandiff = mean(diff_sm),
#               meandiff_adj = mean(diff_adj))

    tbl <- tibble(pixel = 1:ncell(rd), 
                  zone = values(zonerast),
                  rd = values(rd), 
                  wd = values(wd), 
                  adjust_amount = values(sr),
                  rsm = values(runoff_smooth),
                  wsm = values(wateruse_smooth),
                  rsm_wsm = rsm-wsm,
                  diff_sm = values(difference_smooth),
                  diff_orig = values(difference)) %>% 
        mutate(adjust_closer = adjust_amount < 0,
               rd_wd = rsm >= wsm, 
               #both_same = rd_pos == wd_pos,
               adjust_r = ifelse(adjust_closer,
                                 case_when(rd_wd ~ rsm-adjust_amount*0.5,
                                           !rd_wd ~ rsm+adjust_amount*0.5),
                                 case_when(rd_wd ~ rsm-adjust_amount*0.5,
                                           !rd_wd ~ rsm+adjust_amount*0.5)),
               adjust_w = ifelse(adjust_closer,
                                 case_when(rd_wd ~ wsm+adjust_amount*0.5,
                                           !rd_wd ~ wsm-adjust_amount*0.5),
                                 case_when(rd_wd ~ wsm+adjust_amount*0.5,
                                           !rd_wd ~ wsm-adjust_amount*0.5)),
               new_r = rsm + adjust_r,
               new_w = wsm + adjust_w) %>% 
        print(n=30)
    
    tbl %>%
        group_by(zone) %>%
        summarise(orig_diff = mean(diff_orig),
                  sm_diff = mean(diff_sm),
                  adjust_amount = mean(adjust_amount),
                  adjusted_diff = mean(adjust_r-adjust_w),
                  test = round(orig_diff,10) == round(adjusted_diff,10))
    
        # tbl <- tibble(pixel = 1:ncell(rd), 
        #               zone = values(zonerast),
        #               rd = values(rd), 
        #               wd = values(wd), 
        #               #rdwd = rd-wd,
        #               adjust_amount = values(sr),
        #               rsm = values(runoff_smooth),
        #               wsm = values(wateruse_smooth),
        #               rsm_wsm_diff = rsm-wsm,
        #               diff_orig = values(difference)) %>% 
        #     group_by(zone) %>%
        #     mutate(r_mult = 1+(rd/sum(rd)),
        #            w_mult = 1+(wd/sum(wd)),
        #            r_adj = adjust_amount*r_mult,
        #            w_adj = adjust_amount*w_mult,
        #            new_r = rsm+r_adj,
        #            new_w = wsm+r_adj) %>%
        #     select(-pixel,-zone) %>% 
        #     summarise(meandiff = mean(diff_orig),
        #               meandiff_sm = mean(rsm-wsm),
        #               meanadjust = mean(adjust_amount)) %>%
        #     print(n=20)

trun <- tuse <- runoff
values(trun) <- tbl$adjust_r
values(tuse) <- tbl$adjust_w

diff <- trun-tuse

temp <- tibble(zone = values(zonerast), 
       ajd_sm = values(diff), 
       diff = values(difference)) %>%
    group_by(zone) %>%
    summarise_all(mean) %>%
    mutate(d= ajd_sm-diff) %>%
    print


# adj_rh srh rdh wdh +-
# r T T F +  
# r T T T -
# r T F T +
# r T F F -
#     
# r F T F -
# r F T T +
# r F F T -
# r F F F +


# checking
adjzonal <- zonal(addLayer(trun,tuse), zonerast) %>%
    as_tibble %>%
    left_join(zones) %>%
    mutate(simdiff = sim1.1-sim1.2)


zones <- zonerast
interpol <- function(r1, 
                     r2, 
                     r1_nb = matrix(1,3,3), 
                     r2_nb = matrix(1,3,3), 
                     n=100, 
                     zones) {
    
    # smooth r1 and r2 with the neighbourhood matrices
    r1_sm <- raster::focal(r1, w=r1_nb, fun=mean, pad=TRUE, na.rm=TRUE)
    r2_sm <- raster::focal(r2, w=r2_nb, fun=mean, pad=TRUE, na.rm=TRUE)
    
    # calculate differences
    diff_r1r2 <- r1-r2
    diff_r1r2_sm <- r1_sm-r2_sm
    diff_r1 <- r1_sm-r1
    diff_r2 <- r2_sm-r2
    
    # differences in given zones
    zonal_diff <- raster::zonal(diff_r1r2, zones)
    zonal_diff_sm <- raster::zonal(diff_r1r2_sm, zones)
    zd <- cbind(zonal_diff, mean_sm = zonal_diff_sm[,2])
    zd <- cbind(zd, diff = zd[,3]-zd[,2])
    # make raster out of zonal differences
    zonal_diff <- raster::reclassify(zones, zd[,c(1,4)])
    
    # pycnophylactic interpolation for difference in original and smoothed
    # differences. Iterate n times.
    zonal_diff_i <- zonal_diff
    for(i in 1:n) {
        # smooth
        zonal_diff_i <- focal(zonal_diff_i, 
                            w=matrix(1,3,3), 
                            fun=mean, 
                            pad=TRUE, 
                            na.rm=TRUE)
        
        # adjust
        temp <- tibble(zone = values(zones), 
                     smooth_deviation = values(zonal_diff_i),
                     orig_deviation = values(zonal_diff)) %>%
            group_by(zone) %>%
            mutate(mean_sm = mean(smooth_deviation),
                   adjust = ifelse(orig_deviation == 0, 
                                   1, 
                                   mean_sm/orig_deviation),
                   adjusted = smooth_deviation/adjust,
                   mean_adjusted = mean(adjusted),
                   equal = round(orig_deviation,10) == round(mean_adjusted,10)) 
        
        values(zonal_diff_i) <- pull(temp, adjusted)
    } 
    
    # adjust the smoothed r1 and r2 so that within given zones, their difference
    # remains equal
    tbl <- tibble(pixel = 1:ncell(rd), 
                  zone = values(zones),
                  rd = values(diff_r1), 
                  wd = values(diff_r2), 
                  adjust_amount = values(zonal_diff_i),
                  r1sm = values(r1_sm),
                  r2sm = values(r2_sm),
                  rdiff = r1sm-r2sm,
                  diff_sm = values(diff_r1r2_sm),
                  diff_orig = values(diff_r1r2)) %>% 
        mutate(adjust_closer = adjust_amount < 0,
               r1_r2 = r1sm >= r2sm, 
               rdiff_pos = rdiff >= 0,
               adjust_amount = ifelse(rdiff_pos, adjust_amount, -adjust_amount),
               adjusted_r1 = ifelse(adjust_closer,
                                 case_when(r1_r2 ~ r1sm-adjust_amount*0.5,
                                           !r1_r2 ~ r1sm+adjust_amount*0.5),
                                 case_when(r1_r2 ~ r1sm-adjust_amount*0.5,
                                           !r1_r2 ~ r1sm+adjust_amount*0.5)),
               adjusted_r2 = ifelse(adjust_closer,
                                 case_when(r1_r2 ~ r2sm+adjust_amount*0.5,
                                           !r1_r2 ~ r2sm-adjust_amount*0.5),
                                 case_when(r1_r2 ~ r2sm+adjust_amount*0.5,
                                           !r1_r2 ~ r2sm-adjust_amount*0.5))) %>% 
               group_by(zone) %>%
               mutate(orig_diff = mean(diff_orig),
                      sm_diff = mean(diff_sm),
                      adjust_amount = mean(adjust_amount),
                      adjusted_diff = mean(adjusted_r1-adjusted_r2),
                      test = round(orig_diff,10) == round(adjusted_diff,10))
               # 
    r1_out <- r1
    values(r1_out) <- tbl$adjusted_r1
    r2_out <- r2
    values(r2_out) <- tbl$adjusted_r2
    
    out <- addLayer(r1_out, r2_out)
    names(out) <- c("r1_smoothed", "r2_smoothed")
    
    return(out)
    
}

smoothed <- interpol(runoff, wateruse, 
                     r1_nb = matrix(1,9,9), 
                     r2_nb = matrix(1,3,3),
                     n=100,
                     zones)

## test everythings alright
original_difference <- raster::zonal(runoff-wateruse, zones)
smoothed_difference <- raster::zonal(smoothed[[1]]-smoothed[[2]], zones)

all.equal(original_difference, smoothed_difference)

differences <- cbind(original_difference, 
                     mean_sm = smoothed_difference[,2],
                     diff = original_difference[,2]-smoothed_difference[,2])






# smooth r1 and r2 with the neighbourhood matrices
r1_sm <- raster::focal(r1, w=r1_nb, fun=mean, pad=TRUE, na.rm=TRUE)
diff_r1 <- r1_sm-r1

r_zonal_mean <- raster::zonal(r1, zones)
rsm_zonal_mean <- raster::zonal(r1_sm, zones)
zonal_diff <- cbind(r_zonal_mean, sm_mean = rsm_zonal_mean[,2])
zonal_diff <- cbind(zonal_diff, diff = zonal_diff[,3]-zonal_diff[,2])
#zonal_diff <- tibble::as_tibble(zonal_diff)

rzd <- raster::reclassify(zones, zonal_diff[,c(1,4)])

rzd_i <- as.matrix(rzd)
ncol <- ncol(rzd)
nrow <- nrow(rzd)
# boundary condition
rleft <- rzd_i[,1]
rright <- rzd_i[,ncol]
rtop <- c(NA, rzd_i[1,], NA)
rbottom <- c(NA, rzd_i[nrow,], NA)
for(i in 1:n) {
    # add boundary condition
    rzd_i <- cbind(rleft, rzd_i, rright)
    rzd_i <- rbind(rtop, rzd_i, rbottom)
    
    rzd_i <- raster::raster(rzd_i)
    rzd_i <- raster::focal(rzd_i, 
                          w=matrix(1,3,3), 
                          fun=mean, 
                          pad=TRUE, 
                          na.rm=TRUE)
    rzd_i <- as.matrix(rzd_i)[c(-1,-(nrow+1)), c(-1,-(ncol+1))]
    
    # adjust
    temp <- tibble::tibble(zone = values(zones), 
                   smooth_mean = as.vector(t(rzd_i)),
                   orig_mean = values(rzd)) %>%
        dplyr::group_by(zone) %>%
        dplyr::mutate(mean_sm = mean(smooth_mean),
               adjust = ifelse(orig_mean == 0, 
                               1, 
                               mean_sm/orig_mean),
               adjusted = smooth_mean/adjust,
               mean_adjusted = mean(adjusted),
               equal = round(orig_mean,10) == round(mean_adjusted,10)) 
    
    rzd_i <- matrix(dplyr::pull(temp, adjusted), nrow=nrow, ncol=ncol,
                    byrow=TRUE)
} 
values(rzd) <- rzd_i 

temp <- raster::zonal(rzd, zones)

# adjust the smoothed r1 and r2 so that within given zones, their difference
# remains equal
tbl <- dplyr::tibble(zone = values(zones),
                     r = values(r1),
                     sm = values(r1_sm),
                     diff = sm-r,
                     adjust_amount = values(rzd),
                     adjusted = sm - adjust_amount) #%>%
    # dplyr::group_by(zone) %>%
    # mutate(mean_r = mean(r),
    #        mean_adj = mean(adjusted),
    #        test = round(mean_r,10) == round(mean_adj,10))


sm_out <- r1
raster::values(sm_out) <- tbl$adjusted


out <- raster::addLayer(r1_out, r2_out)
names(out) <- c("r1_smoothed", "r2_smoothed")

return(out)



pycnophylactic_smoothing <- function(r1, r2, n=50, zones) {

    diff_r1 <- r2-r1
    
    r_zonal_mean <- raster::zonal(r1, zones)
    rsm_zonal_mean <- raster::zonal(r2, zones)
    zonal_diff <- cbind(r_zonal_mean, sm_mean = rsm_zonal_mean[,2])
    zonal_diff <- cbind(zonal_diff, diff = zonal_diff[,3]-zonal_diff[,2])
    #zonal_diff <- tibble::as_tibble(zonal_diff)
    
    rzd <- raster::reclassify(zones, zonal_diff[,c(1,4)])
    
    rzd_i <- as.matrix(rzd)
    ncol <- ncol(rzd)
    nrow <- nrow(rzd)
    # boundary condition
    rleft <- rzd_i[,1]
    rright <- rzd_i[,ncol]
    rtop <- c(NA, rzd_i[1,], NA)
    rbottom <- c(NA, rzd_i[nrow,], NA)
    for(i in 1:n) {
        # add boundary condition
        rzd_i <- cbind(rleft, rzd_i, rright)
        rzd_i <- rbind(rtop, rzd_i, rbottom)
        
        rzd_i <- raster::raster(rzd_i)
        rzd_i <- raster::focal(rzd_i, 
                               w=matrix(1,3,3), 
                               fun=mean, 
                               pad=TRUE, 
                               na.rm=TRUE)
        rzd_i <- as.matrix(rzd_i)[c(-1,-(nrow+1)), c(-1,-(ncol+1))]
        
        # adjust
        temp <- tibble::tibble(zone = values(zones), 
                               smooth_mean = as.vector(t(rzd_i)),
                               orig_mean = values(rzd)) %>%
            dplyr::group_by(zone) %>%
            dplyr::mutate(mean_sm = mean(smooth_mean),
                          adjust = ifelse(orig_mean == 0, 
                                          1, 
                                          mean_sm/orig_mean),
                          adjusted = smooth_mean/adjust,
                          mean_adjusted = mean(adjusted),
                          equal = round(orig_mean,10) == round(mean_adjusted,10)) 
        
        rzd_i <- matrix(dplyr::pull(temp, adjusted), nrow=nrow, ncol=ncol,
                        byrow=TRUE)
    } 
    values(rzd) <- rzd_i 
    
    
    sm_out <- r2 - rzd
    return(sm_out)
}
    
 
## testing
test_runoff <- pycnophylactic_smoothing(r1, r1_sm, n=50, zones)
test_wateruse <- pycnophylactic_smoothing(r2, r2_sm, n=50, zones)

zd <- zonal(r1-r2, zones)
zd2 <- zonal(test_runoff - test_wateruse, zones)
all.equal(zd,zd2)




test_iter <- list()
ns <- c(1, 5, seq(25,250,by=25))
for(i in seq_along(ns)) {
    n <- ns[i]
    test2 <- pycnophylactic_smoothing(runoff, 
                                      runoff_smooth, 
                                      zones, 
                                      n=n,
                                      smoothing_matrix = matrix(1,3,3),
                                      return_error=TRUE)
    test_iter[[i]] <- test2$Error
}
test_iter <- do.call(stack, test_iter)
#names(test_iter) <- paste("n_iter_",ns)
# levelplot(test_iter, 
#           names.attr=paste("n_iter_",ns))

nl <- nlayers(test_iter)
m <- matrix(1:nl, nrow=3)
#themes <- list(RdBuTheme(), BTCTheme(), GrTheme(), PuOrTheme())
for (i in 1:nl){
    p <- levelplot(test_iter, layers=i,
                   #par.settings=themes[[i]],
                   margin=FALSE,
                   main=paste("n_iter_",ns[i]))
    print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
}



test_matrix <- list()
side <- rep(seq(3,33,by=6),each=2)
ns <- rep(c(25,100), 6)
for(i in seq_along(side)) {
    s <- side[i]
    n <- ns[i]
    test2 <- pycnophylactic_smoothing(runoff, 
                                      runoff_smooth, 
                                      zones, 
                                      n=n,
                                      smoothing_matrix = matrix(1,s,s),
                                      return_error=TRUE)
    test_matrix[[i]] <- test2$Error
}
test_matrix <- do.call(stack, test_matrix)
#names(test_matrix) <- paste("dist_",side)
# levelplot(test_matrix, 
#           names.attr=paste("dist_",side))

    

nl <- nlayers(test_matrix)
m <- matrix(1:nl, nrow=3)
#themes <- list(RdBuTheme(), BTCTheme(), GrTheme(), PuOrTheme())
for (i in 1:nl){
    p <- levelplot(test_matrix, layers=i,
                   #par.settings=themes[[i]],
                   margin=FALSE,
                   main=paste("dist_",(side[i]-1)/2,", iter=",ns[i]))
    print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
}
