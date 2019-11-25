---
title: "Smoodjustment"
author: "Marko Kallio"
date: "25 November 2019"
output:
      html_document:
        keep_md: true
---


## A quick demo of the function. Data not distributed 


The purpose of the function is to adjust smoothed surface so that it retains the original mass of the unsmoothed surface within specified zones.


the following shows how to smooth a layer while keeping the mean of the original layer within specified zones.



```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(raster)
```

```
## Loading required package: sp
```

```
## 
## Attaching package: 'raster'
```

```
## The following object is masked from 'package:dplyr':
## 
##     select
```

```r
library(rasterVis)
```

```
## Loading required package: lattice
```

```
## Loading required package: latticeExtra
```

```
## Loading required package: RColorBrewer
```

```r
library(knitr)
library(kableExtra) # nice html tables
```

```
## 
## Attaching package: 'kableExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     group_rows
```


### Define function

The function adjusts a smoothed raster based on the difference in means between the original raster and the smoothed raster). Pycnophylactic interpolation (PP) is applied to the difference, resulting in a smooth surface where of the errors, but where the mean error within the zone equals the original mean error. The smoothed input raster is then adjusted with the new surface. In practice, smoothing this way is a combination of the two surfaces: the original smoothed raster, and the smoothed error obtained from PP.
The function does not allow negative values in the smoothed raster. If they are present, they are set to zero. The function also keeps cells with value 0 as 0, and therefore any negative value are set to zero and kept there.

There is a current caveat in the function that it does not allow missing values within the specified zones (yet), or the function will result in error.


the inputs are:

* r1 = unsmoothed *RasterLayer*
* r2 = smoothed *RasterLayer*
* zones = rasterlayer defining zones for *raster::zonal()* function
* adjust_threshold = How much of the values in r2 is allowed to change in the adjustment. A value of 0.5 will allow 50% adjustment in r2 in order to restore balance with r1. The function checks whether this is possible - with small values, the zonal difference between r1 and r2 can be so high that it cannot be adjusted to match! 
* n = number of iterations for pycnophylactic interpolation
* smoothing_matrix = neighbourhood used in PP. Experimental, not recommended!
* intensive = TRUE (default) if the raster contains intensive (densitied, e.g. runoff depth), or extensive (counts, e.g. runoff volume) values.
* return_error =*TRUE* or *FALSE* whether to return the estimated error surface. If true, the function will return a *RasterStack* with the adjusted smooth layer, the original smooth layer, and the error layer.
* verbose = print progress indication?

function defitinion not shown - see pp_function.R for that.




### load wateruse data, and remove all values outside the zones. wateruse is resampled to the zones raster.

```r
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
```


### test water use, intensive

The function retains water balance (the volume of water within the zone), rather than the wateruse depth. This is because, if the input rasters cells are not rectangular, and in an equal area projection, the volume of water will change. 


```r
system.time({
    test <- pycnophylactic_smoothing(w_rs, wsm_rs, 
                                     zones1, n=5, adjust_threshold = 0.5,
                                     return_error = TRUE, verbose=TRUE)
})
```

```
## Preparing..
```

```
## Pycnophylactic interpolation of errors..
```

```
##   |                                                                         |                                                                 |   0%  |                                                                         |=============                                                    |  20%  |                                                                         |==========================                                       |  40%  |                                                                         |=======================================                          |  60%  |                                                                         |====================================================             |  80%  |                                                                         |=================================================================| 100%
```

```
## Preparing output..
```

```
##    user  system elapsed 
##  114.58   16.53  147.80
```

### test whether the wateruse depth is preserved - NO

```r
zm_runoff <- raster::zonal(w_rs, zones1)
zm_runoff_smooth <- raster::zonal(test[[1]], zones1)
all.equal(zm_runoff, zm_runoff_smooth)
```

```
## [1] "Mean relative difference: 0.001085118"
```

### test whether the volume of water use is preserved - YES

```r
zm_runoff <- raster::zonal(w_rs*area(w_rs), zones1)
zm_runoff_smooth <- raster::zonal(test[[1]]*area(test[[1]]), zones1)
all.equal(zm_runoff, zm_runoff_smooth)
```

```
## [1] TRUE
```

### test also for any negative values in the raster

```r
any(round(values(test[[1]]),10) < 0, na.rm=TRUE)
```

```
## [1] FALSE
```


### test water use, intensive

if the input water use was expressed in volume, we don't need to do extra work converting depth to volume to preserve water balance. Here, depth (as the imaginary volume) is preserved when intensive is set to FALSE.


```r
system.time({
    test <- pycnophylactic_smoothing(w_rs, wsm_rs, zones1, 
                                     intensive = FALSE, n=5, 
                                     adjust_threshold = 0.5,
                                     return_error = TRUE, verbose=TRUE)
})
```

```
## Preparing..
```

```
## Pycnophylactic interpolation of errors..
```

```
##   |                                                                         |                                                                 |   0%  |                                                                         |=============                                                    |  20%  |                                                                         |==========================                                       |  40%  |                                                                         |=======================================                          |  60%  |                                                                         |====================================================             |  80%  |                                                                         |=================================================================| 100%
```

```
## Preparing output..
```

```
##    user  system elapsed 
##   95.18   12.58  117.69
```

### Test whether depth is preserved - YES

```r
zm_runoff <- raster::zonal(w_rs, zones1)
zm_runoff_smooth <- raster::zonal(test[[1]], zones1)
all.equal(zm_runoff, zm_runoff_smooth)
```

```
## [1] TRUE
```

### Test whether volume (water balance) is preserved - NO

```r
zm_runoff <- raster::zonal(w_rs*area(w_rs), zones1)
zm_runoff_smooth <- raster::zonal(test[[1]]*area(test[[1]]), zones1)
all.equal(zm_runoff, zm_runoff_smooth)
```

```
## [1] "Mean relative difference: 0.0009633507"
```

## test for any negative values

```r
any(round(values(test[[1]]),10) < 0, na.rm=TRUE)
```

```
## [1] FALSE
```

## load runoff data and test that water balance is preserved - yes it is

```r
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
```


```r
system.time({
    test <- pycnophylactic_smoothing(runoff_rs, runoff_sm_rs, 
                                     zones1, n=5, 
                                     adjust_threshold = 0.5,
                                     return_error = TRUE, verbose=TRUE)
})
```

```
## Preparing..
```

```
## Pycnophylactic interpolation of errors..
```

```
##   |                                                                         |                                                                 |   0%  |                                                                         |=============                                                    |  20%  |                                                                         |==========================                                       |  40%  |                                                                         |=======================================                          |  60%  |                                                                         |====================================================             |  80%  |                                                                         |=================================================================| 100%
```

```
## Preparing output..
```

```
##    user  system elapsed 
##  106.87   13.14  127.22
```

```r
zm_runoff <- raster::zonal(runoff_rs*area(runoff_rs), zones1)
zm_runoff_smooth <- raster::zonal(test[[1]]*area(test[[1]]), zones1)
all.equal(zm_runoff, zm_runoff_smooth)
```

```
## [1] TRUE
```


