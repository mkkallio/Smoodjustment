---
title: "Smoodjustment 0.0.9999"
author: "Marko Kallio"
date: "5 December 2019"
output:
      html_document:
        keep_md: true
---

Smoothing a raster layer results in a change in zonal counts over. For some applications, this is not desirable. Smoodjustment adjusts a smoothed raster so that, within specified zones the mass of the original unsmoothed raster is preserved. 

In the most simple case, the adjustment can be made by increasing and lowering the values within a zone so that the two rasters' zonal mass match. However, this results in potential steps at the zone boundaries - a feature that smoothing is often used to eliminate. Smoodjustment smooths the adjustment layer using Tobler's pycnophylactic interpolation to create a smooth adjustment layer.

**Note: this is a very early version of the package designed for a specific purpose. Unexpected errors may occur with other uses.**




### Create some rasters 
We create a random raster to be smoothed and adjusted


```r
set.seed(5122019)
side <- 100 # size of one side
points <- data.frame(x = rep(1:side, each = side),
                     y = rep(1:side, times = side))
variogram <- vgm(psill=1, range=50, model='Exp')
spatial_model <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=variogram, nmax=20)
prediction <- predict(spatial_model, newdata=points, nsim=1)
```

```
## [using unconditional Gaussian simulation]
```

```r
gridded(prediction) = ~x+y
r <- raster(prediction) +5
crs(r) <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
levelplot(r)
```

![](README_files/figure-html/unnamed-chunk-1-1.png)<!-- -->
And lets create some zones


```r
zones <- aggregate(r, 20)
raster::values(zones) <- 1:ncell(zones)
zones <- disaggregate(zones, 20)
levelplot(zones)
```

![](README_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

### smooth raster and adjust

```r
r_smooth <- focal(r, w = matrix(1,9,9), fun='mean', pad=TRUE, na.rm=TRUE)
levelplot(r_smooth)
```

![](README_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
adjusted <- pycnophylactic_adjustment(r, r_smooth, zones, adjust_threshold = 0.5, n =50, intensive=FALSE, return_error = TRUE, verbose=FALSE)
```

```
## Preparing output..
```

```r
levelplot(adjusted)
```

![](README_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
levelplot(adjusted$Adjustment.raster)
```

![](README_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

### Check that the masses hold

```r
original <- raster::zonal(r, zones)
adj_smooth <- raster::zonal(adjusted$Adjusted.r2, zones)
all.equal(original, adj_smooth)
```

```
## [1] TRUE
```


```r
differences <- cbind(original, 
                     adj_smooth = adj_smooth[,2],
                     difference = original[,2]-adj_smooth[,2])
colnames(differences)[2] <- "original"

differences  %>%
    kable() %>%
    kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> zone </th>
   <th style="text-align:right;"> original </th>
   <th style="text-align:right;"> adj_smooth </th>
   <th style="text-align:right;"> difference </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6.202929 </td>
   <td style="text-align:right;"> 6.202929 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 6.355647 </td>
   <td style="text-align:right;"> 6.355647 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 6.493551 </td>
   <td style="text-align:right;"> 6.493551 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 6.300888 </td>
   <td style="text-align:right;"> 6.300888 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6.209866 </td>
   <td style="text-align:right;"> 6.209866 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6.767272 </td>
   <td style="text-align:right;"> 6.767272 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 6.612633 </td>
   <td style="text-align:right;"> 6.612633 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 5.903989 </td>
   <td style="text-align:right;"> 5.903989 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 5.478813 </td>
   <td style="text-align:right;"> 5.478813 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 5.353948 </td>
   <td style="text-align:right;"> 5.353948 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 5.783103 </td>
   <td style="text-align:right;"> 5.783103 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 5.898303 </td>
   <td style="text-align:right;"> 5.898303 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 5.778704 </td>
   <td style="text-align:right;"> 5.778704 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 4.948187 </td>
   <td style="text-align:right;"> 4.948187 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 5.097298 </td>
   <td style="text-align:right;"> 5.097298 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 5.053852 </td>
   <td style="text-align:right;"> 5.053852 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 4.888858 </td>
   <td style="text-align:right;"> 4.888858 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 4.810975 </td>
   <td style="text-align:right;"> 4.810975 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 3.708400 </td>
   <td style="text-align:right;"> 3.708400 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 3.755489 </td>
   <td style="text-align:right;"> 3.755489 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;"> 4.732661 </td>
   <td style="text-align:right;"> 4.732661 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 4.805980 </td>
   <td style="text-align:right;"> 4.805980 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 5.218134 </td>
   <td style="text-align:right;"> 5.218134 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 4.084446 </td>
   <td style="text-align:right;"> 4.084446 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 3.932795 </td>
   <td style="text-align:right;"> 3.932795 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>


