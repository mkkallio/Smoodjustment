% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pycnophylactic_adjustment.R
\name{pycnophylactic_adjustment}
\alias{pycnophylactic_adjustment}
\title{Adjust a smoothed raster to preserve mass of unsmoothed raster}
\usage{
pycnophylactic_adjustment(r1, r2, zones, adjust_threshold = 0.5, n = 5,
  smoothing_matrix = matrix(1, 3, 3), intensive = TRUE,
  return_error = FALSE, verbose = FALSE)
}
\arguments{
\item{r1}{unsmoothed \code{RasterLayer}, from which target zonal mass is
extracted.}

\item{r2}{smoothed \code{RasterLayer}, which' zonal mass we adjust}

\item{zones}{a \code{RasterLayer} with zones in which adjustment takes place}

\item{adjust_threshold}{How much the values of r2 are allowed to change in 
the adjustment. A value of 0.5 will allow a maximum adjustment of 50% in 
r2 in order to restore balance with r1. Within range [0,1].}

\item{n}{number of iterations for pycnophylactic adjustment}

\item{smoothing_matrix}{neighbourhood for pycnophylactic interpolation. 
Experimental, default 3x3 neighbourhood suggested.}

\item{intensive}{Are values r1 and r2 densities (\code{TRUE}) or counts 
(\code{FALSE}).}

\item{return_error}{If \code{TRUE}, returns the pycnophylactic adjustment 
layer in addition to the adjusted.}

\item{verbose}{Print progress indication?}
}
\value{
Volume corrected r2 as single \code{RasterLayer}, or , if 
  \code{return_error == TRUE}, a four-layer raster with corrected r2, 
  uncorrected r2, r1, and the adjustment raster.
}
\description{
The function adjusts a smoothed raster based on the difference in means 
between the original raster and the smoothed raster). Pycnophylactic 
interpolation (PP) is applied to the difference, resulting in a smooth 
surface, but where the mean difference within the zone 
equals the original mean error. The smoothed input raster is then adjusted 
with the new surface. In practice, smoothing this way is a combination of 
the two surfaces: the original smoothed raster, and the smoothed difference
obtained with PP. The function sets all negative values in the 
smoothed raster as 0. The function also keeps cells with value 0 as 0, and 
therefore any negative value are set to zero and kept as such.
}
