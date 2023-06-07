normalizeBRDF <- function(ref, band_names,
                          sun_zenith, sun_azimuth, view_zenith, view_azimuth,
                          BRDF_coeffs, 
                          sun_zenith_norm, lat_norm,
                          degrees=FALSE,
                          ...){
  ##  BRDF normalization for Landsat TM/ETM+ 
  ##  
  ##  Arguments:
  ##    ref:              SpatRaster. Reflectance values.
  ##    band_names:       character. Band names. Should be of length nlyr(ref). Optional if band names are provided in names(ref)
  ##    sun_zenith:       SpatRaster
  ##    sun_azimuth:      SpatRaster
  ##    view_zenith:      SpatRaster
  ##    view_azimuth:     SpatRaster
  ##    BRDF_coeffs:      character or list. Defines the BRDF correction coefficients
  ##                        If character, possible values are: 
  ##                          "Roy" (default): set of global parameters obtained from MODIS  (Roy, D. P., H. K. Zhang, J. Ju, J. L. Gomez-Dans, P. E. Lewis, C. B. Schaaf, Q. Sun, J. Li, H. Huang, and V. Kovalskyy. 2016. "A General Method to Normalize Landsat Reflectance Data to Nadir BRDF Adjusted Reflectance." Remote Sensing of Environment 176 (April): 255-71. https://doi.org/10.1016/j.rse.2016.01.023.)
  ##                          "Flood": set of parameters derived for Australia (Flood, N., T. Danaher, T. Gill, and S. Gillingham. 2013. "An Operational Scheme for Deriving Standardised Surface Reflectance from Landsat TM/ETM+ and SPOT HRG Imagery for Eastern Australia." Remote Sensing 5 (1): 83-109. https://doi.org/10.3390/rs5010083.)
  ##                          "Vandoninck": set of parameters derived for Amazonian forests (Van doninck, J., and H. Tuomisto. 2017. "Evaluation of Directional Normalization Methods for Landsat TM/ETM+ over Primary Amazonian Lowland Forests." International Journal of Applied Earth Observation and Geoinformation 58 (June): 249-63. https://doi.org/10.1016/j.jag.2017.01.017.)
  ##                        If list, should be in form list(f_vol_prime=c(...), f_geo_prime(...))
  ##    sun_zenith_norm:  numeric. Solar zenith angle to normalize to. Either this parameter or lat_norm should be provided.
  ##    lat_norm:         numeric. Latitude from which to derive solar zenith angle to normalize to. Either this parameter or sun_zenith_norm should be provided.
  ##    degrees:          logical. Are angles given in degrees instead of radians? 
  ##
  ##  Details:
  ##    Parameter "BRDF_coeffs" . 
  ##    Possible values: 
  ##
  
  library(terra)
  if(missing(ref) | missing(sun_zenith) | missing(sun_azimuth) | missing(view_zenith) | missing(view_azimuth)) stop("missing input data")
  
  ##  Assign band names to reflectance raster
  if(!missing(band_names)) {
    tmp_band_names <- names(ref)
    names(ref) <- band_names
  }
  
  ##  Normalized solar zenith angle
  if(!missing(lat_norm) & missing(sun_zenith_norm)){
    sun_zenith_norm <- 31.0076 - 
      0.1272*lat_norm + 0.01187*lat_norm^2 + 2.4*10^(-5)*lat_norm^3  - 
      9.48*10^(-7)*lat_norm^4 - 1.95*10^(-9)*lat_norm^5 + 6.15*10^(-11)*lat_norm^6
    sun_zenith_norm <- sun_zenith_norm*pi/180
  } else if(!missing(sun_zenith_norm)){
    if(degrees) sun_zenith_norm <- sun_zenith_norm*pi/180
  } else {
    stop("Must provide either sun_zenith_norm or lat_norm")
  }

  ##  Possible values of BRDF corrections coefficients
  if(missing(BRDF_coeffs)) BRDF_coeffs <- "roy"
  if (class(BRDF_coeffs)=="character"){
    BRDF_coeffs <- switch(tolower(BRDF_coeffs),
                          "roy"        = list(f_vol_prime=c(0.0372, 0.0580, 0.0574, 0.1535, 0.1154, 0.0639)/c(0.0774, 0.1306, 0.1690, 0.3093, 0.3430, 0.2658),
                                              f_geo_prime=c(0.0079, 0.0178, 0.0227, 0.0330, 0.0453, 0.0387)/c(0.0774, 0.1306, 0.1690, 0.3093, 0.3430, 0.2658)),
                          "flood"      = list(f_vol_prime=c(0.93125413991, 0.687401438519, 0.645033011917, 0.704036740665, 0.360201003097, 0.290061903555),
                                              f_geo_prime=c(0.260953557124, 0.213872135374, 0.180032152925, 0.093518142066, 0.162796996525, 0.147723009593)),  
                          "vandoninck" = list(f_vol_prime=c(1.2933, 0.5550, 0.9937, 0.3489, 1.0129, 1.2742),
                                              f_geo_prime=c(0.1902, 0.2456, 0.0000, 0.2755, 0.2348, 0.2426)),
                          stop('Incorrect value for "BRDF_coeffs"'))
  } else if (class(BRDF_coeffs)=="list") {
    #check whether list has correct format
  } else {
    stop('Incorrect value for "BRDF_coeffs"')
  }

  ##  Sun/view angles in degrees or radians?
  if(degrees){
    sun_zenith <- sun_zenith*pi/180
    sun_azimuth <- sun_azimuth*pi/180
    view_zenith <- view_zenith*pi/180
    view_azimuth <- view_azimuth*pi/180
  }
  
  kernels <- function(sunza, satza, relaz){
    ##  Kvol - RossThick kernel for volume scattering (Roujean et al., 1992)
    cos_pa <- cos(sunza)*cos(satza) + sin(sunza)*sin(satza)*cos(relaz)
    pa <- acos(max(-1, min(cos_pa,1)))
    kvol <- ((0.5*pi-pa)*cos(pa)+sin(pa))/(cos(sunza)+cos(satza))-(pi/4)
    
    ##  Kgeo - LiSparse-R kernel for geometric-optical surface scattering (Wanner et al., 1995)
    b <- 1
    #r <- 1
    h <- 2
    
    #cos_pa <- cos(sunza)*cos(satza) + sin(sunza)*sin(satza)*cos(relaz)
    #pa <- acos(max(-1,min(1,cos_pa)))
    
    D <- tan(sunza)^2 + tan(satza)^2 - 2*tan(sunza)*tan(satza)*cos(relaz)
    D <- sqrt(max(D,0))
    
    cos_t <- h/b*(sqrt(D^2 + (tan(sunza)*tan(satza)*sin(relaz))^2)) / (cos(sunza)^(-1)+cos(satza)^(-1))
    t <- acos(max(-1,min(1,cos_t)))
    
    O <- (1/pi)*(t-sin(t)*cos(t))*(cos(sunza)^(-1)+cos(satza)^(-1))
    O <- max(O,0)
    
    kgeo <- O - cos(sunza)^(-1)-cos(satza)^(-1) + 0.5*(1+cos(pa))*cos(sunza)^(-1)*cos(satza)^(-1)
    
    return(list(vol=kvol, geo=kgeo))
  }
  
  ##  Calculate the kernels for initial and normalized angles
  K_init <- kernels(sun_zenith, view_zenith, sun_azimuth-view_azimuth)
  K_norm <- kernels(sun_zenith_norm, 0,0)

  ##  Normalize all layers in ref
  normBand <- function(r){
    band <- names(r)
    b <- switch(tolower(band),
                "blue" = 1,
                "green" = 2,
                "red" = 3,
                "nir" = 4,
                "swir1" = 5,
                "swir2" = 6)
    f_vol_p <- BRDF_coeffs$f_vol_prime[b]
    f_geo_p <- BRDF_coeffs$f_geo_prime[b]
    
    gamma <- (1 + f_vol_p*Kvol_norm + f_geo_p*Kgeo_norm)/(1 + f_vol_p*K_init$vol + f_geo_p*K_init$geo)
    
    r_norm <- round(gamma*r)
  }
  ref_norm <- rast(lapply(1:nlyr(ref), function(x) normBand(subset(ref,x))))
  
  #Re-assign original layer names
  if(!missing(band_names)) {
    names(ref_norm) <- tmp_band_names
  }
  
  return(ref_norm)
}

