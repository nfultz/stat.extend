requireVGAM <- function(){
  #Check for required package
  if (!requireNamespace('VGAM', quietly = TRUE)) {
    stop('Package \'VGAM\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
}

#' @rdname HDR
HDR.benini <- function(cover.prob, shape, y0, scale = y0,
                       gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireVGAM()
  
  #Check inputs
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  if (shape < 0)            { stop('Error: shape is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale <= 0)           { stop('Error: scale should be positive'); }
  
  #Set text for distribution
  DIST <- paste0('Benini distribution with shape = ', shape, 
                 ' and scale = ', scale);
  
  #Compute HDR
  HDR <- hdr(cover.prob, modality = unimodal, Q = VGAM::qbenini, f = VGAM::dbenini,
             distribution = DIST,
             y0=scale, shape=shape,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }

#' @rdname HDR
HDR.frechet <- function(cover.prob, shape, scale = 1, location = 0,
                        gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireVGAM()
  
  #Check inputs
  if (!is.numeric(shape))    { stop('Error: shape should be numeric') }
  if (length(shape) != 1)    { stop('Error: shape should be a single value'); }
  if (shape < 0)             { stop('Error: shape is negative'); }
  if (!is.numeric(scale))    { stop('Error: scale should be numeric') }
  if (length(scale) != 1)    { stop('Error: scale should be a single value'); }
  if (scale < 0)             { stop('Error: scale is negative'); }
  if (!is.numeric(location)) { stop('Error: location should be numeric') }
  if (length(location) != 1) { stop('Error: location should be a single value'); }
  
  #Set text for distribution
  DIST <- paste0('Frechet distribution with shape = ', shape,
                 ', scale = ', scale, ' and location = ', location);
  
  
  
  HDR <- hdr(cover.prob, modality = unimodal, Q = VGAM::qfrechet, f = VGAM::dfrechet,
             distribution = DIST, 
             shape = shape, scale = scale, location = location,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);  
  
  
  
  HDR; }

#' @rdname HDR
HDR.gengamma <- function(cover.prob, d, k, shape1 = d, shape2 = k,
                         rate = 1, scale = 1/rate,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireVGAM()
  
  #Check inputs
  if (!is.numeric(shape1))  { stop('Error: shape1 should be numeric') }
  if (length(shape1) != 1)  { stop('Error: shape1 should be a single value'); }
  if (shape1 < 0)           { stop('Error: shape1 is negative'); }
  if (!is.numeric(shape2))  { stop('Error: shape2 should be numeric') }
  if (length(shape2) != 1)  { stop('Error: shape2 should be a single value'); }
  if (shape2 < 0)           { stop('Error: shape2 is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale < 0)            { stop('Error: scale is negative'); }
  if (!is.numeric(rate))    { stop('Error: rate should be numeric') }
  if (length(rate) != 1)    { stop('Error: rate should be a single value'); }
  if (rate < 0)             { stop('Error: rate is negative'); }
  if ((!missing(rate) && !missing(scale))) {
    if (abs(rate*scale - 1) < 1e-15) 
      warning('specify rate or scale but not both') else
        stop('Error: specify rate or scale but not both'); }
  
  #Set text for distribution
  DIST <- paste0('generalised gamma (Stacy) distribution with shape1 = ', shape1,
                 ', shape2 = ', shape2, ' and scale = ', scale);
  
  #Compute HDR in monotone case
  if (shape1 <= 1) {
    modality <- monotone
  } else if (shape1 > 1) {
    modality <- unimodal 
  }
  
  HDR <- hdr(cover.prob, modality = modality, Q = VGAM::qgengamma.stacy, f = VGAM::dgengamma.stacy, decreasing=TRUE,
             distribution = DIST, 
             scale = scale, d = shape1, k = shape2,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);  
  
  
  HDR; }

#' @rdname HDR
HDR.gumbelII <- function(cover.prob, shape, scale = 1,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireVGAM()
  
  #Check inputs
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  if (shape < 0)            { stop('Error: shape is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale < 0)            { stop('Error: scale is negative'); }
  
  #Set text for distribution
  DIST <- paste0('Gumbel (Type II) distribution with shape = ', shape,
                 ' and scale = ', scale);
  
  
  HDR <- hdr(cover.prob, modality = unimodal, Q = VGAM::qgumbelII, f = VGAM::dgumbelII,
             distribution = DIST, 
             shape = shape, scale = scale,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);  
  
  
  HDR; }

#' @rdname HDR
HDR.lgamma <- function(cover.prob, shape = 1, scale = 1, location = 0,
                       gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireVGAM()
  
  #Check inputs
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  if (shape < 0)            { stop('Error: shape is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale < 0)            { stop('Error: scale is negative'); }
  if (!is.numeric(location)) { stop('Error: location should be numeric') }
  if (length(location) != 1) { stop('Error: location should be a single value'); }
  
  #Set text for distribution
  DIST <- paste0('log-gamma distribution with shape = ', shape,
                 ', scale = ', scale, ' and location = ', location);
  
  HDR <- hdr(cover.prob, modality = unimodal, Q = VGAM::qlgamma, f = VGAM::dlgamma,
             distribution = DIST, 
             shape = shape, scale = scale, location = location,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);  
  
  
  
  HDR; }
