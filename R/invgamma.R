requireInvgamma <- function() {
  #Check for required package
  if (!requireNamespace('invgamma', quietly = TRUE)) {
    stop('Package \'invgamma\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
}


#' @rdname HDR
HDR.invchisq <- function(cover.prob, df, ncp = 0,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireInvgamma()
  
  #Check inputs
  if (!is.numeric(df))      { stop('Error: df should be numeric') }
  if (length(df) != 1)      { stop('Error: df should be a single value'); }
  if (df < 0)               { stop('Error: df is negative'); }
  if (!is.numeric(ncp))     { stop('Error: ncp should be numeric') }
  if (length(ncp) != 1)     { stop('Error: ncp should be a single value'); }
  if (ncp < 0)              { stop('Error: ncp is negative'); }
  
  #Set text for distribution
  DIST <- ifelse(ncp == 0,
                 paste0('inverse chi-squared distribution with ', df, ' degrees-of-freedom'),
                 paste0('inverse chi-squared distribution with ', df, 
                        ' degrees-of-freedom and non-centrality parameter = ', ncp));
  
  #Compute HDR
  HDR <- hdr(cover.prob, modality=unimodal, Q = invgamma::qinvchisq, f = invgamma::dinvchisq,
             distribution = DIST,
             df=df, ncp=ncp,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }

#' @rdname HDR
HDR.invexp <- function(cover.prob, rate = 1,
                       gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireInvgamma()
  
  #Check inputs
  if (!is.numeric(rate))    { stop('Error: rate should be numeric') }
  if (length(rate) != 1)    { stop('Error: rate should be a single value'); }
  if (rate < 0)             { stop('Error: rate is negative'); }
  
  #Set text for distribution
  DIST <- paste0('inverse exponential distribution with rate = ', rate);
  
  #Compute HDR
  HDR <- hdr(cover.prob, modality = unimodal, Q = invgamma::qinvexp, f = invgamma::dinvexp,
             distribution = DIST,
             rate=rate,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }



#' @rdname HDR
HDR.invgamma <- function(cover.prob, shape, rate = 1, scale = 1/rate,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireInvgamma()
  
  #Check inputs
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  if (shape < 0)            { stop('Error: shape is negative'); }
  if (!is.numeric(rate))    { stop('Error: rate should be numeric') }
  if (length(rate) != 1)    { stop('Error: rate should be a single value'); }
  if (rate < 0)             { stop('Error: rate is negative'); }
  if ((!missing(rate) && !missing(scale))) {
    if (abs(rate*scale - 1) < 1e-15) 
      warning('specify rate or scale but not both') else
        stop('Error: specify rate or scale but not both'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale < 0)            { stop('Error: scale is negative'); }
  
  #Set text for distribution
  DIST <- paste0('inverse gamma distribution with shape = ', shape,
                 ' and scale = ', scale);
  
  #Compute HDR
  HDR <- hdr(cover.prob, modality = unimodal, Q = invgamma::qinvgamma, f = invgamma::dinvgamma,
             distribution = DIST,
             shape = shape, scale = scale,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }

