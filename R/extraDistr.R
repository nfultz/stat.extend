
requireExtraDistr <- function() {
  #Check for required package
  if (!requireNamespace('extraDistr', quietly = TRUE)) {
    stop('Package \'extraDistr\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
}



#' @rdname HDR
HDR.betapr <- function(cover.prob, shape1, shape2, scale = 1,
                       gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireExtraDistr()
  
  #Check inputs
  if (!is.numeric(shape1))  { stop('Error: shape1 should be numeric') }
  if (length(shape1) != 1)  { stop('Error: shape1 should be a single value'); }
  if (shape1 < 0)           { stop('Error: shape1 is negative'); }
  if (!is.numeric(shape2))  { stop('Error: shape2 should be numeric') }
  if (length(shape2) != 1)  { stop('Error: shape2 should be a single value'); }
  if (shape2 < 0)           { stop('Error: shape2 is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale <= 0)           { stop('Error: scale should be positive'); }
  
  #Set text for distribution
  DIST <- ifelse(scale == 1, 
                 paste0('beta prime distribution with shape1 = ', shape1, 
                        ' and shape2 = ', shape2),
                 paste0('generalised beta prime distribution with shape1 = ', shape1, 
                        ', shape2 = ', shape2, ' and scale = ', scale));
  
  #Compute HDR modality
  modality <- if(shape1 <= 1) monotone else unimodal;

  #Compute HDR
  HDR <- hdr(cover.prob, modality = modality, Q = extraDistr::qbetapr, f = extraDistr::dbetapr,
             distribution = DIST,
             shape1=shape1, shape2=shape2, scale=scale,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }

#' @rdname HDR
HDR.fatigue <- function(cover.prob, alpha, beta = 1, mu = 0, 
                        gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireExtraDistr()
  
  #Check inputs
  if (!is.numeric(alpha))   { stop('Error: alpha should be numeric') }
  if (length(alpha) != 1)   { stop('Error: alpha should be a single value'); }
  if (alpha < 0)            { stop('Error: alpha is negative'); }
  if (!is.numeric(beta))    { stop('Error: beta should be numeric') }
  if (length(beta) != 1)    { stop('Error: beta should be a single value'); }
  if (beta < 0)             { stop('Error: beta is negative'); }
  if (!is.numeric(mu))      { stop('Error: mu should be numeric') }
  if (length(mu) != 1)      { stop('Error: mu should be a single value'); }
  
  #Set text for distribution
  DIST <- ifelse(((mu == 0)&(beta == 1)), 
                 paste0('standard fatigue-life (Birnbaum-Saunders) distribution with shape = ', alpha),
                 paste0('fatigue-life (Birnbaum-Saunders) distribution with shape = ', alpha, ', scale = ', beta, ' and location = ', mu));
  
  #Compute HDR
  HDR <- hdr(cover.prob, modality = unimodal, Q = extraDistr::qfatigue, f = extraDistr::dfatigue,
             distribution = DIST,
             alpha = alpha, beta = beta, mu = mu,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  
  HDR; }


#' @rdname HDR
HDR.gompertz <- function(cover.prob, shape = 1, scale = 1, 
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireExtraDistr()
  
  #Check inputs
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  if (shape < 0)            { stop('Error: shape is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale < 0)            { stop('Error: scale is negative'); }
  
  #Set text for distribution
  DIST <- paste0('Gompertz distribution with shape = ', shape,
                 ' and scale = ', scale);
  
  
  HDR <- hdr(cover.prob, modality = unimodal, Q = extraDistr::qgompertz, f = extraDistr::dgompertz,
             distribution = DIST, decreasing = shape >= -1,
             b = scale, a = shape,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);  
  
  
  HDR; }


#' @rdname HDR
HDR.gpd <- function(cover.prob, mu = 0, sigma = 1, xi = 0,
                    location = mu, scale = sigma, shape = xi, 
                    gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireExtraDistr()
  
  #Check inputs
  if (!is.numeric(location))      { stop('Error: location should be numeric') }
  if (length(location) != 1) { stop('Error: location should be a single value'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale <= 0)           { stop('Error: scale should be positive'); }
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  
  #Set text for distribution
  DIST <- paste0('generalised Pareto distribution with location = ', location,
                 ', scale = ', scale, ' and shape = ', shape);
  
  
  HDR <- hdr(cover.prob, modality = unimodal, Q = extraDistr::qgpd, f = extraDistr::dgpd,
             distribution = DIST, decreasing = shape >= -1,
             mu = location, sigma = scale, xi= shape,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);  
  
  HDR; }

#' @rdname HDR
HDR.huber <- function(cover.prob, mu, sigma, epsilon, 
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireExtraDistr()
  
  #Check inputs
  if (!is.numeric(mu))      { stop('Error: mu should be numeric') }
  if (length(mu) != 1)      { stop('Error: mu should be a single value'); }
  if (!is.numeric(sigma))   { stop('Error: sigma should be numeric') }
  if (length(sigma) != 1)   { stop('Error: sigma should be a single value'); }
  if (sigma <= 0)           { stop('Error: sigma must be positive'); }
  if (!is.numeric(epsilon)) { stop('Error: epsilon should be numeric') }
  if (length(epsilon) != 1) { stop('Error: epsilon should be a single value'); }
  if (epsilon <= 0)         { stop('Error: epsilon must be positive'); }
  
  
  #Set text for distribution
  if (epsilon == Inf) {
    DIST <- ifelse(((mu == 0)&(sigma == 1)),
                   'standard normal distribution',
                   paste0('normal distribution with mean = ', mu, 
                          ' and standard deviation = ', sigma)); }
  if (epsilon <  Inf) {
    DIST <- paste0('Huber distribution with mean = ', mu,
                   ', scale = ', sigma, ' and cut-point = ', epsilon); }
  
  
  HDR <- hdr(cover.prob, modality = unimodal, Q = extraDistr::qhuber, f = extraDistr::dhuber,
             distribution = DIST,
             mu = mu, sigma = sigma, epsilon = epsilon,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);  
  
  
  HDR; }

#' @rdname HDR
HDR.kumar <- function(cover.prob, a = 1, b = 1, shape1 = a, shape2 = b, 
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireExtraDistr()
  
  #Check inputs
  if (!is.numeric(shape1))  { stop('Error: shape1 should be numeric') }
  if (length(shape1) != 1)  { stop('Error: shape1 should be a single value'); }
  if (shape1 <= 0)           { stop('Error: shape1 is negative'); }
  if (!is.numeric(shape2))  { stop('Error: shape2 should be numeric') }
  if (length(shape2) != 1)  { stop('Error: shape2 should be a single value'); }
  if (shape2 <= 0)           { stop('Error: shape2 is negative'); }
  
  #Set text for distribution
  DIST <- ifelse(((shape1 == 1)&(shape2 == 1)), 
                 'uniform distribution',
                 paste0('Kumaraswamy distribution with shape1 = ', shape1, 
                        ' and shape2 = ', shape2));
  
  decreasing <- FALSE;
  
  #Compute HDR in monotone decreasing case
  if ((shape1 <= 1) && (shape2 >  1)) {
    modality <- monotone;
    decreasing <- TRUE;
  }

  #Compute HDR in monotone increasing case
  if ((shape1 >  1) && (shape2 <= 1)) {
    modality <- monotone;
  }

  #Compute HDR in uniform case
  if ((shape1 == 1) && (shape2 == 1)) {
    modality <- unimodal
  }

  #Compute HDR in unimodal case
  if ((shape1 > 1) && (shape2 > 1)) {
    modality <- unimodal
  }
  
  #Compute HDR in bimodal case
  if ((shape1 < 1) && (shape2 < 1)) {
    modality <- bimodal
  }
  
  HDR <- hdr(cover.prob, modality = modality, Q = extraDistr::qkumar, f = extraDistr::dkumar,
             distribution = DIST, decreasing = decreasing,
             a = shape1, b = shape2,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);  
  

  HDR; }




#' @rdname HDR
HDR.tnorm <- function(cover.prob, mean = 0, sd = 1, 
                      a = -Inf, b = Inf, min = a, max = b,
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  requireExtraDistr()
  
  #Check inputs
  if (!is.numeric(mean))    { stop('Error: mean should be numeric') }
  if (length(mean) != 1)    { stop('Error: mean should be a single value'); }
  if (!is.numeric(sd))      { stop('Error: sd should be numeric') }
  if (length(sd) != 1)      { stop('Error: sd should be a single value'); }
  if (sd < 0)               { stop('Error: sd is negative'); }
  if (!is.numeric(min))     { stop('Error: min should be numeric') }
  if (length(min) != 1)     { stop('Error: min should be a single value'); }
  if (!is.numeric(max))     { stop('Error: max should be numeric') }
  if (length(max) != 1)     { stop('Error: max should be a single value'); }
  if (min > max)            { stop('Error: min cannot be larger than max '); }
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { extraDistr::qtnorm(L, mean = mean, sd = sd, 
                                         a = min, b = max); }
  DD <- function(L) { extraDistr::dtnorm(L, mean = mean, sd = sd, 
                                         a = min, b = max); }
  
  #Set text for distribution
  if ((min == -Inf)&(max == Inf))  {
    DIST <- ifelse(((mean == 0)&(sd == 1)),
                   paste0('standard normal distribution'),
                   paste0('normal distribution with mean = ', mean,
                          ' and standard deviation = ', sd)); }
  if ((min == -Inf)&(max <  Inf)) {
    DIST <- ifelse(((mean == 0)&(sd == 1)),
                   paste0('right-truncated standard normal distribution with max = ', max),
                   paste0('right-truncated normal distribution with mean = ', mean,
                          ' and standard deviation = ', sd, ' with max = ', max)); }
  if ((min >  -Inf)&(max == Inf)) {
    DIST <- ifelse(((mean == 0)&(sd == 1)),
                   paste0('left-truncated standard normal distribution with min = ', min),
                   paste0('left-truncated normal distribution with mean = ', mean,
                          ' and standard deviation = ', sd, ' with min = ', min)); }
  if ((min >  -Inf)&(max <  Inf)) {
    DIST <- ifelse(((mean == 0)&(sd == 1)),
                   paste0('truncated standard normal distribution with min = ', 
                          min, ' and max = ', max),
                   paste0('truncated normal distribution with mean = ', mean,
                          ' and standard deviation = ', sd, ' with min = ', 
                          min, ' and max = ', max)); }
  
  
  HDR <- hdr(cover.prob, modality = unimodal, Q = extraDistr::qtnorm, f = extraDistr::dtnorm,
             distribution = DIST, 
             mean = mean, sd = sd, a = min, b = max,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);  
  
  
  HDR; }
