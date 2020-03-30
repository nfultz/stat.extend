#' @rdname HDR
HDR.arcsine <- function(cover.prob, min = 0, max = 1,
                        gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check inputs
  if (!is.numeric(min))     { stop('Error: min should be numeric') }
  if (length(min) != 1)     { stop('Error: min should be a single value'); }
  if (!is.numeric(max))     { stop('Error: max should be numeric') }
  if (length(max) != 1)     { stop('Error: max should be a single value'); }
  if (min > max)            { stop('Error: min is larger than max'); }
  

  #Set text for distribution
  DIST <- ifelse(((min == 0)&(max == 1)),
                 'standard arcsine distribution',
                 paste0('arcsine distribution with minimum = ', min, 
                        ' and maximum = ', max));
  
  #Compute HDR using beta(1/2, 1/2) distribution
  
  betaHDR <- hdr(cover.prob, modality=bimodal, Q = qbeta, f = dbeta, distribution = DIST,
                 shape1 = 1/2, shape2 = 1/2,
                 gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
    
  HDR <- (max-min)*betaHDR + min; 
  HDR}




HDR.invchisq <- function(cover.prob, df, ncp = 0,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('invgamma', quietly = TRUE)) {
    stop('Package \'invgamma\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
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
##############################################3


HDR.invgamma <- function(cover.prob, shape, rate = 1, scale = 1/rate,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('invgamma', quietly = TRUE)) {
    stop('Package \'invgamma\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
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


HDR.invexp <- function(cover.prob, rate = 1,
                       gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('invgamma', quietly = TRUE)) {
    stop('Package \'invgamma\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
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


HDR.betapr <- function(cover.prob, shape1, shape2, scale = 1,
                       gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('extraDistr', quietly = TRUE)) {
    stop('Package \'extraDistr\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
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

###################################################3

HDR.benini <- function(cover.prob, shape, y0, scale = y0,
                       gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('VGAM', quietly = TRUE)) {
    stop('Package \'VGAM\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
  #Check inputs
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  if (shape < 0)            { stop('Error: shape is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale <= 0)           { stop('Error: scale should be positive'); }
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { VGAM::qbenini(L, y0 = scale, shape = shape); }
  DD <- function(L) { VGAM::dbenini(L, y0 = scale, shape = shape); }
  
  #Set text for distribution
  DIST <- paste0('Benini distribution with shape = ', shape, 
                 ' and scale = ', scale);
  
  #Compute HDR
  HDR <- hdr(cover.prob, modality = unimodal, Q = VGAM::qbenini, f = VGAM::dbenini,
             distribution = DIST,
             y0=scale, shape=shape,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }

#################################

HDR.fatigue <- function(cover.prob, alpha, beta = 1, mu = 0, 
                        gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('extraDistr', quietly = TRUE)) {
    stop('Package \'extraDistr\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
  #Check inputs
  if (!is.numeric(alpha))   { stop('Error: alpha should be numeric') }
  if (length(alpha) != 1)   { stop('Error: alpha should be a single value'); }
  if (alpha < 0)            { stop('Error: alpha is negative'); }
  if (!is.numeric(beta))    { stop('Error: beta should be numeric') }
  if (length(beta) != 1)    { stop('Error: beta should be a single value'); }
  if (beta < 0)             { stop('Error: beta is negative'); }
  if (!is.numeric(mu))      { stop('Error: mu should be numeric') }
  if (length(mu) != 1)      { stop('Error: mu should be a single value'); }
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { extraDistr::qfatigue(L, alpha, beta, mu); }
  DD <- function(L) { extraDistr::dfatigue(L, alpha, beta, mu); }
  
  #Set text for distribution
  DIST <- ifelse(((location == 0)&(scale == 1)), 
                 paste0('standard fatigue-life (Birnbaum-Saunders) distribution with shape = ', alpha),
                 paste0('fatigue-life (Birnbaum-Saunders) distribution with shape = ', alpha, ', scale = ', beta, ' and location = ', mu));
  
  #Compute HDR
  HDR <- HDR.unimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                      gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }


HDR.gpd <- function(cover.prob, mu = 0, sigma = 1, xi = 0,
                    location = mu, scale = sigma, shape = xi, 
                    gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('extraDistr', quietly = TRUE)) {
    stop('Package \'extraDistr\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
  #Check inputs
  if (!is.numeric(location))      { stop('Error: location should be numeric') }
  if (length(location) != 1) { stop('Error: location should be a single value'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale <= 0)           { stop('Error: scale should be positive'); }
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { extraDistr::qgpd(L, mu = location, 
                                       sigma = scale, xi = shape); }
  DD <- function(L) { extraDistr::dgpd(L, mu = location, 
                                       sigma = scale, xi = shape); }
  
  #Set text for distribution
  DIST <- paste0('generalised Pareto distribution with location = ', location,
                 ', scale = ', scale, ' and shape = ', shape);
  
  #Compute HDR in monotone case
  if (shape >= -1) {
    HDR <- HDR.monotone(cover.prob, Q = QQ, f = DD, distribution = DIST); }
  
  #Compute HDR in monotone case
  if (shape < -1) {
    HDR <- HDR.monotone(cover.prob, Q = QQ, f = DD, distribution = DIST,
                        decreasing = FALSE); }
  
  HDR; }


HDR.gompertz <- function(cover.prob, shape = 1, scale = 1, 
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('extraDistr', quietly = TRUE)) {
    stop('Package \'extraDistr\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
  #Check inputs
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  if (shape < 0)            { stop('Error: shape is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale < 0)            { stop('Error: scale is negative'); }
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { extraDistr::qgompertz(L, a = shape, b = scale); }
  DD <- function(L) { extraDistr::dgompertz(L, a = shape, b = scale); }
  
  #Set text for distribution
  DIST <- paste0('Gompertz distribution with shape = ', shape,
                 ' and scale = ', scale);
  
  #Compute HDR in monotonic case
  if (shape >= 1) {
    HDR <- HDR.monotone(cover.prob, Q = QQ, f = DD, distribution = DIST); }
  
  #Compute HDR in unimodal case
  if (shape <  1) {
    HDR <- HDR.unimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                        gradtol = gradtol, steptol = steptol, iterlim = iterlim); }
  
  HDR; }


HDR.huber <- function(cover.prob, mu, sigma, epsilon, 
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('extraDistr', quietly = TRUE)) {
    stop('Package \'extraDistr\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
  #Check inputs
  if (!is.numeric(mu))      { stop('Error: mu should be numeric') }
  if (length(mu) != 1)      { stop('Error: mu should be a single value'); }
  if (!is.numeric(sigma))   { stop('Error: sigma should be numeric') }
  if (length(sigma) != 1)   { stop('Error: sigma should be a single value'); }
  if (sigma <= 0)           { stop('Error: sigma must be positive'); }
  if (!is.numeric(epsilon)) { stop('Error: epsilon should be numeric') }
  if (length(epsilon) != 1) { stop('Error: epsilon should be a single value'); }
  if (epsilon <= 0)         { stop('Error: epsilon must be positive'); }
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { extraDistr::qhuber(L, mu = mu, sigma = sigma, 
                                         epsilon = epsilon); }
  DD <- function(L) { extraDistr::dhuber(L, mu = mu, sigma = sigma, 
                                         epsilon = epsilon); }
  
  #Set text for distribution
  if (epsilon == Inf) {
    DIST <- ifelse(((mu == 0)&(sigma == 1)),
                   'standard normal distribution',
                   paste0('normal distribution with mean = ', mu, 
                          ' and standard deviation = ', sigma)); }
  if (epsilon <  Inf) {
    DIST <- paste0('Huber distribution with mean = ', mu,
                   ', scale = ', sigma, ' and cut-point = ', epsilon); }
  
  #Compute HDR in unimodal case
  HDR <- HDR.unimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                      gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }


HDR.kumar <- function(cover.prob, a = 1, b = 1, shape1 = a, shape2 = b, 
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('extraDistr', quietly = TRUE)) {
    stop('Package \'extraDistr\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
  #Check inputs
  if (!is.numeric(shape1))  { stop('Error: shape1 should be numeric') }
  if (length(shape1) != 1)  { stop('Error: shape1 should be a single value'); }
  if (shape1 < 0)           { stop('Error: shape1 is negative'); }
  if (!is.numeric(shape2))  { stop('Error: shape2 should be numeric') }
  if (length(shape2) != 1)  { stop('Error: shape2 should be a single value'); }
  if (shape2 < 0)           { stop('Error: shape2 is negative'); }
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { extraDistr::qkumar(L, a = shape1, b = shape2); }
  DD <- function(L) { extraDistr::dkumar(L, a = shape1, b = shape2); }
  
  #Set text for distribution
  DIST <- ifelse(((shape1 == 1)&(shape2 == 1)), 
                 'uniform distribution',
                 paste0('Kumaraswamy distribution with shape1 = ', shape1, 
                        ' and shape2 = ', shape2));
  
  #Compute HDR in monotone decreasing case
  if ((shape1 <= 1) && (shape2 >  1)) {
    HDR <- HDR.monotone(cover.prob, Q = QQ, f = DD, distribution = DIST,
                        decreasing = TRUE); }
  
  #Compute HDR in monotone increasing case
  if ((shape1 >  1) && (shape2 <= 1)) {
    HDR <- HDR.monotone(cover.prob, Q = QQ, f = DD, distribution = DIST,
                        decreasing = FALSE); }
  
  #Compute HDR in uniform case
  if ((shape1 == 1) && (shape2 == 1)) {
    HDR <- HDR.unimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                        gradtol = gradtol, steptol = steptol, iterlim = iterlim); }
  
  #Compute HDR in unimodal case
  if ((shape1 > 1) && (shape2 > 1)) {
    HDR <- HDR.unimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                        gradtol = gradtol, steptol = steptol, iterlim = iterlim); }
  
  #Compute HDR in bimodal case
  if ((shape1 < 1) && (shape2 < 1)) {
    HDR <- HDR.bimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                       gradtol = gradtol, steptol = steptol, iterlim = iterlim); }
  
  HDR; }


HDR.frechet <- function(cover.prob, shape, scale = 1, location = 0,
                        gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('VGAM', quietly = TRUE)) {
    stop('Package \'VGAM\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
  #Check inputs
  if (!is.numeric(shape))    { stop('Error: shape should be numeric') }
  if (length(shape) != 1)    { stop('Error: shape should be a single value'); }
  if (shape < 0)             { stop('Error: shape is negative'); }
  if (!is.numeric(scale))    { stop('Error: scale should be numeric') }
  if (length(scale) != 1)    { stop('Error: scale should be a single value'); }
  if (scale < 0)             { stop('Error: scale is negative'); }
  if (!is.numeric(location)) { stop('Error: location should be numeric') }
  if (length(location) != 1) { stop('Error: location should be a single value'); }
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { VGAM::qfrechet(L, shape = shape, scale = scale, location = location); }
  DD <- function(L) { VGAM::dfrechet(L, shape = shape, scale = scale, location = location); }
  
  #Set text for distribution
  DIST <- paste0('Fréchet distribution with shape = ', shape,
                 ', scale = ', scale, ' and location = ', location);
  
  #Compute HDR
  HDR <- HDR.unimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                      gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }


HDR.gumbelII <- function(cover.prob, shape, scale = 1,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('VGAM', quietly = TRUE)) {
    stop('Package \'VGAM\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
  #Check inputs
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  if (shape < 0)            { stop('Error: shape is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale < 0)            { stop('Error: scale is negative'); }
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { VGAM::qgumbelII(L, shape = shape, scale = scale); }
  DD <- function(L) { VGAM::dgumbelII(L, shape = shape, scale = scale); }
  
  #Set text for distribution
  DIST <- paste0('Gumbel (Type II) distribution with shape = ', shape,
                 ' and scale = ', scale);
  
  #Compute HDR
  HDR <- HDR.unimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                      gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }


HDR.lgamma <- function(cover.prob, shape = 1, scale = 1, location = 0,
                       gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('VGAM', quietly = TRUE)) {
    stop('Package \'VGAM\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
  #Check inputs
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  if (shape < 0)            { stop('Error: shape is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale < 0)            { stop('Error: scale is negative'); }
  if (!is.numeric(location)) { stop('Error: location should be numeric') }
  if (length(location) != 1) { stop('Error: location should be a single value'); }
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { VGAM::qlgamma(L, shape = shape, scale = scale, location = location); }
  DD <- function(L) { VGAM::dlgamma(L, shape = shape, scale = scale, location = location); }
  
  #Set text for distribution
  DIST <- paste0('log-gamma distribution with shape = ', shape,
                 ', scale = ', scale, ' and location = ', location);
  
  #Compute HDR
  HDR <- HDR.unimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                      gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }


HDR.gengamma <- function(cover.prob, d, shape1 = d, k, shape2 = k,
                         rate = 1, scale = 1/rate,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('VGAM', quietly = TRUE)) {
    stop('Package \'VGAM\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
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
  
  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { VGAM::qgengamma.stacy(L, scale = scale, d = shape1, k = shape2); }
  DD <- function(L) { VGAM::dgengamma.stacy(L, scale = scale, d = shape1, k = shape2); }
  
  #Set text for distribution
  DIST <- paste0('generalised gamma (Stacy) distribution with shape1 = ', shape1,
                 ', shape2 = ', shape2, ' and scale = ', scale);
  
  #Compute HDR in monotone case
  if (shape1 <= 1) {
    HDR <- HDR.monotone(cover.prob, Q = QQ, f = DD, distribution = DIST,
                        gradtol = gradtol, steptol = steptol, iterlim = iterlim); }
  
  #Compute HDR in unimodal case
  if (shape1 > 1) {
    HDR <- HDR.unimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                        gradtol = gradtol, steptol = steptol, iterlim = iterlim); }
  
  HDR; }


HDR.tnorm <- function(cover.prob, mean = 0, sd = 1, 
                      a = -Inf, b = Inf, min = a, max = b,
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check for required package
  if (!requireNamespace('extraDistr', quietly = TRUE)) {
    stop('Package \'extraDistr\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
  
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
  
  #Compute HDR
  HDR <- HDR.unimodal(cover.prob, Q = QQ, f = DD, distribution = DIST,
                      gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  
  HDR; }
