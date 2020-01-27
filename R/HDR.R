#' Highest density region (HDR)
#'
#' \code{HDR.xxxx} returns the highest density region (HDR) for a chosen distribution.
#'
#' This function computes the highest density region (HDR) for a univariate distribution in the \code{stats} package.  The functions for
#' the HDR for different distributions are named in the form \code{HDR.xxxx} where the \code{xxxx} refers to the distribution
#' (e.g., \code{HDR.chisq}, \code{HDR.gamma}, \code{HDR.norm}, etc.).  The user can use any univariate distribution in the \code{stats} package,
#' and the function accepts parameters from the specified distribution (see table below).  The output of the function is an interval of classes
#' \code{hdr} and \code{interval} giving the highest density region and some related information pertaining to the distribution and the
#' computation of the HDR (for information on intervals, see the \code{sets} package). If the input distribution is continuous then the
#' HDR is a real interval, and if the input distribution discrete then the HDR is a discrete interval.  For non-trivial cases the computation
#' is done by optimisation using the \code{nlm} function.
#'
#' \tabular{lccc}{
#'   HDR.beta    \tab shape1   \tab shape2 \tab ncp   \cr
#'   HDR.cauchy  \tab location \tab scale  \tab       \cr
#'   HDR.chisq   \tab df       \tab ncp    \tab       \cr
#'   HDR.exp     \tab rate     \tab        \tab       \cr
#'   HDR.f       \tab df1      \tab df2    \tab ncp   \cr
#'   HDR.gamma   \tab shape    \tab rate   \tab scale \cr
#'   HDR.lnorm   \tab meanlog  \tab sdlog  \tab       \cr
#'   HDR.norm    \tab mean     \tab sd     \tab       \cr
#'   HDR.t       \tab df       \tab ncp    \tab       \cr
#'   HDR.unif    \tab min      \tab max    \tab       \cr
#'   HDR.weibull \tab shape    \tab scale  \tab       \cr
#'   HDR.binom   \tab size     \tab prob   \tab       \cr
#'   HDR.geom    \tab prob     \tab        \tab       \cr
#'   HDR.hyper   \tab m        \tab n      \tab k     \cr
#'   HDR.nbinom  \tab size     \tab prob   \tab mu    \cr
#'   HDR.pois    \tab lambda   \tab        \tab       \cr }
#'
#' The table above shows the parameters in each of the distributions.  Some have default values, but most need to be specified.  (For the gamma
#' distribution you should specify either the \code{rate} or \code{scale} but not both.)
#'
#' @param alpha The significance level for the HDR (scalar between zero and one).  The probability coverage for the HDR is \code{1-alpha}.
#' @param shape1,shape2,ncp,location,scale,df,rate,df1,df2,meanlog,sdlog,mean,sd,min,max,shape,size,prob,m,n,k,mu,lambda Distribution parameters.
#' @param gradtol Parameter for the nlm optimisation - a positive scalar giving the tolerance at which the scaled gradient is considered close enough to zero to terminate the algorithm (see [\code{nlm} doccumentation](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/nlm.html)).
#' @param steptol Parameter for the nlm optimisation - a positive scalar providing the minimum allowable relative step length (see [\code{nlm} doccumentation](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/nlm.html)).
#' @param iterlim Parameter for the nlm optimisation - a positive integer specifying the maximum number of iterations to be performed before the program is terminated (see [\code{nlm} doccumentation](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/nlm.html)).
#' @return An interval object with classes \code{hdr} and \code{interval} containing the highest density region and related information.
#' @name HDR
NULL

#########################
# Continous Distributions
#########################

#' @examples
#' HDR.norm(.05)
#' @rdname HDR
HDR.norm <- function(alpha, mean = 0, sd = 1,
                     gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(mean))    { stop('Error: mean should be numeric') }
  if (length(mean) != 1)    { stop('Error: mean should be a single value'); }
  if (!is.numeric(sd))      { stop('Error: sd should be numeric') }
  if (length(sd) != 1)      { stop('Error: sd should be a single value'); }
  if (sd < 0)               { stop('Error: sd is negative'); }

  #Set text for distribution
  DIST <- ifelse(((mean == 0)&(sd == 1)),
                 paste0('standard normal distribution'),
                 paste0('normal distribution with mean = ', mean,
                        ' and standard deviation = ', sd));

  hdr(alpha, unimodal, Q = qnorm, f = dnorm, distribution = DIST,
               mean = mean, sd = sd,
               gradtol = gradtol, steptol = steptol, iterlim = iterlim); }

#' @rdname HDR
HDR.lnorm <- function(alpha, meanlog = 0, sdlog = 1,
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(meanlog)) { stop('Error: meanlog should be numeric') }
  if (length(meanlog) != 1) { stop('Error: meanlog should be a single value'); }
  if (!is.numeric(sdlog))   { stop('Error: sdlog should be numeric') }
  if (length(sdlog) != 1)   { stop('Error: sdlog should be a single value'); }
  if (sdlog < 0)            { stop('Error: sdlog is negative'); }

  #Set text for distribution
  DIST <- ifelse(((meanlog == 0)&(sdlog == 1)),
                 paste0('standard log-normal distribution'),
                 paste0('log-normal distribution with log-mean = ', meanlog,
                        ' and log-standard deviation = ', sdlog));

  hdr(alpha, unimodal, Q = qlnorm, f = dlnorm, distribution = DIST,
           meanlog = meanlog, sdlog = sdlog,
           gradtol = gradtol, steptol = steptol, iterlim = iterlim); }

#' @rdname HDR
HDR.t <- function(alpha, df, ncp = 0,
                  gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(df))      { stop('Error: df should be numeric') }
  if (length(df) != 1)      { stop('Error: df should be a single value'); }
  if (df < 0)               { stop('Error: df is negative'); }
  if (!is.numeric(ncp))     { stop('Error: ncp should be numeric') }
  if (length(ncp) != 1)     { stop('Error: ncp should be a single value'); }
  if (ncp < 0)              { stop('Error: ncp is negative'); }

  #Set text for distribution
  DIST <- ifelse(ncp == 0,
                 paste0('Student\'s T distribution with ', df, ' degrees-of-freedom'),
                 paste0('Student\'s T distribution with ', df,
                        ' degrees-of-freedom and non-centrality parameter = ', ncp));

  hdr(alpha, unimodal, Q = qt, f = dt, distribution = DIST,
               df = df, ncp = ncp,
               gradtol = gradtol, steptol = steptol, iterlim = iterlim); }


#' @rdname HDR
HDR.cauchy <- function(alpha, location = 0, scale = 1,
                       gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(location)) { stop('Error: location should be numeric') }
  if (length(location) != 1) { stop('Error: location should be a single value'); }
  if (!is.numeric(scale))    { stop('Error: scale should be numeric') }
  if (length(scale) != 1)    { stop('Error: scale should be a single value'); }
  if (scale < 0)             { stop('Error: scale is negative'); }

  #Set text for distribution
  DIST <- ifelse(((location == 0)&(scale == 1)),
                 paste0('standard Cauchy distribution'),
                 paste0('Cauchy distribution with location = ', location,
                        ' and scale = ', scale));

  hdr(alpha, unimodal, Q = qcauchy, f = dcauchy, distribution = DIST,
           location = location, scale = scale,
           gradtol = gradtol, steptol = steptol, iterlim = iterlim); }


#' @rdname HDR
HDR.f <- function(alpha, df1, df2, ncp = 0,
                  gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(df1))     { stop('Error: df1 should be numeric') }
  if (length(df1) != 1)     { stop('Error: df1 should be a single value'); }
  if (df1 < 0)              { stop('Error: df1 is negative'); }
  if (!is.numeric(df2))     { stop('Error: df2 should be numeric') }
  if (length(df2) != 1)     { stop('Error: df2 should be a single value'); }
  if (df2 < 0)              { stop('Error: df2 is negative'); }
  if (!is.numeric(ncp))     { stop('Error: ncp should be numeric') }
  if (length(ncp) != 1)     { stop('Error: ncp should be a single value'); }
  if (ncp < 0)              { stop('Error: ncp is negative'); }

  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { qf(L, df1, df2, ncp); }
  DD <- function(L) { df(L, df1, df2, ncp); }

  #Set text for distribution
  DIST <- ifelse(ncp == 0,
                 paste0('F distribution with ', df1,
                        ' numerator degrees-of-freedom and ', df2,
                        ' denominator degrees-of-freedom'),
                 paste0('F distribution with ', df1,
                        ' numerator degrees-of-freedom and ', df2,
                        ' denominator degrees-of-freedom and non-centrality parameter = ', ncp));

  modality <- if(df1 <= 2) monotone else unimodal;

  HDR <- hdr(alpha, modality, qf, df, distribution = DIST,
                  df1 = df1, df2 = df2, ncp = ncp,
                  decreasing = TRUE,
                  gradtol = gradtol, steptol = steptol, iterlim = iterlim);

  HDR; }


#' @rdname HDR
HDR.beta <- function(alpha, shape1, shape2, ncp = 0,
                     gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(shape1))  { stop('Error: shape1 should be numeric') }
  if (length(shape1) != 1)  { stop('Error: shape1 should be a single value'); }
  if (shape1 < 0)           { stop('Error: shape1 is negative'); }
  if (!is.numeric(shape2))  { stop('Error: shape2 should be numeric') }
  if (length(shape2) != 1)  { stop('Error: shape2 should be a single value'); }
  if (shape2 < 0)           { stop('Error: shape2 is negative'); }
  if (!is.numeric(ncp))     { stop('Error: ncp should be numeric') }
  if (length(ncp) != 1)     { stop('Error: ncp should be a single value'); }
  if (ncp < 0)              { stop('Error: ncp is negative'); }

  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { qbeta(L, shape1, shape2, ncp); }
  DD <- function(L) { dbeta(L, shape1, shape2, ncp); }

  #Set text for distribution
  DIST <- ifelse(ncp == 0,
                 ifelse(((shape1 == 1) && (shape2 == 1)),
                        'standard uniform distribution',
                        paste0('beta distribution with shape1 = ', shape1,
                               ' and shape2 = ', shape2)),
                 paste0('beta distribution with shape1 = ', shape1,
                        ' and shape2 = ', shape2,
                        ' and non-centrality parameter = ', ncp));

  decreasing <- FALSE;

  #Compute HDR in monotone decreasing case
  if ((shape1 <= 1) && (shape2 >  1)) {
    modality <- monotone;
    decreasing <- TRUE;}


  #Compute HDR in monotone increasing case
  if ((shape1 >  1) && (shape2 <= 1)) {
    modality <- monotone;}

  #Compute HDR in uniform case
  if ((shape1 == 1) && (shape2 == 1)) {
    modality <- unimodal;}

  #Compute HDR in unimodal case
  if ((shape1 > 1) && (shape2 > 1)) {
    modality <- unimodal;}

  #Compute HDR in bimodal case
  if ((shape1 < 1) && (shape2 < 1)) {
    modality <- bimodal;}

  HDR <- hdr(alpha, modality, Q = qbeta, f = dbeta, distribution = DIST,
                  shape1 = shape1, shape2 = shape2, ncp = ncp,
                  decreasing = decreasing,
                  gradtol = gradtol, steptol = steptol, iterlim = iterlim);

  HDR; }


#' @rdname HDR
HDR.chisq <- function(alpha, df, ncp = 0,
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(df))      { stop('Error: df should be numeric') }
  if (length(df) != 1)      { stop('Error: df should be a single value'); }
  if (df < 0)               { stop('Error: df is negative'); }
  if (!is.numeric(ncp))     { stop('Error: ncp should be numeric') }
  if (length(ncp) != 1)     { stop('Error: ncp should be a single value'); }
  if (ncp < 0)              { stop('Error: ncp is negative'); }


  #Set text for distribution
  DIST <- ifelse(ncp == 0,
                 paste0('chi-squared distribution with ', df, ' degrees-of-freedom'),
                 paste0('chi-squared distribution with ', df,
                        ' degrees-of-freedom and non-centrality parameter = ', ncp));

  #Compute HDR in monotone case;
  if (df <= 2) {
    modality = monotone; }

  #Compute HDR in unimodal case;
  if (df > 2) {
    modality <- unimodal; }

  HDR <- hdr(alpha, modality, Q = qchisq, f = dchisq, distribution = DIST,
                  df = df, ncp = ncp,
                  gradtol = gradtol, steptol = steptol, iterlim = iterlim);

  HDR; }


#' @rdname HDR
HDR.gamma <- function(alpha, shape, rate = 1, scale = 1/rate,
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

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
  DIST <- ifelse((shape == 1),
                 paste0('exponential distribution with scale = ', scale),
                 paste0('gamma distribution with shape = ', shape,
                        ' and scale = ', scale));

  #Compute HDR in monotone case;
  if (shape <= 1) { modality <- monotone; }

  #Compute HDR in unimodal case;
  if (shape > 1) { modality <- unimodal; }

  HDR <- hdr(alpha, modality, Q = qgamma, f = dgamma, distribution = DIST,
                  shape = shape,  scale = scale,
                  gradtol = gradtol, steptol = steptol, iterlim = iterlim);

  HDR; }


#' @rdname HDR
HDR.weibull <- function(alpha, shape, scale = 1,
                        gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(shape))   { stop('Error: shape should be numeric') }
  if (length(shape) != 1)   { stop('Error: shape should be a single value'); }
  if (shape < 0)            { stop('Error: shape is negative'); }
  if (!is.numeric(scale))   { stop('Error: scale should be numeric') }
  if (length(scale) != 1)   { stop('Error: scale should be a single value'); }
  if (scale < 0)            { stop('Error: scale is negative'); }

  #Set text for distribution
  DIST <- ifelse((shape == 1),
                 paste0('exponential distribution with scale = ', scale),
                 paste0('Weibull distribution with shape = ', shape,
                        ' and scale = ', scale));

  #Compute HDR in monotone case;
  if (shape <= 1) { modality <- monotone; }

  #Compute HDR in unimodal case;
  if (shape > 1) { modality <- unimodal; }


  HDR <- hdr(alpha, modality, Q = qweibull, f = dweibull, distribution = DIST,
                      shape = shape, scale = scale,
                      decreasing = TRUE,
                      gradtol = gradtol, steptol = steptol, iterlim = iterlim);

  HDR; }


#' @rdname HDR
HDR.exp <- function(alpha, rate,
                    gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(rate))    { stop('Error: rate should be numeric') }
  if (length(rate) != 1)    { stop('Error: rate should be a single value'); }
  if (rate < 0)             { stop('Error: rate is negative'); }

  #Set text for distribution
  DIST <- paste0('exponential distribution with scale = ', 1/rate);

  #Compute the HDR
  hdr(alpha, monotone, Q = qexp, f = dexp, distribution = DIST,
           rate = rate,
           decreasing = TRUE);}


#' @rdname HDR
HDR.unif <- function(alpha, min = 0, max = 1,
                     gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(min))   { stop('Error: min should be numeric') }
  if (length(min) != 1)   { stop('Error: min should be a single value'); }
  if (!is.numeric(max))   { stop('Error: max should be numeric') }
  if (length(max) != 1)   { stop('Error: max should be a single value'); }
  if (min > max)          { stop('Error: min is greater than max'); }

  #Set text for distribution
  DIST <- ifelse(((min == 1)&(max == 1)),
                 paste0('standard continuous uniform distribution'),
                 paste0('continuous uniform distribution with minimum = ', min,
                        ' and maximum = ', max));

  hdr(alpha, unimodal, Q = qunif, f = dunif, distribution = DIST,
           min = min, max = max,
           gradtol = gradtol, steptol = steptol, iterlim = iterlim); }




#' @rdname HDR
HDR.hyper <- function(alpha, m, n, k,
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(m))  { stop('Error: m should be numeric') }
  if (length(m) != 1)  { stop('Error: m should be a single value'); }
  if (m < 1)           { stop('Error: m should be at least one'); }
  if (!is.numeric(n))  { stop('Error: n should be numeric') }
  if (length(n) != 1)  { stop('Error: n should be a single value'); }
  if (n < 1)           { stop('Error: n should be at least one'); }
  if (!is.numeric(k))  { stop('Error: k should be numeric') }
  if (length(k) != 1)  { stop('Error: k should be a single value'); }
  if (k < 1)           { stop('Error: k should be at least one'); }
  if (k > n + m)       { stop('Error: k is greater than n + m'); }

  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { qhyper(L, m, n, k); }
  FF <- function(L) { phyper(L, m, n, k); }

  #Set text for distribution
  DIST <- paste0('hypergeometric distribution with ', m, ' white balls, ', n,
                 ' black balls, and ', k, ' balls drawn');

  hdr(alpha, modality=discrete.unimodal, Q = qhyper, F = phyper, distribution = DIST,
           m = m, n = n, k = k,
           gradtol = gradtol, steptol = steptol, iterlim = iterlim); }


#' @rdname HDR
HDR.geom <- function(alpha, prob,
                     gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(prob))  { stop('Error: prob should be numeric') }
  if (length(prob) != 1)  { stop('Error: prob should be a single value'); }
  if (prob < 0)           { stop('Error: prob should be between zero and one'); }
  if (prob > 1)           { stop('Error: prob should be between zero and one'); }

  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { qgeom(L, prob); }
  FF <- function(L) { pgeom(L, prob); }

  #Set text for distribution
  DIST <- paste0('geometric distribution with probability = ', prob);

  hdr(alpha, discrete.unimodal, Q = qgeom, F = pgeom, distribution = DIST,
                        prob = prob,
                        gradtol = gradtol, steptol = steptol, iterlim = iterlim); }


#' @rdname HDR
HDR.binom <- function(alpha, size, prob,
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(size))  { stop('Error: size should be numeric') }
  if (length(size) != 1)  { stop('Error: size should be a single value'); }
  if (size < 0)           { stop('Error: size is negative'); }
  if (!is.numeric(prob))  { stop('Error: prob should be numeric') }
  if (length(prob) != 1)  { stop('Error: prob should be a single value'); }
  if (prob < 0)           { stop('Error: prob should be between zero and one'); }
  if (prob > 1)           { stop('Error: prob should be between zero and one'); }

  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { qbinom(L, size, prob); }
  FF <- function(L) { pbinom(L, size, prob); }

  #Set text for distribution
  DIST <- paste0('binomial distribution with size = ', size,
                 ' and probability = ', prob);

  hdr(alpha, discrete.unimodal, Q = qbinom, F = pbinom, distribution = DIST,
           size = size, prob = prob,
                        gradtol = gradtol, steptol = steptol, iterlim = iterlim); }


#' @rdname HDR
HDR.pois <- function(alpha, lambda,
                     gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(lambda)) { stop('Error: lambda should be numeric') }
  if (length(lambda) != 1) { stop('Error: lambda should be a single value'); }
  if (lambda < 0)          { stop('Error: lambda is negative'); }

  #Simplify probability functions (with stipulated parameters)
  QQ <- function(L) { qpois(L, lambda); }
  FF <- function(L) { ppois(L, lambda); }

  #Set text for distribution
  DIST <- paste0('Poisson distribution with rate = ', lambda);

  hdr(alpha, discrete.unimodal, Q = qpois, F = ppois, distribution = DIST,
                        lambda = lambda,
                        gradtol = gradtol, steptol = steptol, iterlim = iterlim); }


#' @rdname HDR
HDR.nbinom <- function(alpha, size, prob, mu,
                       gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  if (!is.numeric(size)) { stop('Error: size should be numeric') }
  if (length(size) != 1) { stop('Error: size should be a single value'); }
  if (size < 0)          { stop('Error: size is negative'); }
  if ((missing(prob))&(missing(mu)))
  { stop('Error: Must specify either prob or mu'); }
  if (!(missing(prob)))  {
    if (!is.numeric(prob)) { stop('Error: prob should be numeric') }
    if (length(prob) != 1) { stop('Error: prob should be a single value'); }
    if (prob < 0)          { stop('Error: prob should be between zero and one'); }
    if (prob > 1)          { stop('Error: prob should be between zero and one'); }}
  if ((!missing(prob))&(!missing(mu)))
  { stop('Error: Both prob and mu are specified'); }
  if (!(missing(mu)))    {
    if (!is.numeric(mu))   { stop('Error: mu should be numeric') }
    if (length(mu) != 1)   { stop('Error: mu should be a single value'); }
    if (mu < 0)            { stop('Error: mu is negative'); }}

  #Simplify probability functions (with stipulated parameters)
  # NOTE we are doing this here because p/qnbinom require mu to be *missing* (not NULL)
  if (!(missing(prob)))  {
    QQ <- function(L, ...) { qnbinom(L, size, prob = prob); }
    FF <- function(L, ...) { pnbinom(L, size, prob = prob); } }
  if (missing(prob)   )  {
    QQ <- function(L, ...) { qnbinom(L, size, mu = mu); }
    FF <- function(L, ...) { pnbinom(L, size, mu = mu); } }

  #Set text for distribution
  DIST <- ifelse(!missing(prob),
                 paste0('negative binomial distribution with size = ',
                        size, ' and probability = ', prob),
                 paste0('negative binomial distribution with size = ',
                        size, ' and mean = ', mu));


  hdr(alpha, discrete.unimodal, 
           #Q = qnbinom, F = pbinom, 
           Q = QQ, F = FF, 
           distribution = DIST,
           #size = size, prob = prob, mu = mu,
           gradtol = gradtol, steptol = steptol, iterlim = iterlim); }

#' @export
print.hdr <- function(x, ...) {

  #Print description of HDR
  cat('\n        Highest Density Region (HDR) \n \n');
  cat(paste0(sprintf(100*attributes(x)$probability, fmt = '%#.2f'), '%'),
      'HDR for', attributes(x)$distribution, '\n');

  #Print method
  if (!is.na(attributes(x)$method)) {
    cat(attributes(x)$method, '\n'); }

  #Print HDR interval
  cat('\n');
  writeLines(as.character(c(x)))
  cat('\n');
  invisible(c(x)); }

