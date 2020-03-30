
#The functions below are the CONF.mean and CONF.var functions, which give optimised confidence intervals for the mean and variance of a population, taken from a set of sample data.  These functions allow specification of the population size, and in cases where this is finite, the confidence intervals impose a “finite population correction”.  These functions generate an interval object with additional class ci, which has its own custom print method.


CONF.mean <- function(alpha, x = NULL, sample.mean = mean(x), 
                      sample.variance = var(x), n = length(x),
                      N = Inf, kurt = 3, unsampled = FALSE,
                      gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check input alpha
  if (!is.numeric(alpha))   { stop('Error: alpha should be numeric') }
  if (length(alpha) != 1)   { stop('Error: alpha should be a single value'); }
  if (alpha < 0)            { stop('Error: alpha is negative'); }
  if (alpha > 1)            { stop('Error: alpha is greater than one'); }
  
  #Check congruence of data inputs
  if ((!missing(x) && !missing(sample.mean))) {
    if (abs(sample.mean - mean(x)) < 1e-15) 
      warning('specify data or sample mean but not both') else
        stop('Error: specify data or sample mean but not both'); }
  if ((!missing(x) && !missing(sample.variance))) {
    if (abs(sample.variance - var(x)) < 1e-15) 
      warning('specify data or sample variance but not both') else
        stop('Error: specify data or sample variance but not both'); }
  if ((!missing(x) && !missing(n))) {
    if (n != length(x)) 
      stop('Error: specify data or n but not both'); }
  
  #Check data inputs
  if (!is.null(x)) {
    if (!is.numeric(x))     { stop('Error: x should be numeric') } }
  if (!is.numeric(sample.mean))    { stop('Error: mean should be numeric') }
  if (length(sample.mean) != 1)    { stop('Error: mean should be a single value'); }
  if (!is.numeric(sample.variance)) { stop('Error: variance should be numeric') }
  if (length(sample.variance) != 1) { stop('Error: variance should be a single value'); }
  if (sample.variance < 0)          { stop('Error: variance is negative'); }
  if (!is.numeric(n))       { stop('Error: n should be numeric') }
  if (as.integer(n) != n)   { stop('Error: n should be an integer') }
  if (length(n) != 1)       { stop('Error: n should be a single value'); }
  if (n < 3)                { stop('Error: Method requires n > 2'); }
  
  #Check inputs N, kurt and unsampled
  if (!is.numeric(N))       { stop('Error: N should be numeric') }
  if (N != Inf) { 
    if (as.integer(N) != N) { stop('Error: N should be an integer') } }
  if (length(N) != 1)       { stop('Error: N should be a single value'); }
  if (N < n)                { stop('Error: N cannot be smaller than n'); }
  if (!is.numeric(kurt))    { stop('Error: kurt should be numeric') }
  if (length(kurt) != 1)    { stop('Error: kurt should be a single value'); }
  if (kurt < 1)             { stop('Error: kurt is less than one'); }
  if (!is.logical(unsampled)) { stop('Error: unsampled should be TRUE/FALSE') }
  if (length(unsampled) != 1) { stop('Error: unsampled should be a single value'); }
  
  #Check inputs for nlm
  if (!is.numeric(gradtol)) { stop('Error: gradtol should be numeric') }
  if (length(gradtol) != 1) { stop('Error: gradtol should be a single value'); }
  if (gradtol <= 0)         { stop('Error: gradtol should be positive'); }
  if (!is.numeric(steptol)) { stop('Error: steptol should be numeric') }
  if (length(steptol) != 1) { stop('Error: steptol should be a single value'); }
  if (steptol <= 0)         { stop('Error: steptol should be positive'); }
  if (!is.numeric(iterlim)) { stop('Error: iterlim should be numeric') }
  if (length(iterlim) != 1) { stop('Error: iterlim should be a single value'); }
  if (iterlim <= 0)         { stop('Error: iterlim should be positive'); }
  
  #############
  
  #Compute CONF
  df  <- 2*n/(kurt - (n-3)/(n-1));
  TT  <- abs(qt(alpha/2, df, ncp = 0));
  UU  <- if (N == Inf) { 1 } else { (N-n)/N };
  ADJ <- if (!unsampled) { sqrt(UU) } else { 1/sqrt(UU) }
  L   <- sample.mean - sample.variance*ADJ*TT/sqrt(n);
  U   <- sample.mean + sample.variance*ADJ*TT/sqrt(n);
  CONF <- sets::interval(l = L, r = U, bounds = 'closed');
  attr(CONF, 'method') <- as.character(NA);
  
  #Add class and attributes
  class(CONF) <- c('ci', 'interval');
  attr(CONF, 'probability') <- 1 - alpha;
  attr(CONF, 'parameter')   <- if (N == Inf) {
    paste0('mean parameter for infinite population') } else {
      if (unsampled) {
        paste0('mean for unsampled population of size ', N-n) } else {
          paste0('mean for population of size ', N) } }
  
  CONF; }


CONF.variance <- function(alpha, x = NULL, 
                          sample.variance = var(x), n = length(x), 
                          N = Inf, kurt = 3, unsampled = FALSE, 
                          gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  
  #Check input alpha
  if (!is.numeric(alpha))   { stop('Error: alpha should be numeric') }
  if (length(alpha) != 1)   { stop('Error: alpha should be a single value'); }
  if (alpha < 0)            { stop('Error: alpha is negative'); }
  if (alpha > 1)            { stop('Error: alpha is greater than one'); }
  
  #Check congruence of data inputs
  if ((!missing(x) && !missing(sample.variance))) {
    if (abs(sample.variance - var(x)) < 1e-15) 
      warning('specify data or sample variance but not both') else
        stop('Error: specify data or sample variance but not both'); }
  if ((!missing(x) && !missing(n))) {
    if (n != length(x)) 
      stop('Error: specify data or n but not both'); }
  
  #Check data inputs
  if (!is.null(x)) {
    if (!is.numeric(x))     { stop('Error: x should be numeric') } }
  if (!is.numeric(sample.variance)) { stop('Error: variance should be numeric') }
  if (length(sample.variance) != 1) { stop('Error: variance should be a single value'); }
  if (sample.variance < 0)          { stop('Error: variance is negative'); }
  if (!is.numeric(n))       { stop('Error: n should be numeric') }
  if (as.integer(n) != n)   { stop('Error: n should be an integer') }
  if (length(n) != 1)       { stop('Error: n should be a single value'); }
  if (n < 3)                { stop('Error: Method requires n > 2'); }
  
  #Check inputs N, kurt and unsampled
  if (!is.numeric(N))       { stop('Error: N should be numeric') }
  if (N != Inf) { 
    if (as.integer(N) != N) { stop('Error: N should be an integer') } }
  if (length(N) != 1)       { stop('Error: N should be a single value'); }
  if (N <= n)               { stop('Error: N should be larger than n'); }
  if (!is.numeric(kurt))    { stop('Error: kurt should be numeric') }
  if (length(kurt) != 1)    { stop('Error: kurt should be a single value'); }
  if (kurt < 1)             { stop('Error: kurt is less than one'); }
  if (!is.logical(unsampled)) { stop('Error: unsampled should be TRUE/FALSE') }
  if (length(unsampled) != 1) { stop('Error: unsampled should be a single value'); }
  
  #Check inputs for nlm
  if (!is.numeric(gradtol)) { stop('Error: gradtol should be numeric') }
  if (length(gradtol) != 1) { stop('Error: gradtol should be a single value'); }
  if (gradtol <= 0)         { stop('Error: gradtol should be positive'); }
  if (!is.numeric(steptol)) { stop('Error: steptol should be numeric') }
  if (length(steptol) != 1) { stop('Error: steptol should be a single value'); }
  if (steptol <= 0)         { stop('Error: steptol should be positive'); }
  if (!is.numeric(iterlim)) { stop('Error: iterlim should be numeric') }
  if (length(iterlim) != 1) { stop('Error: iterlim should be a single value'); }
  if (iterlim <= 0)         { stop('Error: iterlim should be positive'); }
  
  #############
  
  #Simplify probability functions (with stipulated parameters)
  df1 <- 2*n/(kurt - (n-3)/(n-1));
  df2 <- 2*(N-n)/(2 + (kurt-3)*(1-2/N-1/(N*n)));
  Q <- function(L) { qf(L, df1, df2); }
  f <- function(L) { df(L, df1, df2); }
  
  #Set objective function
  WW <- function(phi) { 
    
    #Set parameter functions
    T0 <- alpha/(1+exp(-phi));
    T1 <- T0/(1+exp(phi));
    T2 <- T1*((1-exp(2*phi))/(1+2*exp(phi)+exp(2*phi)));
    
    #Set interval bounds and objective
    L  <- Q(T0);
    U  <- Q(T0 + 1 - alpha);
    W0 <- 1/L - 1/U;
    
    #Set gradient of objective
    if (!is.null(f)) { 
      attr(W0, 'gradient') <- T1*(1/(f(U)*U^2) - 1/(f(L)*L^2)); }
    
    W0; }
  
  #Compute the HDR
  #The starting value for the parameter phi is set to zero
  #This is the exact optima in the case of a symmetric distribution
  OPT <- nlm(WW, p = 0, 
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  TT <- alpha/(1+exp(-OPT$estimate));
  A <- (n-1)/(N-1);
  B <- if (N == Inf) { 1 } else { (N-n)/(N-1) }
  L   <- A + B/Q(TT + 1 - alpha);
  U   <- A + B/Q(TT);
  CONF <- sample.variance*sets::interval(l = L, r = U, bounds = 'closed');
  
  #Add the description of the method
  METHOD <- ifelse((OPT$iterations == 1), 
                   paste0('Computed using nlm optimisation with ', 
                          OPT$iterations, ' iteration (code = ', OPT$code, ')'),
                   paste0('Computed using nlm optimisation with ', 
                          OPT$iterations, ' iterations (code = ', OPT$code, ')'));
  attr(CONF, 'method') <- METHOD;
  
  #Add class and attributes
  class(CONF) <- c('CI', 'interval');
  attr(CONF, 'confidence')  <- 1 - alpha;
  attr(CONF, 'parameter')   <- if (N == Inf) {
    paste0('variance parameter for infinite population') } else {
      if (unsampled) {
        paste0('variance for unsampled population of size ', N-n) } else {
          paste0(' variance for population of size ', N) } }
  
  CONF; }


print.ci <- function(object) {
  
  #Print description of confidence interval
  cat('\n        Confidence Interval (CI) \n \n');
  cat(paste0(sprintf(100*attributes(object)$confidence, fmt = '%#.2f'), '%'),
      'CI for', attributes(object)$parameter, '\n');
  
  #Print method
  if (!is.na(attributes(object)$method)) {
    cat(attributes(object)$method, '\n'); }
  
  #Print confidence interval
  cat('\n');
  writeLines(as.character(c(object)))
  invisible(c(object))
  cat('\n'); }

