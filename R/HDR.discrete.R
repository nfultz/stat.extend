#' Highest density region (HDR) for an arbitrary discrete distribution
#'
#' This function computes the highest density region (HDR) with support on the integers.  The distribution can be any discrete
#' distribution concentrated on the integers --- it does not have to have any shape properties for the function to work.  The user
#' must give the density function ```f``` for the distribution.  To improve the search properties of the algorithm, the user can
#' also give lower and upper bounds for the support of the distribution if these are available.  (Warning: If the user specifies
#' incorrect bounds on the support, that do not contain the full support of the distribution, then the algorithm may continue to
#' search without end, in which case the function will not terminate.  Similarly, if the user specifies a sequence function E that
#' is not a proper bijection to the integers then the algorithm may continue to search without end, in which case the function will
#' not terminate.) The output of the function is a \code{'hdr'} object containing the HDR for the discrete distribution.
#'
#' @param cover.prob The minimum coverage probability for the region
#' @param f The density (mass) function for the distribution
#' @param supp.min A minimum bound for the support of the distribution
#' @param supp.max A maximum bound for the support of the distribution
#' @param E A bijective function mapping the natural numbers (1,2,3,...) to a set covering the support of the distribution (optional);
#' if included, the algorithm will search the support of the distribution in the order specified by this function; if not included,
#' the algorithm will search the integers in a default order.
#' @param ... additional parameters of f
#' @param distribution a label
#' @return If all inputs are correctly specified (i.e., arguments and parameters are in allowable range)
#' then the output will be a list of class ```hdr``` containing the HDR and related information.

HDR.discrete <- function(cover.prob, f, supp.min = -Inf, supp.max = Inf, E = NULL,
                         ...,
                         distribution = 'an unspecified input distribution') {
  
  #Check inputs
  if (!is.numeric(cover.prob))   { stop('Error: cover.prob should be numeric') }
  if (length(cover.prob) != 1)   { stop('Error: cover.prob should be a single value'); }
  if (cover.prob < 0)            { stop('Error: cover.prob is negative'); }
  if (cover.prob > 1)            { stop('Error: cover.prob is greater than one'); }
  if (!is.numeric(supp.min))   { stop('Error: supp.min should be numeric') }
  if (length(supp.min) != 1)   { stop('Error: supp.min should be a single value'); }
  if (!is.numeric(supp.max))   { stop('Error: supp.max should be numeric') }
  if (length(supp.max) != 1)   { stop('Error: supp.max should be a single value'); }
  supp.min <- ceiling(supp.min);
  supp.max <- floor(supp.max);
  if (supp.max < supp.min)       { stop('Error: specified supp.min and supp.max give empty support'); }
  
  #############
  
  f <- partial(f, ...);
  
  
  #Compute the HDR in trivial cases where cover.prob is 0 or 1
  #When cover.prob = 0 the HDR is the empty region
  if (cover.prob == 0) {
    HDR <- sets::interval();
    attr(HDR, 'probability') <- 0; }
  #When cover.prob = 1 the HDR is the support of the distribution
  if (cover.prob == 1) {
    if ((supp.min >  -Inf)&(supp.max <  Inf)) {
      VAL <- (supp.min:supp.max);
      for (i in 1:length(VAL)) { if (f(VAL[i]) == 0) { VAL[i] <- NA; } }
      VAL <- VAL[!is.na(VAL)];
      HDR <- sets::set(VAL);
      HDR <- sets::as.interval(HDR);
      attr(HDR, 'probability') <- 1; } else {
        MSG <- paste0('Attempting to find the support of discrete distribution ', deparse(substitute(f)),
                      ' with no finite bound on support given by user \n',
                      '--- Output region is set to be the full range from supp.min to supp.max \n',
                      '--- This might not be the HDR for the distribution');
        warning(MSG);
        HDR <- sets::set(c(supp.min, supp.max));
        HDR <- sets::as.interval(HDR);
        attr(HDR, 'probability') <- 1; } }
  
  #############
  
  #Compute the HDR in non-trivial cases where 0 < cover.prob < 1
  #Computation uses the iterative one-at-a-time method
  
  if ((cover.prob > 0) & (cover.prob < 1)) {
    
    #Set element sequence for distribution (if not specified by user)
    if (is.null(E)) {
      #Finite support (take values from left to right)
      if ((supp.min >  -Inf)&(supp.max <  Inf)) {
        ELEM <- function(i) { supp.min - 1 + i; } }
      #Left-bounded support (take values from left to right)
      if ((supp.min >  -Inf)&(supp.max == Inf)) {
        ELEM <- function(i) { supp.min - 1 + i; } }
      #Right-bounded support (take values from right to left)
      if ((supp.min == -Inf)&(supp.max <  Inf)) {
        ELEM <- function(i) { supp.max + 1 - i; } }
      #Unbounded support (take values oscillating out from zero)
      if ((supp.min == -Inf)&(supp.max == Inf)) {
        ELEM <- function(i) { ifelse(i%%2 == 1, floor(i/2), -i/2); } } } else { ELEM <- E }
    
    #Start by computing a minimum covering set
    #The complex vector HHH keeps track of the highest density values and their elements
    #These are represented as complex numbers where the real part is the density and the imaginary part is the element
    INDEX <- 0;
    PROB  <- 0;
    HHH <- rep(NA_complex_, 128)    #Pre-allocating to avoid multiple resizes
    while(PROB < cover.prob) {
      INDEX   <- INDEX + 1;
      if (INDEX >= length(HHH)) { length(HHH) <- 2*length(HHH) }
      ITEM    <- ELEM(INDEX);
      NEWPROB <- f(ITEM);
      HHH[INDEX] <- NEWPROB + ITEM*1i;
      PROB  <- PROB + NEWPROB; }
    OUTPROB <- 1 - PROB;
    HHH     <- sort(HHH, decreasing = FALSE);                      #Sorting step
    while (PROB - Re(HHH[1]) >= cover.prob) {
      PROB  <- PROB - Re(HHH[1]);
      HHH   <- HHH[-1]; }                                          #Trim step
    MINPROB <- Re(HHH[1]);
    
    #Iteratively update the matrix HHH until it contains the values in the HDR
    while (MINPROB < OUTPROB) {
      INDEX   <- INDEX + 1;
      ITEM    <- ELEM(INDEX);
      NEWPROB <- f(ITEM);
      OUTPROB <- OUTPROB - NEWPROB;
      if (NEWPROB > MINPROB) {
        IND <- findInterval(NEWPROB, Re(HHH))+1;
        hhh <- complex(length.out = length(HHH)+1);
        hhh[-IND] <- HHH;
        hhh[IND]  <- NEWPROB + ITEM*1i;
        PROB  <- PROB + NEWPROB;
        while (PROB - Re(hhh[1]) >= cover.prob) {
          PROB  <- PROB - Re(hhh[1]);
          hhh   <- hhh[-1]; }                                      #Trim step
        HHH <- hhh;
        MINPROB <- Re(HHH[1]); } }
    
    #Set the HDR
    HDR <- sets::as.interval(sets::set(Im(HHH)));
    attr(HDR, 'probability') <- sum(Re(HHH)); }
  
  #Add method attribute
  attr(HDR, 'method') <- paste0('Computed using discrete optimisation with minimum coverage probability = ', sprintf(100*cover.prob, fmt = '%#.2f'), '%');
  
  #Add class and attributes
  class(HDR) <- c('hdr', class(HDR));
  attr(HDR, 'distribution') <- distribution;
  
  HDR; }