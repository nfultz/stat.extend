
hdr <- function(cover.prob, modality, Q, distribution, ...) {
  if (!is.numeric(cover.prob))   { stop('Error: cover.prob should be numeric') }
  if (length(cover.prob) != 1)   { stop('Error: cover.prob should be a single value'); }
  if (cover.prob < 0)            { stop('Error: cover.prob is negative'); }
  if (cover.prob > 1)            { stop('Error: cover.prob is greater than one'); }

  #Compute the HDR in trivial cases where cover.prob is 0 or 1
  #When cover.prob = 1 the HDR is the support of the distribution
  if (cover.prob == 1) {
    Q <- partial(Q, ...);
    HDR <- structure(sets::interval(l = Q(0), r = Q(1), bounds = 'closed'),
                     method = NA_character_);}

  #When cover.prob = 0 the HDR is the empty set
  if (cover.prob == 0) {
    HDR <- structure(sets::interval(), method=NA_character_);}

  #Compute the HDR in non-trivial cases where 0 < (1-cover.prob) < 1
  if (((1-cover.prob) > 0) && ((1-cover.prob) < 1)) {
    HDR <- modality(cover.prob=cover.prob, Q=Q, ...);}

  HDR <- structure(HDR,
                   class = c('hdr',class(HDR)),
                   probability = attr(HDR, "probability") %||% cover.prob,
                   distribution = distribution);

  HDR;
}

# f is not used below, but needed to capture for distributions that switch modalities by param 
monotone <- function(cover.prob, Q, f=NULL, decreasing = TRUE, ...) {
  Q <- partial(Q, ...);

  L <- if (decreasing) Q(0)       else Q((1-cover.prob));
  U <- if (decreasing) Q(1-(1-cover.prob)) else Q(1);

  HDR <- structure(sets::interval(l = L, r = U, bounds = 'closed'),
                   method= 'Computed using monotone optimisation');

  HDR; }


unimodal <- function(cover.prob, Q, f = NULL, u = NULL, ...,
                     gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  checkIterArgs(gradtol, steptol, iterlim)

  # Capture distribution params
  Q <- partial(Q, ...);
  if(is.function(f)) f <- partial(f, ...);
  if(is.function(u)) u <- partial(u, ...);

  #Computation is done using nonlinear optimisation using nlm


  #Set objective function
  WW <- function(phi) {

    #Set parameter functions
    T0 <- (1-cover.prob)/(1+exp(-phi));
    T1 <- T0/(1+exp(phi));
    T2 <- T1*((1-exp(2*phi))/(1+2*exp(phi)+exp(2*phi)));

    #Set interval bounds and objective
    L  <- Q(T0);
    U  <- Q(T0 + 1 - (1-cover.prob));
    W0 <- U - L;

    #Set gradient and Hessian of objective (if able)
    if (is.function(f)) {
      attr(W0, 'gradient') <- T1*(1/f(U) - 1/f(L)); }
    if (is.function(f) & is.function(u)) {
      attr(W0, 'hessian')  <- T2*(1/f(U) - 1/f(L)) +
        T1^2*(u(L)/(f(L)^2) - u(U)/(f(U)^2)); }

    W0; }

  #Compute the HDR
  #The starting value for the parameter phi is set to zero
  #This is the exact optima in the case of a symmetric distribution
  OPT <- nlm(WW, p = 0,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  TT <- (1-cover.prob)/(1+exp(-OPT$estimate));
  L   <- Q(TT);
  U   <- Q(TT + 1 - (1-cover.prob));

  #Add the description of the method
  METHOD <- ifelse((OPT$iterations == 1),
                   paste0('Computed using nlm optimisation with ',
                          OPT$iterations, ' iteration (code = ', OPT$code, ')'),
                   paste0('Computed using nlm optimisation with ',
                          OPT$iterations, ' iterations (code = ', OPT$code, ')'));

  HDR <- structure(sets::interval(l = L, r = U, bounds = 'closed'),
                   method = METHOD);

  HDR; }


bimodal <- function(cover.prob, Q, f = NULL, u = NULL, ...,
                    distribution = 'an unspecified input distribution',
                    gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  checkIterArgs(gradtol, steptol, iterlim);

  # Capture distribution params
  Q <- partial(Q, ...);
  if(is.function(f)) f <- partial(f, ...);
  if(is.function(u)) u <- partial(u, ...);


  #Computation is done using nonlinear optimisation using nlm

  #Set objective function
  WW <- function(phi) {

    #Set parameter functions
    T0 <- (1-(1-cover.prob))/(1+exp(-phi));
    T1 <- T0/(1+exp(phi));
    T2 <- T1*((1-exp(2*phi))/(1+2*exp(phi)+exp(2*phi)));

    #Set interval bounds and objective
    L  <- Q(T0);
    U  <- Q(T0 + (1-cover.prob));
    W0 <- 1 - U + L;

    #Set gradient of objective (if able)
    if (is.function(f)) {
      attr(W0, 'gradient') <- - T1*(1/f(U) - 1/f(L)); }

    #Set Hessian of objective (if able)
    if (is.function(f) & is.function(u)) {
      attr(W0, 'hessian')  <- - T2*(1/f(U) - 1/f(L)) -
        T1^2*(u(L)/(f(L)^2) - u(U)/(f(U)^2)); }

    W0; }

  #Compute the HDR
  #The starting value for the parameter phi is set to zero
  #This is the exact optima in the case of a symmetric distribution
  OPT <- nlm(f = WW, p = 0,
             gradtol = gradtol, steptol = steptol, iterlim = iterlim);
  TT <- (1-cover.prob)/(1+exp(-OPT$estimate));
  L  <- Q(TT*(1-(1-cover.prob))/(1-cover.prob));
  U  <- Q(TT*(1-(1-cover.prob))/(1-cover.prob) + (1-cover.prob));
  HDR1 <- sets::interval(l = Q(0), r = L, bounds = 'closed');
  HDR2 <- sets::interval(l = U, r = Q(1), bounds = 'closed');

  #Add the description of the method
  METHOD <- ifelse((OPT$iterations == 1),
                   paste0('Computed using nlm optimisation with ',
                          OPT$iterations, ' iteration (code = ', OPT$code, ')'),
                   paste0('Computed using nlm optimisation with ',
                          OPT$iterations, ' iterations (code = ', OPT$code, ')'));

  HDR  <- structure(sets::interval_union(HDR1, HDR2), method=METHOD);

  HDR; }


discrete.unimodal <- function(cover.prob, Q, F, f = NULL, s = NULL, ...,
                              gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  checkIterArgs(gradtol, steptol, iterlim);

  Q <- partial(Q, ...);
  F <- partial(F, ...);


  #Compute the HDR

  MIN <- Q(0);
  MAX <- Q((1-cover.prob));
  TT  <- F(MIN:MAX);
  W   <- rep(NA, length(TT));
  P   <- rep(NA, length(TT));
  for (L in MIN:MAX) { LP     <- ifelse(L > MIN, F(L-1), 0);
  U      <- Q(LP+1-(1-cover.prob));
  W[L-MIN+1] <- U-L+1;
  P[L-MIN+1] <- F(U) - LP; }
  for (i in 1:length(TT)) { if (W[i] != min(W)) P[i] <- 0; }

  L   <- which.max(P)+MIN-1;
  U   <- L+W[L-MIN+1]-1;
  HDR <- structure(sets::integers(l = L, r = U),
                   probability = F(U) - ifelse(L > 0, F(L-1), 0),
                   method = paste0('Computed using discrete optimisation with minimum coverage probability = ',
                                   sprintf(100*(1-(1-cover.prob)), fmt = '%#.2f'), '%'))

  HDR; }


discrete <- function(cover.prob, Q, f, ...) {

    # Capture distribution params
    Q <- partial(Q, ...);
    f <- partial(f, ...);  
    
    #Compute the HDR
    alpha <- 0.5*(1-cover.prob);
    S     <- 0L;
    PCUT  <- 0L;
    while (S <= 1-PCUT) {
      MIN <- Q(alpha/2);
      MAX <- Q(1- alpha/2);
      VALS <- seq(from = MIN, to = MAX);
      PVEC <- f(VALS);
      S    <- sum(PVEC);
      if (S < cover.prob)  {
        PCUT <- 0L; } else {
          SORT <- sort(PVEC, decreasing = TRUE);
          NN   <- 0L;
          PP   <- 0L;
          while (PP < cover.prob) {
            NN   <- NN+1;
            PCUT <- SORT[NN];
            PP   <- PP+PCUT; } }
      alpha <- 0.5*alpha; }
    HDRVALS <- sort((MIN-1) + order(PVEC, decreasing = TRUE)[1:NN],
                    decreasing = FALSE);
    if (NN == (max(HDRVALS)-min(HDRVALS)+1)) {
      HDR <- sets::as.interval(HDRVALS); } else {
        HDR <- sets::as.set(HDRVALS); }
    attr(HDR, 'probability') <- PP; 
  
  #Add method attribute
  attr(HDR, 'method') <- paste0('Computed using discrete optimisation with minimum coverage probability = ', sprintf(100*cover.prob, fmt = '%#.2f'), '%');
  
  #Add class and attributes
  class(HDR) <- c('hdr', class(HDR));

  HDR; }

# Version 2 of discrete, more efficient implementation
discrete2 <- function(cover.prob, f, ..., supp.min = -Inf, supp.max = Inf, E = NULL,
                         distribution = 'an unspecified input distribution') {
  
  #Check inputs
  if (!is.numeric(supp.min))   { stop('Error: supp.min should be numeric') }
  if (length(supp.min) != 1)   { stop('Error: supp.min should be a single value'); }
  if (!is.numeric(supp.max))   { stop('Error: supp.max should be numeric') }
  if (length(supp.max) != 1)   { stop('Error: supp.max should be a single value'); }
  supp.min <- ceiling(supp.min);
  supp.max <- floor(supp.max);
  if (supp.max < supp.min)       { stop('Error: specified supp.min and supp.max give empty support'); }
  
  #############
  
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
                      '--- This might not be the smallest HDR for the distribution');
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
        ELEM <- function(i) { ifelse(i%%2 == 1, floor(i/2), -i/2); } } }
    
    #Start by computing a minimum covering set
    #The matrix HHH keeps track of the highest density values
    INDEX <- 0;
    PROBS <- 0;
    while (PROBS < cover.prob) {
      INDEX <- INDEX + 1;
      PROBS <- PROBS + f(ELEM(INDEX)); }
    ELEMS <- rep(0, INDEX);
    PROBS <- rep(0, INDEX);
    for (i in 1: INDEX) {
      ELEMS[i] <- ELEM(i);
      PROBS[i] <- f(ELEM(i)); }
    HHH  <- matrix(c(ELEMS[order(PROBS, decreasing = TRUE)], sort(PROBS, decreasing = TRUE)),
                   byrow = FALSE, ncol = 2);
    OUTPROB <- 1 - sum(HHH[,2]);
    LAST <- sum(cumsum(HHH[,2]) < cover.prob) + 1;
    HHH  <- HHH[1:LAST, ];
    MINPROB <- min(HHH[,2]);
    colnames(HHH) <- c('Value', 'Probability');
    
    #Iteratively update the matrix HHH until it contains the values in the HDR
    while (MINPROB < OUTPROB) {
      INDEX   <- INDEX + 1;
      NEWPROB <- f(ELEM(INDEX));
      OUTPROB <- OUTPROB - NEWPROB;
      if (NEWPROB > MINPROB) {
        IND <- 1 + sum(NEWPROB <= HHH[,2]);
        hhh <- matrix(NA, nrow = nrow(HHH)+1, ncol = 2);
        hhh[IND,  ] <- c(ELEM(INDEX), NEWPROB);
        hhh[-IND, ] <- HHH;
        LAST <- sum(cumsum(hhh[,2]) < cover.prob) + 1;
        hhh  <- hhh[1:LAST, ];
        colnames(HHH) <- c('Value', 'Probability');
        HHH <- hhh;
        MINPROB <- min(HHH[,2]); } }
    
    #Set the HDR
    HDR <- sets::as.interval(sets::set(HHH[,1]));
    attr(HDR, 'probability') <- sum(HHH[,2]); }
  
  #Add method attribute
  attr(HDR, 'method') <- paste0('Computed using discrete optimisation with minimum coverage probability = ', sprintf(100*cover.prob, fmt = '%#.2f'), '%');
  
  #Add class and attributes
  class(HDR) <- c('hdr', class(HDR));
  attr(HDR, 'distribution') <- distribution;
  
  HDR; }





