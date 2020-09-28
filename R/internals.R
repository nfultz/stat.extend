
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



