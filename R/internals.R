
hdr <- function(alpha, modality, Q, distribution, ...) {
  if (!is.numeric(alpha))   { stop('Error: alpha should be numeric') }
  if (length(alpha) != 1)   { stop('Error: alpha should be a single value'); }
  if (alpha < 0)            { stop('Error: alpha is negative'); }
  if (alpha > 1)            { stop('Error: alpha is greater than one'); }

  #Compute the HDR in trivial cases where alpha is 0 or 1
  #When alpha = 0 the HDR is the support of the distribution
  if (alpha == 0) {
    Q <- partial(Q, ...);
    HDR <- structure(sets::interval(l = Q(0), r = Q(1), bounds = 'closed'),
                     method = NA_character_);}

  #When alpha = 1 the HDR is the empty region
  if (alpha == 1) {
    HDR <- structure(sets::interval(), method=NA_character_);}

  #Compute the HDR in non-trivial cases where 0 < alpha < 1
  if ((alpha > 0) && (alpha < 1)) {
    HDR <- switch(modality,
                  monotone = monotone(alpha=alpha, Q=Q, ...),
                  unimodal = unimodal(alpha=alpha, Q=Q, ...),
                  bimodal  = bimodal(alpha=alpha, Q=Q, ...),
                  discrete.unimodal = discrete.unimodal(alpha=alpha, Q=Q, ...)); }

  HDR <- structure(HDR,
                   class = c('hdr','interval'),
                   probability = 1 - alpha,
                   distribution = distribution);

  HDR;
}


monotone <- function(alpha, Q, decreasing = TRUE, ...) {
  Q <- partial(Q, ...);

  L <- if (decreasing) Q(0)       else Q(alpha);
  U <- if (decreasing) Q(1-alpha) else Q(1);

  HDR <- structure(sets::interval(l = L, r = U, bounds = 'closed'),
                   method= 'Computed using monotone optimisation');

  HDR; }


unimodal <- function(alpha, Q, f = NULL, u = NULL, ...,
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
    T0 <- alpha/(1+exp(-phi));
    T1 <- T0/(1+exp(phi));
    T2 <- T1*((1-exp(2*phi))/(1+2*exp(phi)+exp(2*phi)));

    #Set interval bounds and objective
    L  <- Q(T0);
    U  <- Q(T0 + 1 - alpha);
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
  TT <- alpha/(1+exp(-OPT$estimate));
  L   <- Q(TT);
  U   <- Q(TT + 1 - alpha);

  #Add the description of the method
  METHOD <- ifelse((OPT$iterations == 1),
                   paste0('Computed using nlm optimisation with ',
                          OPT$iterations, ' iteration (code = ', OPT$code, ')'),
                   paste0('Computed using nlm optimisation with ',
                          OPT$iterations, ' iterations (code = ', OPT$code, ')'));

  HDR <- structure(sets::interval(l = L, r = U, bounds = 'closed'),
                   method = METHOD);

  HDR; }


bimodal <- function(alpha, Q, f = NULL, u = NULL, ...,
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
    T0 <- (1-alpha)/(1+exp(-phi));
    T1 <- T0/(1+exp(phi));
    T2 <- T1*((1-exp(2*phi))/(1+2*exp(phi)+exp(2*phi)));

    #Set interval bounds and objective
    L  <- Q(T0);
    U  <- Q(T0 + alpha);
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
  TT <- alpha/(1+exp(-OPT$estimate));
  L  <- Q(TT*(1-alpha)/alpha);
  U  <- Q(TT*(1-alpha)/alpha + alpha);
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





alphaparam <- function(alpha) {NULL};
iterparams <- function(gradtol, steptol, iterlim) {NULL};

checkIterArgs <- function(gradtol, steptol, iterlim) {
  if (!is.numeric(gradtol)) { stop('Error: gradtol should be numeric') }
  if (length(gradtol) != 1) { stop('Error: gradtol should be a single value'); }
  if (gradtol <= 0)         { stop('Error: gradtol should be positive'); }
  if (!is.numeric(steptol)) { stop('Error: steptol should be numeric') }
  if (length(steptol) != 1) { stop('Error: steptol should be a single value'); }
  if (steptol <= 0)         { stop('Error: steptol should be positive'); }
  if (!is.numeric(iterlim)) { stop('Error: iterlim should be numeric') }
  if (length(iterlim) != 1) { stop('Error: iterlim should be a single value'); }
  if (iterlim <= 0)         { stop('Error: iterlim should be positive'); }
}

partial <- function(f, ...) {
  args <- list(...);
  args <- args[names(args) %in% names(formals(f))];

  QQ <- function(L) { do.call("f", c(L, args)); }

  QQ;}

`%||%` <- function(l, r) {if (is.null(l)) r else l}



discrete.unimodal <- function(alpha, Q, F, f = NULL, s = NULL, ...,
                              gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {

  #Check inputs
  checkIterArgs();

  Q <- partial(Q, ...);
  if(is.function(f)) f <- partial(f, ...);
  if(is.function(u)) u <- partial(u, ...);


  #Compute the HDR

  MIN <- Q(0);
  MAX <- Q(alpha);
  TT  <- F(MIN:MAX);
  W   <- rep(NA, length(TT));
  P   <- rep(NA, length(TT));
  for (L in MIN:MAX) { LP     <- ifelse(L > MIN, F(L-1), 0);
  U      <- Q(LP+1-alpha);
  W[L-MIN+1] <- U-L+1;
  P[L-MIN+1] <- F(U) - LP; }
  for (i in 1:length(TT)) { if (W[i] != min(W)) P[i] <- 0; }

  L   <- which.max(P)+MIN-1;
  U   <- L+W[L-MIN+1]-1;
  HDR <- structure(sets::integers(l = L, r = U),
                   probability = F(U) - ifelse(L > 0, F(L-1), 0),
                   method = paste0('Computed using discrete optimisation with minimum coverage probability = ',
                                   sprintf(100*(1-alpha), fmt = '%#.2f'), '%'))

  HDR; }
