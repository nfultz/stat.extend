#' Reformat HDRs and confidence intervals objects
#'
#' This function reformats HDRs and confidence intervals back and forth between set format and data frame format.  If the object is a 'hdr'
#' object (HDR presented as a set) it is reformatted into a 'hdr.df' object (HDR presented as a data frame) and *vice versa*.  If the object
#' is a 'ci' object (confidence interval as a set) it is reformatted into a 'ci.df' object (confidence interval presented as a data frame)
#' and *vice versa*.  All attributes and information is preserved when changing formats.  If the object is not of a recognised kind (or is
#' of multiple recognised kinds) then it is returned unchanged and the function gives a warning.
#'
#' @param OBJ An object to reformat (either a HDR or a confidence interval)
#' @return Returns the reformatted object

reformat <- function(OBJ) {
  
  HDR     <- ('hdr'    %in% class(OBJ));
  HDR.DF  <- ('hdr.df' %in% class(OBJ));
  CI      <- ('ci'     %in% class(OBJ));
  CI.DF   <- ('ci.df'  %in% class(OBJ));
  
  if (HDR + HDR.DF + CI + CI.DF == 0) {
    warning('object to reformat is not recognised --- returning object unchanged');
    return(OBJ); }
  
  if (HDR + HDR.DF + CI + CI.DF > 1) {
    warning('object to reformat has multiple recognised formats --- returning object unchanged');
    return(OBJ); }
  
  ATTR <- attributes(OBJ);
  
  if (HDR) {
    VEC  <- unlist(OBJ);
    m    <- ceiling(length(VEC)/4);
    OUT  <- data.frame(Lower  = rep(0, m),    Upper  = rep(0, m),
                       LC     = character(m), RC     = character(m),
                       stringsAsFactors = FALSE);
    for (i in 1:m) {
      OUT$Lower[i]  <- VEC[1 + (i-1)*4];
      OUT$Upper[i]  <- VEC[2 + (i-1)*4];
      OUT$LC[i]     <- ifelse(VEC[3 + (i-1)*4] == 1, 'closed', 'open');
      OUT$RC[i]     <- ifelse(VEC[4 + (i-1)*4] == 1, 'closed', 'open'); }
    rownames(OUT)   <- sprintf('Interval[%s]', 1:m);
    class(OUT)      <- c('data.frame', 'hdr.df');
    attr(OUT, 'domain')       <- ATTR$domain;
    attr(OUT, 'method')       <- ATTR$method;
    attr(OUT, 'probability')  <- ATTR$probability;
    attr(OUT, 'distribution') <- ATTR$distribution; }
  
  if (HDR.DF) {
    m <- nrow(OBJ);
    OUT <- sets::interval(l = 0, r = 0, bounds = '()');
    for (i in 1:m) {
      BOUNDS <- ifelse(OBJ[i, 3] == 'closed',
                       ifelse(OBJ[i, 4] == 'closed', '[]', '[)'),
                       ifelse(OBJ[i, 4] == 'closed', '(]', '()'));
      REGION <- sets::interval(l = OBJ[i, 1], r = OBJ[i, 2], bounds = BOUNDS);
      OUT    <- sets::interval_union(OUT, REGION); }
    class(OUT) <- c('hdr', 'interval');
    attr(OUT, 'domain')       <- ATTR$domain;
    attr(OUT, 'method')       <- ATTR$method;
    attr(OUT, 'probability')  <- ATTR$probability;
    attr(OUT, 'distribution') <- ATTR$distribution; }
  
  if (CI) {
    VEC  <- unlist(OBJ);
    m    <- ceiling(length(VEC)/4);
    OUT  <- data.frame(Lower  = rep(0, m),    Upper  = rep(0, m),
                       LC     = character(m), RC     = character(m),
                       stringsAsFactors = FALSE);
    for (i in 1:m) {
      OUT$Lower[i]  <- VEC[1 + (i-1)*4];
      OUT$Upper[i]  <- VEC[2 + (i-1)*4];
      OUT$LC[i]     <- ifelse(VEC[3 + (i-1)*4] == 1, 'closed', 'open');
      OUT$RC[i]     <- ifelse(VEC[4 + (i-1)*4] == 1, 'closed', 'open'); }
    rownames(OUT)   <- sprintf('Interval[%s]', 1:m);
    class(OUT)      <- c('data.frame', 'ci.df');
    attr(OUT, 'domain')       <- ATTR$domain;
    attr(OUT, 'method')       <- ATTR$method;
    attr(OUT, 'data')         <- ATTR$data;
    attr(OUT, 'confidence')   <- ATTR$confidence;
    attr(OUT, 'parameter')    <- ATTR$parameter; }
  
  if (CI.DF) {
    m <- nrow(OBJ);
    OUT <- sets::interval(l = 0, r = 0, bounds = '()');
    for (i in 1:m) {
      BOUNDS <- ifelse(OBJ[i, 3] == 'closed',
                       ifelse(OBJ[i, 4] == 'closed', '[]', '[)'),
                       ifelse(OBJ[i, 4] == 'closed', '(]', '()'));
      REGION <- sets::interval(l = OBJ[i, 1], r = OBJ[i, 2], bounds = BOUNDS);
      OUT    <- sets::interval_union(OUT, REGION); }
    class(OUT) <- c('ci', 'interval');
    attr(OUT, 'domain')       <- ATTR$domain;
    attr(OUT, 'method')       <- ATTR$method;
    attr(OUT, 'data')         <- ATTR$data;
    attr(OUT, 'confidence')   <- ATTR$confidence;
    attr(OUT, 'parameter')    <- ATTR$parameter; }
  
  OUT; }

#' @export
#' @rdname reformat
#' @param x an R object
#' @param ... unused
as.data.frame.hdr <- function(x, ...) reformat(x)

#' @export
#' @rdname reformat
as.data.frame.ci <- function(x, ...) reformat(x)

#' Used to inherit roxygen docs
#'
#' @param gradtol Parameter for the nlm optimisation - a positive scalar giving the tolerance at which the scaled gradient is considered close enough to zero to terminate the algorithm (see [\code{nlm} doccumentation](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/nlm.html)).
#' @param steptol Parameter for the nlm optimisation - a positive scalar providing the minimum allowable relative step length (see [\code{nlm} doccumentation](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/nlm.html)).
#' @param iterlim Parameter for the nlm optimisation - a positive integer specifying the maximum number of iterations to be performed before the program is terminated (see [\code{nlm} doccumentation](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/nlm.html)).
#' @keywords internal
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

partial <- function(FUN, ...) {
  args <- list(...);
  args <- args[names(args) %in% names(formals(FUN))];
  
  QQ <- function(L) { do.call("FUN", c(list(L), args)); }
  
  QQ;}

`%||%` <- function(l, r) {if (is.null(l)) r else l}

