UNSPECIFIED_LABEL='an unspecified input distribution'

#' Custom HDR Functions
#' 
#' 
#' @param cover.prob The probability coverage for the HDR (scalar between zero and one).  The significance level for the HDR i is \code{1-cover.prob}.  
#' @param Q a CDF of a distribution
#' @param f a PDF of a distribution
#' @param u a log-derivative of f
#' @param decreasing Direction of monotone distribution
#' @param distribution a label
#' @param ... Arguments for Q, f and u
#' @inheritParams checkIterArgs 
#' @return An interval object with classes \code{hdr} and \code{interval} containing the highest density region and related information.
#' 
#' @rdname custom
HDR.monotone <- function(cover.prob, Q, decreasing = TRUE, distribution = UNSPECIFIED_LABEL, ...) {
  hdr(cover.prob = cover.prob, modality=monotone, Q=Q, distribution=distribution, decreasing=decreasing, ...)
}



#' @rdname custom
HDR.unimodal <- function(cover.prob, Q, f = NULL, u = NULL, distribution = UNSPECIFIED_LABEL, 
                         ...,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  hdr(cover.prob = cover.prob, modality=unimodal, Q=Q, f=f, u=u, distribution=distribution,
      ...,
      gradtol = gradtol, steptol = steptol, iterlim = iterlim)
}

#' @rdname custom
HDR.bimodal <- function(cover.prob, Q, f = NULL, u = NULL, distribution = UNSPECIFIED_LABEL, 
                         ...,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  hdr(cover.prob = cover.prob, modality=bimodal, Q=Q, f=f, u=u, distribution=distribution,
      ...,
      gradtol = gradtol, steptol = steptol, iterlim = iterlim)
}

#' @rdname custom
HDR.discrete.unimodal <- function(cover.prob, Q, f = NULL, u = NULL, distribution = UNSPECIFIED_LABEL, 
                         ...,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  hdr(cover.prob = discrete.unimodal, modality=monotone, Q=Q, f=f, u=u, distribution=distribution,
      ...,
      gradtol = gradtol, steptol = steptol, iterlim = iterlim)
}

#' @rdname custom
HDR.discrete <- function(cover.prob, Q, f = NULL, u = NULL, distribution = UNSPECIFIED_LABEL, 
                                  ...,
                                  gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  hdr(cover.prob = discrete, modality=discrete, Q=Q, f=f, u=u, distribution=distribution,
      ...,
      gradtol = gradtol, steptol = steptol, iterlim = iterlim)
}




