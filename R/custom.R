UNSPECIFIED_LABEL='an unspecified input distribution'

#' Highest density region (HDR) for an arbitrary distributions
#' 
#' 
#' @param cover.prob The probability coverage for the HDR (scalar between zero and one).  The significance level for the HDR i is \code{1-cover.prob}.  
#' @param Q an inverse CDF of a distribution
#' @param f a PDF of a distribution
#' @param F a CDF of a distribution
#' @param u a log-derivative of f
#' @param decreasing Direction of monotone distribution
#' @param distribution a label
#' @param ... Arguments for Q, f and u
#' @inheritParams checkIterArgs 
#' @return An interval object with classes \code{hdr} and \code{interval} containing the highest density region and related information.
#' 
#' @seealso HDR.discrete
#' 
#' @examples 
#' 
#' HDR.monotone(.95, Q=qexp)
#' 
#' HDR.unimodal(.95, Q=qnorm)
#' 
#' HDR.bimodal(.95, Q=qbeta, shape1=1/2, shape2=1/2)
#' 
#' HDR.discrete.unimodal(.95, Q=qpois, F=ppois, lambda=1)
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
HDR.discrete.unimodal <- function(cover.prob, Q, F, f = NULL, u = NULL, distribution = UNSPECIFIED_LABEL, 
                         ...,
                         gradtol = 1e-10, steptol = 1e-10, iterlim = 100) {
  hdr(cover.prob = cover.prob, modality=discrete.unimodal, Q=Q, F=F, f=f, u=u, distribution=distribution,
      ...,
      gradtol = gradtol, steptol = steptol, iterlim = iterlim)
}




#' Highest density region (HDR) for an arbitrary discrete distribution
#'
#' This function computes the highest density region (HDR) with support on the integers.  The distribution can be any discrete
#' distribution concentrated on the integers --- it does not have to have any shape properties for the function to work.  
#' 
#' The user must give the density function ```f``` for the distribution.  To improve the search properties of the algorithm, they can
#' also give lower and upper bounds for the support of the distribution if these are available.  (Warning: If the user specifies
#' incorrect bounds on the support, that do not contain the full support of the distribution, then the algorithm may continue to
#' search without end, in which case the function will not terminate.  Similarly, if the user specifies a sequence function E that
#' is not a proper bijection to the integers then the algorithm may continue to search without end, in which case the function will
#' not terminate.) The output of the function is a \code{'hdr'} object containing the HDR for the discrete distribution.
#'
#' @usage \code{HDR.discrete(cover.prob, Q, f, E = NULL, distribution = 'an unspecified input distribution')}
#' @param cover.prob The minimum coverage probability for the region
#' @param f The density (mass) function for the distribution
#' @param supp.min A minimum bound for the support of the distribution
#' @param supp.max A maximum bound for the support of the distribution
#' @param E A bijective function mapping the natural numbers (1,2,3,...) to a set covering the support of the distribution (optional);
#' if included, the algorithm will search the support of the distribution in the order specified by this function; if not included,
#' the algorithm will search the integers in a default order.
#' @return If all inputs are correctly specified (i.e., arguments and parameters are in allowable range)
#' then the output will be a list of class ```hdr``` containing the HDR and related information.
HDR.discrete <- function(cover.prob, f, distribution = 'an unspecified input distribution', 
                         ...,
                         supp.min = -Inf, supp.max = Inf, E = NULL) {
  
  if(identical(cover.prob, 1)) {
    # NOTE: If user asks for 1, set it machine epsilon smaller, so that HDR base function does not return 
    # the support interval, but instead filters out anything having f(0) = 0.
    cover.prob <- 1 - .Machine$double.eps/2
    
  }
  
  hdr(cover.prob = cover.prob, modality=discrete, Q=NULL, distribution=distribution, ...,
        f=f,  supp.min = supp.min, supp.max = supp.max, E = E )
}  