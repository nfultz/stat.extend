

stopifnot(all.equal(
  stat.extend::HDR.arcsine(.5, 0, 1), # should match beta(.5, .5)
  structure(list(list(l = 0, r = 0.146446609406726, lc = TRUE, 
                      rc = TRUE), list(l = 0.853553390593274, r = 1, lc = TRUE, 
                                       rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 0 iterations (code = 1)", probability = 0.5, distribution = "standard arcsine distribution")
))



stopifnot(all.equal(
  stat.extend::HDR.invgamma(.5, shape=1),
  structure(list(list(l = 0.224157477631746, r = 1.49184317616973, 
                    lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 11 iterations (code = 1)", probability = 0.5, distribution = "inverse gamma distribution with shape = 1 and scale = 1")
))