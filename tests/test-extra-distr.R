

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


stopifnot(all.equal(
  stat.extend::HDR.betapr(.5, shape1=1, shape2=4),
  structure(list(list(l = 0, r = 0.189207115002721, lc = TRUE, 
                      rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using monotone optimisation", probability = 0.5, distribution = "beta prime distribution with shape1 = 1 and shape2 = 4")
))

stopifnot(all.equal(
  stat.extend::HDR.betapr(.5, shape1=2, shape2=4),
  structure(list(list(l = 0.0649354493196968, r = 0.492937442908933, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 10 iterations (code = 1)", probability = 0.5, distribution = "beta prime distribution with shape1 = 2 and shape2 = 4")
))

stopifnot(all.equal(
  stat.extend::HDR.benini(.5, 1, 1),
  structure(list(list(l = 1.24098584574423, r = 2.43041994178578, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 9 iterations (code = 1)", probability = 0.5, distribution = "Benini distribution with shape = 1 and scale = 1")
))


stopifnot(all.equal(
  stat.extend::HDR.fatigue(.5, alpha=3, beta=2, mu=1),
  structure(list(list(l = 1.01685734320788, r = 3.00239281275403, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 10 iterations (code = 1)", probability = 0.5, distribution = "fatigue-life (Birnbaum-Saunders) distribution with shape = 3, scale = 2 and location = 1")
))
