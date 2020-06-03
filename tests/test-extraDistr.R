

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
  stat.extend::HDR.fatigue(.5, alpha=3, beta=2, mu=1),
  structure(list(list(l = 1.01685734320788, r = 3.00239281275403, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 10 iterations (code = 1)", probability = 0.5, distribution = "fatigue-life (Birnbaum-Saunders) distribution with shape = 3, scale = 2 and location = 1")
))




stopifnot(all.equal(
  stat.extend::HDR.gompertz(.5, shape = 1, scale = 1),
  structure(list(list(l = 1.5731194125e-11, r = 0.526589034157626, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 35 iterations (code = 1)", probability = 0.5, distribution = "Gompertz distribution with shape = 1 and scale = 1")
  ))


stopifnot(all.equal(
  stat.extend::HDR.gpd(.5, mu = 0, sigma = 1, xi = .5),
  structure(list(list(l = 1.48681067457801e-12, r = 0.828427124750396, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 38 iterations (code = 1)", probability = 0.5, distribution = "generalised Pareto distribution with location = 0, scale = 1 and shape = 0.5")
))

stopifnot(all.equal(
  stat.extend::HDR.huber(.5, mu = 1, sigma = 1, epsilon = .1),
  structure(list(list(l = -5.93163780852227, r = 7.93163780852227, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 0 iterations (code = 1)", probability = 0.5, distribution = "Huber distribution with mean = 1, scale = 1 and cut-point = 0.1")
))



stopifnot(all.equal(
  stat.extend::HDR.kumar(.5, 1, 1),
  structure(list(list(l = 0.25, r = 0.75, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                            "interval"), domain = "R", method = "Computed using nlm optimisation with 0 iterations (code = 1)", probability = 0.5, distribution = "uniform distribution")
))



stopifnot(all.equal(
  stat.extend::HDR.tnorm(.5, 1, 1, 0, 1),
  structure(list(list(l = 0.558229453270944, r = 0.999999999961474, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 32 iterations (code = 1)", probability = 0.5, distribution = "truncated normal distribution with mean = 1 and standard deviation = 1 with min = 0 and max = 1")
))
