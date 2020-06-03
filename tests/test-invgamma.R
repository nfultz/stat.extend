
stopifnot(all.equal(
  stat.extend::HDR.invchisq(.5, 10, 0),
  structure(list(list(l = 0.0607447079278582, r = 0.11865708180454, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 8 iterations (code = 1)", probability = 0.5, distribution = "inverse chi-squared distribution with 10 degrees-of-freedom")))



stopifnot(all.equal(
  stat.extend::HDR.invexp(.5, rate=1),
  structure(list(list(l = 0.224157477631746, r = 1.49184317616973, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 11 iterations (code = 1)", probability = 0.5, distribution = "inverse exponential distribution with rate = 1")
))


stopifnot(all.equal(
  stat.extend::HDR.invgamma(.5, shape=1),
  structure(list(list(l = 0.224157477631746, r = 1.49184317616973, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 11 iterations (code = 1)", probability = 0.5, distribution = "inverse gamma distribution with shape = 1 and scale = 1")
))