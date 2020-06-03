
stopifnot(all.equal(
  stat.extend::HDR.benini(.5, 1, 1),
  structure(list(list(l = 1.24098584574423, r = 2.43041994178578, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 9 iterations (code = 1)", probability = 0.5, distribution = "Benini distribution with shape = 1 and scale = 1")
))


stopifnot(all.equal(
  stat.extend::HDR.frechet(.5, 1, 2, 3),
  structure(list(list(l = 3.44831495526452, r = 5.98368635234049, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 7 iterations (code = 1)", probability = 0.5, distribution = "Frechet distribution with shape = 1, scale = 2 and location = 3")
))

stopifnot(all.equal(
  stat.extend::HDR.gengamma(.5, 1, 2),
  structure(list(list(l = 0, r = 1.67834699001666, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                                     "interval"), domain = "R", method = "Computed using monotone optimisation", probability = 0.5, distribution = "generalised gamma (Stacy) distribution with shape1 = 1, shape2 = 2 and scale = 1")
))

stopifnot(all.equal(
  stat.extend::HDR.gengamma(.5, 2, 2),
  structure(list(list(l = 0.912572951370049, r = 1.56638534772855, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 5 iterations (code = 1)", probability = 0.5, distribution = "generalised gamma (Stacy) distribution with shape1 = 2, shape2 = 2 and scale = 1")
))


stopifnot(all.equal(
  stat.extend::HDR.gumbelII(.5, 2, 2),
  structure(list(list(l = 1.1496471220439, r = 2.58072199114017, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 9 iterations (code = 1)", probability = 0.5, distribution = "Gumbel (Type II) distribution with shape = 2 and scale = 2")
))

stopifnot(all.equal(
  stat.extend::HDR.lgamma(.5, 2, 2, 0),
  structure(list(list(l = 0.302904714056459, r = 2.30417608071809, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 6 iterations (code = 1)", probability = 0.5, distribution = "log-gamma distribution with shape = 2, scale = 2 and location = 0")
))
