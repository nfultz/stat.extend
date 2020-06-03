
stopifnot(all.equal(
  stat.extend::HDR.arcsine(.5, 0, 1), # should match beta(.5, .5)
  structure(list(list(l = 0, r = 0.146446609406726, lc = TRUE, 
                      rc = TRUE), list(l = 0.853553390593274, r = 1, lc = TRUE, 
                                       rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 0 iterations (code = 1)", probability = 0.5, distribution = "standard arcsine distribution")
))

stopifnot(all.equal(
  stat.extend::HDR.cauchy(1-.05, 1, 2),
  structure(list(list(l = -24.4124094723494, r = 26.4124094723493, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 0 iterations (code = 1)", probability = 0.95, distribution = "Cauchy distribution with location = 1 and scale = 2")
  
))


stopifnot(all.equal(
  stat.extend::HDR.lnorm(1-.05, 1, .3),
  structure(list(list(l = 1.34512113868324, r = 4.58832908980926, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 7 iterations (code = 1)", probability = 0.95, distribution = "log-normal distribution with log-mean = 1 and log-standard deviation = 0.3")
))


stopifnot(all.equal(
  stat.extend::HDR.norm(1-.05,1,3),
  structure(list(list(l = -4.87989195362016, r = 6.87989195362016, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 0 iterations (code = 1)", probability = 0.95, distribution = "normal distribution with mean = 1 and standard deviation = 3")
))


stopifnot(all.equal(
  stat.extend::HDR.t(1-.05, 29, .3),
  structure(list(list(l = -1.73385540922536, r = 2.35964477149638, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 4 iterations (code = 1)", probability = 0.95, distribution = "Student's T distribution with 29 degrees-of-freedom and non-centrality parameter = 0.3")
))




# monotone
stopifnot(all.equal(
  stat.extend::HDR.f(1-.05, 1, 2, 3),
  structure(list(list(l = 0, r = 76.2364236939022, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                                     "interval"), domain = "R", method = "Computed using monotone optimisation", probability = 0.95, distribution = "F distribution with 1 numerator degrees-of-freedom and 2 denominator degrees-of-freedom and non-centrality parameter = 3")
))

# unimodal
stopifnot(all.equal(
  stat.extend::HDR.f(1-.05, 3, 2, 3),
  structure(list(list(l = 6.68324314416018e-13, r = 38.4926454487213, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 1 iteration (code = 1)", probability = 0.95, distribution = "F distribution with 3 numerator degrees-of-freedom and 2 denominator degrees-of-freedom and non-centrality parameter = 3")  
))

#monotone down
stopifnot(all.equal(
  stat.extend::HDR.beta(1-.05, 1, 2),
  structure(list(list(l = 0, r = 0.776393202250021, lc = TRUE, 
                      rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using monotone optimisation", probability = 0.95, distribution = "beta distribution with shape1 = 1 and shape2 = 2")
))

# monotone up
stopifnot(all.equal(
  stat.extend::HDR.beta(1-.05, 2, 1),
  structure(list(list(l = 0.223606797749979, r = 1, lc = TRUE, 
                      rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using monotone optimisation", probability = 0.95, distribution = "beta distribution with shape1 = 2 and shape2 = 1")
))

# beta(1,1) => uniform
stopifnot(all.equal(
  stat.extend::HDR.beta(1-.05, 1, 1),
  structure(list(list(l = 0.025, r = 0.975, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                              "interval"), domain = "R", method = "Computed using nlm optimisation with 0 iterations (code = 1)", probability = 0.95, distribution = "standard uniform distribution")  
))

# unimodal
stopifnot(all.equal(
stat.extend::HDR.beta(1-.05, 2, 2),
structure(list(list(l = 0.0942993240502461, r = 0.905700675949753, 
                    lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 0 iterations (code = 1)", probability = 0.95, distribution = "beta distribution with shape1 = 2 and shape2 = 2")
))

# bimodal
stopifnot(all.equal(
  stat.extend::HDR.beta(1-.05, .1, .1),
  structure(list(list(l = 0, r = 0.361776246451287, lc = TRUE,
                      rc = TRUE), list(l = 0.638223753548714, r = 1, lc = TRUE, 
                                       rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 0 iterations (code = 1)", probability = 0.95, distribution = "beta distribution with shape1 = 0.1 and shape2 = 0.1")
))

# Monotone 
stopifnot(all.equal(
  stat.extend::HDR.chisq(1-.05, 2, .1),
  structure(list(list(l = 0, r = 6.28726305532095, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                                     "interval"), domain = "R", method = "Computed using monotone optimisation", probability = 0.95, distribution = "chi-squared distribution with 2 degrees-of-freedom and non-centrality parameter = 0.1")
))

# Unimodal
stopifnot(all.equal(
  stat.extend::HDR.chisq(1-.05, 200, .1),
  structure(list(list(l = 161.567164972462, r = 239.761750052931, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 4 iterations (code = 1)", probability = 0.95, distribution = "chi-squared distribution with 200 degrees-of-freedom and non-centrality parameter = 0.1")
))



# Gamma(1, x) => Exponential(x)
stopifnot(all.equal(
  stat.extend::HDR.gamma(1-.05, shape = 1, scale=5),
  structure(list(list(l = 0, r = 14.97866136777, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                                   "interval"), domain = "R", method = "Computed using monotone optimisation", probability = 0.95, distribution = "exponential distribution with scale = 5")  
))


# Unimodal
stopifnot(all.equal(
  stat.extend::HDR.gamma(1-.05, shape = 2, scale=5),
  structure(list(list(l = 0.211816667149859, r = 23.8258412369454, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 10 iterations (code = 1)", probability = 0.95, distribution = "gamma distribution with shape = 2 and scale = 5")
))

# Weibull(1, X) => Exponential(x)
stopifnot(all.equal(
  stat.extend::HDR.weibull(1-.05, shape = 1, scale=5),
  structure(list(list(l = 0, r = 14.97866136777, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                                   "interval"), domain = "R", method = "Computed using monotone optimisation", probability = 0.95, distribution = "exponential distribution with scale = 5")
))

# unimodal
stopifnot(all.equal(
  stat.extend::HDR.weibull(1-.05, shape = 2, scale=5),
  structure(list(list(l = 0.390576258175476, r = 8.83948980233535, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 8 iterations (code = 1)", probability = 0.95, distribution = "Weibull distribution with shape = 2 and scale = 5")
))

# Monotone down
stopifnot(all.equal(
  stat.extend::HDR.exp(1-.05, rate=5),
  structure(list(list(l = 0, r = 0.599146454710798, lc = TRUE, 
                      rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using monotone optimisation", probability = 0.95, distribution = "exponential distribution with scale = 0.2")
))

# Unimodal
stopifnot(all.equal(
  stat.extend::HDR.unif(1-.05),
  structure(list(list(l = 0.025, r = 0.975, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                              "interval"), domain = "R", method = "Computed using nlm optimisation with 0 iterations (code = 1)", probability = 0.95, distribution = "continuous uniform distribution with minimum = 0 and maximum = 1")
))

#Discrete unimodal
stopifnot(all.equal(
  stat.extend::HDR.hyper(1-.05, 50, 5, 5),
  structure(list(list(l = 3, r = 5, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                      "interval"), domain = "Z", probability = 0.996406479203372, method = "Computed using discrete optimisation with minimum coverage probability = 95.00%", distribution = "hypergeometric distribution with 50 white balls, 5 black balls, and 5 balls drawn")
))

# Discrete unimodal
stopifnot(all.equal(
  stat.extend::HDR.geom(1-.05, .50),
  structure(list(list(l = 0, r = 4, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                      "interval"), domain = "Z", probability = 0.96875, method = "Computed using discrete optimisation with minimum coverage probability = 95.00%", distribution = "geometric distribution with probability = 0.5")
))

# Discrete unimodal
stopifnot(all.equal(
  stat.extend::HDR.binom(1-.05, 20, .50),
  structure(list(list(l = 6, r = 14, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                       "interval"), domain = "Z", probability = 0.958610534667969, method = "Computed using discrete optimisation with minimum coverage probability = 95.00%", distribution = "binomial distribution with size = 20 and probability = 0.5")
))

# discrete unimodal
stopifnot(all.equal(
  stat.extend::HDR.pois(1-.05, 20,),
  structure(list(list(l = 12, r = 29, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                        "interval"), domain = "Z", probability = 0.956794960887162, method = "Computed using discrete optimisation with minimum coverage probability = 95.00%", distribution = "Poisson distribution with rate = 20")
))

# Discrete unimodal
stopifnot(all.equal(
  stat.extend::HDR.nbinom(1-.05, 10, .20),
  structure(list(list(l = 15, r = 68, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                        "interval"), domain = "Z", probability = 0.951352378174766, method = "Computed using discrete optimisation with minimum coverage probability = 95.00%", distribution = "negative binomial distribution with size = 10 and probability = 0.2")
))

# nbinom "mu mode"
stopifnot(all.equal(
  stat.extend::HDR.nbinom(1-.05, 4, mu=3),
  structure(list(list(l = 0, r = 7, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                    "interval"), domain = "Z", probability = 0.95479730230908, method = "Computed using discrete optimisation with minimum coverage probability = 95.00%", distribution = "negative binomial distribution with size = 4 and mean = 3")
))

# short circuit when alpha = 1
stopifnot(all.equal(
  stat.extend::HDR.norm(1-0, 10, .20),
  structure(list(list(l = -Inf, r = Inf, lc = TRUE, rc = TRUE)), class = c("hdr", 
                                                                           "interval"), domain = "R", method = NA_character_, probability = 1, distribution = "normal distribution with mean = 10 and standard deviation = 0.2")
))

# short circuit when alpha = 0
stopifnot(all.equal(
  stat.extend::HDR.norm(1-1, 10, .20),
  structure(list(), class = c("hdr", "interval"), domain = NA, method = NA_character_, probability = 0, distribution = "normal distribution with mean = 10 and standard deviation = 0.2")
))



