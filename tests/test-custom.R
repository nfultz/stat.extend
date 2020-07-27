
stopifnot(all.equal(
  stat.extend::HDR.monotone(.95, Q=qexp, 1)[[1]],
  stat.extend::HDR.exp(.95, 1)[[1]]
))

stopifnot(all.equal(
  stat.extend::HDR.unimodal(.95, Q=qnorm, f=dnorm, mean=0, sd=1)[[1]],
  stat.extend::HDR.norm(.95)[[1]]
))

stopifnot(all.equal(
  stat.extend::HDR.bimodal(.95, Q=qbeta, f=dbeta, shape1=.5, shape2=.5)[1:2],
  stat.extend::HDR.beta(.95, .5, .5)[1:2]
))


stopifnot(all.equal(
  stat.extend::HDR.discrete(.95, Q=qpois, f=dpois, lambda=1)[1:2],
  stat.extend::HDR.pois(.95, 1)[1:2]
))

stopifnot(all.equal(
  stat.extend::HDR.discrete.unimodal(.95, Q=qpois, F=ppois, lambda=1)[1:2],
  stat.extend::HDR.pois(.95, 1)[1:2]
))

# Below has 0 prob at 2; shall not throw an error.
Q = function(p) ifelse(p <= .08, 1, ifelse(p <= .14, 3, ifelse(p <= .4, 4, 5)))
f = function(x) c(.08, 0, .06, .26, .6)[floor(x)]

stopifnot(inherits(  
  stat.extend:::HDR.discrete(.95, Q, f),
  "set"
))


