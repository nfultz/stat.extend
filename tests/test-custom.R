
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
