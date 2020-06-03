
stopifnot(all.equal(
  stat.extend::HDR.monotone(.95, Q=qexp, 1)[[1]],
  stat.extend::HDR.exp(.95, 1)[[1]]
))
