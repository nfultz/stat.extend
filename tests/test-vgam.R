
stopifnot(all.equal(
  stat.extend::HDR.benini(.5, 1, 1),
  structure(list(list(l = 1.24098584574423, r = 2.43041994178578, 
                      lc = TRUE, rc = TRUE)), class = c("hdr", "interval"), domain = "R", method = "Computed using nlm optimisation with 9 iterations (code = 1)", probability = 0.5, distribution = "Benini distribution with shape = 1 and scale = 1")
))

