pkgname <- "stat.extend"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('stat.extend')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("HDR")
### * HDR

flush(stderr()); flush(stdout())

### Name: HDR
### Title: Highest density region (HDR)
### Aliases: HDR HDR.norm HDR.lnorm HDR.t HDR.cauchy HDR.f HDR.beta
###   HDR.chisq HDR.gamma HDR.weibull HDR.exp HDR.unif HDR.hyper HDR.geom
###   HDR.binom HDR.pois HDR.nbinom

### ** Examples

HDR.norm(.05)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
