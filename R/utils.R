requireVGAM <- function(){
  #Check for required package
  if (!requireNamespace('VGAM', quietly = TRUE)) {
    stop('Package \'VGAM\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
}

requireExtraDist <- function() {
  #Check for required package
  if (!requireNamespace('extraDistr', quietly = TRUE)) {
    stop('Package \'extraDistr\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
}


requireInvgamma <- function() {
  #Check for required package
  if (!requireNamespace('invgamma', quietly = TRUE)) {
    stop('Package \'invgamma\' is required for this function; if you install that package you can run this function.', call. = FALSE) }
}