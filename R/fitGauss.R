# Assessing whether a peak is a Gaussian peak
# Barry Song
# 250707

#' @title fitGauss
#' @description
#' Assessing whether a peak is a Gaussian peak
#'
#' @param int `numeric`, intensity of a peak.
#'
#' @returns `numeric(1)`, r2 value
#' @export
#'
#' @examples
#' fitGauss(int)
fitGauss <- function(int){
  md <- max(int)
  d1 <- int / md
  lri <- which(d1 > 0.5)
  li <- lri[1]
  ri <- tail(lri, 1)
  fwhm <- ri - li
  c <- fwhm / 2.355 # fwhm = 2 * sqrt(2Ln(2)) = 2.355 * c
  x <- 1:length(d1)
  # Gaussian fitting
  fit <- try(stats::nls(d1 ~ a*exp(-((x-b)^2)/(2*c^2)),
                        start = list(a = max(d1),
                                     b = which.max(d1),
                                     c = c)), silent = TRUE)
  if(class(fit) != "try-error"){
    d1_fit <- stats::predict(fit, x = x)
    r2_fit <- round(1 - sum((d1 - d1_fit)^2) / sum((d1 - mean(d1_fit))^2), 2)
  }else{
    r2_fit <- 0
  }
  return(r2_fit)
}
