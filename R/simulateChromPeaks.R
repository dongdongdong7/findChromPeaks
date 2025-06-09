# Simulate chromatographic peaks
# Barry Song
# 250608

#' @rdname simulateChromPeaks
#' @title Simulate chromatographic peaks
#' @description
#' Simulate three types of chromatographic peaks:
#' gaussian symmetrical peak, tailing peak and leading peak
#' @param rtime `numeric()`, retention time vector of chromatogram
#' @param peakRt `numeric(1)`, retention time of peak apex
#' @param peakWidth `numeric(1)`, peak width
#' @param peakHeight `numeric(1)`, peak height
#'
#' @returns `numeric()`, intensity vector of chromatogram
#' @export
#'
#' @examples
#' rtime <- seq(0, 1000, 0.5)
#' peak1 <- gaussian_peak(rtime, peakRt = 240, peakWidth = 20, peakHeight = 10000)
#' peak2 <- tailing_peak(rtime, peakRt = 480, peakWidth = 25, peakHeight = 1300, tau = 10)
#' peak3 <- leading_peak(rtime, peakRt = 720, peakWidth = 18, peakHeight = 5000, tau = -8)
#' combined <- peak1 + peak2 + peak3
#' baseline <- 100 + 1 * rtime + 100 * sin(rtime/100)
#' final_chrom <- combined + baseline  + abs(rnorm(length(rtime), mean = 0, sd = 25))
#' plot(rtime, final_chrom, type = "l", xlab = "Retention Time (sec)", ylab = "Intensity")
gaussian_peak <- function(rtime, peakRt, peakWidth, peakHeight){
  # peakWidth 现在表示基线峰宽(W)，转换为标准差σ = W/4
  sigma <- peakWidth / 4  # 对于高斯峰，W ≈ 4σ

  peakHeight * exp(-((rtime - peakRt)^2) / (2 * sigma^2))
}

#' @rdname simulateChromPeaks
#' @param tau `numeric(1)`, the degree of tailing or leading of chromatographic peak;
#' for tailing peas, greater than 0; for leading less than 0
tailing_peak <- function(rtime, peakRt, peakWidth, peakHeight, tau = 10){
  sigma <- peakWidth / 4
  a <- (sigma^2 - tau*(rtime - peakRt)) / (sqrt(2)*sigma*tau)
  term1 <- peakHeight * sigma / tau * sqrt(pi/2)
  term2 <- exp(0.5 * (sigma/tau)^2 - (rtime - peakRt)/tau)
  term3 <- 1 - pracma::erf(a)

  term1 * term2 * term3
}

#' @rdname simulateChromPeaks
leading_peak <- function(rtime, peakRt, peakWidth, peakHeight, tau = -10){
  sigma <- peakWidth / 4

  a <- (sigma^2 - tau*(rtime - peakRt)) / (sqrt(2)*sigma*tau)
  term1 <- peakHeight * sigma / abs(tau) * sqrt(pi/2)
  term2 <- exp(0.5 * (sigma/tau)^2 - (rtime - peakRt)/tau)
  term3 <- 1 + pracma::erf(a)

  term1 * term2 * term3
}

