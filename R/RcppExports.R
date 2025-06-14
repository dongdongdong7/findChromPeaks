# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

C_noiseEs <- function(intensity, mag = 3.0) {
    .Call(`_findChromPeaks_C_noiseEs`, intensity, mag)
}

C_localMaxima <- function(y, halfWindowSize) {
    .Call(`_findChromPeaks_C_localMaxima`, y, halfWindowSize)
}

C_DescendMinTol <- function(d, startpos, maxDescOutlier) {
    .Call(`_findChromPeaks_C_DescendMinTol`, d, startpos, maxDescOutlier)
}

C_RectUnique <- function(m, order, nrow, ncol, xdiff, ydiff) {
    .Call(`_findChromPeaks_C_RectUnique`, m, order, nrow, ncol, xdiff, ydiff)
}

