# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

C_find_insert_indices <- function(peaks_mz, peaks_rt, ref_mz, ref_rt, tol_mz, rt_tol) {
    .Call(`_findChromPeaks_C_find_insert_indices`, peaks_mz, peaks_rt, ref_mz, ref_rt, tol_mz, rt_tol)
}

C_match_peaks <- function(sample_mz, sample_rt, ref_mz, ref_rt, a, b, rt_diff_tol, ppm) {
    .Call(`_findChromPeaks_C_match_peaks`, sample_mz, sample_rt, ref_mz, ref_rt, a, b, rt_diff_tol, ppm)
}

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

