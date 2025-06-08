rtime <- seq(0, 1000, 0.5)

peak1 <- .gaussian_peak(rtime, peakRt = 240, peakWidth = 20, peakHeight = 10000)
peak2 <- .emg_peak(rtime, peakRt = 480, peakWidth = 25, peakHeight = 1000, tau = 10)
peak3 <- .leading_peak(rtime, peakRt = 720, peakWidth = 18, peakHeight = 5000, tau = -8)

combined <- peak1 + peak2 + peak3 + rnorm(length(rtime), sd = 0.8)

baseline <- 3 + 0.003 * rtime + 2 * sin(rtime/150)
final_chrom <- combined + baseline
plot(rtime, final_chrom, type = "l")

findChromPeaks_CWT(int = final_chrom, rt = rtime, peakwidth = c(5, 30), snthresh = 3, noise = 100)
xcms::peaksWithCentWave(int = final_chrom, rt = rtime, snthresh = 1)
