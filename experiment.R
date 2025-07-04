# chromData1
chromData <- readRDS(file = "./inst/demo/chromData1.rds")
plot(x = chromData$rtime, y = chromData$intensity, type = "o")
findChromPeaks_CWT(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(5, 20), snthresh = 3)
xcms::peaksWithCentWave(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(5, 20), snthresh = 3)

# chromData2
chromData <- readRDS(file = "./inst/demo/chromData2.rds")
plot(x = chromData$rtime, chromData$intensity, type = "o")
findChromPeaks_CWT(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(5, 20), snthresh = 3)
findChromPeaks_CWT(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(5, 20), snthresh = 3, r2thresh = 0.8)

# chromData3
chromData <- readRDS(file = "./inst/demo/chromData3.rds")
plot(x = chromData$rtime, chromData$intensity, type = "o")
findChromPeaks_CWT(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(5, 20), snthresh = 3, extendLengthMSW = TRUE)

# chromData4
chromData <- readRDS(file = "./inst/demo/chromData4.rds")
plot(x = chromData$rtime, chromData$intensity, type = "o")
findChromPeaks_CWT(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(5, 20), snthresh = 3)

# chromData5
# Change minpack.lm::nlsLM back to nls
chromData <- readRDS(file = "./inst/demo/chromData5.rds")
plot(x = chromData$rtime, chromData$intensity, type = "o")
findChromPeaks_CWT(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(5, 20), snthresh = 3)

# chromData6
# Add peak points number judgment
chromData <- readRDS(file = "./inst/demo/chromData6.rds")
plot(x = chromData$rtime, chromData$intensity, type = "o")
findChromPeaks_CWT(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(5, 20), snthresh = 3)

# chromData7
chromData <- readRDS(file = "./inst/demo/chromData7.rds")
plot(x = chromData$rtime, chromData$intensity, type = "o")
findChromPeaks_CWT(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(1, 10), snthresh = 3)


system.time({findChromPeaks_CWT(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(5, 20), snthresh = 3)})
system.time({xcms::peaksWithCentWave(int = chromData$intensity, rt = chromData$rtime, peakwidth = c(5, 20), snthresh = 3, fitgauss = TRUE)})

rtime <- seq(0, 1000, 0.5)
peak1 <- gaussian_peak(rtime, peakRt = 240, peakWidth = 20, peakHeight = 10000)
peak2 <- tailing_peak(rtime, peakRt = 480, peakWidth = 25, peakHeight = 1300, tau = 10)
peak3 <- leading_peak(rtime, peakRt = 720, peakWidth = 18, peakHeight = 5000, tau = -8)
combined <- peak1 + peak2 + peak3
baseline <- 100 + 1 * rtime + 100 * sin(rtime/100)
final_chrom <- combined + baseline  + abs(rnorm(length(rtime), mean = 0, sd = 25))
plot(rtime, final_chrom, type = "l", xlab = "Retention Time (sec)", ylab = "Intensity")
system.time({
  findChromPeaks_CWT(int = final_chrom, rt = rtime, peakwidth = c(5, 30), snthresh = 2, noise = 100)
})
system.time({
  xcms::peaksWithCentWave(int = final_chrom, rt = rtime, peakwidth = c(5, 30), snthresh = 2, noise = 100, fitgauss = TRUE)
})
