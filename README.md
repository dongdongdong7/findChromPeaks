# findChromPeaks

In the preprocessing of metabolomics data, a basic requirement is to find metabolite peaks in chromatograms.

## Simulate chromatographic peaks

```R
rtime <- seq(0, 1000, 0.5)
peak1 <- gaussian_peak(rtime, peakRt = 240, peakWidth = 20, peakHeight = 10000)
peak2 <- tailing_peak(rtime, peakRt = 480, peakWidth = 25, peakHeight = 1300, tau = 10)
peak3 <- leading_peak(rtime, peakRt = 720, peakWidth = 18, peakHeight = 5000, tau = -8)
combined <- peak1 + peak2 + peak3
baseline <- 100 + 1 * rtime + 100 * sin(rtime/100)
final_chrom <- combined + baseline  + abs(rnorm(length(rtime), mean = 0, sd = 25))
plot(rtime, final_chrom, type = "l", xlab = "Retention Time (sec)", ylab = "Intensity")
```

<img src=".\assets\image-20250609084604328.png" alt="image-20250609084604328" style="zoom:80%;" />

## Find chromatographic peaks

```R
findChromPeaks_CWT(int = final_chrom, rt = rtime, peakwidth = c(5, 30), snthresh = 10, noise = 100)
#       rt rtmin rtmax sci scimin scimax     into baseline     maxo sn f scale scpos scmin scmax lmin lmax
# [1,] 240   222 258.5 481    445    518 141563.8      430 10416.98 24 1     5   481   476   486   23   98
#      r2 cs
# [1,]  1  1
findChromPeaks_CWT(int = final_chrom, rt = rtime, peakwidth = c(5, 30), snthresh = 3, noise = 100)
#         rt rtmin rtmax  sci scimin scimax      into baseline      maxo sn f scale scpos scmin scmax lmin
# [1,] 240.0   222 258.5  481    445    518 141563.84      430 10416.980 24 1     5   481   476   486   23
# [2,] 713.5   685 731.5 1428   1371   1464  97857.16      939  4058.996  4 2    11  1432  1421  1443    1
#      lmax   r2 cs
# [1,]   98 1.00  1
# [2,]   94 0.82  1
findChromPeaks_CWT(int = final_chrom, rt = rtime, peakwidth = c(5, 30), snthresh = 2, noise = 100)
#         rt rtmin rtmax  sci scimin scimax      into baseline      maxo sn  f scale scpos scmin scmax lmin
# [1,] 240.0   222 258.5  481    445    518 141563.84      428 10416.980 24 25     5   481   476   486   23
# [2,] 487.5   459 488.0  976    919    977  44596.60      546  1387.139  3 43    11   972   961   983    8
# [3,] 713.5   685 731.5 1428   1371   1464  97857.16      916  4058.996  4 58    11  1432  1421  1443    1
#      lmax   r2   cs
# [1,]   98 1.00 1.00
# [2,]  108 0.90 0.44
# [3,]   94 0.82 1.00
```



