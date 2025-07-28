# Find chromatogram peaks
# These are modified from xcms peaksWithCentWave function
# Barry Song
# 250507

# Noise estimation of chromatogram such as EIC or mass track
# int: intensisy of chromatogram,  `numeric()`
# mag: `integer()`
.noiseEs <- function(int, mag = 3){
  C_noiseEs(intensity = int, mag = mag)
}

# Find local maxima
# y: numeric
# halfWindowSize: half window size;
# returns a logical vector of the same length as y
.localMaxima <- function(y, halfWindowSize){
  C_localMaxima(y = y, halfWindowSize = halfWindowSize)
}

# Find ROI
# int: intensity of mass_track or chromatogram
# rt: retention time of mass_track or chromatogram
# peakwidth: minimum and maximum values of peak width
# noise: noise threshold
# prefilter: c(k, I), only regions of interest with at least k centroids with signal >= I are returned.
# return a matrix with scmin scmax and sccent
# .getRtROI(int = c(100, 200, 300, 1000, 500, 200, 50), rt = c(100, 102, 105, 108, 110, 115, 120))
.getRtROI <- function(int, rt, peakwidth = c(20, 50), noise = 0,
                      prefilter = c(3, 100)) {
  peakwidth <- range(peakwidth)
  if (length(prefilter) != 2)
    stop("'prefilter' has to be a 'numeric' of length 2")
  int_len <- length(int)
  if (int_len != length(rt))
    stop("lengths of 'int' and 'rt' have to match")
  ## rt to halfWindowSize:
  rt_step <- mean(diff(rt), na.rm = TRUE)
  up_bound <- ceiling(peakwidth[2] / rt_step)
  pk_idx <- which(.localMaxima(int, floor(peakwidth[1] / rt_step)))
  ## First filter: int > noise
  pk_idx <- pk_idx[int[pk_idx] >= noise]
  if (!length(pk_idx))
    return(matrix(ncol = 2, nrow = 0))
  ## Define the ROIs
  scmin <- sapply(pk_idx - up_bound, max, y = 1)
  scmax <- sapply(pk_idx + up_bound, min, y = int_len)
  ## Second filter: at least k values larger I
  roi_idxs <- mapply(scmin, scmax, FUN = seq, SIMPLIFY = FALSE)
  ok <- vapply(roi_idxs,
               FUN = function(x, k, I) {
                 sum(int[x] >= I) >= k
               },
               FUN.VALUE = logical(1),
               k = prefilter[1], I = prefilter[2],
               USE.NAMES = FALSE)
  if (any(ok))
    cbind(scmin = scmin[ok], scmax = scmax[ok], sccent = pk_idx[ok])
  else
    matrix(ncol = 3, nrow = 0)
}

# Finding peak boundaries from the left and right
# d: intensity values
# startpos: `integer(2)`, left and right start position for searching
# maxDescOutlier: maximum number of outliers allowed
.descendMinTol <- function (d, startpos, maxDescOutlier){
  C_DescendMinTol(d = d, startpos = startpos, maxDescOutlier = maxDescOutlier)
}

# Data points greater than thresh are retained for narrowing the range of peak
# lm: `integer(2)`, start position and end position of one peak on d
# d: intensity values
# thresh: noise
.narrow_rt_boundaries <- function (lm, d, thresh = 1){
  lm_seq <- lm[1]:lm[2]
  above_thresh <- d[lm_seq] >= thresh
  if (any(above_thresh)) {
    above_thresh <- above_thresh | c(above_thresh[-1], FALSE) |
      c(FALSE, above_thresh[-length(above_thresh)])
    lm <- range(lm_seq[above_thresh], na.rm = TRUE)
  }
  lm
}

# Delete rows in matrix where mz and rt overlap
# m: a matrix of four column mzmin, mzmax, rtmin and rtmax
# order: priority of rows in the matrix
.rectUnique <- function(m, order = seq(length = nrow(m)), xdiff = 0, ydiff = 0){
  nr <- nrow(m)
  nc <- ncol(m)
  if(is.null(nr) || nr == 0) return(m)
  if(!is.double(m)) m <- as.double(m)
  as.logical(C_RectUnique(m = m, order = as.integer(order), nrow = nr, ncol = nc, xdiff = as.double(xdiff), ydiff = as.double(ydiff)))
}

#' @title findChromPeaks_CWT
#' @description
#' Find chromatographic peaks using CentWave.
#'
#' @param int `numeric`, intensity of chromatogram.
#' @param rt `numeric`, retention time of chromatogram.
#' @param peakwidth `numeric(2)` with the lower and upper boun of the expected peak width.
#' @param snthresh `numeric(1)` defining the signal to noise ratio cutoff.
#' @param minPs `integer(1)`, the ROI region requires a minimum of minPs of signals greater than the noise.
#' @param noise `numeric(1)`, noise of chromatogram.
#' @param estimateNoise `logical(1)`, whether to estimate noise
#' @param extendLengthMSW `logical(1)` please see: [xcms::centWave]
#' @param r2thresh `numeric(1)` threshold of peak shape.
#' @param csthresh `numeric(1)` threshold of cs.
#'
#' @returns `matrix`
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
#' findChromPeaks_CWT(int = final_chrom, rt = rtime, peakwidth = c(5, 30), snthresh = 2, noise = 100)
#' xcms::peaksWithCentWave(int = final_chrom, rt = rtime, peakwidth = c(5, 30), snthresh = 2, noise = 100)
findChromPeaks_CWT <- function(int, rt,
                               peakwidth = c(5, 20),
                               snthresh = 10,
                               minPs = 3,
                               noise = 100,
                               estimateNoise = TRUE,
                               extendLengthMSW = TRUE,
                               r2thresh = 0.6,
                               csthresh = 0.2){
  if(length(peakwidth) != 2) stop("'peakwidth' has to be a numeric of length 2")

  int[is.na(int)] <- 0 # avoid NAs
  # scanIndex <- 1:length(int)

  dt <- round(mean(diff(rt)), 1) # dt: mean of diff rt
  rt_uniform <- seq(min(rt), max(rt), by = dt) # new rt
  interp_linear <- approx(rt, int, xout = rt_uniform) # Get the int interpolation on rt_uniform
  int_o <- int;rt_o <- rt # store original int and rt
  int <- interp_linear$y;rt <- interp_linear$x # new int and rt

  # estimate noise
  if(estimateNoise){
    mag <- min(snthresh, 3)
    noise_es <- .noiseEs(int = int_o, mag = mag)
    if(length(which(int_o > noise_es)) > minPs) noise <- noise_es # avoid noise misestimation
  }
  # Get ROIs
  rois <- .getRtROI(int = int, rt = rt, peakwidth = peakwidth, noise = noise, prefilter = c(minPs, noise))

  # column names
  basenames <- c("mz", "mzmin", "mzmax",
                 "rt", "rtmin", "rtmax",
                 "sci", "scimin", "scimax",
                 "into", "intb", "baseline", "maxo", "sn", "ps")
  # f: index of rois; scale: best scale for a peak;
  # scpos: centre of roi; scmin: beginning of roi; scmax: end of roi
  # lmin: initial start of peak on roi
  # lmax: initial end of peak on roi
  # r2: peak shape r2 fitted using a gaussian function
  # cs: cs value same with asari
  verbosenames <- c("f", "scale", "scpos", "scmin", "scmax", "lmin", "lmax", "r2", "cs")
  peaks_names <- c(basenames, verbosenames)
  peaks_ncols <- length(peaks_names)
  peakinfo_names <- c("scale", "scaleNr", "scpos", "scmin", "scmax")

  # nopeaks <- matrix(nrow = 0, ncol = peaks_ncols,
  #                   dimnames = list(character(), peaks_names))
  # nopeaks <- nopeaks[, -(1:3), drop = FALSE]

  ## Peak width: seconds to scales
  # scale is the number of data points in the observation peaks, 2 is an empirical constant
  scalerange <- round((peakwidth / dt) / 2)
  # scalerange <- round((peakwidth / mean(diff(rt))) / 2)
  scalerange[1] <- max(3, scalerange[1]) # scales must be greater than 3, because a peak cannot be smaller than 3 data points
  scales <- seq(from = scalerange[1], to = scalerange[2], by = 2)
  minPeakWidth <-  scales[1]
  noiserange <- ceiling(minPeakWidth * 3)
  maxDescOutlier <- floor(minPeakWidth / 2)
  scanrange <- c(1, length(rt))
  sci <- scanrange[1]:scanrange[length(scanrange)]
  Nscantime <- length(int)
  ## mzdiff <- -0.001
  mzdiff <- 0

  peaklist <- NULL

  for(i in seq_len(nrow(rois))){
    scmin <- rois[i, "scmin"]
    scmax <- rois[i, "scmax"]

    N <- scmax - scmin + 1 # length of roi
    peaks <- matrix(ncol = peaks_ncols, nrow = 0,
                    dimnames = list(character(), peaks_names))
    peakinfo <- matrix(ncol = 5, nrow = 0,
                       dimnames = list(character(), peakinfo_names))
    ## Could also return the "correct one..."
    sccenter <- scmin + floor(N/2) - 1 # ROI centre
    ## sccenter <- rois[i, "sccent"]
    scrange <- c(scmin, scmax) # # start and end of roi

    ## the ROI
    otd <- scmin:scmax # ROI index
    od <- int[otd] # ROI intensity

    wCoefs <- MSW.cwt(od, scales = scales, wavelet = 'mexh',
                      extendLengthMSW = extendLengthMSW)
    if(is.null(dim(wCoefs))) next
    if (otd[length(otd)] == Nscantime) ## workaround, localMax fails otherwise
      wCoefs[nrow(wCoefs),] <- wCoefs[nrow(wCoefs) - 1, ] * 0.99
    localMax <- MSW.getLocalMaximumCWT(wCoefs)
    rL <- MSW.getRidge(localMax)
    wpeaks <- sapply(rL, function(x) {
      w <- min(1:length(x), ncol(wCoefs))
      any(wCoefs[x,w] >= 1)
    })
    if(any(wpeaks)){
      wpeaksidx <- which(wpeaks)
      for(p in 1:length(wpeaksidx)){
        opp <- rL[[wpeaksidx[p]]]
        pp <- unique(opp)   ## NOTE: this is not ordered!
        if(length(pp) >= 1){
          dv <- otd[pp]
          if(any(dv)){
            if(any(od[pp] > noise)){
              inti <- numeric(length(opp))
              irange <- rep(ceiling(scales[1]/2), length(opp))
              for (k in 1:length(opp)) {
                kpos <- opp[k]
                r1_ <- ifelse(kpos - irange[k] > 1,
                              kpos - irange[k], 1)
                r2_ <- ifelse(kpos + irange[k] < length(od),
                              kpos + irange[k], length(od))
                inti[k] <- sum(od[r1_:r2_])
              }
              maxpi <- which.max(inti)
              if (length(maxpi) > 1) {
                m <- wCoefs[opp[maxpi], maxpi]
                bestcol <- which(m == max(m),
                                 arr.ind = TRUE)[2]
                best.scale.nr <- maxpi[bestcol]
              } else {
                best.scale.nr <- maxpi
              }

              best.scale <- scales[best.scale.nr]
              best.scale.pos <- opp[best.scale.nr]

              pprange <- min(pp):max(pp)
              lwpos <- max(1, best.scale.pos - best.scale)
              rwpos <- min(best.scale.pos + best.scale, length(otd))
              p1 <- match(otd[lwpos], otd)[1]
              p2 <- match(otd[rwpos], otd)
              p2 <- p2[length(p2)]
              if (is.na(p1)) p1 <- 1
              if (is.na(p2)) p2 <- N
              maxint <- max(od[p1:p2])

              peaks <- rbind(
                peaks,
                c(1, 1, 1,    # mz, mzmin, mzmax,
                  NA, NA, NA, # rt, rtmin, rtmax,
                  NA, NA, NA, # sci, scimin, scimax,
                  NA,         # intensity (sum)
                  NA,         # intensity (-bl)
                  NA,         # baseline
                  maxint,     # max intensity
                  NA, # S/N Ratio
                  NA, # ps
                  i,        # ROI Position
                  best.scale, # Scale
                  otd[best.scale.pos],
                  otd[lwpos],
                  otd[rwpos], # Peak positions guessed from the wavelet's (scan nr)
                  NA, NA, # Peak limits (scan nr)
                  NA, # gauss peak shape r2
                  NA)) # cSelectivity
              peakinfo <- rbind(
                peakinfo,
                c(best.scale, best.scale.nr,
                  best.scale.pos, lwpos, rwpos))
              ## Peak positions guessed from the wavelet's
            }
          }
        }
      }
    }

    # postprocessing
    for(p in seq_len(nrow(peaks))){
      ## find minima, assign rt and intensity values
      lm <- .descendMinTol(od, startpos = c(peakinfo[p, "scmin"], peakinfo[p, "scmax"]), maxDescOutlier)
      lm <- .narrow_rt_boundaries(lm, od)
      lm_range <- lm[1]:lm[2]
      pd <- od[lm_range]

      peakrange <- otd[lm]
      peaks[p, "rtmin"] <- rt[peakrange[1]]
      peaks[p, "rtmax"] <- rt[peakrange[2]]
      peaks[p, "scimin"] <- sci[peakrange[1]]
      peaks[p, "scimax"] <- sci[peakrange[2]]
      peaks[p, "maxo"] <- max(pd)
      pwid <- (rt[peakrange[2]] - rt[peakrange[1]]) /
        (peakrange[2] - peakrange[1])
      if(is.na(pwid)) pwid <- 1
      peaks[p, "into"] <- pwid * sum(pd)
      peaks[p, "lmin"] <- lm[1]
      peaks[p, "lmax"] <- lm[2]
      peaks[p, "rt"] <- rt[peaks[p, "scpos"]]
      peaks[p, "sci"] <- sci[peaks[p, "scpos"]]
    }

    if(nrow(peaks) == 0) next

    if(nrow(peaks) > 1){
      uorder <- order(peaks[, "scale"], decreasing = FALSE) # retain peaks with smaller scales
      uindex <- .rectUnique(peaks[,c("mzmin", "mzmax", "scmin", "scmax")], order = uorder, xdiff = mzdiff, ydiff = -0.00001)
      peaks <- peaks[uindex, ]
      peakinfo <- peakinfo[uindex, ]
    }

    if (!is.null(peaks)) peaklist[[length(peaklist) + 1]] <- peaks
  }

  if(length(peaklist) == 0){
    warning("No peaks found!")
    nopeaks <- matrix(nrow = 0, ncol = peaks_ncols,
                      dimnames = list(character(), peaks_names))
    return(nopeaks[, -(1:3), drop = FALSE])
  }

  p <- do.call(rbind, peaklist)

  uorder <- order(p[, "into"], decreasing = TRUE)
  pm <- p[, c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE]
  uindex <- .rectUnique(pm, uorder, mzdiff, ydiff = -0.00001) # allow adjacent peaks
  p <- p[uindex, -(1:3), drop = FALSE]

  # calculate sn
  peakrange <- unique(purrr::list_c(lapply(1:nrow(p), function(i) {
    p[i, "scimin"]:p[i, "scimax"]
  })))
  for(i in 1:nrow(p)){
    tdrange <- c(p[i, "scimin"], p[i, "scimax"])
    td <- tdrange[1]:tdrange[2]
    ntd_left <- 1:p[i, "scimin"]
    ntd_left <- ntd_left[!ntd_left %in% peakrange]
    ntd_left <- ntd_left[max(1, length(ntd_left) - noiserange):length(ntd_left)]
    ntd_right <- p[i, "scimax"]:Nscantime
    ntd_right <- ntd_right[!ntd_right %in% peakrange]
    ntd_right <- ntd_right[1:min(1+noiserange, length(ntd_right))]
    ntd1 <- c(ntd_left, ntd_right)
    nd <- int[ntd1]
    nd <- nd[!is.na(nd)]
    if(length(nd) == 0) baseline <- 1
    else baseline <- mean(nd, na.rm = TRUE)
    if(baseline < 1) baseline <- 1
    p[i, "baseline"] <- round(baseline)
    pwid <- (rt[tdrange[2]] - rt[tdrange[1]]) / (tdrange[2] - tdrange[1])
    p[i, "intb"] <- pwid * sum(int[td] - p[i, "baseline"])
    p[i, "sn"] <- round(p[i, "maxo"] / baseline)
  }
  p <- p[which(p[, "sn"] > snthresh), , drop = FALSE]

  # Conversion of interpolated intensity to orignal intensity
  if(nrow(p) != 0){
    int_ap <- int;rt_ap <- rt
    int <- int_o;rt <- rt_o
    # peaks
    for(i in 1:nrow(p)){
      p[i, "sci"] <- which.min(abs(p[i, "rt"] - rt))
      p[i, "scimin"] <- which.min(abs(p[i, "rtmin"] - rt))
      p[i, "scimax"] <- which.min(abs(p[i, "rtmax"] - rt))
      p[i, "rt"] <- rt[p[i, "sci"]]
      p[i, "rtmin"] <- rt[p[i, "scimin"]]
      p[i, "rtmax"] <- rt[p[i, "scimax"]]
      # rt_scpos <- rt_ap[p[i, "scpos"]]
      # rt_scmin <- rt_ap[p[i, "scmin"]]
      # rt_scmax <- rt_ap[p[i, "scmax"]]
      # p[i, "scpos"] <- which.min(abs(rt_scpos - rt))
      # p[i, "scmin"] <- which.min(abs(rt_scmin - rt))
      # p[i, "scmax"] <- which.min(abs(rt_scmax - rt))
    }
    # rois
    for(i in 1:nrow(rois)){
      rt_scmin <- rt_ap[rois[i, "scmin"]]
      rt_scmax <- rt_ap[rois[i, "scmax"]]
      rt_sccent <- rt_ap[rois[i, "sccent"]]
      rois[i, "scmin"] <- which.min(abs(rt - rt_scmin))
      rois[i, "scmax"] <- which.min(abs(rt - rt_scmax))
      rois[i, "sccent"] <- which.min(abs(rt - rt_sccent))
    }
  }

  # fitgauss
  # otd_lm <- otd[lm_range]
  if(nrow(p) != 0){
    for(i in 1:nrow(p)){
      ptd <- p[i, "scimin"]:p[i, "scimax"]
      pd <- int[ptd]
      md <- max(pd)
      d1 <- pd / md
      lri <- which(d1 > 0.5)
      li <- lri[1]
      ri <- lri[length(lri)]
      fwhm <- ri - li
      c <- fwhm / 2.355 # fwhm = 2*sqrt(2*ln(2)) = 2.355 * c
      x <- 1:length(d1)
      # Gaussian fitting
      fit <- try(nls(d1 ~ a*exp(-((x-b)^2)/(2*c^2)),
                     start = list(a = max(d1),
                                  b = which.max(d1),
                                  c = c),
                     control = nls.control(maxiter = 50, tol = 1e-5),
                     algorithm = "port"), silent = TRUE)
      # fit <- try(minpack.lm::nlsLM(d1 ~ a*exp(-((x-b)^2)/(2*c^2)),
      #                              start = c(a = max(d1),
      #                                        b = which.max(d1),
      #                                        c = c)), silent = TRUE)
      if (class(fit) != "try-error"){
        d1_fit <- predict(fit, x = x)
        peak_pos <- which.max(d1_fit)
        if(peak_pos == 1 | peak_pos == length(d1_fit)) next
        deriv1 <- diff(d1_fit) / diff(x)
        x_deriv1 <- 1:length(deriv1)
        inflection_left_pos <- which.max(deriv1)
        inflection_right_pos <- which.min(deriv1)
        # The first order derivatives close to 0 are the boundaries of the peaks
        zero_pos <- which(abs(deriv1)< 0.001)
        left_pos <- zero_pos[which(zero_pos < inflection_left_pos)]
        if(length(left_pos) == 0) left_pos <- 1
        else left_pos <- max(1, left_pos[length(left_pos)])
        right_pos <- zero_pos[which(zero_pos > inflection_right_pos)]
        if(length(right_pos) == 0) right_pos <- length(d1_fit)
        else right_pos <- min(right_pos[1], length(x_deriv1))
        p[i, "scimin"] <- ptd[left_pos];p[i, "scimax"] <- ptd[right_pos]
        p[i, "rtmin"] <- rt[ptd[left_pos]];p[i, "rtmax"] <- rt[ptd[right_pos]]
        p[i, "sci"] <- ptd[which.max(d1_fit)];p[i, "rt"] <- rt[p[i, "sci"]]
        peakrange <- c(p[i, "scimin"], p[i, "scimax"])
        pwid <- (rt[peakrange[2]] - rt[peakrange[1]]) / (peakrange[2] - peakrange[1])
        if(is.na(pwid)) pwid <- 1
        p[i, "into"] <- pwid * sum(pd)
        p[i, "intb"] <- pwid * sum(pd - p[i, "baseline"])
        r2_fit <- round(1 - sum((d1 - d1_fit)^2) / sum((d1 - mean(d1_fit))^2), 2)
      }
      else{
        r2_fit <- 0
      }
      p[i, "r2"] <- r2_fit
    }
    p <- p[which(p[, "r2"] > r2thresh), , drop = FALSE]
  }

  # Peak points
  if(nrow(p) != 0){
    for(i in 1:nrow(p)){
      ptd <- p[i, "scimin"]:p[i, "scimax"]
      pd <- int[ptd]
      p[i, "ps"] <- length(which(pd > noise))
    }
  }
  p <- p[which(p[, "ps"] > minPs), , drop = FALSE]

  # calculate cSelectivity
  if(nrow(p) != 0){
    peakrange <- unique(purrr::list_c(lapply(1:nrow(p), function(i) {
      p[i, "scimin"]:p[i, "scimax"]
    })))
    for(i in 1:nrow(p)){
      # roi
      roi_ <- rois[p[i, "f"], ]
      roi_range <- roi_["scmin"]:roi_["scmax"]
      peak_dps <- int[peakrange[peakrange %in% roi_range]]
      roi_dps <- int[roi_range]
      pd <- int[p[i, "scimin"]:p[i, "scimax"]]
      md <- max(pd)
      half_md <- md / 2
      peak_level_size <- length(which(peak_dps > half_md))
      bg_level_size <- length(which(roi_dps > half_md))
      if(bg_level_size >= peak_level_size & bg_level_size > 0){
        cs_ <- round(peak_level_size / bg_level_size, 2)
      }else{
        cs_ <- 0
      }
      p[i, "cs"] <- cs_
    }
    p <- p[which(p[, "cs"] > csthresh), , drop = FALSE]
  }

  return(p)

}
