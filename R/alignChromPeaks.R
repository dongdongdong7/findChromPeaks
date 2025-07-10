# Align chromatogram peaks
# Barry Song
# 250709

.find_insert_indices <- function(peaks_mz, peaks_rt, ref_mz, ref_rt, ppm = 10, rt_diff_tol = 10){
  tol_mz <- MsCoreUtils::ppm(peaks_mz, ppm = ppm)
  C_find_insert_indices(peaks_mz = peaks_mz, peaks_rt = peaks_rt, ref_mz = ref_mz, ref_rt = ref_rt, tol_mz = tol_mz, rt_tol = rt_diff_tol)
}

.match_peaks <- function(sample_mz, sample_rt, ref_mz, ref_rt, a, b, rt_diff_tol, ppm) {
  C_match_peaks(sample_mz = sample_mz, sample_rt = sample_rt, ref_mz = ref_mz, ref_rt = ref_rt,
                a = a, b = b, rt_diff_tol = rt_diff_tol, ppm = ppm)
}

#' @rdname alignChromPeaks
#' @title Align chromatographic peaks
#' @description
#' The algorithm of peak alignment here is based on the idea of Joint Aligner implemented in MZmine and MS-DIAL.
#' @param chromPeaksDT `data.table()`, with mz, mzmin, mzmax, rt, rtmin, rtmax, into, maxo, sn and sample
#' @param ppm `numeric()` defining the maximal tolerated m/z deviation between two peaks from one compound.
#' @param method merge methods are ppm or range
#' @export
#' @examples
#' chromPeaksDT <- readRDS(file = "./inst/demo/chromPeaksDT1.rds")
#' chromPeaksDT <- mergePeaks(chromPeaksDT = chromPeaksDT)
mergePeaks <- function(chromPeaksDT, ppm = 10, method = c("ppm", "range")){
  method <- match.arg(method)
  sample_id <- unique(chromPeaksDT$sample)
  pb <- progress::progress_bar$new(
    format = "[:bar] :elapsed | progress: :percent",
    total = length(sample_id),
    clear = FALSE,
    width = 60
  )
  chromPeaksDT_new <- data.table::rbindlist(lapply(sample_id, function(i) {
    pb$tick()
    chromPeaksDT_i <- chromPeaksDT[chromPeaksDT$sample == i, ]
    if(method == "ppm"){
      mz_tol_vec <- MsCoreUtils::ppm(chromPeaksDT_i$mz, ppm = ppm) * 1e3
      mz_ranges <- IRanges::IRanges(start = chromPeaksDT_i$mz * 1e3 - mz_tol_vec, end = chromPeaksDT_i$mz * 1e3 + mz_tol_vec)
      rt_ranges <- IRanges::IRanges(start = chromPeaksDT_i$rtmin * 1e3, end = chromPeaksDT_i$rtmax * 1e3)
    }else{
      mz_ranges <- IRanges::IRanges(start = chromPeaksDT_i$mzmin * 1e3, end = chromPeaksDT_i$mzmax * 1e3)
      rt_ranges <- IRanges::IRanges(start = chromPeaksDT_i$rtmin * 1e3, end = chromPeaksDT_i$rtmax * 1e3)
    }
    mz_overlaps <- IRanges::findOverlaps(mz_ranges, mz_ranges, type = "any", minoverlap = 2)
    mz_groups <- data.table::as.data.table(mz_overlaps)[, .(group = min(subjectHits)), by = queryHits]
    chromPeaksDT_i[, mz_group := mz_groups$group[.I]]

    rt_overlaps <- IRanges::findOverlaps(rt_ranges, rt_ranges, type = "any", minoverlap = 2)
    rt_groups <- data.table::as.data.table(rt_overlaps)[, .(group = min(subjectHits)), by = queryHits]
    chromPeaksDT_i[, rt_group := rt_groups$group[.I]]

    chromPeaksDT_i[, group_id := .GRP, by = .(mz_group, rt_group)]

    chromPeaksDT_i_new <- chromPeaksDT_i[, .(
      mz = mean(mz),
      mzmin = min(mzmin),
      mzmax = max(mzmax),
      rt = mean(rt),
      rtmin = min(rtmin),
      rtmax = max(rtmax),
      into = sum(into),
      intb = sum(intb),
      maxo = max(maxo),
      sn = min(sn),
      sample = sample[1]
    ), by = group_id]

    chromPeaksDT_i_new[, .(mz, mzmin, mzmax, rt, rtmin, rtmax, into, intb, maxo, sn, sample)]
  }))
  chromPeaksDT_new[, cpid := paste0("CP", formatC(seq_len(.N), flag = "0", width = ceiling(log10(.N + 1))))]
  data.table::setcolorder(chromPeaksDT_new, neworder = c("cpid", setdiff(colnames(chromPeaksDT_new), "cpid")))
  return(chromPeaksDT_new)
}

#' @rdname alignChromPeaks
#' @param `integer(1)` idx of reference sample
#' @param rt_diff_tol tolerance value for rt deviation
#' @export
#' @examples
#' chromPeaksDT: data.table with at least cpid, mz, rt, sample
#' refFeaturesDT <- makingRefFeatures(chromPeaksDT)
makingRefFeatures <- function(chromPeaksDT, reference = 1, ppm = 10, rt_diff_tol = 10){
  chromPeaksDT_reference <- chromPeaksDT[chromPeaksDT$sample == reference, ]

  featuresDT_reference <- chromPeaksDT_reference[, .(mz, mzmin, mzmax, rt, rtmin, rtmax)]

  sample_idx_noref <- setdiff(unique(chromPeaksDT$sample), reference)
  pb <- progress::progress_bar$new(
    format = "[:bar] :elapsed | progress: :percent",
    total = length(sample_idx_noref),
    clear = FALSE,
    width = 60
  )
  for(i in sample_idx_noref){
    pb$tick()
    chromPeaksDT_i <- chromPeaksDT[sample == i, ]
    insert_idx <- .find_insert_indices(peaks_mz = chromPeaksDT_i$mz, peaks_rt = chromPeaksDT_i$rt,
                                       ref_mz = featuresDT_reference$mz, ref_rt = featuresDT_reference$rt,
                                       ppm = ppm, rt_diff_tol = rt_diff_tol)
    # insert_idx <- which(sapply(1:nrow(chromPeaksDT_i), function(j) {
    #   tol_mz_j <- MsCoreUtils::ppm(chromPeaksDT_i[j, mz], ppm = ppm)
    #   all(abs(chromPeaksDT_i[j, mz] - featuresDT_reference$mz) > tol_mz_j |
    #         abs(chromPeaksDT_i[j, rt] - featuresDT_reference$rt) > rt_diff_tol)
    # }))
    featuresDT_reference <- rbind(featuresDT_reference, chromPeaksDT_i[insert_idx, .(mz, mzmin, mzmax, rt, rtmin, rtmax)])
  }
  featuresDT_reference
}

#' @rdname alignChromPeaks
#' @param refFeaturesDT `data.table()`, reference features data.table
#' @param a `numeric(1)`, factor for rt score
#' @param b `numeric(1)`, factor for m/z score
#' @export
#' @examples
#' featuresDT <- fittingFeatures(chromPeaksDT = chromPeaksDT, refFeaturesDT = refFeaturesDT, ppm = 10, rt_diff_tol = 10)
fittingFeatures <- function(chromPeaksDT, refFeaturesDT, ppm = 10, rt_diff_tol = 10, a = 0.5, b = 0.5){
  sample_idx <- unique(chromPeaksDT$sample)
  pb <- progress::progress_bar$new(
    format = "[:bar] :elapsed | progress: :percent",
    total = length(sample_idx),
    clear = FALSE,
    width = 60
  )
  matched_peaks_idx_list <- lapply(sample_idx, function(i) {
    pb$tick()
    peaks_idx_i <- which(chromPeaksDT$sample == i)
    ref_mz <- refFeaturesDT$mz
    ref_rt <- refFeaturesDT$rt
    sample_mz <- chromPeaksDT[peaks_idx_i, mz]
    sample_rt <- chromPeaksDT[peaks_idx_i, rt]
    matched_idx_i <- .match_peaks(sample_mz = sample_mz, sample_rt = sample_rt, ref_mz = ref_mz, ref_rt = ref_rt,
                                  a = a, b = b, rt_diff_tol = rt_diff_tol, ppm = ppm)
    peaks_idx_i[matched_idx_i]
  })
  refFeaturesDT[, (paste0("sample", seq_along(matched_peaks_idx_list))) := matched_peaks_idx_list]
  return(refFeaturesDT)
}

#' @rdname alignChromPeaks
#' @param minFraction `numeric(1)` defining the minimum fraction of samples in at least one sample group in which the peaks have to be present to be considered as a peak group (feature).
#' @param minPeaks `numeric(1)` with the minimum number of samples in at least one sample group in which the peaks have to be detected to be considered a peak group (feature).
#' @export
#' @examples
#' featuresDT <- filteringFeatures(featuresDT = featuresDT)
filteringFeatures <- function(featuresDT, minFraction = 0.5, minPeaks = 1){
  sample_names <- setdiff(colnames(featuresDT), c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax"))
  n_samples <- length(sample_names)
  rows_notNA_num <- rowSums(!is.na(featuresDT[, ..sample_names]))
  featuresDT <- featuresDT[rows_notNA_num > round(n_samples * minFraction) & rows_notNA_num > minPeaks, ]
  return(featuresDT)
}

#' @rdname alignChromPeaks
#' @param output values output type
#' @export
#' @examples
#' featuresValueDT <- getFeaturesValue(featuresDT = featuresDT, chromPeaksDT = chromPeaksDT)
getFeaturesValue <- function(featuresDT, chromPeaksDT, output = c("into", "intb", "maxo")){
  output <- match.arg(output)
  sample_names <- setdiff(colnames(featuresDT), c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax"))
  sample_idx <- unique(chromPeaksDT$sample)
  featuresValueDT <- featuresDT[, c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax")]
  for(i in seq_along(sample_names)){
    matched_idx <- featuresDT[[sample_names[i]]]
    featuresValueDT[, sample_names[i] := chromPeaksDT[[output]][matched_idx]]
  }
  return(featuresValueDT)
}

# .mz_paired_mapping(mzVec1 = chromPeaksDT$mz[chromPeaksDT$sample == 1], mzVec2 = chromPeaksDT$mz[chromPeaksDT$sample == 2], ppm = 20)
.mz_paired_mapping <- function(mzVec1, mzVec2, ppm = 5){
  map_l <- lapply(1:length(mzVec1), function(i) {
    x <- mzVec1[i]
    j <- which.min(abs(mzVec2 - x))
    y <- mzVec2[j]
    if(abs(x - y) < MsCoreUtils::ppm(max(x, y), ppm = ppm)){return(c(i, j))}
    else return(NULL)
  })
  map_l <- map_l[!sapply(map_l, is.null)]
  mz1_i <- sapply(map_l, function(x) {x[1]})
  mz2_j <- sapply(map_l, function(x) {x[2]})
  dt_tmp <- data.table::data.table(
    mz1 = mzVec1[mz1_i],
    mz2 = mzVec2[mz2_j],
    mz1_i = mz1_i,
    mz2_j = mz2_j
  )
  dp1 <- which(duplicated(dt_tmp$mz1_i, fromLast = TRUE) | duplicated(dt_tmp$mz1_i, fromLast = FALSE))
  dp2 <- which(duplicated(dt_tmp$mz2_j, fromLast = TRUE) | duplicated(dt_tmp$mz2_j, fromLast = FALSE))
  dp_index <- c(dp1, dp2)
  if(length(dp_index) > 0){
    dt_tmp <- dt_tmp[-dp_index, ]
  }
  dt_tmp$ratio_deltas <- abs(dt_tmp$mz1 - dt_tmp$mz2) / dt_tmp$mz1
  return(dt_tmp)
}



