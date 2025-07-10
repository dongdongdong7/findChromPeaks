// C++ functions used in alignChromPeaks.R

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

// Helper function to calculate ppm tolerance
inline double C_calculate_ppm(double mz, double ppm) {
  return mz * ppm * 1e-6;
}

// [[Rcpp::export]]
Rcpp::IntegerVector C_find_insert_indices(Rcpp::NumericVector peaks_mz, Rcpp::NumericVector peaks_rt,
                                          Rcpp::NumericVector ref_mz, Rcpp::NumericVector ref_rt,
                                          Rcpp::NumericVector tol_mz, double rt_tol) {
  int n = peaks_mz.size();
  int m = ref_mz.size();
  Rcpp::IntegerVector res;

  for(int i = 0; i < n; i++) {
    bool insert = true;
    for(int j = 0; j < m; j++) {
      if(std::abs(peaks_mz[i] - ref_mz[j]) <= tol_mz[i] &&
         std::abs(peaks_rt[i] - ref_rt[j]) <= rt_tol) {
        insert = false;
        break;
      }
    }
    if(insert) res.push_back(i + 1); // R uses 1-based indexing
  }

  return res;
}

// [[Rcpp::export]]
Rcpp::IntegerVector C_match_peaks(
    Rcpp::NumericVector sample_mz,
    Rcpp::NumericVector sample_rt,
    Rcpp::NumericVector ref_mz,
    Rcpp::NumericVector ref_rt,
    double a,
    double b,
    double rt_diff_tol,
    double ppm) {

  int n_ref = ref_mz.size();
  int n_peaks = sample_mz.size();

  // Pre-allocate integer vector for results (1-based indices)
  Rcpp::IntegerVector result(n_ref, NA_INTEGER);

  // Track which sample peaks have been matched
  std::vector<bool> matched_peaks(n_peaks, false);

  for (int j = 0; j < n_ref; ++j) {
    double current_ref_mz = ref_mz[j];
    double current_ref_rt = ref_rt[j];
    double tol_mz = C_calculate_ppm(current_ref_mz, ppm);

    double max_score = -1.0;
    int best_idx = -1;  // Will be converted to 1-based

    for (int i = 0; i < n_peaks; ++i) {
      // Skip already matched peaks
      if (matched_peaks[i]) {
        continue;
      }

      double rt_diff = sample_rt[i] - current_ref_rt;
      double mz_diff = sample_mz[i] - current_ref_mz;

      if(fabs(mz_diff) > tol_mz || fabs(rt_diff) > rt_diff_tol){
        continue;
      }

      double score = a * std::exp(-0.5 * std::pow(rt_diff / rt_diff_tol, 2)) +
        b * std::exp(-0.5 * std::pow(mz_diff / tol_mz, 2));

      if (score > max_score) {
        max_score = score;
        best_idx = i;
      }
    }

    if (best_idx != -1) {
      matched_peaks[best_idx] = true;
      result[j] = best_idx + 1; // Convert to 1-based index
    }
  }

  return result;
}
