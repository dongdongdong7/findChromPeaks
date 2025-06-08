// C++ functions used in findChromPeaks.R
#include <Rcpp.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::export]]
double C_noiseEs(Rcpp::NumericVector intensity, double mag = 3.0) {
  // Step 1: Filter positive values and sort
  intensity = intensity[intensity > 0];
  std::sort(intensity.begin(), intensity.end());

  const int n = intensity.size();
  if (n == 1) return intensity[0];

  // Step 2: Iterate to find noise threshold
  for (int i = 1; i < n - 1; ++i) {
    const double current_int = intensity[i + 1];
    Rcpp::NumericVector sub_vec = intensity[Rcpp::seq(0, i)];
    const double noiEsti = Rcpp::mean(sub_vec) + mag * Rcpp::sd(sub_vec);
    if (current_int > noiEsti) return intensity[i];
  }

  // Step 3: Return last element if no early stop
  return intensity[n - 1];
}

// 辅助函数：在区间 [start, end] 内返回最大值索引
int C_windowMaxIdx(const Rcpp::NumericVector& x, int start, int end) {
  int m = start;
  for (int i = start + 1; i <= end; ++i) {
    if (x[m] < x[i]) m = i;
  }
  return m;
}

// 主函数：检测局部最大值
// [[Rcpp::export]]
Rcpp::LogicalVector C_localMaxima(const Rcpp::NumericVector& y, int halfWindowSize) {
  int n = y.size();
  int windowSize = 2 * halfWindowSize;

  // 扩展 y：前后各填充 halfWindowSize 个 0
  Rcpp::NumericVector y_extended(n + windowSize, 0.0);
  std::copy(y.begin(), y.end(), y_extended.begin() + halfWindowSize);

  // 在扩展后的数据上计算局部最大值
  Rcpp::LogicalVector output_extended(y_extended.size(), false);
  int m = C_windowMaxIdx(y_extended, 0, windowSize);
  if (m == halfWindowSize) output_extended[m] = true;

  for (int i = windowSize + 1, l = i - windowSize, mid = (l + i) / 2;
       i < y_extended.size();
       ++i, ++l, ++mid) {
    if (m < l) {
      m = C_windowMaxIdx(y_extended, l, i);
    } else if (y_extended[i] > y_extended[m]) {
      m = i;
    }
    if (m == mid) output_extended[m] = true;
  }

  // 截取有效部分（去掉扩展的边界）
  Rcpp::LogicalVector output(n, false);
  std::copy(
    output_extended.begin() + halfWindowSize,
    output_extended.end() - halfWindowSize,
    output.begin()
  );

  return output;
}

// [[Rcpp::export]]
Rcpp::IntegerVector C_DescendMinTol(Rcpp::NumericVector d, Rcpp::IntegerVector startpos, int maxDescOutlier) {
  int l = startpos[0] - 1;  // C++ 用 0-based 索引
  int r = startpos[1] - 1;
  int N = d.size();
  int outl = 0, opos = 0, vpos = 0;

  // 向左搜索
  while ((l > 0) && (d[l] > 0) && (outl <= maxDescOutlier)) {
    vpos = (outl > 0) ? opos : l;
    if (d[l - 1] > d[vpos]) {
      outl++;
    } else {
      outl = 0;
    }
    if (outl == 1) opos = l;
    l--;
  }
  if (outl > 0) l += outl;

  // 向右搜索
  outl = 0;
  while ((r < N - 1) && (d[r] > 0) && (outl <= maxDescOutlier)) {
    vpos = (outl > 0) ? opos : r;
    if (d[r + 1] > d[vpos]) {
      outl++;
    } else {
      outl = 0;
    }
    if (outl == 1) opos = r;
    r++;
  }
  if (outl > 0) r -= outl;

  return Rcpp::IntegerVector::create(l + 1, r + 1);  // 返回 1-based 索引
}

// [[Rcpp::export]]
Rcpp::IntegerVector C_RectUnique(Rcpp::NumericVector m, Rcpp::IntegerVector order,
                                 int nrow, int ncol,
                                 double xdiff, double ydiff) {

  int i, j, io, jo;
  int x1 = 0, x2 = nrow, y1 = nrow*2, y2 = nrow*3;
  Rcpp::IntegerVector keep(nrow, 1); // 初始化为全1

  for (i = 0; i < nrow; i++) {
    io = order[i] - 1; // R的索引从1开始，C++从0开始
    keep[io] = 1;

    for (j = 0; j < i; j++) {
      jo = order[j] - 1; // 同样调整索引

      if (keep[jo] &&
          !(m[x1+io] - m[x2+jo] > xdiff ||
          m[x1+jo] - m[x2+io] > xdiff ||
          m[y1+io] - m[y2+jo] > ydiff ||
          m[y1+jo] - m[y2+io] > ydiff)) {
        keep[io] = 0;
        break;
      }
    }
  }

  return keep;
}
