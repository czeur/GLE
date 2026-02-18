// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Trim matrix: W - min(c, c_max) * ones
static mat trim_matrix_cpp(const mat& W, double c) {
  int d = W.n_rows;
  mat Winv;
  double c_use;
  if (inv_sympd(Winv, W)) {
    double c_max = 1.0 / accu(Winv);
    c_use = std::min(c, c_max);
  } else {
    // W is singular — c_max is infinity, so c_use = c
    c_use = c;
  }
  return W - c_use * ones(d, d);
}

// [[Rcpp::export]]
List est_W_cpp(arma::mat data, double q_threshold) {
  int n = data.n_rows;
  int d = data.n_cols;

  // Rank transformation: mdata[i,j] = 1 / (1 - rank(data[,j]) / (n+1))
  mat mdata(n, d);
  const double inv_np1 = 1.0 / (n + 1.0);
  for (int j = 0; j < d; j++) {
    uvec order = sort_index(data.col(j));
    for (int k = 0; k < n; ++k) {
      const double rx = (k + 1) * inv_np1;
      mdata(order(k), j) = 1.0 / (1.0 - rx);
    }
  }

  mat cov_mat(d, d, fill::zeros);
  List cov_k_list(d);
  double trace_sum = 0.0;

  double threshold_val = 1.0 / (1.0 - q_threshold);

  for (int k = 0; k < d; k++) {
    // Find rows where mdata[,k] > threshold
    uvec keep = find(mdata.col(k) > threshold_val);
    int nk = keep.n_elem;

    if (nk < 2) {
      // Not enough data points, store zero matrix
      cov_k_list[k] = mat(d - 1, d - 1, fill::zeros);
      continue;
    }

    // Extract submatrix: all columns except k, only kept rows
    uvec notk(d - 1);
    int pos = 0;
    for (int j = 0; j < d; j++) {
      if (j != k) notk(pos++) = j;
    }

    mat data_trunc_notk = mdata(keep, notk);
    vec data_trunc_k = mdata(keep, uvec({(uword)k}));

    // w = log(data_trunc[, -k]) - log(data_trunc[, k]) * ones_row
    mat w = log(data_trunc_notk);
    w.each_col() -= log(data_trunc_k);

    // Compute covariance of w
    mat cov_k = cov(w);

    // Accumulate
    // cov_mat[-k, -k] += cov_k / d
    cov_mat(notk, notk) += cov_k / d;
    cov_k_list[k] = cov_k;
    trace_sum += accu(cov_k) / (d * d * d);
  }

  mat trimmed = trim_matrix_cpp(cov_mat, trace_sum);

  // Convert to correlation matrix: cor = D^{-1/2} * cov * D^{-1/2}
  vec sd_inv = 1.0 / sqrt(trimmed.diag());
  mat result = trimmed % (sd_inv * sd_inv.t());

  return List::create(Named("cov") = result,
                      Named("subcovlist") = cov_k_list);
}
