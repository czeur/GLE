// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Gram-based coordinate descent LASSO solver
// Solves: min 0.5 * alpha^T * G * alpha + h^T * alpha + lambda * ||alpha||_1
// G is (d-1) x (d-1) symmetric positive definite (Gram matrix)
// h is (d-1) linear term
// alpha_init is warm-start initialization (zeros if empty)
static vec lasso_gram_cd(const mat& G, const vec& h, double lambda,
                         const vec& alpha_init,
                         int max_iter = 1000, double tol = 1e-7) {
  int p = G.n_cols;
  vec alpha = alpha_init.n_elem == (uword)p ? alpha_init : vec(p, fill::zeros);

  // Precompute G * alpha for residual tracking
  vec Galpha = G * alpha;

  for (int iter = 0; iter < max_iter; iter++) {
    double max_change = 0.0;

    for (int j = 0; j < p; j++) {
      double alpha_old = alpha(j);

      // Partial gradient excluding j-th diagonal: -(h_j + sum_{k!=j} G_{jk} alpha_k)
      double rho_j = -(h(j) + Galpha(j) - G(j, j) * alpha_old);

      // Soft thresholding
      double alpha_new;
      if (rho_j > lambda) {
        alpha_new = (rho_j - lambda) / G(j, j);
      } else if (rho_j < -lambda) {
        alpha_new = (rho_j + lambda) / G(j, j);
      } else {
        alpha_new = 0.0;
      }

      if (alpha_new != alpha_old) {
        double diff = alpha_new - alpha_old;
        // Update Galpha incrementally: O(p) instead of O(p^2)
        Galpha += G.col(j) * diff;
        alpha(j) = alpha_new;
        double change = std::abs(diff);
        if (change > max_change) max_change = change;
      }
    }

    if (max_change < tol) break;
  }

  return alpha;
}

// [[Rcpp::export]]
List glasso_c_cpp(arma::mat S, double lambda, double c = 0.0,
                  int iter_max = 200) {
  int d = S.n_rows;

  // Compute convergence tolerance from scale of S
  mat S_offdiag = S;
  S_offdiag.diag().zeros();
  double tol = 1e-2 * mean(mean(abs(S_offdiag)));

  // Initial guess
  mat W = S;

  // Sparsity pattern from column solutions (preserves exact LASSO zeros)
  umat graph(d, d, fill::zeros);

  // Theta assembled from column solutions (fallback when inv_sympd fails)
  mat Theta_direct(d, d, fill::zeros);

  // Warm-start storage: previous alpha for each column
  std::vector<vec> alpha_prev(d);
  for (int j = 0; j < d; j++) {
    alpha_prev[j] = vec(d - 1, fill::zeros);
  }

  // Iteration
  double delta = 1.0;
  int iteration = 0;

  while (delta > tol && iteration < iter_max) {
    iteration++;
    mat W_old = W;

    for (int j = 0; j < d; j++) {
      // Build index vector for "not j"
      uvec notj(d - 1);
      int pos = 0;
      for (int k = 0; k < d; k++) {
        if (k != j) notj(pos++) = k;
      }

      mat W11 = W(notj, notj);
      uvec jj = {(uword)j};
      vec w12 = vec(W(notj, jj));
      double w22 = W(j, j);
      double w22star = S(j, j);
      vec s12 = vec(S(notj, jj));

      mat Theta11Inv = W11 - w12 * w12.t() / w22;

      // Linear term for the Gram-based lasso:
      // min 0.5 * alpha^T * Theta11Inv * alpha + h^T * alpha + lambda * ||alpha||_1
      vec h = s12 + (c * w22star) * (Theta11Inv * ones<vec>(d - 1));

      // Solve lasso with warm start from previous iteration
      vec alpha = lasso_gram_cd(Theta11Inv, h, lambda, alpha_prev[j]);
      alpha_prev[j] = alpha;

      // theta12 = alpha / w22star + c
      vec theta12 = alpha / w22star + c;

      // Update W
      vec tmp_vec = Theta11Inv * theta12;
      W11 = Theta11Inv + w22star * tmp_vec * tmp_vec.t();
      w12 = -w22star * tmp_vec;
      w22 = w22star;

      W(notj, notj) = W11;
      W(notj, jj) = w12;
      W(jj, notj) = w12.t();
      W(j, j) = w22;
      W = (W + W.t()) / 2.0;

      // Track sparsity from exact LASSO zeros: alpha==0 means theta12==c
      for (int i = 0; i < d - 1; i++) {
        graph(notj(i), j) = (alpha(i) == 0.0) ? 1 : 0;
      }

      // Store column solution in Theta_direct (preserves exact zeros)
      for (int i = 0; i < d - 1; i++) {
        Theta_direct(notj(i), j) = theta12(i);
      }
      Theta_direct(j, j) = 1.0 / w22star + dot(theta12, tmp_vec);
    }

    delta = mean(mean(abs(W - W_old)));
  }

  mat Theta;
  if (!inv_sympd(Theta, W)) {
    // W not positive definite — fall back to column-assembled Theta
    Theta = (Theta_direct + Theta_direct.t()) / 2.0;
  }

  // Symmetrize graph: absent only if absent from both column solutions
  graph = graph % graph.t();
  graph.diag().zeros();

  return List::create(Named("Theta") = Theta,
                      Named("Sigma") = W,
                      Named("graph") = graph);
}

// [[Rcpp::export]]
List glasso_c_reest_cpp(arma::mat S, double lambda, double c = 0.0,
                        int iter_max = 200) {
  int d = S.n_rows;

  if (c == 0.0) {
    vec eigvals = eig_sym(S);
    c = 1.0 / (d * eigvals(1));  // smallest positive eigenvalue (eigvals[0] ≈ 0)
  }

  // Modify S
  S = S + (1.0 / c / (d * d)) * ones(d, d);

  List res = glasso_c_cpp(S, lambda, c, iter_max);
  mat Theta = as<mat>(res["Theta"]);
  umat select = as<umat>(res["graph"]);

  mat Theta_hat = Theta - c;

  // Zero out entries identified as absent by LASSO sparsity
  Theta_hat.elem(find(select)).zeros();

  return List::create(Named("Theta_hat") = Theta_hat,
                      Named("graph") = select);
}
