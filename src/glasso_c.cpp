// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Coordinate descent LASSO solver
// Solves: min 0.5 * ||A * beta - b||^2 + lambda * ||beta||_1
// A is (d-1) x (d-1), b is (d-1), lambda is scalar
static vec lasso_cd(const mat& A, const vec& b, double lambda,
                    int max_iter = 1000, double tol = 1e-7) {
  int p = A.n_cols;
  int n = A.n_rows;
  vec beta(p, fill::zeros);
  vec residual = -b;  // residual = A*beta - b, starts as -b since beta=0

  // Precompute A^T A columns norms and A^T b
  // For coordinate descent: update_j = (A_j^T (b - A_{-j} beta_{-j})) / (A_j^T A_j)
  // We use the residual form for efficiency
  vec AtA_diag(p);
  for (int j = 0; j < p; j++) {
    AtA_diag(j) = dot(A.col(j), A.col(j));
  }

  for (int iter = 0; iter < max_iter; iter++) {
    double max_change = 0.0;

    for (int j = 0; j < p; j++) {
      double beta_old = beta(j);

      // partial residual: rho_j = A_j^T (b - A_{-j} beta_{-j})
      //                         = A_j^T b - A_j^T A_{-j} beta_{-j}
      //                         = A_j^T (b - A*beta + A_j * beta_j)
      //                         = -A_j^T residual + AtA_diag(j) * beta_j
      double rho_j = -dot(A.col(j), residual) + AtA_diag(j) * beta_old;

      // Soft thresholding
      double beta_new;
      if (rho_j > lambda) {
        beta_new = (rho_j - lambda) / AtA_diag(j);
      } else if (rho_j < -lambda) {
        beta_new = (rho_j + lambda) / AtA_diag(j);
      } else {
        beta_new = 0.0;
      }

      if (beta_new != beta_old) {
        // Update residual: residual += A_j * (beta_new - beta_old)
        residual += A.col(j) * (beta_new - beta_old);
        double change = std::abs(beta_new - beta_old);
        if (change > max_change) max_change = change;
        beta(j) = beta_new;
      }
    }

    if (max_change < tol) break;
  }

  return beta;
}

// [[Rcpp::export]]
List glasso_c_cpp(arma::mat S, double lambda, double c = 0.0,
                  int iter_max = 1000) {
  int d = S.n_rows;

  // Compute convergence tolerance from scale of S
  mat S_offdiag = S;
  S_offdiag.diag().zeros();
  double tol = 1e-2 * mean(mean(abs(S_offdiag)));

  // Initial guess
  mat W = S;

  // Sparsity pattern from column solutions (preserves exact LASSO zeros)
  umat graph(d, d, fill::zeros);

  // Iteration
  double delta = 1.0;
  int iteration = 1;

  while (delta > tol && iteration <= iter_max) {
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

      // Eigen decomposition
      vec eigval;
      mat eigvec;
      eig_sym(eigval, eigvec, Theta11Inv);

      // Compute A = Q * diag(sqrt(eigenvalues)) * Q^T
      vec sqrt_eigval = sqrt(clamp(eigval, 1e-12, datum::inf));
      vec inv_sqrt_eigval = 1.0 / sqrt_eigval;

      mat A = eigvec * diagmat(sqrt_eigval) * eigvec.t();
      mat AInv = eigvec * diagmat(inv_sqrt_eigval) * eigvec.t();

      // b = AInv * (-s12 - Theta11Inv * c * w22star * ones)
      vec b = AInv * (-s12 - (c * w22star) * sum(Theta11Inv, 1));

      // glmnet minimizes (1/(2n))||Ax-b||^2 + lam*||x||_1 with n=d-1, lam=lambda/(d-1)
      // which is equivalent to minimizing (1/2)||Ax-b||^2 + lambda*||x||_1
      // Our lasso_cd minimizes (1/2)||Ax-b||^2 + lam_cd*||x||_1, so lam_cd = lambda
      vec alpha = lasso_cd(A, b, lambda);

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

      // Track sparsity from exact LASSO zeros: alpha==0 means theta12==c
      for (int i = 0; i < d - 1; i++) {
        graph(notj(i), j) = (alpha(i) == 0.0) ? 1 : 0;
      }
    }

    // Symmetrize (once per outer iteration for numerical safety)
    W = (W + W.t()) / 2.0;

    delta = mean(mean(abs(W - W_old)));
  }

  mat Theta;
  if (!inv_sympd(Theta, W)) {
    // W not positive definite — add small ridge and retry
    Theta = inv_sympd(W + 1e-10 * eye<mat>(d, d));
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
                        int iter_max = 1000) {
  int d = S.n_rows;

  if (c == 0.0) {
    vec eigvals = eig_sym(S);
    c = 1.0 / (d * eigvals(d - 2));  // second largest eigenvalue (0-indexed, sorted ascending)
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
