#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Observed-data log-likelihood, with Omega0 = A(beta) Sigma A(beta)' + Omega.
static double observed_loglik(const arma::mat& gtg, int p, double n1, double n2,
                              const arma::mat& R, const arma::vec& beta,
                              double sigma2, const arma::vec& tau2,
                              double sigmax2, const arma::vec& sigma) {
  const int m = beta.n_elem, d = m + 1;
  arma::mat Omega0(d, d, arma::fill::zeros);
  Omega0(0, 0) = sigma2 + sigmax2 / n1;
  Omega0(arma::span(1, m), 0) = sigma2 * beta;
  Omega0(0, arma::span(1, m)) = sigma2 * beta.t();
  Omega0(arma::span(1, m), arma::span(1, m)) =
    sigma2 * beta * beta.t() + arma::diagmat(tau2) + R % (sigma * sigma.t()) / n2;
  double log_det, sign;
  arma::log_det(log_det, sign, Omega0);
  return -0.5 * p * d * std::log(2.0 * M_PI) - 0.5 * p * log_det -
         0.5 * arma::trace(arma::solve(Omega0, gtg, arma::solve_opts::likely_sympd));
}

// PX-EM algorithm described in the Supplementary Information. The causal
// effects whose zero-based indices are in `fixed` are held at zero.
// [[Rcpp::export]]
Rcpp::List mrmoss_pxem(const arma::vec& gamma_hat, const arma::mat& Gamma_hat,
                       const arma::mat& R, double n1, double n2,
                       const arma::vec& theta0, const arma::uvec& fixed,
                       int maxiter, double tol) {
  const int p = Gamma_hat.n_rows, m = Gamma_hat.n_cols, d = m + 1;
  const arma::mat g = arma::join_horiz(gamma_hat, Gamma_hat);
  const arma::mat gtg = g.t() * g;
  const arma::mat R_inv = arma::inv_sympd(R);
  arma::uvec is_fixed(m, arma::fill::zeros);
  is_fixed.elem(fixed).ones();

  arma::vec beta = theta0.subvec(0, m - 1);
  beta.elem(fixed).zeros();
  double sigma2 = theta0(m);
  arma::vec tau2 = theta0.subvec(m + 1, 2 * m);
  double sigmax2 = theta0(2 * m + 1);
  arma::vec sigma = theta0.subvec(2 * m + 2, 3 * m + 1);

  double loglik = observed_loglik(gtg, p, n1, n2, R, beta, sigma2, tau2, sigmax2, sigma);
  int iter = 0;
  while (iter < maxiter) {
    ++iter;
    // E-step: q_j(r_j) = N(mu_j, Sigma_tilde), with mu_j = M g_j.
    arma::mat A = arma::eye(d, d);
    A(arma::span(1, m), 0) = beta;
    arma::mat Omega_inv(d, d, arma::fill::zeros);
    Omega_inv(0, 0) = n1 / sigmax2;
    Omega_inv(arma::span(1, m), arma::span(1, m)) = n2 * R_inv / (sigma * sigma.t());
    arma::vec prior(d);
    prior(0) = sigma2;
    prior.tail(m) = tau2;
    const arma::mat Sigma_tilde =
      arma::inv_sympd(A.t() * Omega_inv * A + arma::diagmat(1.0 / prior));
    const arma::mat M = Sigma_tilde * A.t() * Omega_inv;
    const arma::mat mu_g = M * gtg;        // sum_j mu_j g_j'
    const arma::mat mu_mu = mu_g * M.t();  // sum_j mu_j mu_j'

    // M-step.
    const double lambda = mu_g(0, 0) / (p * Sigma_tilde(0, 0) + mu_mu(0, 0));
    const double h = mu_mu(0, 0) + p * Sigma_tilde(0, 0);
    for (int i = 0; i < m; ++i) {
      if (is_fixed(i)) continue;
      double numerator = 0.0;
      for (int k = 0; k < m; ++k) {
        double c_k = mu_g(0, k + 1) - mu_mu(0, k + 1) - p * Sigma_tilde(k + 1, 0);
        if (k != i) c_k -= h * beta(k);
        numerator += Omega_inv(i + 1, k + 1) * c_k;
      }
      beta(i) = numerator / (Omega_inv(i + 1, i + 1) * h);
    }
    for (int i = 0; i < m; ++i) {
      tau2(i) = Sigma_tilde(i + 1, i + 1) + mu_mu(i + 1, i + 1) / p;
    }
    arma::mat A_lambda = arma::eye(d, d);
    A_lambda(0, 0) = lambda;
    A_lambda(arma::span(1, m), 0) = beta;
    const arma::mat B = arma::eye(d, d) - A_lambda * M;
    const arma::mat G = p * A_lambda * Sigma_tilde * A_lambda.t() + B * gtg * B.t();  // sum_j G_j
    for (int i = 0; i < m; ++i) {
      double U = 0.0;
      for (int k = 0; k < m; ++k) {
        if (k != i) U += G(k + 1, i + 1) * R_inv(k, i) / sigma(k);
      }
      U *= n2;
      const double W = n2 * G(i + 1, i + 1) * R_inv(i, i);
      sigma(i) = 2.0 * W / (std::sqrt(U * U + 4.0 * p * W) - U);
    }
    sigma2 = Sigma_tilde(0, 0) + mu_mu(0, 0) / p;
    sigmax2 = n1 * lambda * lambda * Sigma_tilde(0, 0) +
              n1 / p * (gtg(0, 0) - 2.0 * lambda * mu_g(0, 0) + lambda * lambda * mu_mu(0, 0));

    // Reduction step.
    beta /= lambda;
    sigma2 *= lambda * lambda;

    const double loglik_new = observed_loglik(gtg, p, n1, n2, R, beta, sigma2, tau2, sigmax2, sigma);
    const double change = std::abs(loglik_new - loglik) / std::abs(loglik_new);
    loglik = loglik_new;
    if (change < tol) break;
  }
  arma::vec theta = arma::join_vert(arma::join_vert(beta, arma::vec({sigma2}), tau2),
                                    arma::join_vert(arma::vec({sigmax2}), sigma));
  return Rcpp::List::create(Rcpp::_["theta"] = theta, Rcpp::_["loglik"] = loglik,
                            Rcpp::_["iteration"] = iter);
}
