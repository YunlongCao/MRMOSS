#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

static arma::vec pack_theta_xp(const arma::vec& beta,
                               double sigma2,
                               const arma::vec& tau2,
                               double sigmax2,
                               const arma::vec& sigma) {
  return arma::join_vert(
    arma::join_vert(
      arma::join_vert(
        arma::join_vert(beta, arma::vec({sigma2})),
        tau2
      ),
      arma::vec({sigmax2})
    ),
    sigma
  );
}

static double quad_xp(const arma::mat& GtG, const arma::vec& v) {
  return arma::as_scalar(v.t() * GtG * v);
}

static double loglikelihood_xp_core(const arma::mat& GtG,
                                    const arma::mat& R,
                                    int p,
                                    int n1,
                                    int n2,
                                    const arma::vec& theta) {
  const int m = R.n_rows;
  arma::vec beta = theta.subvec(0, m - 1);
  double sigma2 = theta(m);
  arma::vec tau2 = theta.subvec(m + 1, 2 * m);
  double sigmax2 = theta(2 * m + 1);
  arma::vec sigma = theta.subvec(2 * m + 2, 3 * m + 1);

  arma::mat Omega0 = arma::zeros<arma::mat>(m + 1, m + 1);
  Omega0(0, 0) = sigma2 + sigmax2 / n1;
  Omega0.submat(1, 0, m, 0) = sigma2 * beta;
  Omega0.submat(0, 1, 0, m) = (sigma2 * beta).t();
  arma::mat lower = sigma2 * (beta * beta.t());
  lower.diag() += tau2;
  lower += (R % (sigma * sigma.t())) / n2;
  Omega0.submat(1, 1, m, m) = lower;

  double log_det_val;
  double sign;
  arma::log_det(log_det_val, sign, Omega0);

  return -p * (m + 1) * std::log(2.0 * M_PI) / 2.0
         -p * log_det_val / 2.0
         -arma::trace(arma::inv(Omega0) * GtG) / 2.0;
}

static arma::mat residual_cross_xp(const arma::mat& GtG,
                                   const arma::mat& C,
                                   const arma::mat& mu0) {
  arma::mat B = arma::eye<arma::mat>(GtG.n_rows, GtG.n_cols) -
                (C * mu0).t();
  return B.t() * GtG * B;
}

static arma::mat make_W_xp(const arma::mat& R_inv,
                           const arma::vec& sigma_old,
                           int n2) {
  arma::mat W = R_inv;
  W.each_col() /= sigma_old;
  W.each_row() /= sigma_old.t();
  W *= n2;
  return W;
}

static void estep_from_W_xp(const arma::mat& W,
                            double omega_x,
                            const arma::vec& beta,
                            double sigma2,
                            const arma::vec& tau2,
                            arma::mat& Sigma_tilde,
                            arma::mat& mu0) {
  const int m = beta.n_elem;
  const int d = m + 1;
  arma::vec Wbeta = W * beta;

  arma::mat K = arma::zeros<arma::mat>(d, d);
  K(0, 0) = 1.0 / sigma2 + omega_x + arma::dot(beta, Wbeta);
  K.submat(1, 0, m, 0) = Wbeta;
  K.submat(0, 1, 0, m) = Wbeta.t();
  K.submat(1, 1, m, m) = W;
  for (int i = 0; i < m; ++i) {
    K(i + 1, i + 1) += 1.0 / tau2(i);
  }

  Sigma_tilde = arma::inv_sympd(K);

  arma::mat H = arma::zeros<arma::mat>(d, d);
  H(0, 0) = omega_x;
  H.submat(0, 1, 0, m) = Wbeta.t();
  H.submat(1, 1, m, m) = W;
  mu0 = Sigma_tilde * H;
}

static void update_sigma_xp(const arma::mat& G_const,
                            const arma::mat& residual_cross,
                            const arma::mat& R_inv,
                            const arma::vec& sigma_old,
                            int p,
                            int n2,
                            arma::vec& sigma_new,
                            double rd) {
  const int m = sigma_old.n_elem;
  arma::mat S_out = p * G_const.submat(1, 1, m, m) +
                    residual_cross.submat(1, 1, m, m);

  for (int i = 0; i < m; ++i) {
    arma::vec weights = R_inv.row(i).t() / sigma_old;
    double Sii = S_out(i, i);
    double Ui = n2 * (arma::dot(S_out.col(i), weights) -
                      Sii * R_inv(i, i) / sigma_old(i));
    double Vi = n2 * Sii * R_inv(i, i);
    sigma_new(i) = 2.0 * Vi / (std::sqrt(Ui * Ui + 4.0 * p * Vi) - Ui);
    sigma_new(i) *= rd;
  }
}

static double beta_num_xp(const arma::mat& GtG,
                          const arma::mat& Omega_inv,
                          const arma::mat& C,
                          const arma::mat& Sigma_tilde,
                          const arma::mat& mu0,
                          int p,
                          int r) {
  arma::vec u0 = mu0.row(0).t();
  arma::mat B = arma::eye<arma::mat>(GtG.n_rows, GtG.n_cols) -
                (C * mu0).t();
  arma::rowvec weighted_residual = u0.t() * GtG * B;

  double out = -p * arma::as_scalar(Omega_inv.row(r) * C *
                                    Sigma_tilde.col(0));
  out += arma::as_scalar(weighted_residual * Omega_inv.row(r).t());
  return out;
}

static arma::vec beta_num_all_xp(const arma::mat& GtG,
                                 const arma::mat& W,
                                 const arma::mat& Sigma_tilde,
                                 const arma::mat& mu0,
                                 const arma::vec& beta_old,
                                 int p) {
  const int m = beta_old.n_elem;
  arma::vec Wdiag = W.diag();

  double s00 = Sigma_tilde(0, 0);
  arma::vec s_out = Sigma_tilde.submat(1, 0, m, 0);
  arma::vec z1 = beta_old * s00 + s_out;
  arma::vec term1 = -static_cast<double>(p) *
                    (W * z1 - Wdiag % beta_old * s00);

  arma::vec u0 = mu0.row(0).t();
  arma::vec a = GtG * u0;
  arma::vec termA = W * a.subvec(1, m);

  arma::vec h = mu0 * a;
  double h0 = h(0);
  arma::vec z2 = beta_old * h0 + h.subvec(1, m);
  arma::vec second = W * z2 - Wdiag % beta_old * h0;

  return term1 + termA - second;
}

// [[Rcpp::export]]
List MRMOSS_PX_xprod_joint_cpp(arma::vec gamma_hat,
                               arma::mat Gamma_hat,
                               arma::mat R,
                               int n1,
                               int n2,
                               arma::vec theta0,
                               int maxiter,
                               double rd,
                               double tol = 1e-8) {
  const int p = Gamma_hat.n_rows;
  const int m = Gamma_hat.n_cols;
  const int d = m + 1;
  arma::mat g = arma::join_horiz(gamma_hat, Gamma_hat);
  arma::mat GtG = g.t() * g;
  arma::mat R_inv = arma::inv_sympd(R);
  arma::vec e0 = arma::zeros<arma::vec>(d);
  e0(0) = 1.0;

  arma::vec theta_new;
  arma::vec beta_new = theta0.subvec(0, m - 1);
  double sigma2_new = theta0(m);
  arma::vec tau2_new = theta0.subvec(m + 1, 2 * m);
  double sigmax2_new = theta0(2 * m + 1);
  arma::vec sigma_new = theta0.subvec(2 * m + 2, 3 * m + 1);

  double loglikeli = 0.0;
  double dif = 1.0;
  int iter = 0;
  arma::vec loglikeli_vec(maxiter + 2);
  loglikeli_vec(0) = 0.0;

  while ((iter <= maxiter) && (dif >= tol)) {
    arma::vec beta_old = beta_new;
    double sigma2_old = sigma2_new;
    arma::vec tau2_old = tau2_new;
    double sigmax2_old = sigmax2_new;
    arma::vec sigma_old = sigma_new;

    arma::mat A = arma::eye<arma::mat>(d, d);
    A.col(0).subvec(1, m) = beta_old;

    arma::mat W = make_W_xp(R_inv, sigma_old, n2);
    arma::mat Sigma_tilde;
    arma::mat mu0;
    estep_from_W_xp(W, n1 / sigmax2_old, beta_old, sigma2_old, tau2_old,
                    Sigma_tilde, mu0);
    arma::vec u0 = mu0.row(0).t();
    double mu0_sq_sum = quad_xp(GtG, u0);
    double gamma_mu0_sum = arma::as_scalar(e0.t() * GtG * u0);

    double lambda = gamma_mu0_sum /
                    (p * Sigma_tilde(0, 0) + mu0_sq_sum);

    arma::mat A_lambda = A;
    A_lambda(0, 0) = lambda;
    arma::mat G_lambda = A_lambda * Sigma_tilde * A_lambda.t();
    arma::mat residual_cross = residual_cross_xp(GtG, A_lambda, mu0);

    double beta_den_common = p * Sigma_tilde(0, 0) + mu0_sq_sum;
    arma::vec beta_num = beta_num_all_xp(
      GtG, W, Sigma_tilde, mu0, beta_old, p
    );
    arma::vec beta_den = W.diag() * beta_den_common;
    for (int i = 0; i < m; ++i) {
      int r = i + 1;
      beta_new(i) = beta_num(i) / beta_den(i);

      arma::vec ur = mu0.row(r).t();
      tau2_new(i) = Sigma_tilde(r, r) + quad_xp(GtG, ur) / p;
    }

    update_sigma_xp(G_lambda, residual_cross, R_inv, sigma_old, p, n2,
                    sigma_new, rd);

    sigma2_new = Sigma_tilde(0, 0) + mu0_sq_sum / p;
    arma::vec vx = e0 - lambda * u0;
    sigmax2_new = n1 * lambda * Sigma_tilde(0, 0) +
                  n1 * quad_xp(GtG, vx) / p;

    beta_new = beta_new / lambda;
    sigma2_new = lambda * lambda * sigma2_new;

    theta_new = pack_theta_xp(beta_new, sigma2_new, tau2_new, sigmax2_new,
                              sigma_new);

    loglikeli = loglikelihood_xp_core(GtG, R, p, n1, n2, theta_new);
    loglikeli_vec(iter + 1) = loglikeli;
    dif = std::abs(loglikeli_vec(iter + 1) - loglikeli_vec(iter));
    dif = std::min(dif, dif / (1.0 + std::abs(loglikeli_vec(iter + 1))));
    iter++;
  }

  return Rcpp::List::create(
    Rcpp::Named("beta") = theta_new.subvec(0, m - 1).t(),
    Rcpp::Named("theta") = theta_new.t(),
    Rcpp::Named("loglikeli") = loglikeli_vec.subvec(1, iter),
    Rcpp::Named("iteration") = iter - 1
  );
}

// [[Rcpp::export]]
List MRMOSS_PX_xprod_lrt_diag_cpp(arma::vec gamma_hat,
                                  arma::mat Gamma_hat,
                                  arma::mat R,
                                  int n1,
                                  int n2,
                                  arma::vec theta0,
                                  arma::uvec test,
                                  int maxiter,
                                  double rd,
                                  double tol = 1e-5,
                                  int miniter = 20) {
  const int p = Gamma_hat.n_rows;
  const int m = Gamma_hat.n_cols;
  const int d = m + 1;
  arma::mat g = arma::join_horiz(gamma_hat, Gamma_hat);
  arma::mat GtG = g.t() * g;
  arma::mat R_inv = arma::inv_sympd(R);
  arma::vec e0 = arma::zeros<arma::vec>(d);
  e0(0) = 1.0;

  arma::vec theta_new;
  arma::vec beta_new = theta0.subvec(0, m - 1);
  double sigma2_new = theta0(m);
  arma::vec tau2_new = theta0.subvec(m + 1, 2 * m);
  double sigmax2_new = theta0(2 * m + 1);
  arma::vec sigma_new = theta0.subvec(2 * m + 2, 3 * m + 1);

  double loglikeli = 0.0;
  double dif = 1.0;
  int iter = 0;
  arma::vec loglikeli_vec(maxiter + 2);
  loglikeli_vec(0) = 0.0;

  while ((iter <= maxiter) && ((iter < miniter) || (dif >= tol))) {
    arma::vec beta_old = beta_new;
    double sigma2_old = sigma2_new;
    arma::vec tau2_old = tau2_new;
    double sigmax2_old = sigmax2_new;
    arma::vec sigma_old = sigma_new;

    arma::mat A = arma::eye<arma::mat>(d, d);
    A.col(0).subvec(1, m) = beta_old;

    arma::mat W = make_W_xp(R_inv, sigma_old, n2);
    arma::mat Sigma_tilde;
    arma::mat mu0;
    estep_from_W_xp(W, n1 / sigmax2_old, beta_old, sigma2_old, tau2_old,
                    Sigma_tilde, mu0);
    arma::vec u0 = mu0.row(0).t();
    double mu0_sq_sum = quad_xp(GtG, u0);
    double gamma_mu0_sum = arma::as_scalar(e0.t() * GtG * u0);

    double lambda = gamma_mu0_sum /
                    (p * Sigma_tilde(0, 0) + mu0_sq_sum);

    arma::mat A_lambda = A;
    A_lambda(0, 0) = lambda;
    arma::mat G_lambda = A_lambda * Sigma_tilde * A_lambda.t();
    arma::mat residual_cross = residual_cross_xp(GtG, A_lambda, mu0);

    double beta_den_common = p * Sigma_tilde(0, 0) + mu0_sq_sum;
    arma::vec beta_num = beta_num_all_xp(
      GtG, W, Sigma_tilde, mu0, beta_old, p
    );
    arma::vec beta_den = W.diag() * beta_den_common;
    for (int i = 0; i < m; ++i) {
      int r = i + 1;
      beta_new(i) = beta_num(i) / beta_den(i);

      arma::vec ur = mu0.row(r).t();
      tau2_new(i) = Sigma_tilde(r, r) + quad_xp(GtG, ur) / p;
    }

    update_sigma_xp(G_lambda, residual_cross, R_inv, sigma_old, p, n2,
                    sigma_new, rd);

    sigma2_new = Sigma_tilde(0, 0) + mu0_sq_sum / p;
    arma::vec vx = e0 - lambda * u0;
    sigmax2_new = n1 * lambda * Sigma_tilde(0, 0) +
                  n1 * quad_xp(GtG, vx) / p;

    beta_new = beta_new / lambda;
    sigma2_new = lambda * lambda * sigma2_new;

    theta_new = pack_theta_xp(beta_new, sigma2_new, tau2_new, sigmax2_new,
                              sigma_new);

    loglikeli = loglikelihood_xp_core(GtG, R, p, n1, n2, theta_new);
    loglikeli_vec(iter + 1) = loglikeli;
    dif = std::abs(loglikeli_vec(iter + 1) - loglikeli_vec(iter));
    iter++;
  }

  const arma::vec beta_full = beta_new;
  const double sigma2_full = sigma2_new;
  const arma::vec tau2_full = tau2_new;
  const double sigmax2_full = sigmax2_new;
  const arma::vec sigma_full = sigma_new;

  arma::mat loglikeli_null = arma::zeros<arma::mat>(test.n_elem, maxiter + 2);
  arma::mat theta_null = arma::zeros<arma::mat>(test.n_elem, theta0.n_elem);
  arma::vec LRT = arma::zeros<arma::vec>(test.n_elem);
  arma::vec pvalue = arma::zeros<arma::vec>(test.n_elem);
  arma::vec log_pvalue = arma::zeros<arma::vec>(test.n_elem);
  arma::ivec iteration_null = arma::zeros<arma::ivec>(test.n_elem);

  for (arma::uword t_idx = 0; t_idx < test.n_elem; ++t_idx) {
    arma::uword t = test(t_idx) - 1;
    double dif0 = 1.0;
    int iter0 = 0;
    arma::vec beta_null = beta_full;
    double sigma2_null = sigma2_full;
    arma::vec tau2_null = tau2_full;
    double sigmax2_null = sigmax2_full;
    arma::vec sigma_null = sigma_full;

    while ((iter0 <= maxiter) && ((iter0 < miniter) || (dif0 >= tol))) {
      arma::vec beta_old = beta_null;
      beta_old(t) = 0.0;
      double sigma2_old = sigma2_null;
      arma::vec tau2_old = tau2_null;
      double sigmax2_old = sigmax2_null;
      arma::vec sigma_old = sigma_null;

      arma::mat A = arma::eye<arma::mat>(d, d);
      A.col(0).subvec(1, m) = beta_old;

      arma::mat W = make_W_xp(R_inv, sigma_old, n2);
      arma::mat Sigma_tilde;
      arma::mat mu0;
      estep_from_W_xp(W, n1 / sigmax2_old, beta_old, sigma2_old, tau2_old,
                      Sigma_tilde, mu0);
      arma::mat G = A * Sigma_tilde * A.t();
      arma::mat residual_cross = residual_cross_xp(GtG, A, mu0);
      arma::vec u0 = mu0.row(0).t();
      double mu0_sq_sum = quad_xp(GtG, u0);
      double beta_den_common = p * Sigma_tilde(0, 0) + mu0_sq_sum;
      arma::vec beta_num = beta_num_all_xp(
        GtG, W, Sigma_tilde, mu0, beta_old, p
      );
      arma::vec beta_den = W.diag() * beta_den_common;

      for (arma::uword i = 0; i < static_cast<arma::uword>(m); ++i) {
        int r = i + 1;
        if (i == t) {
          beta_null(i) = 0.0;
        } else {
          beta_null(i) = beta_num(i) / beta_den(i);
        }
        arma::vec ur = mu0.row(r).t();
        tau2_null(i) = Sigma_tilde(r, r) + quad_xp(GtG, ur) / p;
      }

      update_sigma_xp(G, residual_cross, R_inv, sigma_old, p, n2,
                      sigma_null, rd);

      sigma2_null = Sigma_tilde(0, 0) + mu0_sq_sum / p;
      arma::vec vx = e0 - u0;
      sigmax2_null = n1 * Sigma_tilde(0, 0) +
                    n1 * quad_xp(GtG, vx) / p;

      arma::vec theta_current = pack_theta_xp(beta_null, sigma2_null, tau2_null,
                                              sigmax2_null, sigma_null);
      theta_null.row(t_idx) = theta_current.t();
      loglikeli_null(t_idx, iter0 + 1) =
        loglikelihood_xp_core(GtG, R, p, n1, n2, theta_current);

      ++iter0;
      dif0 = std::abs(loglikeli_null(t_idx, iter0) -
                      loglikeli_null(t_idx, iter0 - 1));
      dif0 = std::min(dif0,
                      dif0 / (1.0 + std::abs(loglikeli_null(t_idx, iter0))));
    }

    LRT(t_idx) = 2.0 * (loglikeli_vec(iter) - loglikeli_null(t_idx, iter0));
    if (!R_finite(LRT(t_idx))) {
      pvalue(t_idx) = NA_REAL;
      log_pvalue(t_idx) = NA_REAL;
    } else {
      if (LRT(t_idx) < 0.0) LRT(t_idx) = 0.0;
      pvalue(t_idx) = R::pchisq(LRT(t_idx), 1.0, false, false);
      log_pvalue(t_idx) = R::pchisq(LRT(t_idx), 1.0, false, true);
    }
    iteration_null(t_idx) = iter0 - 1;
  }

  arma::mat loglikeli_overall = arma::zeros<arma::mat>(1, maxiter + 2);
  arma::mat theta_overall = arma::zeros<arma::mat>(1, theta0.n_elem);
  double LRT_overall = 0.0;
  double pvalue_overall = 0.0;
  double log_pvalue_overall = 0.0;

  double dif_overall = 1.0;
  int iter_overall = 0;
  double sigma2_overall = sigma2_full;
  arma::vec tau2_overall = tau2_full;
  double sigmax2_overall = sigmax2_full;
  arma::vec sigma_overall = sigma_full;

  while ((iter_overall <= maxiter) &&
         ((iter_overall < miniter) || (dif_overall >= tol))) {
    arma::vec beta_new_overall = arma::zeros<arma::vec>(m);
    double sigma2_old = sigma2_overall;
    arma::vec tau2_old = tau2_overall;
    double sigmax2_old = sigmax2_overall;
    arma::vec sigma_old = sigma_overall;

    arma::mat A = arma::eye<arma::mat>(d, d);

    arma::mat W = make_W_xp(R_inv, sigma_old, n2);
    arma::mat Sigma_tilde;
    arma::mat mu0;
    estep_from_W_xp(W, n1 / sigmax2_old, beta_new_overall, sigma2_old,
                    tau2_old, Sigma_tilde, mu0);
    arma::mat G = A * Sigma_tilde * A.t();
    arma::mat residual_cross = residual_cross_xp(GtG, A, mu0);

    for (int i = 0; i < m; ++i) {
      int r = i + 1;
      arma::vec ur = mu0.row(r).t();
      tau2_overall(i) = Sigma_tilde(r, r) + quad_xp(GtG, ur) / p;
    }

    update_sigma_xp(G, residual_cross, R_inv, sigma_old, p, n2,
                    sigma_overall, rd);

    arma::vec u0 = mu0.row(0).t();
    sigma2_overall = Sigma_tilde(0, 0) + quad_xp(GtG, u0) / p;
    arma::vec vx = e0 - u0;
    sigmax2_overall = n1 * Sigma_tilde(0, 0) +
                  n1 * quad_xp(GtG, vx) / p;

    arma::vec theta_current = pack_theta_xp(beta_new_overall, sigma2_overall,
                                            tau2_overall, sigmax2_overall,
                                            sigma_overall);
    theta_overall = theta_current.t();
    loglikeli_overall(iter_overall + 1) =
      loglikelihood_xp_core(GtG, R, p, n1, n2, theta_current);

    ++iter_overall;
    dif_overall = std::abs(loglikeli_overall(iter_overall) -
                           loglikeli_overall(iter_overall - 1));
    dif_overall = std::min(
      dif_overall,
      dif_overall / (1.0 + std::abs(loglikeli_overall(iter_overall)))
    );
  }

  LRT_overall = 2.0 * (loglikeli_vec(iter) -
                       loglikeli_overall(iter_overall));
  if (!R_finite(LRT_overall)) {
    pvalue_overall = NA_REAL;
    log_pvalue_overall = NA_REAL;
  } else {
    if (LRT_overall < 0.0) LRT_overall = 0.0;
    pvalue_overall = R::pchisq(LRT_overall, m, false, false);
    log_pvalue_overall = R::pchisq(LRT_overall, m, false, true);
  }

  return Rcpp::List::create(
    Rcpp::Named("beta") = theta_new.subvec(0, m - 1).t(),
    Rcpp::Named("theta") = theta_new.t(),
    Rcpp::Named("theta_null") = theta_null,
    Rcpp::Named("loglikeli") = loglikeli_vec.subvec(1, iter),
    Rcpp::Named("loglikeli_null") =
      loglikeli_null.cols(1, loglikeli_null.n_cols - 1),
    Rcpp::Named("iteration") = iter - 1,
    Rcpp::Named("iteration_null") = iteration_null,
    Rcpp::Named("iteration_overall") = iter_overall - 1,
    Rcpp::Named("LRT") = LRT,
    Rcpp::Named("LRT_overall") = LRT_overall,
    Rcpp::Named("pvalue") = pvalue,
    Rcpp::Named("log_pvalue") = log_pvalue,
    Rcpp::Named("theta_overall") = theta_overall,
    Rcpp::Named("pvalue_overall") = pvalue_overall,
    Rcpp::Named("log_pvalue_overall") = log_pvalue_overall
  );
}

// [[Rcpp::export]]
List MRMOSS_PX_xprod_subset_lrt_cpp(arma::vec gamma_hat,
                                    arma::mat Gamma_hat,
                                    arma::mat R,
                                    int n1,
                                    int n2,
                                    arma::vec theta0,
                                    List subsets,
                                    CharacterVector subset_names,
                                    int maxiter,
                                    double rd,
                                    double tol = 1e-5,
                                    int miniter = 20) {
  const int p = Gamma_hat.n_rows;
  const int m = Gamma_hat.n_cols;
  const int d = m + 1;
  arma::mat g = arma::join_horiz(gamma_hat, Gamma_hat);
  arma::mat GtG = g.t() * g;
  arma::mat R_inv = arma::inv_sympd(R);
  arma::vec e0 = arma::zeros<arma::vec>(d);
  e0(0) = 1.0;

  arma::vec theta_new;
  arma::vec beta_new = theta0.subvec(0, m - 1);
  double sigma2_new = theta0(m);
  arma::vec tau2_new = theta0.subvec(m + 1, 2 * m);
  double sigmax2_new = theta0(2 * m + 1);
  arma::vec sigma_new = theta0.subvec(2 * m + 2, 3 * m + 1);

  double loglikeli = 0.0;
  double dif = 1.0;
  int iter = 0;
  arma::vec loglikeli_vec(maxiter + 2);
  loglikeli_vec(0) = 0.0;

  while ((iter <= maxiter) && ((iter < miniter) || (dif >= tol))) {
    arma::vec beta_old = beta_new;
    double sigma2_old = sigma2_new;
    arma::vec tau2_old = tau2_new;
    double sigmax2_old = sigmax2_new;
    arma::vec sigma_old = sigma_new;

    arma::mat A = arma::eye<arma::mat>(d, d);
    A.col(0).subvec(1, m) = beta_old;

    arma::mat W = make_W_xp(R_inv, sigma_old, n2);
    arma::mat Sigma_tilde;
    arma::mat mu0;
    estep_from_W_xp(W, n1 / sigmax2_old, beta_old, sigma2_old, tau2_old,
                    Sigma_tilde, mu0);
    arma::vec u0 = mu0.row(0).t();
    double mu0_sq_sum = quad_xp(GtG, u0);
    double gamma_mu0_sum = arma::as_scalar(e0.t() * GtG * u0);

    double lambda = gamma_mu0_sum /
                    (p * Sigma_tilde(0, 0) + mu0_sq_sum);

    arma::mat A_lambda = A;
    A_lambda(0, 0) = lambda;
    arma::mat G_lambda = A_lambda * Sigma_tilde * A_lambda.t();
    arma::mat residual_cross = residual_cross_xp(GtG, A_lambda, mu0);

    double beta_den_common = p * Sigma_tilde(0, 0) + mu0_sq_sum;
    arma::vec beta_num = beta_num_all_xp(
      GtG, W, Sigma_tilde, mu0, beta_old, p
    );
    arma::vec beta_den = W.diag() * beta_den_common;
    for (int i = 0; i < m; ++i) {
      int r = i + 1;
      beta_new(i) = beta_num(i) / beta_den(i);

      arma::vec ur = mu0.row(r).t();
      tau2_new(i) = Sigma_tilde(r, r) + quad_xp(GtG, ur) / p;
    }

    update_sigma_xp(G_lambda, residual_cross, R_inv, sigma_old, p, n2,
                    sigma_new, rd);

    sigma2_new = Sigma_tilde(0, 0) + mu0_sq_sum / p;
    arma::vec vx = e0 - lambda * u0;
    sigmax2_new = n1 * lambda * Sigma_tilde(0, 0) +
                  n1 * quad_xp(GtG, vx) / p;

    beta_new = beta_new / lambda;
    sigma2_new = lambda * lambda * sigma2_new;

    theta_new = pack_theta_xp(beta_new, sigma2_new, tau2_new, sigmax2_new,
                              sigma_new);

    loglikeli = loglikelihood_xp_core(GtG, R, p, n1, n2, theta_new);
    loglikeli_vec(iter + 1) = loglikeli;
    dif = std::abs(loglikeli_vec(iter + 1) - loglikeli_vec(iter));
    iter++;
  }

  const arma::vec beta_full = beta_new;
  const double sigma2_full = sigma2_new;
  const arma::vec tau2_full = tau2_new;
  const double sigmax2_full = sigmax2_new;
  const arma::vec sigma_full = sigma_new;
  const double loglik_full = loglikeli_vec(iter);
  const int n_subsets = subsets.size();

  NumericVector LRT(n_subsets);
  NumericVector pvalue(n_subsets);
  NumericVector log_pvalue(n_subsets);
  NumericVector neglog10p(n_subsets);
  NumericVector loglik_null(n_subsets);
  IntegerVector iteration_null(n_subsets);
  IntegerVector df(n_subsets);

  for (int s = 0; s < n_subsets; ++s) {
    IntegerVector idx = subsets[s];
    std::vector<int> fixed(m, 0);
    int n_fixed = 0;
    for (int j = 0; j < idx.size(); ++j) {
      int z = idx[j] - 1;
      if (z >= 0 && z < m && fixed[z] == 0) {
        fixed[z] = 1;
        n_fixed++;
      }
    }
    if (n_fixed == 0) {
      LRT[s] = NA_REAL;
      pvalue[s] = NA_REAL;
      log_pvalue[s] = NA_REAL;
      neglog10p[s] = NA_REAL;
      loglik_null[s] = NA_REAL;
      iteration_null[s] = NA_INTEGER;
      df[s] = 0;
      continue;
    }
    df[s] = n_fixed;

    double dif0 = 1.0;
    int iter0 = 0;
    arma::vec beta_null = beta_full;
    for (int i = 0; i < m; ++i) {
      if (fixed[i]) beta_null(i) = 0.0;
    }
    double sigma2_null = sigma2_full;
    arma::vec tau2_null = tau2_full;
    double sigmax2_null = sigmax2_full;
    arma::vec sigma_null = sigma_full;
    arma::vec loglikeli_null_vec(maxiter + 2);
    loglikeli_null_vec(0) = 0.0;

    while ((iter0 <= maxiter) && ((iter0 < miniter) || (dif0 >= tol))) {
      arma::vec beta_old = beta_null;
      for (int i = 0; i < m; ++i) {
        if (fixed[i]) beta_old(i) = 0.0;
      }
      double sigma2_old = sigma2_null;
      arma::vec tau2_old = tau2_null;
      double sigmax2_old = sigmax2_null;
      arma::vec sigma_old = sigma_null;

      arma::mat A = arma::eye<arma::mat>(d, d);
      A.col(0).subvec(1, m) = beta_old;

      arma::mat W = make_W_xp(R_inv, sigma_old, n2);
      arma::mat Sigma_tilde;
      arma::mat mu0;
      estep_from_W_xp(W, n1 / sigmax2_old, beta_old, sigma2_old, tau2_old,
                      Sigma_tilde, mu0);
      arma::mat G = A * Sigma_tilde * A.t();
      arma::mat residual_cross = residual_cross_xp(GtG, A, mu0);
      arma::vec u0 = mu0.row(0).t();
      double mu0_sq_sum = quad_xp(GtG, u0);
      double beta_den_common = p * Sigma_tilde(0, 0) + mu0_sq_sum;
      arma::vec beta_num = beta_num_all_xp(
        GtG, W, Sigma_tilde, mu0, beta_old, p
      );
      arma::vec beta_den = W.diag() * beta_den_common;

      for (int i = 0; i < m; ++i) {
        int r = i + 1;
        if (fixed[i]) {
          beta_null(i) = 0.0;
        } else {
          beta_null(i) = beta_num(i) / beta_den(i);
        }
        arma::vec ur = mu0.row(r).t();
        tau2_null(i) = Sigma_tilde(r, r) + quad_xp(GtG, ur) / p;
      }

      update_sigma_xp(G, residual_cross, R_inv, sigma_old, p, n2,
                      sigma_null, rd);

      sigma2_null = Sigma_tilde(0, 0) + mu0_sq_sum / p;
      arma::vec vx = e0 - u0;
      sigmax2_null = n1 * Sigma_tilde(0, 0) +
                    n1 * quad_xp(GtG, vx) / p;

      arma::vec theta_current = pack_theta_xp(beta_null, sigma2_null, tau2_null,
                                              sigmax2_null, sigma_null);
      loglikeli_null_vec(iter0 + 1) =
        loglikelihood_xp_core(GtG, R, p, n1, n2, theta_current);

      ++iter0;
      dif0 = std::abs(loglikeli_null_vec(iter0) -
                      loglikeli_null_vec(iter0 - 1));
      dif0 = std::min(dif0,
                      dif0 / (1.0 + std::abs(loglikeli_null_vec(iter0))));
    }

    loglik_null[s] = loglikeli_null_vec(iter0);
    LRT[s] = 2.0 * (loglik_full - loglik_null[s]);
    if (R_IsNA(LRT[s]) || !R_finite(LRT[s]) || LRT[s] < 0.0) {
      pvalue[s] = NA_REAL;
      log_pvalue[s] = NA_REAL;
      neglog10p[s] = NA_REAL;
    } else {
      pvalue[s] = R::pchisq(LRT[s], static_cast<double>(n_fixed), false, false);
      log_pvalue[s] = R::pchisq(LRT[s], static_cast<double>(n_fixed), false, true);
      neglog10p[s] = -log_pvalue[s] / std::log(10.0);
    }
    iteration_null[s] = iter0 - 1;
  }

  return Rcpp::List::create(
    Rcpp::Named("subset_names") = subset_names,
    Rcpp::Named("beta") = beta_new.t(),
    Rcpp::Named("theta") = theta_new.t(),
    Rcpp::Named("loglikeli") = loglikeli_vec.subvec(1, iter),
    Rcpp::Named("loglik_full") = loglik_full,
    Rcpp::Named("iteration") = iter - 1,
    Rcpp::Named("LRT_subset") = LRT,
    Rcpp::Named("df_subset") = df,
    Rcpp::Named("pvalue_subset") = pvalue,
    Rcpp::Named("log_pvalue_subset") = log_pvalue,
    Rcpp::Named("neglog10p_subset") = neglog10p,
    Rcpp::Named("loglik_null_subset") = loglik_null,
    Rcpp::Named("iteration_null_subset") = iteration_null
  );
}
