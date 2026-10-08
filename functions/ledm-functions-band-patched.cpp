#include <Rcpp.h>
using namespace Rcpp;

inline int idx3(int i, int j, int k, int d1, int d2) {
  // R column-major for array(dim = c(d1, d2, d3))
  return i + d1 * j + d1 * d2 * k;
}

inline double logsumexp_vec(const NumericVector& x) {
  int n = x.size();
  
  double m = R_NegInf;
  for (int i = 0; i <= n - 1; ++i) {
    if (R_finite(x[i]) && x[i] > m) m = x[i];
  }
  
  if (!R_finite(m)) return R_NegInf;
  
  double s = 0.0;
  for (int i = 0; i <= n - 1; ++i) {
    if (R_finite(x[i])) s += std::exp(x[i] - m);
  }
  
  return m + std::log(s);
}

// [[Rcpp::export]]
NumericVector get_log_prob_emissions_cpp(IntegerVector x,
                                         IntegerVector z,
                                         int B) {
  int p = z.size();
  int l = x.size();

  NumericVector out(p * (l + 1) * 4, R_NegInf);
  
  const double log_eps3 = std::log(1.0 / 3.0);
  const double log_quarter = std::log(0.25);
  
  int stride_j = p;
  int stride_k = p * (l + 1);
  
  // k offsets
  int K = 0 * stride_k;
  int D = 1 * stride_k;
  int S = 2 * stride_k;
  int I = 3 * stride_k;
  
  for (int j = 0; j <= p - 1; ++j) {
    int t_min = std::max(0, j - B);
    int t_max = std::min(l, j + B);
    
    // ---- D ----
    int idxD = j + D + t_min * stride_j;
    for (int t = t_min; t <= t_max; ++t, idxD += stride_j) {
      out[idxD] = 0.0;
    }
    
    // ---- K and S ----
    int idxK = j + K + t_min * stride_j;
    int idxS = j + S + t_min * stride_j;
    
    int t_max_m1 = std::min(t_max, l - 1);
    for (int t = t_min; t <= t_max_m1; ++t, idxK += stride_j, idxS += stride_j) {
      bool eq = (z[j] == x[t]);
      out[idxK] = eq ? 0.0 : R_NegInf;
      out[idxS] = eq ? R_NegInf : log_eps3;
    }
    
    // ---- I ----
    int idxI = j + I + t_min * stride_j;
    idxK = j + K + t_min * stride_j;
    
    int t_max_m2 = std::min(t_max, l - 2);
    for (int t = t_min; t <= t_max_m2; ++t, idxI += stride_j, idxK += stride_j) {
      out[idxI] = out[idxK] + log_quarter;
    }
  }
  
  out.attr("dim") = IntegerVector::create(p, l + 1, 4);
  return out;
}

// [[Rcpp::export]]
NumericMatrix get_log_alpha_cpp(NumericVector log_theta,
                                NumericVector log_prob_emissions,
                                int B) {
  IntegerVector dims = log_prob_emissions.attr("dim");
  int p = dims[0];
  int l = dims[1] - 1;

  double logK = log_theta[0];
  double logD = log_theta[1];
  double logS = log_theta[2];
  double logI = log_theta[3];
  
  NumericMatrix log_alpha(p + 1, l + 1);
  std::fill(log_alpha.begin(), log_alpha.end(), R_NegInf);
  log_alpha(0, 0) = 0.0;
  NumericVector vals(4);
  
  for (int j = 1; j <= p; ++j) {
    int t_min = std::max(0, j - B);
    int t_max = std::min(l, j + B);
    
    for (int t = t_min; t <= t_max; ++t) {
      vals[0] = logD + log_alpha(j - 1, t) +
        log_prob_emissions[idx3(j - 1, t, 1, p, l + 1)];
      
      vals[1] = R_NegInf;
      if (t >= 1) {
        vals[1] = logK + log_alpha(j - 1, t - 1) +
          log_prob_emissions[idx3(j - 1, t - 1, 0, p, l + 1)];
      }
      
      vals[2] = R_NegInf;
      if (t >= 1) {
        vals[2] = logS + log_alpha(j - 1, t - 1) +
          log_prob_emissions[idx3(j - 1, t - 1, 2, p, l + 1)];
      }
      
      vals[3] = R_NegInf;
      if (t >= 2) {
        vals[3] = logI + log_alpha(j - 1, t - 2) +
          log_prob_emissions[idx3(j - 1, t - 2, 3, p, l + 1)];
      }
      
      log_alpha(j, t) = logsumexp_vec(vals);
    }
  }
  return log_alpha;
}

// [[Rcpp::export]]
NumericMatrix get_log_beta_cpp(NumericVector log_theta,
                               NumericVector log_prob_emissions,
                               int B) {
  IntegerVector dims = log_prob_emissions.attr("dim");
  int p = dims[0];
  int l = dims[1] - 1;

  double logK = log_theta[0];
  double logD = log_theta[1];
  double logS = log_theta[2];
  double logI = log_theta[3];
  
  NumericMatrix log_beta(p + 1, l + 1);
  std::fill(log_beta.begin(), log_beta.end(), R_NegInf);
  log_beta(p, l) = 0.0;
  NumericVector vals(4);
  
  for (int j = p - 1; j >= 0; --j) {
    int t_min = std::max(0, j - B);
    int t_max = std::min(l, j + B);
    
    for (int t = t_min; t <= t_max; ++t) {
      vals[0] = logD + log_beta(j + 1, t) +
        log_prob_emissions[idx3(j, t, 1, p, l + 1)];
      
      vals[1] = R_NegInf;
      if (t <= l - 1) {
        vals[1] = logK + log_beta(j + 1, t + 1) +
          log_prob_emissions[idx3(j, t, 0, p, l + 1)];
      }
      
      vals[2] = R_NegInf;
      if (t <= l - 1) {
        vals[2] = logS + log_beta(j + 1, t + 1) +
          log_prob_emissions[idx3(j, t, 2, p, l + 1)];
      }
      
      vals[3] = R_NegInf;
      if (t <= l - 2) {
        vals[3] = logI + log_beta(j + 1, t + 2) +
          log_prob_emissions[idx3(j, t, 3, p, l + 1)];
      }
      
      log_beta(j, t) = logsumexp_vec(vals);
    }
  }
  return log_beta;
}

// [[Rcpp::export]]
NumericVector get_log_gamma_num_cpp(NumericVector log_theta,
                                    NumericMatrix log_alpha,
                                    NumericMatrix log_beta,
                                    NumericVector log_prob_emissions,
                                    int B) {
  IntegerVector dims = log_prob_emissions.attr("dim");
  int p = dims[0];
  int l = dims[1] - 1;

  double logK = log_theta[0];
  double logD = log_theta[1];
  double logS = log_theta[2];
  double logI = log_theta[3];
  
  NumericVector log_gamma(p * (l + 1) * 4, R_NegInf);
  log_gamma.attr("dim") = IntegerVector::create(p, l + 1, 4);
  
  for (int j = 0; j <= p - 1; ++j) {
    int t_min = std::max(0, j - B);
    int t_max = std::min(l, j + B);
    
    for (int t = t_min; t <= t_max; ++t) {
      // D
      log_gamma[idx3(j, t, 1, p, l + 1)] =
        log_alpha(j, t) + logD +
        log_prob_emissions[idx3(j, t, 1, p, l + 1)] +
        log_beta(j + 1, t);
    }
    
    int t_max_m1 = std::min(t_max, l - 1);
    for (int t = t_min; t <= t_max_m1; ++t) {
      
      // K
      log_gamma[idx3(j, t, 0, p, l + 1)] =
        log_alpha(j, t) + logK +
        log_prob_emissions[idx3(j, t, 0, p, l + 1)] +
        log_beta(j + 1, t + 1);
      
      // S
      log_gamma[idx3(j, t, 2, p, l + 1)] =
        log_alpha(j, t) + logS +
        log_prob_emissions[idx3(j, t, 2, p, l + 1)] +
        log_beta(j + 1, t + 1);
    }
    
    int t_max_m2 = std::min(t_max, l - 2);
    for (int t = t_min; t <= t_max_m2; ++t) {
        
      // I
      log_gamma[idx3(j, t, 3, p, l + 1)] =
        log_alpha(j, t) + logI +
        log_prob_emissions[idx3(j, t, 3, p, l + 1)] +
        log_beta(j + 1, t + 2);
    }
  }
  
  return log_gamma;
}


// local emission at one site only
// state coding: 0=K, 1=D, 2=S, 3=I
inline double local_log_emission_one_site(const IntegerVector& x,
                                          int l,
                                          int t,
                                          int zcand,
                                          int state) {
  const double log_eps3 = std::log(1.0 / 3.0);
  const double log_quarter = std::log(0.25);
  
  if (state == 1) {
    // D
    return 0.0;
  }
  
  if (state == 0) {
    // K
    if (t <= l - 1 && x[t] == zcand) return 0.0;
    return R_NegInf;
  }
  
  if (state == 2) {
    // S, strict substitution
    if (t <= l - 1 && x[t] != zcand) return log_eps3;
    return R_NegInf;
  }
  
  // I: emit template symbol, then inserted symbol
  if (t <= l - 2 && x[t] == zcand) return log_quarter;
  return R_NegInf;
}


inline void update_log_alpha_next_row_one_site(const IntegerVector& x,
                                               int l,
                                               int j,
                                               int z_j,
                                               const NumericVector& log_theta,
                                               NumericMatrix& log_alpha,
                                               int B) {
  int row = j + 1;
  int t_min = std::max(0, row - B);
  int t_max = std::min(l, row + B);
  NumericVector vals(4);

  // Clear the row first. Some entries may have been valid under the old Z/alpha.
  for (int t = 0; t <= l; ++t) {
    log_alpha(row, t) = R_NegInf;
  }

  for (int t = t_min; t <= t_max; ++t) {
    // D: template deletion, read position unchanged
    vals[0] = log_theta[1] + log_alpha(j, t) +
      local_log_emission_one_site(x, l, t, z_j, 1);

    // K: match, consume one read symbol
    vals[1] = R_NegInf;
    if (t >= 1) {
      double le = local_log_emission_one_site(x, l, t - 1, z_j, 0);
      if (R_finite(le)) {
        vals[1] = log_theta[0] + log_alpha(j, t - 1) + le;
      }
    }

    // S: substitution, consume one read symbol
    vals[2] = R_NegInf;
    if (t >= 1) {
      double le = local_log_emission_one_site(x, l, t - 1, z_j, 2);
      if (R_finite(le)) {
        vals[2] = log_theta[2] + log_alpha(j, t - 1) + le;
      }
    }

    // I: template symbol plus inserted read symbol, consume two read symbols
    vals[3] = R_NegInf;
    if (t >= 2) {
      double le = local_log_emission_one_site(x, l, t - 2, z_j, 3);
      if (R_finite(le)) {
        vals[3] = log_theta[3] + log_alpha(j, t - 2) + le;
      }
    }

    log_alpha(row, t) = logsumexp_vec(vals);
  }
}

// exact local score for one read i, one site j, one candidate zcand
inline double score_candidate_one_read(const IntegerVector& x,
                                       int l,
                                       int j,
                                       int zcand,
                                       const NumericVector& log_theta,
                                       const NumericMatrix& log_alpha,
                                       const NumericMatrix& log_beta,
                                       int B) {
  NumericVector vals(4 * (l + 1), R_NegInf);
  int c = 0;

  // K
  int t_min = std::max(0, j - B);
  int t_max = std::min(l, j + B);
  
  for (int t = t_min; t <= t_max; ++t) {
    if (t <= l - 1) {
      double le = local_log_emission_one_site(x, l, t, zcand, 0);
      if (R_finite(le)) {
        vals[c] = log_alpha(j, t) + log_theta[0] + le + log_beta(j + 1, t + 1);
      }
    }
    c += 1;
  }
  
  // D
  for (int t = t_min; t <= t_max; ++t) {
    double le = local_log_emission_one_site(x, l, t, zcand, 1);
    vals[c] = log_alpha(j, t) + log_theta[1] + le + log_beta(j + 1, t);
    c += 1;
  }
  
  // S
  for (int t = t_min; t <= t_max; ++t) {
    if (t <= l - 1) {
      double le = local_log_emission_one_site(x, l, t, zcand, 2);
      if (R_finite(le)) {
        vals[c] = log_alpha(j, t) + log_theta[2] + le + log_beta(j + 1, t + 1);
      }
    }
    c += 1;
  }
  
  // I
  for (int t = t_min; t <= t_max; ++t) {
    if (t <= l - 2) {
      double le = local_log_emission_one_site(x, l, t, zcand, 3);
      if (R_finite(le)) {
        vals[c] = log_alpha(j, t) + log_theta[3] + le + log_beta(j + 1, t + 2);
      }
    }
    c += 1;
  }
  
  return logsumexp_vec(vals);
}

// [[Rcpp::export]]
IntegerVector update_Z_cpp_fast(List X,
                                IntegerVector Z,
                                NumericVector log_theta,
                                NumericVector tau,
                                IntegerVector B_vec) {
  int n = X.size();
  int p = Z.size();
  
  Z = clone(Z);
  
  // Cache reads and the dynamic-programming tables for the starting Z.
  // During the left-to-right sweep, beta can stay fixed: beta(j + 1, .) only
  // depends on the suffix sites j + 1:p, which have not yet been changed.
  // Alpha is updated one row at a time after each accepted site update.
  std::vector<IntegerVector> xs(n);
  std::vector<int> ls(n);
  std::vector<NumericMatrix> alphas(n);
  std::vector<NumericMatrix> betas(n);
  
  for (int i = 0; i <= n - 1; ++i) {
    xs[i] = as<IntegerVector>(X[i]);
    ls[i] = xs[i].size();
    NumericVector log_prob_emissions = get_log_prob_emissions_cpp(xs[i], Z, B_vec[i]);
    alphas[i] = get_log_alpha_cpp(log_theta, log_prob_emissions, B_vec[i]);
    betas[i]  = get_log_beta_cpp(log_theta, log_prob_emissions, B_vec[i]);
  }
  
  for (int j = 0; j <= p - 1; ++j) {
    NumericVector scores(4, 0.0);
    
    for (int h = 0; h <= 3; ++h) {
      int zcand = h + 1;
      
      for (int i = 0; i <= n - 1; ++i) {
        scores[h] += tau[i] * score_candidate_one_read(
          xs[i], ls[i], j, zcand, log_theta, alphas[i], betas[i], B_vec[i]
        );
      }
    }
    
    double log_denom = logsumexp_vec(scores);
    NumericVector probs(4);
    for (int h = 0; h <= 3; ++h) {
      probs[h] = std::exp(scores[h] - log_denom);
    }
    
    double u = R::runif(0.0, 1.0);
    double cum = 0.0;
    int chosen = 3;
    for (int h = 0; h <= 3; ++h) {
      cum += probs[h];
      if (u <= cum) {
        chosen = h;
        break;
      }
    }
    
    Z[j] = chosen + 1;
    
    // Keep the prefix cache exact for the next site. A full alpha/beta refresh is
    // unnecessary in a left-to-right sweep; only alpha(j + 1, .) is needed next.
    for (int i = 0; i <= n - 1; ++i) {
      update_log_alpha_next_row_one_site(
        xs[i], ls[i], j, Z[j], log_theta, alphas[i], B_vec[i]
      );
    }
  }
  
  return Z;
}
