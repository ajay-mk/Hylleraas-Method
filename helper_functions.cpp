//
// Created by Ajay Melekamburath on 5/8/23.
//

#include "helper_functions.h"

// Function Definitions

int factorial(const int n) {
  assert(n >= 0);
  if (n == 0)
    return 1;
  else
    return (n * factorial(n - 1));
}

int binomial_coeff(const int n, const int k) {
  std::vector<int> solutions(k);
  solutions[0] = n - k + 1;

  for (auto i = 1; i < k; ++i) {
    solutions[i] = solutions[i - 1] * (n - k + i + 1) / (i + 1);
  }
  return solutions[k - 1];
}

double K_nlm(const int n, const int l, const int m, const double alpha,
             const double beta, const double gamma) {

  auto pre_fac = 16 * pow(M_PI, 2) * factorial(n + 1) * factorial(l + 1) *
                 factorial(m + 1);
  double value = 0.0;

  for (auto a = 0; a <= n + 1; ++a) {
    for (auto b = 0; b <= l + 1; ++b) {
      for (auto c = 0; c <= m + 1; ++c) {
        value += binomial_coeff(l + 1 - b + a, a) *
                 binomial_coeff(m + 1 - c + b, b) *
                 binomial_coeff(n + 1 - c + c, c) /
                 (pow(alpha + beta, l - b + a + 2) *
                  pow(alpha + gamma, n - a + c + 2) *
                  pow(beta + gamma, m - c + b + 2));
      }
    }
  }
  const auto result = pre_fac * value;
  return result;
}

double eval_S(const int ni, const int li, const int mi, const int nj,
              const int lj, const int mj, const double alpha, const double beta,
              const double gamma) {
  // bra and ket have the same exponential parameters alpha, beta and gamma
  return K_nlm(ni + nj, li + lj, mi + mj, alpha, beta, gamma);
}

double eval_nuc_attr(const double Z, const int ni, const int li, const int mi,
                     const int nj, const int lj, const int mj,
                     const double alpha, const double beta,
                     const double gamma) {
  // bra and ket have the same exponential parameters alpha, beta and gamma
  // term 1
  const auto t1 = K_nlm(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma);
  // term 2
  const auto t2 = K_nlm(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma);

  return (-Z * (t1 + t2));
}

double eval_elec_rep(const int ni, const int li, const int mi, const int nj,
                     const int lj, const int mj, const double alpha,
                     const double beta, const double gamma) {
  // bra and ket have the same exponential parameters alpha, beta and gamma
  return K_nlm(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma);
}

double eval_T(int ni, int li, int mi, int nj, int lj, int mj, double alpha,
              double beta, double gamma) {
  ;
}
