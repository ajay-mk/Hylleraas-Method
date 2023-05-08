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

  if (k > 0)
    return 0;

  if (k == 0 || k == n)
    return 1;

  return binomial_coeff(n - 1, k - 1) + binomial_coeff(n - 1, k);
}

double K_nlm(const int n, const int l, const int m, const double alpha,
             const double beta, const double gamma) {
  assert(n >= 0 && l >= 0 && m >= 0);

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

double eval_V_nuc(const double Z, const int ni, const int li, const int mi,
                     const int nj, const int lj, const int mj,
                     const double alpha, const double beta,
                     const double gamma) {
  assert(Z > 0);
  // bra and ket have the same exponential parameters alpha, beta and gamma
  // term 1
  const auto t1 = K_nlm(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma);
  // term 2
  const auto t2 = K_nlm(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma);

  return (-Z * (t1 + t2));
}

double eval_V_elec(const int ni, const int li, const int mi, const int nj,
                     const int lj, const int mj, const double alpha,
                     const double beta, const double gamma) {
  // bra and ket have the same exponential parameters alpha, beta and gamma
  return K_nlm(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma);
}

double eval_T(int ni, int li, int mi, int nj, int lj, int mj, double alpha,
              double beta, double gamma) {
  double value = 0.0;
  // first line
  value += (-1.0 / 8.0) * (pow(alpha, 2) + pow(beta, 2) + 2 * pow(gamma, 2)) *
           eval_S(ni, li, mi, nj, lj, mj, alpha, beta, gamma);

  // second line
  value += (nj * alpha / 2.0) *
           K_nlm(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma);
  value += (lj * beta / 2.0) *
           K_nlm(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma);
  value +=
      (mj * gamma) * K_nlm(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma);

  // third line
  value += (-nj * (nj - 1) / 2.0) *
           K_nlm(ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma);
  value += (-lj * (lj - 1) / 2.0) *
           K_nlm(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma);

  // fourth line
  value += (mj * (mj - 1)) *
           K_nlm(ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma);

  // fifth line
  value +=
      (alpha / 2.0) * K_nlm(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma);
  value +=
      (beta / 2.0) * K_nlm(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma);
  value += gamma * K_nlm(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma);

  // sixth amd seventh lines
  value += -nj * K_nlm(ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma);
  value += -lj * K_nlm(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma);
  value += (-mj * 2) * K_nlm(ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma);

  // eighth line
  value += (-alpha * gamma / 8.0) *
           (K_nlm(ni + nj - 1, li + lj, mi + mj + 1, alpha, beta, gamma) +
            K_nlm(ni + nj + 1, li + lj, mi + mj - 1, alpha, beta, gamma) -
            K_nlm(ni + nj - 1, li + lj + 2, mi + mj - 1, alpha, beta, gamma));

  // ninth line
  value += (-beta * gamma / 8.0) *
           (K_nlm(ni + nj, li + lj - 1, mi + mj + 1, alpha, beta, gamma) +
            K_nlm(ni + nj, li + lj + 1, mi + mj - 1, alpha, beta, gamma) -
            K_nlm(ni + nj + 2, li + lj - 1, mi + mj - 1, alpha, beta, gamma));

  // tenth line
  value += (nj * gamma / 4.0) *
           (K_nlm(ni + nj - 2, li + lj, mi + mj + 1, alpha, beta, gamma) +
            K_nlm(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma) -
            K_nlm(ni + nj - 2, lj + lj + 2, mi + mj - 1, alpha, beta, gamma));

  // eleventh line
  value += (mj * alpha / 4.0) *
           (K_nlm(ni + nj - 1, li + lj, mi + mj, alpha, beta, gamma) +
            K_nlm(ni + nj + 1, li + lj, mi + mj - 2, alpha, beta, gamma) -
            K_nlm(ni + nj - 1, li + lj + 2, mi + mj - 2, alpha, beta, gamma));

  // twelfth line
  value += (-nj * mj / 2.0) *
           (K_nlm(ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma) +
            K_nlm(ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma) -
            K_nlm(ni + nj - 2, li + lj + 2, mi + mj - 2, alpha, beta, gamma));

  // thirteenth line
  value += (lj * gamma / 4.0) *
           (K_nlm(ni + nj, li + lj - 2, mi + mj + 1, alpha, beta, gamma) +
            K_nlm(ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma) -
            K_nlm(ni + nj + 2, li + lj - 2, mi + mj - 1, alpha, beta, gamma));

  // fourteenth line
  value += (mj * beta / 4.0) *
           (K_nlm(ni + nj, li + lj - 1, mi + mj, alpha, beta, gamma) +
            K_nlm(ni + nj, li + lj + 1, mi + mj - 2, alpha, beta, gamma) -
            K_nlm(ni + nj + 2, li + lj - 1, mi + mj - 1, alpha, beta, gamma));

  // last (fifteenth) line
  value += (-lj * mj / 2.0) *
           (K_nlm(ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma) +
            K_nlm(ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma) -
            K_nlm(ni + nj + 2, li + lj - 2, mi + mj - 2, alpha, beta, gamma));

  return value;
}

Eigen::MatrixXd compute_overlap(const BasisFn &basis) {
  using namespace Eigen;
  const int n = basis.size();
  auto result = MatrixXd(n, n);
  result.fill(0.0);

  for (auto i = 0; i < n; ++i) {
    for (auto j = 0; j < n; ++j) {
      // fetch quantum numbers from basis
      auto qn_i = basis[i].first;
      auto qn_j = basis[j].first;
      // fetch exponents from basis
      auto exps = basis[i].second; // same exponents for bra and ket, so reuse
      result(i, j) = eval_S(qn_i[0], qn_i[1], qn_i[2], qn_j[0], qn_j[1],
                            qn_j[2], 2 * exps[0], 2 * exps[1], 2 * exps[2]);
    }
  }
  return result;
}
