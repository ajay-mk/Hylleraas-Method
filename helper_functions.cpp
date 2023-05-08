//
// Created by Ajay Melekamburath on 5/8/23.
//

#include "helper_functions.h"

// Function Definitions

int factorial(const int n) {
  assert(n >= 0);
  return (n == 1 || n == 0) ? 1 : n * factorial(n - 1);
}

int binomial_coeff(const int n, const int k) {
  if (k > n)
    return 0;
  if (k == 0 || k == n)
    return 1;
  return binomial_coeff(n - 1, k - 1) + binomial_coeff(n - 1, k);
}

double K_nlm(const int n, const int l, const int m, const double alpha,
             const double beta, const double gamma) {
  // in this project, alpha = beta and gamma = 0, but trying to keep it general

  auto pre_fac = 16 * pow(M_PI, 2) * factorial(n + 1) * factorial(l + 1) *
                 factorial(m + 1);

  double value = 0.0;
  for (std::size_t a = 0; a != n + 2; ++a) {
    for (std::size_t b = 0; b != l + 2; ++b) {
      for (std::size_t c = 0; c != m + 2; ++c) {
        value += binomial_coeff(l + 1 - b + a, a) *
                 binomial_coeff(m + 1 - c + b, b) *
                 binomial_coeff(n + 1 - a + c, c) /
                 (pow(alpha + beta, l - b + a + 2) *
                  pow(alpha + gamma, n - a + c + 2) *
                  pow(beta + gamma, m - c + b + 2));
      }
    }
  }
  auto result = pre_fac * value;
  return result;
}

double eval_S(const int ni, const int li, const int mi, const int nj,
              const int lj, const int mj, const double alpha, const double beta,
              const double gamma) {
  // bra and ket have the same exponential parameters alpha, beta and gamma
  return K_nlm(ni + nj, li + lj, mi + mj, alpha, beta, gamma);
}

double eval_V_nuc(const double Z, const int ni, const int li, const int mi,
                  const int nj, const int lj, const int mj, const double alpha,
                  const double beta, const double gamma) {
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

double eval_T_terms(const double pre_fac, const int n, const int l, const int m,
                    const double alpha, const double beta, const double gamma) {
  if (pre_fac == 0)
    return 0;
  return pre_fac * K_nlm(n, l, m, alpha, beta, gamma);
}

double eval_T(int ni, int li, int mi, int nj, int lj, int mj, double alpha,
              double beta, double gamma) {
  double value = 0.0;
  // first line
  value += (-1.0 / 8.0) * (pow(alpha, 2) + pow(beta, 2) + 2 * pow(gamma, 2)) *
           eval_S(ni, li, mi, nj, lj, mj, alpha, beta, gamma);

  // second line
  value += eval_T_terms(nj * alpha / 2.0, ni + nj - 1, li + lj, mi + mj, alpha,
                        beta, gamma);
  value += eval_T_terms(lj * beta / 2.0, ni + nj, li + lj - 1, mi + mj, alpha,
                        beta, gamma);
  value += eval_T_terms(mj * gamma, ni + nj, li + lj, mi + mj - 1, alpha, beta,
                        gamma);

  // third line
  value += eval_T_terms((-nj * (nj - 1) / 2.0), ni + nj - 2, li + lj, mi + mj,
                        alpha, beta, gamma);
  value += eval_T_terms(-lj * (lj - 1) / 2, ni + nj, li + lj - 2, mi + mj,
                        alpha, beta, gamma);

  // fourth line
  value += eval_T_terms(mj * (mj - 1), ni + nj, li + lj, mi + mj - 2, alpha,
                        beta, gamma);

  // fifth line
  value += eval_T_terms(alpha / 2.0, ni + nj - 1, li + lj, mi + mj, alpha, beta,
                        gamma);
  value += eval_T_terms(beta / 2.0, ni + nj, li + lj - 1, mi + mj, alpha, beta,
                        gamma);
  value +=
      eval_T_terms(gamma, ni + nj, li + lj, mi + mj - 1, alpha, beta, gamma);

  // sixth amd seventh lines
  value += eval_T_terms(-nj, ni + nj - 2, li + lj, mi + mj, alpha, beta, gamma);
  value += eval_T_terms(-lj, ni + nj, li + lj - 2, mi + mj, alpha, beta, gamma);
  value +=
      eval_T_terms(-mj * 2, ni + nj, li + lj, mi + mj - 2, alpha, beta, gamma);

  // eighth line
  value += eval_T_terms(-alpha * gamma / 8.0, ni + nj - 1, li + lj, mi + mj + 1,
                        alpha, beta, gamma) +
           eval_T_terms(-alpha * gamma / 8.0, ni + nj + 1, li + lj, mi + mj - 1,
                        alpha, beta, gamma) -
           eval_T_terms(-alpha * gamma / 8.0, ni + nj - 1, li + lj + 2,
                        mi + mj - 1, alpha, beta, gamma);

  // ninth line
  value += eval_T_terms(-beta * gamma / 8.0, ni + nj, li + lj - 1, mi + mj + 1,
                        alpha, beta, gamma) +
           eval_T_terms(-beta * gamma / 8.0, ni + nj, li + lj + 1, mi + mj - 1,
                        alpha, beta, gamma) -
           eval_T_terms(-beta * gamma / 8.0, ni + nj + 2, li + lj - 1,
                        mi + mj - 1, alpha, beta, gamma);

  // tenth line
  value += eval_T_terms(nj * gamma / 4.0, ni + nj - 2, li + lj, mi + mj + 1,
                        alpha, beta, gamma) +
           eval_T_terms(nj * gamma / 4.0, ni + nj, li + lj, mi + mj - 1, alpha,
                        beta, gamma) -
           eval_T_terms(nj * gamma / 4.0, ni + nj - 2, lj + lj + 2, mi + mj - 1,
                        alpha, beta, gamma);

  // eleventh line
  value += eval_T_terms(mj * alpha / 4.0, ni + nj - 1, li + lj, mi + mj, alpha,
                        beta, gamma) +
           eval_T_terms(mj * alpha / 4.0, ni + nj + 1, li + lj, mi + mj - 2,
                        alpha, beta, gamma) -
           eval_T_terms(mj * alpha / 4.0, ni + nj - 1, li + lj + 2, mi + mj - 2,
                        alpha, beta, gamma);

  // twelfth line
  value += eval_T_terms(-nj * mj / 2.0, ni + nj - 2, li + lj, mi + mj, alpha,
                        beta, gamma) +
           eval_T_terms(-nj * mj / 2.0, ni + nj, li + lj, mi + mj - 2, alpha,
                        beta, gamma) -
           eval_T_terms(-nj * mj / 2.0, ni + nj - 2, li + lj + 2, mi + mj - 2,
                        alpha, beta, gamma);

  // thirteenth line
  value += eval_T_terms(lj * gamma / 4.0, ni + nj, li + lj - 2, mi + mj + 1,
                        alpha, beta, gamma) +
           eval_T_terms(lj * gamma / 4.0, ni + nj, li + lj, mi + mj - 1, alpha,
                        beta, gamma) -
           eval_T_terms(lj * gamma / 4.0, ni + nj + 2, li + lj - 2, mi + mj - 1,
                        alpha, beta, gamma);

  // fourteenth line
  value += eval_T_terms(mj * beta / 4.0, ni + nj, li + lj - 1, mi + mj, alpha,
                        beta, gamma) +
           eval_T_terms(mj * beta / 4.0, ni + nj, li + lj + 1, mi + mj - 2,
                        alpha, beta, gamma) -
           eval_T_terms(mj * beta / 4.0, ni + nj + 2, li + lj - 1, mi + mj - 1,
                        alpha, beta, gamma);

  // last (fifteenth) line
  value += eval_T_terms(-lj * mj / 2.0, ni + nj, li + lj - 2, mi + mj, alpha,
                        beta, gamma) +
           eval_T_terms(-lj * mj / 2.0, ni + nj, li + lj, mi + mj - 2, alpha,
                        beta, gamma) -
           eval_T_terms(-lj * mj / 2.0, ni + nj + 2, li + lj - 2, mi + mj - 2,
                        alpha, beta, gamma);

  return value;
}

Eigen::MatrixXd compute_overlap(const BasisFn &basis) {
  const int n = basis.size();
  auto result = Eigen::MatrixXd(n, n);
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

Eigen::MatrixXd compute_hamiltonian(const BasisFn &basis, const double Z) {
  const int n = basis.size();
  auto result = Eigen::MatrixXd(n, n);
  result.fill(0.0);

  for (auto i = 0; i < n; ++i) {
    for (auto j = 0; j < n; ++j) {
      // fetch quantum numbers from basis
      auto qn_i = basis[i].first;
      auto qn_j = basis[j].first;
      // fetch exponents from basis
      auto exps = basis[i].second; // same exponents for bra and ket, so reuse

      // H_ij = T_ij + V_ij(ne) + V_ij(ee)
      auto t = eval_T(qn_i[0], qn_i[1], qn_i[2], qn_j[0], qn_j[1], qn_j[2],
                      2 * exps[0], 2 * exps[1], 2 * exps[2]);
      auto v_ne = eval_V_nuc(Z, qn_i[0], qn_i[1], qn_i[2], qn_j[0], qn_j[1],
                             qn_j[2], 2 * exps[0], 2 * exps[1], 2 * exps[2]);
      auto v_ee = eval_V_elec(qn_i[0], qn_i[1], qn_i[2], qn_j[0], qn_j[1],
                              qn_j[2], 2 * exps[0], 2 * exps[1], 2 * exps[2]);
      result(i, j) = t + v_ee + v_ne;
    }
  }
  return result;
}
