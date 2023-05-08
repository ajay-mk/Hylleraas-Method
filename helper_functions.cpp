//
// Created by Ajay Melekamburath on 5/8/23.
//

#include "helper_functions.h"

// Function Definitions

double K_nlm(const int n, const double l, const double m, const double alpha,
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

int factorial(const int n){
  assert(n >= 0);
  if (n == 0)
    return 1;
  else
    return (n * factorial(n-1));
}

int binomial_coeff(const int n, const int k){
  std::vector<int> solutions(k);
  solutions[0] = n - k + 1;

  for (auto i = 1 ; i < k; ++i){
    solutions[i] = solutions[i-1] * (n - k + i + 1)/(i + 1);
  }
  return solutions[k-1];
}
