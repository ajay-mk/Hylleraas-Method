//
// Created by Ajay Melekamburath on 5/8/23.
//

#include "helper_functions.h"

double K_nlm(const int n, const double l, const double m, const double alpha,
             const double beta, const double gamma) {
  double value = 0.0;

  return value;
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
