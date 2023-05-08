//
// Created by Ajay Melekamburath on 5/8/23.
//

#include <iostream>
#include <cmath>
#include <vector>
#include "Eigen/Eigen"

#ifndef HYLLERAAS_HELPER_FUNCTIONS_H
#define HYLLERAAS_HELPER_FUNCTIONS_H

// Function Declarations
/// Computes factorial of a number
int factorial(const int n);

/// Computes binomial coefficient for nCk
int binomial_coeff(const int n, const int k);


/// Equation 33
double K_nlm(const int n, const double l, const double m, const double alpha,
             const double beta, const double gamma);

#endif // HYLLERAAS_HELPER_FUNCTIONS_H
