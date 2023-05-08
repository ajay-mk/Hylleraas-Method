//
// Created by Ajay Melekamburath on 5/8/23.
//

#include "Eigen/Eigen"
#include <cmath>
#include <iostream>
#include <vector>

#ifndef HYLLERAAS_HELPER_FUNCTIONS_H
#define HYLLERAAS_HELPER_FUNCTIONS_H

// TypeDefs
typedef std::vector<std::pair<std::vector<int>, std::vector<double>>> BasisFn;
// example: {{{0,0,0},{1.6875, 1.6875, 0.0}}}
// general: {{{n,l,m},{alpha/2, beta/2, gamma/2}}}

// Function Declarations
/// Computes factorial of a number
int factorial(int n);

/// Computes binomial coefficient for nCk
int binomial_coeff(int n, int k);

/// Evaluates the K integral [Equation 33]
double K_nlm(int n, int l, int m, double alpha, double beta, double gamma);

/// Evaluates overlap S_ij [Equation 39]
double eval_S(int ni, int li, int mi, int nj, int lj, int mj, double alpha,
              double beta, double gamma);

/// Evaluates nuclear attaction V_ij(ne) [Equation 41]
double eval_V_nuc(double Z, int ni, int li, int mi, int nj, int lj, int mj,
                     double alpha, double beta, double gamma);

/// Evaluates electron-electron repulsion V_ij(ee) [Equation 43]
double eval_V_elec(int ni, int li, int mi, int nj, int lj, int mj,
                     double alpha, double beta, double gamma);

/// Evaluates kinetic energy T_ij [Equation 45]
double eval_T(int ni, int li, int mi, int nj, int lj, int mj, double alpha,
              double beta, double gamma);

/// Evaluates overlap matrix S
Eigen::MatrixXd compute_overlap(const BasisFn &basis);

#endif // HYLLERAAS_HELPER_FUNCTIONS_H
