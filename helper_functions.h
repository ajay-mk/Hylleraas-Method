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

// Structs
struct hylleraas_results {
  Eigen::MatrixXd evals;
  Eigen::MatrixXd evecs;
  Eigen::MatrixXd H;
  Eigen::MatrixXd S;
};

// Function Declarations
/// Computes factorial of a number
int factorial(int n);

/// Computes binomial coefficient for C(n, k)
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
double eval_V_elec(int ni, int li, int mi, int nj, int lj, int mj, double alpha,
                   double beta, double gamma);

/// Evaluates kinetic energy T_ij [Equation 45]
double eval_T(int ni, int li, int mi, int nj, int lj, int mj, double alpha,
              double beta, double gamma);

/// evaluates terms in T expression, returns zero if pre-factor is zero
/// otherwise returns pre-factor * integral
double eval_T_terms(double pre_fac, int n, int l, int m, double alpha,
                    double beta, double gamma);

/// Computes overlap matrix S
Eigen::MatrixXd compute_overlap(const BasisFn &basis);

/// Computes Hamiltonian matrix H
Eigen::MatrixXd compute_hamiltonian(const BasisFn &basis, double Z);

/// Solves the secular equations
hylleraas_results do_hylleraas_simple(const BasisFn &basis, const double Z);

/// Constructs a basis based on parameters alpha, beta and N
BasisFn construct_basis(int N, double alpha, double gamma);

/// Solves the secular equation, more general
/// constructs basis based on N, alpha and gamma
hylleraas_results do_hylleraas(double Z, int N, double alpha, double gamma);

#endif // HYLLERAAS_HELPER_FUNCTIONS_H
