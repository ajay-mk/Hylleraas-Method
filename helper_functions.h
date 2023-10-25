//
// Created by Ajay Melekamburath on 5/8/23.
//

#ifndef HYLLERAAS_HELPER_FUNCTIONS_H
#define HYLLERAAS_HELPER_FUNCTIONS_H

#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include <thread>

#ifdef USE_MULTIPRECISION
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_MULTIPRECISION
using namespace boost::multiprecision;
using float_type = number<boost::multiprecision::cpp_dec_float<20>>;
using int_type = cpp_int;
#else
using float_type = double;
using int_type = long long int;
#endif

typedef std::vector<std::pair<std::vector<int>, std::vector<double>>> BasisFn;
// example: {{{0,0,0},{1.6875, 1.6875, 0.0}}}
// general: {{{n,l,m},{alpha/2, beta/2, gamma/2}}}

typedef Eigen::Matrix<float_type, Eigen::Dynamic, Eigen::Dynamic,
                      Eigen::RowMajor>
    Matrix;

// Structs
struct hylleraas_results {
  Matrix evals;
  Matrix evecs;
  Matrix H;
  Matrix S;
};

// Function Declarations
/// Computes compute_factorial of a number
template <typename T> T compute_factorial(int n);

/// Computes binomial coefficient for C(n, k)
template <typename T> T compute_binomial_coeff(int n, int k);

/// Evaluates the K integral [Equation 33]
float_type K_nlm(int n, int l, int m, double alpha, double beta, double gamma);

/// Evaluates overlap S_ij [Equation 39]
float_type eval_S(int ni, int li, int mi, int nj, int lj, int mj, double alpha,
                  double beta, double gamma);

/// Evaluates nuclear attraction V_ij(ne) [Equation 41]
float_type eval_V_nuc(double Z, int ni, int li, int mi, int nj, int lj, int mj,
                      double alpha, double beta, double gamma);

/// Evaluates electron-electron repulsion V_ij(ee) [Equation 43]
float_type eval_V_elec(int ni, int li, int mi, int nj, int lj, int mj,
                       double alpha, double beta, double gamma);

/// Evaluates kinetic energy T_ij [Equation 45]
float_type eval_T(int ni, int li, int mi, int nj, int lj, int mj, double alpha,
                  double beta, double gamma);

/// evaluates terms in T expression, returns zero if pre-factor is zero
/// otherwise returns pre-factor * integral
float_type eval_T_terms(double pre_fac, int n, int l, int m, double alpha,
                        double beta, double gamma);

/// Computes overlap matrix S
Matrix compute_overlap(const BasisFn &basis);

/// Computes Hamiltonian matrix H
Matrix compute_hamiltonian(const BasisFn &basis, double Z);

/// Solves the secular equations
hylleraas_results do_hylleraas_simple(const BasisFn &basis, double Z);

/// Constructs a basis based on parameters alpha, beta and N
BasisFn construct_basis(int N, double alpha, double gamma);

/// Solves the secular equation, more general
/// constructs basis based on N, alpha and gamma
hylleraas_results do_hylleraas(int N, double alpha, double gamma, double Z);

#endif // HYLLERAAS_HELPER_FUNCTIONS_H
