//
// Created by Ajay Melekamburath on 5/8/23.
//

// Hylleraas Method Implementation Project for CHEM6664, Spring 23

#include "helper_functions.h"
#include <Eigen/Eigen>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
  std::cout << std::setprecision(10);
  std::cout << "Hylleraas Method\n" << std::endl;

  BasisFn bs_single = {{{0, 0, 0}, {1.6875, 1.6875, 0.0}}};

  auto S_00 = compute_overlap(bs_single);
  auto H_00 = compute_hamiltonian(bs_single, 2);

  std::cout << "Overlap using single element basis: " << S_00 << std::endl;
  std::cout << "Hamiltonian using single element basis: " << H_00 << std::endl;

  std::cout << "Energy using single element basis: " << H_00(0) / S_00(0)
            << std::endl;

  BasisFn bs_3elem = {{{0, 0, 0}, {1.8, 1.8, 0.0}},
                      {{1, 1, 0}, {1.8, 1.8, 0.0}},
                      {{0, 0, 1}, {1.8, 1.8, 0.0}}};

  auto results_3elem = do_hylleraas_simple(bs_3elem, 2);
  std::cout << "\nOverlap using three element basis: \n"
            << results_3elem.S << std::endl;
  std::cout << "\nHamiltonian using three element basis: \n"
            << results_3elem.H << std::endl;

  std::cout << "\nOrbital energies using three element basis: \n"
            << results_3elem.evals << std::endl;

  return 0;
}