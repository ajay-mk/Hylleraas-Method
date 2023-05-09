//
// Created by Ajay Melekamburath on 5/8/23.
//

// Hylleraas Method Implementation Project for CHEM6664, Spring 23

#include "helper_functions.h"
#include <iomanip>
#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << std::setprecision(12);
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
            << results_3elem.evals.transpose()
            << std::endl; // transpose only for pretty printing

  std::cout << "\nGround state energy using three element basis: "
            << results_3elem.evals(0) << " Eh" << std::endl;

  // Calculation of ground state energy of Helium atom
  std::cout << "\nCalculation of ground state energy of Helium atom\n";
  // vary N from 0 to 15

  std::cout << "N"
            << "\t"
            << "Energy (Eh)" << std::endl;
  for (auto n = 0; n <= 10; ++n) {
    // construct basis with alpha = beta = 1.8, gamma = 0;
    // set Z = 2 for Helium atom

    std::cout << "Basis Size: " << construct_basis(n, 1.8, 0).size()
              << std::endl;
    auto result = do_hylleraas(n, 1.8, 0.0, 2.0);

    std::cout << n << "\t" << result.evals(0) << std::endl;
  }

  return 0;
}