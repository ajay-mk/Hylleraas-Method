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
  std::cout << "Hylleraas method for correlation energy\n" << std::endl;

  BasisFn test_basis = {{{0, 0, 0}, {1.6875, 1.6875, 0.0}}}; // see typedef

  std::cout << "Overlap using single element basis: "
            << compute_overlap(test_basis) << std::endl;
  std::cout << "Hamiltonian using single element basis: "
            << compute_hamiltonian(test_basis, 2) << std::endl;

  return 0;
}