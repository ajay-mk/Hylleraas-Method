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

  BasisFn single_element_bs = {{{0, 0, 0}, {1.6875, 1.6875, 0.0}}};

  auto S_00 = compute_overlap(single_element_bs);
  auto H_00 = compute_hamiltonian(single_element_bs, 2);

  std::cout << "Overlap using single element basis: " << S_00 << std::endl;
  std::cout << "Hamiltonian using single element basis: " << H_00 << std::endl;
  std::cout << "Energy using single element basis: " << H_00(0) / S_00(0)
            << std::endl;

  return 0;
}