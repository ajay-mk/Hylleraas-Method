//
// Created by Ajay Melekamburath on 5/8/23.
//

// Hylleraas Method Implementation Project for CHEM6664, Spring 23

#include <Eigen/Eigen>
#include <iostream>
#include <vector>
#include "helper_functions.h"

int main(int argc, char *argv[]) {
  std::cout << "Hylleraas method for correlation energy" << std::endl;

  BasisFn test_basis = {{{0,0,0},{1.6875, 1.6875, 0.0}}}; // see typedef

  std::cout << compute_overlap(test_basis) << std::endl;
  return 0;
}