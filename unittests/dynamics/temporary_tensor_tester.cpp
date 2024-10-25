/*******************************************************************************
 * Copyright (c) 2022 - 2024 NVIDIA Corporation & Affiliates.                  *
 * All rights reserved.                                                        *
 *                                                                             *
 * This source code and the accompanying materials are made available under    *
 * the terms of the Apache License 2.0 which accompanies this distribution.    *
 ******************************************************************************/

#include "cudaq/utils/tensor.h"
#include "cudaq/operators.h"
#include <gtest/gtest.h>


// cudaq::tensor _id_matrix(int size) {
//   auto mat = cudaq::tensor(size, size);
//   for (int i = 0; i < size; i++)
//     mat(i, i) = 1.0 + 0.0j;
//   return mat;
// }

// cudaq::tensor _annihilate_matrix(int size) {
//   auto mat = cudaq::tensor(size, size);
//   for (std::size_t i = 0; i + 1 < size; i++)
//     mat(i, i + 1) = std::sqrt(static_cast<double>(i + 1)) + 0.0 * 'j';
//   return mat;
// }

// cudaq::tensor _create_matrix(int size) {
//   auto mat = cudaq::tensor(size, size);
//   for (std::size_t i = 0; i + 1 < size; i++)
//     mat(i + 1, i) = std::sqrt(static_cast<double>(i + 1)) + 0.0 * 'j';
//   return mat;
// }


/// REMOVEME: This file is temporary and is just working out operator specific
/// routines using the new tensor class.

TEST(ExpressionTester, checkTensor) {

  {
    auto mat = cudaq::tensor({2, 2});
    cudaq::tensor t({1, 2, 1});
    mat.dump();
  }

}