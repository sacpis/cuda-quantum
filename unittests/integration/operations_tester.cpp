/*******************************************************************************
 * Copyright (c) 2022 - 2024 NVIDIA Corporation & Affiliates.                  *
 * All rights reserved.                                                        *
 *                                                                             *
 * This source code and the accompanying materials are made available under    *
 * the terms of the Apache License 2.0 which accompanies this distribution.    *
 ******************************************************************************/

#include <gtest/gtest.h>

#include "cudaq/expressions.h"

TEST(ExpressionTester, checkPreBuiltElementaryOps) {
  auto op_id = cudaq::ElementaryOperator::identity(2);
  auto op_zero = cudaq::ElementaryOperator::zero(2);

  op_id.to_matrix({},{});
  op_zero.to_matrix({},{});

  // Trying different sizes and just visually confirming for now.
  cudaq::ElementaryOperator::identity(3);
  cudaq::ElementaryOperator::zero(3);

  cudaq::ElementaryOperator::identity(4);
  cudaq::ElementaryOperator::zero(4);
}

TEST(ExpressionTester, checkCustomElementaryOps) {
  // pass
}