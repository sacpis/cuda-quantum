/*******************************************************************************
 * Copyright (c) 2022 - 2024 NVIDIA Corporation & Affiliates.                  *
 * All rights reserved.                                                        *
 *                                                                             *
 * This source code and the accompanying materials are made available under    *
 * the terms of the Apache License 2.0 which accompanies this distribution.    *
 ******************************************************************************/

#include "cudaq/matrix.h"
#include "cudaq/operators.h"
#include <gtest/gtest.h>

/// NOTE: Not yet testing any of the matrix conversions. Just testing
/// the attributes of the output data type coming from the arithmetic.
/// These tests should be built upon to do actual numeric checks once
/// the implementations are complete.

TEST(ExpressionTester, checkPreBuiltElementaryOpsScalars) {

  auto function = [](std::map<std::string, std::complex<double>> parameters) {
    return parameters["value"];
  };

  // `elementary_operator + scalar_operator`
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::scalar_operator(1.0);

    auto sum = self + other;
    auto reverse = other + self;

    ASSERT_TRUE(sum.get_terms().size() == 2);
    ASSERT_TRUE(reverse.get_terms().size() == 2);
  }

   // `elementary_operator + scalar_operator`
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::scalar_operator(function);

    auto sum = self + other;
    auto reverse = other + self;

    ASSERT_TRUE(sum.get_terms().size() == 2);
    ASSERT_TRUE(reverse.get_terms().size() == 2);
  }

   // `elementary_operator - scalar_operator`
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::scalar_operator(1.0);

    auto sum = self - other;
    auto reverse = other - self;

    ASSERT_TRUE(sum.get_terms().size() == 2);
    ASSERT_TRUE(reverse.get_terms().size() == 2);
  }

   // `elementary_operator - scalar_operator`
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::scalar_operator(function);

    auto sum = self - other;
    auto reverse = other - self;

    ASSERT_TRUE(sum.get_terms().size() == 2);
    ASSERT_TRUE(reverse.get_terms().size() == 2);
  }

   // `elementary_operator * scalar_operator`
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::scalar_operator(1.0);

    auto product = self * other;
    auto reverse = other * self;

    ASSERT_TRUE(product.get_terms().size() == 2);
    ASSERT_TRUE(reverse.get_terms().size() == 2);
  }

   // `elementary_operator * scalar_operator`
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::scalar_operator(function);

    auto product = self * other;
    auto reverse = other * self;

    ASSERT_TRUE(product.get_terms().size() == 2);
    ASSERT_TRUE(reverse.get_terms().size() == 2);
  }

   // `elementary_operator / scalar_operator`
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::scalar_operator(1.0);

    auto product = self / other;
    /// FIXME:
    // auto reverse = other / self;

    ASSERT_TRUE(product.get_terms().size() == 2);
    // ASSERT_TRUE(reverse.get_terms().size() == 2);
  }

   // `elementary_operator / scalar_operator`
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::scalar_operator(function);

    auto product = self / other;
    /// FIXME:
    // auto reverse = other / self;

    ASSERT_TRUE(product.get_terms().size() == 2);
    // ASSERT_TRUE(reverse.get_terms().size() == 2);
  }
}

/// Prebuilt elementary ops against one another.
TEST(ExpressionTester, checkPreBuiltElementaryOpsSelf) {

  /// TODO: Check the output degrees attribute.

  // Addition, same DOF.
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::elementary_operator::create(0);

    auto sum = self + other;
    ASSERT_TRUE(sum.get_terms().size() == 2);
  }

  // Addition, different DOF's.
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::elementary_operator::create(1);

    auto sum = self + other;
    ASSERT_TRUE(sum.get_terms().size() == 2);
  }

  // Subtraction, same DOF.
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::elementary_operator::create(0);

    auto sum = self - other;
    ASSERT_TRUE(sum.get_terms().size() == 2);
  }

  // Subtraction, different DOF's.
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::elementary_operator::create(1);

    auto sum = self - other;
    ASSERT_TRUE(sum.get_terms().size() == 2);
  }

  // Multiplication, same DOF.
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::elementary_operator::create(0);

    auto product = self * other;
    ASSERT_TRUE(product.get_terms().size() == 2);
  }

  // Multiplication, different DOF's.
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    auto other = cudaq::elementary_operator::create(1);

    auto product = self * other;
    ASSERT_TRUE(product.get_terms().size() == 2);
  }
}

/// Testing arithmetic between elementary operators and operator
/// sums.
TEST(ExpressionTester, checkElementaryOpsAgainstOpSum) {

  /// Addition.
  {
    auto self = cudaq::elementary_operator::annihilate(0);
    /// Creating an arbitrary operator sum to work against.
    auto operator_sum = cudaq::elementary_operator::create(0) +
                        cudaq::elementary_operator::identity(1);

    auto got = self + operator_sum;
    std::cout << "term count of original op sum = "
              << operator_sum.get_terms().size() << "\n";
    std::cout << "term count = " << got.get_terms().size() << "\n";
    ASSERT_TRUE(got.get_terms().size() == 3);
  }
}
