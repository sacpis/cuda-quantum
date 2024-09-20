/*******************************************************************************
 * Copyright (c) 2022 - 2024 NVIDIA Corporation & Affiliates.                  *
 * All rights reserved.                                                        *
 *                                                                             *
 * This source code and the accompanying materials are made available under    *
 * the terms of the Apache License 2.0 which accompanies this distribution.    *
 ******************************************************************************/

#include <gtest/gtest.h>

#include "cudaq/expressions.h"

cudaq::complex_matrix zero(int size) {
  auto mat = cudaq::complex_matrix(size, size);
  mat.set_zero();
  return mat;
}

cudaq::complex_matrix id(int size) {
  auto mat = cudaq::complex_matrix(size, size);
  for (int i = 0; i < size; i++)
    mat(i, i) = 1.0 + 0.0j;
  return mat;
}

TEST(ExpressionTester, checkPreBuiltElementaryOps) {
  // Check 2x2 variations.
  {
    auto id_2 = cudaq::ElementaryOperator::identity(2);
    auto zero_2 = cudaq::ElementaryOperator::zero(2);

    auto want_id_2 = id_2.to_matrix({}, {});
    auto want_zero_2 = zero_2.to_matrix({}, {});
    auto got_id_2 = id(2);
    auto got_zero_2 = zero(2);

    ASSERT_TRUE(want_id_2 == got_id_2);
    ASSERT_TRUE(want_zero_2 == got_zero_2);
  }

  // Check 3x3 variations.
  {
    auto id_3 = cudaq::ElementaryOperator::identity(3);
    auto zero_3 = cudaq::ElementaryOperator::zero(3);

    auto want_id_3 = id_3.to_matrix({}, {});
    auto want_zero_3 = zero_3.to_matrix({}, {});
    auto got_id_3 = id(3);
    auto got_zero_3 = zero(3);

    ASSERT_TRUE(want_id_3 == got_id_3);
    ASSERT_TRUE(want_zero_3 == got_zero_3);
  }

  // Check 4x4 variations.
  {
    auto id_4 = cudaq::ElementaryOperator::identity(4);
    auto zero_4 = cudaq::ElementaryOperator::zero(4);

    auto want_id_4 = id_4.to_matrix({}, {});
    auto want_zero_4 = zero_4.to_matrix({}, {});
    auto got_id_4 = id(4);
    auto got_zero_4 = zero(4);

    ASSERT_TRUE(want_id_4 == got_id_4);
    ASSERT_TRUE(want_zero_4 == got_zero_4);
  }
}

TEST(ExpressionTester, checkCustomElementaryOps) {
  // pass
}

TEST(ExpressionTester, checkScalarOps) {

  std::complex<double> value_0 = 0.1 + 0.1;
  std::complex<double> value_1 = 0.1 + 1.0;
  std::complex<double> value_2 = 2.0 + 0.1;
  std::complex<double> value_3 = 2.0 + 1.0;

  // // From concrete values.
  // {
  //   auto operator_0 = cudaq::ScalarOperator(value_0);
  //   auto operator_1 = cudaq::ScalarOperator(value_1);
  //   auto operator_2 = cudaq::ScalarOperator(value_2);
  //   auto operator_3 = cudaq::ScalarOperator(value_3);

  //   auto got_value_0 = operator_0.evaluate({});
  //   auto got_value_1 = operator_1.evaluate({});
  //   auto got_value_2 = operator_2.evaluate({});
  //   auto got_value_3 = operator_3.evaluate({});

  //   EXPECT_NEAR(std::abs(value_0), std::abs(got_value_0), 1e-5);
  //   EXPECT_NEAR(std::abs(value_1), std::abs(got_value_1), 1e-5);
  //   EXPECT_NEAR(std::abs(value_2), std::abs(got_value_2), 1e-5);
  //   EXPECT_NEAR(std::abs(value_3), std::abs(got_value_3), 1e-5);
  // }

  // // From a lambda function.
  // {
  //   auto function = [](std::vector<std::complex<double>> vec) {
  //     return vec[0];
  //   };

  //   auto operator_0 = cudaq::ScalarOperator(function);
  //   auto operator_1 = cudaq::ScalarOperator(function);
  //   auto operator_2 = cudaq::ScalarOperator(function);
  //   auto operator_3 = cudaq::ScalarOperator(function);

  //   auto got_value_0 = operator_0.evaluate({value_0});
  //   auto got_value_1 = operator_1.evaluate({value_1});
  //   auto got_value_2 = operator_2.evaluate({value_2});
  //   auto got_value_3 = operator_3.evaluate({value_3});

  //   EXPECT_NEAR(std::abs(value_0), std::abs(got_value_0), 1e-5);
  //   EXPECT_NEAR(std::abs(value_1), std::abs(got_value_1), 1e-5);
  //   EXPECT_NEAR(std::abs(value_2), std::abs(got_value_2), 1e-5);
  //   EXPECT_NEAR(std::abs(value_3), std::abs(got_value_3), 1e-5);
  // }

  // Arithmetic overloads against complex doubles.
  {
    // + : Constant scalar operator.
    {
      auto scalar_op = cudaq::ScalarOperator(value_0);
      auto new_scalar_op = value_1 + scalar_op;

      auto scalar_op_1 = cudaq::ScalarOperator(value_0);
      auto reverse_order_op = scalar_op_1 + value_1;
      auto got_value = new_scalar_op.evaluate({});
      auto got_value_1 = reverse_order_op.evaluate({});
      auto want_value = value_1 + value_0;

      EXPECT_NEAR(std::abs(got_value), std::abs(want_value), 1e-5);
      EXPECT_NEAR(std::abs(got_value_1), std::abs(want_value), 1e-5);
    }

    // + : Scalar operator from lambda.
    {
      // Should try to use some local variables in this function so
      // I can test if it's working properly.
      auto function = [&](std::vector<std::complex<double>> vec) {
        return vec[0];
      };
      auto scalar_op = cudaq::ScalarOperator(function);
      auto new_scalar_op = value_0 + scalar_op;

      auto scalar_op_1 = cudaq::ScalarOperator(function);
      auto reverse_order_op = scalar_op_1 + value_0;

      auto got_value = new_scalar_op.evaluate({value_1});
      auto got_value_1 = reverse_order_op.evaluate({value_1});

      EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_1),
                  1e-5);
      EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 + value_0),
                  1e-5);
    }

    // // - : Constant scalar operator.
    // {
    //   auto scalar_op = cudaq::ScalarOperator(value_1);
    //   auto new_scalar_op = value_3 - scalar_op;

    //   auto scalar_op_1 = cudaq::ScalarOperator(value_1);
    //   auto reverse_order_op = scalar_op_1 - value_3;

    //   auto got_value = new_scalar_op.evaluate({});
    //   auto got_value_1 = reverse_order_op.evaluate({});

    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_3 - value_1),
    //               1e-5);
    //   EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 - value_3),
    //               1e-5);
    // }

    // // - : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   auto new_scalar_op = value_2 - scalar_op;

    //   auto scalar_op_1 = cudaq::ScalarOperator(function);
    //   auto reverse_order_op = scalar_op_1 - value_2;

    //   auto got_value = new_scalar_op.evaluate({value_1});
    //   auto got_value_1 = reverse_order_op.evaluate({value_1});

    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_2 - value_1),
    //               1e-5);
    //   EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 - value_2),
    //               1e-5);
    // }

    // // * : Constant scalar operator.
    // {
    //   auto scalar_op = cudaq::ScalarOperator(value_2);
    //   // std::cout << "\nscalar op value = " << scalar_op.evaluate({}) <<
    //   "\n\n";
    //   // std::cout << "\nvalue_3 = " << value_3 << "\n\n";
    //   auto new_scalar_op = value_3 * scalar_op;

    //   // auto scalar_op_1 = cudaq::ScalarOperator(value_2);
    //   // auto reverse_order_op = scalar_op_1 * value_3;

    //   auto got_value = new_scalar_op.evaluate({});
    //   // auto got_value_1 = reverse_order_op.evaluate({});

    //   // std::cout << "got value = " << got_value << "\n";
    //   // std::cout << "want value = " << (value_3 * value_2) << "\n";

    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_3 * value_2),
    //               1e-5);
    //   // EXPECT_NEAR(std::abs(got_value_1), std::abs(value_2 * value_3),
    //   //             1e-5);
    // }

    // // * : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   auto new_scalar_op = value_3 * scalar_op;

    //   auto scalar_op_1 = cudaq::ScalarOperator(function);
    //   auto reverse_order_op = scalar_op_1 * value_3;

    //   auto got_value = new_scalar_op.evaluate({value_1});
    //   auto got_value_1 = reverse_order_op.evaluate({value_1});

    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_3 * value_2),
    //               1e-5);
    //   EXPECT_NEAR(std::abs(got_value_1), std::abs(value_2 * value_3),
    //               1e-5);
    // }

    // // / : Constant scalar operator.
    // {
    //   auto scalar_op = cudaq::ScalarOperator(value_2);
    //   auto new_scalar_op = value_3 / scalar_op;

    //   auto scalar_op_1 = cudaq::ScalarOperator(value_2);
    //   auto reverse_order_op = scalar_op_1 / value_3;

    //   auto got_value = new_scalar_op.evaluate({});
    //   auto got_value_1 = reverse_order_op.evaluate({});

    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_3 / value_2),
    //               1e-5);
    //   EXPECT_NEAR(std::abs(got_value_1), std::abs(value_2 / value_3),
    //               1e-5);
    // }

    // // / : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   auto new_scalar_op = value_3 / scalar_op;

    //   auto scalar_op_1 = cudaq::ScalarOperator(function);
    //   auto reverse_order_op = scalar_op_1 / value_3;

    //   auto got_value = new_scalar_op.evaluate({value_1});
    //   auto got_value_1 = reverse_order_op.evaluate({value_1});

    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_3 / value_2),
    //               1e-5);
    //   EXPECT_NEAR(std::abs(got_value_1), std::abs(value_2 / value_3),
    //               1e-5);
    // }

    // // += : Constant scalar operator.
    // {
    //   auto scalar_op = cudaq::ScalarOperator(value_0);
    //   scalar_op += value_0;

    //   auto got_value = scalar_op.evaluate({});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_0),
    //               1e-5);
    // }

    // // += : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   scalar_op += value_1;

    //   auto got_value = scalar_op.evaluate({value_0});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_1),
    //               1e-5);
    // }

    // // -= : Constant scalar operator.
    // {
    //   auto scalar_op = cudaq::ScalarOperator(value_0);
    //   scalar_op -= value_0;

    //   auto got_value = scalar_op.evaluate({});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_0 - value_0),
    //               1e-5);
    // }

    // // -= : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   scalar_op -= value_1;

    //   auto got_value = scalar_op.evaluate({value_0});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_0 - value_1),
    //               1e-5);
    // }

    // // *= : Constant scalar operator.
    // {
    //   auto scalar_op = cudaq::ScalarOperator(value_2);
    //   scalar_op *= value_3;

    //   auto got_value = scalar_op.evaluate({});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_2 * value_3),
    //               1e-5);
    // }

    // // *= : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   scalar_op *= value_3;

    //   auto got_value = scalar_op.evaluate({value_2});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_2 * value_3),
    //               1e-5);
    // }

    // // /= : Constant scalar operator.
    // {
    //   auto scalar_op = cudaq::ScalarOperator(value_2);
    //   scalar_op /= value_3;

    //   auto got_value = scalar_op.evaluate({});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_2 / value_3),
    //               1e-5);
    // }

    // // /= : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   scalar_op /= value_3;

    //   auto got_value = scalar_op.evaluate({value_2});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_2 / value_3),
    //               1e-5);
    // }
  }
}