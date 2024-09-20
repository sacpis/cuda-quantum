/*******************************************************************************
 * Copyright (c) 2022 - 2024 NVIDIA Corporation & Affiliates.                  *
 * All rights reserved.                                                        *
 *                                                                             *
 * This source code and the accompanying materials are made available under    *
 * the terms of the Apache License 2.0 which accompanies this distribution.    *
 ******************************************************************************/

#include <cmath>
#include <gtest/gtest.h>

#include "cudaq/expressions.h"

cudaq::complex_matrix zero_matrix(int size) {
  auto mat = cudaq::complex_matrix(size, size);
  mat.set_zero();
  return mat;
}

cudaq::complex_matrix id_matrix(int size) {
  auto mat = cudaq::complex_matrix(size, size);
  for (int i = 0; i < size; i++)
    mat(i, i) = 1.0 + 0.0j;
  return mat;
}

cudaq::complex_matrix annihilate_matrix(int size) {
  auto mat = cudaq::complex_matrix(size, size);
  for (std::size_t i = 0; i + 1 < size; i++)
    mat(i, i + 1) = std::sqrt(static_cast<double>(i + 1)) + 0.0 * 'j';
  return mat;
}

cudaq::complex_matrix create_matrix(int size) {
  auto mat = cudaq::complex_matrix(size, size);
  for (std::size_t i = 0; i + 1 < size; i++)
    mat(i + 1, i) = std::sqrt(static_cast<double>(i + 1)) + 0.0 * 'j';
  return mat;
}

cudaq::complex_matrix position_matrix(int size) {
  auto mat = cudaq::complex_matrix(size, size);
  for (std::size_t i = 0; i + 1 < size; i++) {
    mat(i + 1, i) = 0.5 * std::sqrt(static_cast<double>(i + 1)) + 0.0 * 'j';
    mat(i, i + 1) = 0.5 * std::sqrt(static_cast<double>(i + 1)) + 0.0 * 'j';
  }
  return mat;
}

cudaq::complex_matrix momentum_matrix(int size) {
  auto mat = cudaq::complex_matrix(size, size);
  for (std::size_t i = 0; i + 1 < size; i++) {
    mat(i + 1, i) = (0.5j) * std::sqrt(static_cast<double>(i + 1)) + 0.0 * 'j';
    mat(i, i + 1) =
        -1. * (0.5j) * std::sqrt(static_cast<double>(i + 1)) + 0.0 * 'j';
  }
  return mat;
}

cudaq::complex_matrix number_matrix(int size) {
  auto mat = cudaq::complex_matrix(size, size);
  for (int i = 0; i < size; i++)
    mat(i, i) = static_cast<double>(i) + 0.0j;
  return mat;
}

cudaq::complex_matrix parity_matrix(int size) {
  auto mat = cudaq::complex_matrix(size, size);
  for (int i = 0; i < size; i++)
    mat(i, i) = std::pow(-1., static_cast<double>(i)) + 0.0j;
  return mat;
}

cudaq::complex_matrix displace_matrix(int size,
                                      std::complex<double> amplitude) {
  auto mat = cudaq::complex_matrix(size, size);
  for (std::size_t i = 0; i + 1 < size; i++) {
    mat(i + 1, i) =
        amplitude * std::sqrt(static_cast<double>(i + 1)) + 0.0 * 'j';
    mat(i, i + 1) = -1. * std::conj(amplitude) * (0.5 * 'j') *
                        std::sqrt(static_cast<double>(i + 1)) +
                    0.0 * 'j';
  }
  return mat.exp();
}

TEST(ExpressionTester, checkPreBuiltElementaryOps) {
  std::vector<int> sizes = {2, 3, 4, 5};
  // Identity operator.
  {
    for (int size : sizes) {
      auto id = cudaq::ElementaryOperator::identity(size);
      auto got_id = id.to_matrix({}, {});
      auto want_id = id_matrix(size);
      ASSERT_TRUE(want_id == got_id);
    }
  }

  // Zero operator.
  {
    for (int size : sizes) {
      auto zero = cudaq::ElementaryOperator::zero(size);
      auto got_zero = zero.to_matrix({}, {});
      auto want_zero = zero_matrix(size);
      ASSERT_TRUE(want_zero == got_zero);
    }
  }

  // Annihilation operator.
  {
    for (int size : sizes) {
      auto annihilate = cudaq::ElementaryOperator::annihilate(size);
      auto got_annihilate = annihilate.to_matrix({}, {});
      auto want_annihilate = annihilate_matrix(size);
      ASSERT_TRUE(want_annihilate == got_annihilate);
    }
  }

  // Creation operator.
  {
    for (int size : sizes) {
      auto create = cudaq::ElementaryOperator::create(size);
      auto got_create = create.to_matrix({}, {});
      auto want_create = create_matrix(size);
      ASSERT_TRUE(want_create == got_create);
    }
  }

  // Position operator.
  {
    for (int size : sizes) {
      auto position = cudaq::ElementaryOperator::position(size);
      auto got_position = position.to_matrix({}, {});
      auto want_position = position_matrix(size);
      ASSERT_TRUE(want_position == got_position);
    }
  }

  // Momentum operator.
  {
    for (int size : sizes) {
      auto momentum = cudaq::ElementaryOperator::momentum(size);
      auto got_momentum = momentum.to_matrix({}, {});
      auto want_momentum = momentum_matrix(size);
      ASSERT_TRUE(want_momentum == got_momentum);
    }
  }

  // Number operator.
  {
    for (int size : sizes) {
      auto number = cudaq::ElementaryOperator::number(size);
      auto got_number = number.to_matrix({}, {});
      auto want_number = number_matrix(size);
      ASSERT_TRUE(want_number == got_number);
    }
  }

  // Parity operator.
  {
    for (int size : sizes) {
      auto parity = cudaq::ElementaryOperator::parity(size);
      auto got_parity = parity.to_matrix({}, {});
      auto want_parity = parity_matrix(size);
      ASSERT_TRUE(want_parity == got_parity);
    }
  }

  // Displacement operator.
  {
    for (int size : sizes) {
      auto amplitude = 1.0 + 1.0j;
      auto displace = cudaq::ElementaryOperator::displace(size, amplitude);
      auto got_displace = displace.to_matrix({}, {});
      auto want_displace = displace_matrix(size, amplitude);
      ASSERT_TRUE(want_displace == got_displace);
    }
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

  // From concrete values.
  {
    auto operator_0 = cudaq::ScalarOperator(value_0);
    auto operator_1 = cudaq::ScalarOperator(value_1);
    auto operator_2 = cudaq::ScalarOperator(value_2);
    auto operator_3 = cudaq::ScalarOperator(value_3);

    auto got_value_0 = operator_0.evaluate({});
    auto got_value_1 = operator_1.evaluate({});
    auto got_value_2 = operator_2.evaluate({});
    auto got_value_3 = operator_3.evaluate({});

    EXPECT_NEAR(std::abs(value_0), std::abs(got_value_0), 1e-5);
    EXPECT_NEAR(std::abs(value_1), std::abs(got_value_1), 1e-5);
    EXPECT_NEAR(std::abs(value_2), std::abs(got_value_2), 1e-5);
    EXPECT_NEAR(std::abs(value_3), std::abs(got_value_3), 1e-5);
  }

  // From a lambda function.
  {
    auto function = [](std::vector<std::complex<double>> vec) {
      return vec[0];
    };

    auto operator_0 = cudaq::ScalarOperator(function);
    auto operator_1 = cudaq::ScalarOperator(function);
    auto operator_2 = cudaq::ScalarOperator(function);
    auto operator_3 = cudaq::ScalarOperator(function);

    auto got_value_0 = operator_0.evaluate({value_0});
    auto got_value_1 = operator_1.evaluate({value_1});
    auto got_value_2 = operator_2.evaluate({value_2});
    auto got_value_3 = operator_3.evaluate({value_3});

    EXPECT_NEAR(std::abs(value_0), std::abs(got_value_0), 1e-5);
    EXPECT_NEAR(std::abs(value_1), std::abs(got_value_1), 1e-5);
    EXPECT_NEAR(std::abs(value_2), std::abs(got_value_2), 1e-5);
    EXPECT_NEAR(std::abs(value_3), std::abs(got_value_3), 1e-5);
  }

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

    /// BROKEN:
    // + : Scalar operator from lambda.
    // {
    //   // Should try to use some local variables in this function so
    //   // I can test if it's working properly.
    //   auto function = [&](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   auto new_scalar_op = value_0 + scalar_op;

    //   auto scalar_op_1 = cudaq::ScalarOperator(function);
    //   auto reverse_order_op = scalar_op_1 + value_0;

    //   auto got_value = new_scalar_op.evaluate({value_1});
    //   auto got_value_1 = reverse_order_op.evaluate({value_1});

    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_1), 1e-5);
    //   EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 + value_0), 1e-5);
    // }

    // - : Constant scalar operator.
    {
      auto scalar_op = cudaq::ScalarOperator(value_1);
      auto new_scalar_op = value_3 - scalar_op;

      auto scalar_op_1 = cudaq::ScalarOperator(value_1);
      auto reverse_order_op = scalar_op_1 - value_3;

      auto got_value = new_scalar_op.evaluate({});
      auto got_value_1 = reverse_order_op.evaluate({});

      EXPECT_NEAR(std::abs(got_value), std::abs(value_3 - value_1), 1e-5);
      EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 - value_3), 1e-5);
    }

    /// BROKEN:
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

    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_2 - value_1), 1e-5);
    //   EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 - value_2), 1e-5);
    // }

    // * : Constant scalar operator.
    {
      auto scalar_op = cudaq::ScalarOperator(value_2);
      auto new_scalar_op = value_3 * scalar_op;

      auto scalar_op_1 = cudaq::ScalarOperator(value_2);
      auto reverse_order_op = scalar_op_1 * value_3;

      auto got_value = new_scalar_op.evaluate({});
      auto got_value_1 = reverse_order_op.evaluate({});

      EXPECT_NEAR(std::abs(got_value), std::abs(value_3 * value_2), 1e-5);
      EXPECT_NEAR(std::abs(got_value_1), std::abs(value_2 * value_3), 1e-5);
    }

    /// BROKEN:
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

    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_3 * value_2), 1e-5);
    //   EXPECT_NEAR(std::abs(got_value_1), std::abs(value_2 * value_3), 1e-5);
    // }

    // / : Constant scalar operator.
    {
      auto scalar_op = cudaq::ScalarOperator(value_2);
      auto new_scalar_op = value_3 / scalar_op;

      auto scalar_op_1 = cudaq::ScalarOperator(value_2);
      auto reverse_order_op = scalar_op_1 / value_3;

      auto got_value = new_scalar_op.evaluate({});
      auto got_value_1 = reverse_order_op.evaluate({});

      EXPECT_NEAR(std::abs(got_value), std::abs(value_3 / value_2), 1e-5);
      EXPECT_NEAR(std::abs(got_value_1), std::abs(value_2 / value_3), 1e-5);
    }

    /// BROKEN:
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

    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_3 / value_2), 1e-5);
    //   EXPECT_NEAR(std::abs(got_value_1), std::abs(value_2 / value_3), 1e-5);
    // }

    // += : Constant scalar operator.
    {
      auto scalar_op = cudaq::ScalarOperator(value_0);
      scalar_op += value_0;

      auto got_value = scalar_op.evaluate({});
      EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_0), 1e-5);
    }

    /// BROKEN:
    // // += : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   scalar_op += value_1;

    //   auto got_value = scalar_op.evaluate({value_0});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_1), 1e-5);
    // }

    // -= : Constant scalar operator.
    {
      auto scalar_op = cudaq::ScalarOperator(value_0);
      scalar_op -= value_0;

      auto got_value = scalar_op.evaluate({});
      EXPECT_NEAR(std::abs(got_value), std::abs(value_0 - value_0), 1e-5);
    }

    /// BROKEN:
    // // -= : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   scalar_op -= value_1;

    //   auto got_value = scalar_op.evaluate({value_0});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_0 - value_1), 1e-5);
    // }

    // *= : Constant scalar operator.
    {
      auto scalar_op = cudaq::ScalarOperator(value_2);
      scalar_op *= value_3;

      auto got_value = scalar_op.evaluate({});
      EXPECT_NEAR(std::abs(got_value), std::abs(value_2 * value_3), 1e-5);
    }

    /// BROKEN:
    // // *= : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   scalar_op *= value_3;

    //   auto got_value = scalar_op.evaluate({value_2});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_2 * value_3), 1e-5);
    // }

    // /= : Constant scalar operator.
    {
      auto scalar_op = cudaq::ScalarOperator(value_2);
      scalar_op /= value_3;

      auto got_value = scalar_op.evaluate({});
      EXPECT_NEAR(std::abs(got_value), std::abs(value_2 / value_3), 1e-5);
    }

    /// BROKEN:
    // // /= : Scalar operator from lambda.
    // {
    //   auto function = [](std::vector<std::complex<double>> vec) {
    //     return vec[0];
    //   };
    //   auto scalar_op = cudaq::ScalarOperator(function);
    //   scalar_op /= value_3;

    //   auto got_value = scalar_op.evaluate({value_2});
    //   EXPECT_NEAR(std::abs(got_value), std::abs(value_2 / value_3), 1e-5);
    // }
  }
}