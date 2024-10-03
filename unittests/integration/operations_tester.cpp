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

  // Keeping this fixed throughout.
  int degree_index = 0;
  // Dimension map.
  std::map<int, int> dimensions;

  // Identity operator.
  {
    for (int size : sizes) {
      dimensions[degree_index] = size;
      auto id = cudaq::ElementaryOperator::identity(degree_index);
      auto got_id = id.to_matrix(dimensions, {});
      auto want_id = id_matrix(size);
      ASSERT_TRUE(want_id == got_id);
    }
  }

  // Zero operator.
  {
    for (int size : sizes) {
      dimensions[degree_index] = size;
      auto zero = cudaq::ElementaryOperator::zero(degree_index);
      auto got_zero = zero.to_matrix(dimensions, {});
      auto want_zero = zero_matrix(size);
      ASSERT_TRUE(want_zero == got_zero);
    }
  }

  // Annihilation operator.
  {
    for (int size : sizes) {
      dimensions[degree_index] = size;
      auto annihilate = cudaq::ElementaryOperator::annihilate(degree_index);
      auto got_annihilate = annihilate.to_matrix(dimensions, {});
      auto want_annihilate = annihilate_matrix(size);
      ASSERT_TRUE(want_annihilate == got_annihilate);
    }
  }

  // Creation operator.
  {
    for (int size : sizes) {
      dimensions[degree_index] = size;
      auto create = cudaq::ElementaryOperator::create(degree_index);
      auto got_create = create.to_matrix(dimensions, {});
      auto want_create = create_matrix(size);
      ASSERT_TRUE(want_create == got_create);
    }
  }

  // Position operator.
  {
    for (int size : sizes) {
      dimensions[degree_index] = size;
      auto position = cudaq::ElementaryOperator::position(degree_index);
      auto got_position = position.to_matrix(dimensions, {});
      auto want_position = position_matrix(size);
      ASSERT_TRUE(want_position == got_position);
    }
  }

  // Momentum operator.
  {
    for (int size : sizes) {
      dimensions[degree_index] = size;
      auto momentum = cudaq::ElementaryOperator::momentum(degree_index);
      auto got_momentum = momentum.to_matrix(dimensions, {});
      auto want_momentum = momentum_matrix(size);
      ASSERT_TRUE(want_momentum == got_momentum);
    }
  }

  // Number operator.
  {
    for (int size : sizes) {
      dimensions[degree_index] = size;
      auto number = cudaq::ElementaryOperator::number(degree_index);
      auto got_number = number.to_matrix(dimensions, {});
      auto want_number = number_matrix(size);
      ASSERT_TRUE(want_number == got_number);
    }
  }

  // Parity operator.
  {
    for (int size : sizes) {
      dimensions[degree_index] = size;
      auto parity = cudaq::ElementaryOperator::parity(degree_index);
      auto got_parity = parity.to_matrix(dimensions, {});
      auto want_parity = parity_matrix(size);
      ASSERT_TRUE(want_parity == got_parity);
    }
  }

  // // Displacement operator.
  // {
  //   for (int size : sizes) {
  //     dimensions[degree_index] = size;
  //     auto amplitude = 1.0 + 1.0j;
  //     auto displace = cudaq::ElementaryOperator::displace(degree_index,
  //     amplitude); auto got_displace = displace.to_matrix(dimensions, {});
  //     auto want_displace = displace_matrix(size, amplitude);
  //     ASSERT_TRUE(want_displace == got_displace);
  //   }
  // }

  // TODO: Squeeze operator.
}

// TEST(ExpressionTester, checkCustomElementaryOps) {
//   // pass
// }

TEST(ExpressionTester, checkScalarOpsSimple) {

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
    auto function = [](std::map<std::string, std::complex<double>> parameters) {
      return parameters["value"];
    };

    std::map<std::string, std::complex<double>> parameter_map;

    auto operator_0 = cudaq::ScalarOperator(function);
    auto operator_1 = cudaq::ScalarOperator(function);
    auto operator_2 = cudaq::ScalarOperator(function);
    auto operator_3 = cudaq::ScalarOperator(function);

    parameter_map["value"] = value_0;
    auto got_value_0 = operator_0.evaluate(parameter_map);
    parameter_map["value"] = value_1;
    auto got_value_1 = operator_1.evaluate(parameter_map);
    parameter_map["value"] = value_2;
    auto got_value_2 = operator_2.evaluate(parameter_map);
    parameter_map["value"] = value_3;
    auto got_value_3 = operator_3.evaluate(parameter_map);

    EXPECT_NEAR(std::abs(value_0), std::abs(got_value_0), 1e-5);
    EXPECT_NEAR(std::abs(value_1), std::abs(got_value_1), 1e-5);
    EXPECT_NEAR(std::abs(value_2), std::abs(got_value_2), 1e-5);
    EXPECT_NEAR(std::abs(value_3), std::abs(got_value_3), 1e-5);
  }
}

TEST(ExpressionTester, checkScalarOpsArithmeticDoubles) {
  // Arithmetic overloads against complex doubles.
  std::complex<double> value_0 = 0.1 + 0.1;
  std::complex<double> value_1 = 0.1 + 1.0;
  std::complex<double> value_2 = 2.0 + 0.1;
  std::complex<double> value_3 = 2.0 + 1.0;

  auto local_variable = true;
  auto function = [&](std::map<std::string, std::complex<double>> parameters) {
    if (!local_variable)
      throw std::runtime_error("Local variable not detected.");
    return parameters["value"];
  };

  // + : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_0);

    auto new_scalar_op = value_1 + scalar_op;
    auto reverse_order_op = scalar_op + value_1;

    auto got_value = new_scalar_op.evaluate({});
    auto got_value_1 = reverse_order_op.evaluate({});
    auto want_value = value_1 + value_0;

    EXPECT_NEAR(std::abs(got_value), std::abs(want_value), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(want_value), 1e-5);

    // Checking composition of many scalar operators.
    auto third_op = new_scalar_op + reverse_order_op;
    auto got_value_third = third_op.evaluate({});
    EXPECT_NEAR(std::abs(got_value_third), std::abs(want_value + want_value),
                1e-5);
  }

  // + : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);

    auto new_scalar_op = value_0 + scalar_op;
    auto reverse_order_op = scalar_op + value_0;

    auto got_value = new_scalar_op.evaluate({{"value", value_1}});
    auto got_value_1 = reverse_order_op.evaluate({{"value", value_1}});

    EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_1), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 + value_0), 1e-5);

    // Checking composition of many scalar operators.
    auto third_op = new_scalar_op + reverse_order_op;
    auto got_value_third = third_op.evaluate({{"value", value_1}});
    auto want_value = value_0 + value_1 + value_1 + value_0;
    EXPECT_NEAR(std::abs(got_value_third), std::abs(want_value), 1e-5);
  }

  // - : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_1);

    auto new_scalar_op = value_3 - scalar_op;
    auto reverse_order_op = scalar_op - value_3;

    auto got_value = new_scalar_op.evaluate({});
    auto got_value_1 = reverse_order_op.evaluate({});

    EXPECT_NEAR(std::abs(got_value), std::abs(value_3 - value_1), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 - value_3), 1e-5);

    // Checking composition of many scalar operators.
    auto third_op = new_scalar_op - reverse_order_op;
    auto got_value_third = third_op.evaluate({});
    auto want_value = (value_3 - value_1) - (value_1 - value_3);
    EXPECT_NEAR(std::abs(got_value_third), std::abs(want_value), 1e-5);
  }

  // - : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);

    auto new_scalar_op = value_2 - scalar_op;
    auto reverse_order_op = scalar_op - value_2;

    auto got_value = new_scalar_op.evaluate({{"value", value_1}});
    auto got_value_1 = reverse_order_op.evaluate({{"value", value_1}});

    EXPECT_NEAR(std::abs(got_value), std::abs(value_2 - value_1), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 - value_2), 1e-5);

    // Checking composition of many scalar operators.
    auto third_op = new_scalar_op - reverse_order_op;
    auto got_value_third = third_op.evaluate({{"value", value_1}});
    auto want_value = (value_2 - value_1) - (value_1 - value_2);
    EXPECT_NEAR(std::abs(got_value_third), std::abs(want_value), 1e-5);
  }

  // * : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_2);

    auto new_scalar_op = value_3 * scalar_op;
    auto reverse_order_op = scalar_op * value_3;

    auto got_value = new_scalar_op.evaluate({});
    auto got_value_1 = reverse_order_op.evaluate({});

    EXPECT_NEAR(std::abs(got_value), std::abs(value_3 * value_2), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_2 * value_3), 1e-5);

    // Checking composition of many scalar operators.
    auto third_op = new_scalar_op * reverse_order_op;
    auto got_value_third = third_op.evaluate({});
    auto want_value = (value_3 * value_2) * (value_2 * value_3);
    EXPECT_NEAR(std::abs(got_value_third), std::abs(want_value), 1e-5);
  }

  // * : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);

    auto new_scalar_op = value_3 * scalar_op;
    auto reverse_order_op = scalar_op * value_3;

    auto got_value = new_scalar_op.evaluate({{"value", value_2}});
    auto got_value_1 = reverse_order_op.evaluate({{"value", value_2}});

    EXPECT_NEAR(std::abs(got_value), std::abs(value_3 * value_2), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_2 * value_3), 1e-5);

    // Checking composition of many scalar operators.
    auto third_op = new_scalar_op * reverse_order_op;
    auto got_value_third = third_op.evaluate({{"value", value_2}});
    auto want_value = (value_3 * value_2) * (value_2 * value_3);
    EXPECT_NEAR(std::abs(got_value_third), std::abs(want_value), 1e-5);
  }

  // / : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_2);

    auto new_scalar_op = value_3 / scalar_op;
    auto reverse_order_op = scalar_op / value_3;

    auto got_value = new_scalar_op.evaluate({});
    auto got_value_1 = reverse_order_op.evaluate({});

    EXPECT_NEAR(std::abs(got_value), std::abs(value_2 / value_3), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_3 / value_2), 1e-5);

    // Checking composition of many scalar operators.
    auto third_op = new_scalar_op / reverse_order_op;
    auto got_value_third = third_op.evaluate({});
    auto want_value = (value_2 / value_3) / (value_3 / value_2);
    EXPECT_NEAR(std::abs(got_value_third), std::abs(want_value), 1e-5);
  }

  // / : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);

    auto new_scalar_op = value_3 / scalar_op;
    auto reverse_order_op = scalar_op / value_3;

    auto got_value = new_scalar_op.evaluate({{"value", value_1}});
    auto got_value_1 = reverse_order_op.evaluate({{"value", value_1}});

    EXPECT_NEAR(std::abs(got_value), std::abs(value_1 / value_3), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_3 / value_1), 1e-5);

    // Checking composition of many scalar operators.
    auto third_op = new_scalar_op / reverse_order_op;
    auto got_value_third = third_op.evaluate({{"value", value_1}});
    auto want_value = (value_1 / value_3) / (value_3 / value_1);
    EXPECT_NEAR(std::abs(got_value_third), std::abs(want_value), 1e-5);
  }

  // += : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_0);
    scalar_op += value_0;

    auto got_value = scalar_op.evaluate({});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_0), 1e-5);
  }

  // += : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    scalar_op += value_1;

    auto got_value = scalar_op.evaluate({{"value", value_0}});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_1), 1e-5);
  }

  // -= : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_0);
    scalar_op -= value_0;

    auto got_value = scalar_op.evaluate({});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_0 - value_0), 1e-5);
  }

  // -= : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    scalar_op -= value_1;

    auto got_value = scalar_op.evaluate({{"value", value_0}});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_0 - value_1), 1e-5);
  }

  // *= : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_2);
    scalar_op *= value_3;

    auto got_value = scalar_op.evaluate({});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_2 * value_3), 1e-5);
  }

  // *= : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    scalar_op *= value_3;

    auto got_value = scalar_op.evaluate({{"value", value_2}});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_2 * value_3), 1e-5);
  }

  // /= : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_2);
    scalar_op /= value_3;

    auto got_value = scalar_op.evaluate({});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_2 / value_3), 1e-5);
  }

  // /= : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    scalar_op /= value_3;

    auto got_value = scalar_op.evaluate({{"value", value_2}});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_2 / value_3), 1e-5);
  }
}

TEST(ExpressionTester, checkScalarOpsArithmeticScalarOps) {
  // Arithmetic overloads against other scalar ops.
  std::complex<double> value_0 = 0.1 + 0.1;
  std::complex<double> value_1 = 0.1 + 1.0;
  std::complex<double> value_2 = 2.0 + 0.1;
  std::complex<double> value_3 = 2.0 + 1.0;

  auto local_variable = true;
  auto function = [&](std::map<std::string, std::complex<double>> parameters) {
    if (!local_variable)
      throw std::runtime_error("Local variable not detected.");
    return parameters["value"];
  };

  // I use another function here to make sure that local variables
  // that may be unique to each ScalarOp's generators are both kept
  // track of when we merge the generators.
  auto alternative_local_variable = true;
  auto alternative_function =
      [&](std::map<std::string, std::complex<double>> parameters) {
        if (!alternative_local_variable)
          throw std::runtime_error("Local variable not detected.");
        return parameters["other"];
      };

  // + : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_0);
    auto other_scalar_op = cudaq::ScalarOperator(value_1);

    auto new_scalar_op = other_scalar_op + scalar_op;
    auto reverse_order_op = scalar_op + other_scalar_op;

    auto got_value = new_scalar_op.evaluate({});
    auto got_value_1 = reverse_order_op.evaluate({});
    auto want_value = value_1 + value_0;

    EXPECT_NEAR(std::abs(got_value), std::abs(want_value), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(want_value), 1e-5);
  }

  // + : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    auto other_scalar_op = cudaq::ScalarOperator(alternative_function);

    auto new_scalar_op = other_scalar_op + scalar_op;
    auto reverse_order_op = scalar_op + other_scalar_op;

    std::map<std::string, std::complex<double>> parameter_map = {
        {"value", value_1}, {"other", value_0}};

    auto got_value = new_scalar_op.evaluate(parameter_map);
    auto got_value_1 = reverse_order_op.evaluate(parameter_map);

    EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_1), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 + value_0), 1e-5);
  }

  // - : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_2);
    auto other_scalar_op = cudaq::ScalarOperator(value_1);

    auto new_scalar_op = other_scalar_op - scalar_op;
    auto reverse_order_op = scalar_op - other_scalar_op;

    auto got_value = new_scalar_op.evaluate({});
    auto got_value_1 = reverse_order_op.evaluate({});
    auto want_value = value_1 - value_2;

    EXPECT_NEAR(std::abs(got_value), std::abs(want_value), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(want_value), 1e-5);
  }

  // - : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    auto other_scalar_op = cudaq::ScalarOperator(alternative_function);

    auto new_scalar_op = other_scalar_op - scalar_op;
    auto reverse_order_op = scalar_op - other_scalar_op;

    std::map<std::string, std::complex<double>> parameter_map = {
        {"value", value_1}, {"other", value_3}};

    auto got_value = new_scalar_op.evaluate(parameter_map);
    auto got_value_1 = reverse_order_op.evaluate(parameter_map);

    EXPECT_NEAR(std::abs(got_value), std::abs(value_3 - value_1), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 - value_3), 1e-5);
  }

  // * : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_2);
    auto other_scalar_op = cudaq::ScalarOperator(value_3);

    auto new_scalar_op = other_scalar_op * scalar_op;
    auto reverse_order_op = scalar_op * other_scalar_op;

    auto got_value = new_scalar_op.evaluate({});
    auto got_value_1 = reverse_order_op.evaluate({});
    auto want_value = value_3 * value_2;
    auto reverse_want_value = value_2 * value_3;

    EXPECT_NEAR(std::abs(got_value), std::abs(want_value), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(reverse_want_value), 1e-5);
  }

  // * : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    auto other_scalar_op = cudaq::ScalarOperator(alternative_function);

    auto new_scalar_op = other_scalar_op * scalar_op;
    auto reverse_order_op = scalar_op * other_scalar_op;

    std::map<std::string, std::complex<double>> parameter_map = {
        {"value", value_1}, {"other", value_3}};

    auto got_value = new_scalar_op.evaluate(parameter_map);
    auto got_value_1 = reverse_order_op.evaluate(parameter_map);

    EXPECT_NEAR(std::abs(got_value), std::abs(value_3 * value_1), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_1 * value_3), 1e-5);
  }

  // / : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_0);
    auto other_scalar_op = cudaq::ScalarOperator(value_2);

    auto new_scalar_op = other_scalar_op / scalar_op;
    auto reverse_order_op = scalar_op / other_scalar_op;

    auto got_value = new_scalar_op.evaluate({});
    auto got_value_1 = reverse_order_op.evaluate({});
    auto want_value = value_2 / value_0;
    auto reverse_want_value = value_0 / value_2;

    EXPECT_NEAR(std::abs(got_value), std::abs(want_value), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(reverse_want_value), 1e-5);
  }

  // / : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    auto other_scalar_op = cudaq::ScalarOperator(alternative_function);

    auto new_scalar_op = other_scalar_op / scalar_op;
    auto reverse_order_op = scalar_op / other_scalar_op;

    std::map<std::string, std::complex<double>> parameter_map = {
        {"value", value_0}, {"other", value_3}};

    auto got_value = new_scalar_op.evaluate(parameter_map);
    auto got_value_1 = reverse_order_op.evaluate(parameter_map);

    EXPECT_NEAR(std::abs(got_value), std::abs(value_3 / value_0), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_0 / value_3), 1e-5);
  }

  // += : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_0);
    auto other = cudaq::ScalarOperator(value_0);
    scalar_op += other;

    auto got_value = scalar_op.evaluate({});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_0), 1e-5);
  }

  // += : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    auto other = cudaq::ScalarOperator(value_1);
    scalar_op += other;

    auto scalar_op_1 = cudaq::ScalarOperator(function);
    auto other_function = cudaq::ScalarOperator(alternative_function);
    scalar_op_1 += other_function;

    auto got_value = scalar_op.evaluate({{"value", value_0}});
    auto got_value_1 =
        scalar_op_1.evaluate({{"value", value_0}, {"other", value_1}});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_0 + value_1), 1e-5);
    EXPECT_NEAR(std::abs(got_value_1), std::abs(value_0 + value_1), 1e-5);
  }

  // -= : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_0);
    scalar_op -= value_0;

    auto got_value = scalar_op.evaluate({});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_0 - value_0), 1e-5);
  }

  // -= : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    scalar_op -= value_1;

    auto got_value = scalar_op.evaluate({{"value", value_0}});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_0 - value_1), 1e-5);
  }

  // *= : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_2);
    scalar_op *= value_3;

    auto got_value = scalar_op.evaluate({});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_2 * value_3), 1e-5);
  }

  // *= : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    scalar_op *= value_3;

    auto got_value = scalar_op.evaluate({{"value", value_2}});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_2 * value_3), 1e-5);
  }

  // /= : Constant scalar operator.
  {
    auto scalar_op = cudaq::ScalarOperator(value_2);
    scalar_op /= value_3;

    auto got_value = scalar_op.evaluate({});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_2 / value_3), 1e-5);
  }

  // /= : Scalar operator from lambda.
  {
    auto scalar_op = cudaq::ScalarOperator(function);
    scalar_op /= value_3;

    auto got_value = scalar_op.evaluate({{"value", value_2}});
    EXPECT_NEAR(std::abs(got_value), std::abs(value_2 / value_3), 1e-5);
  }
}

TEST(ExpressionTester, checkScalarAgainstElementary) {

  std::complex<double> value_0 = 0.1 + 0.1;
  std::complex<double> value_1 = 0.1 + 1.0;
  std::complex<double> value_2 = 2.0 + 0.1;
  std::complex<double> value_3 = 2.0 + 1.0;

  auto local_variable = true;
  auto function = [&](std::map<std::string, std::complex<double>> parameters) {
    if (!local_variable)
      throw std::runtime_error("Local variable not detected.");
    return parameters["value"];
  };

  // Identity against constant.
  {
    auto id_op = cudaq::ElementaryOperator::identity(0);
    auto scalar_op = cudaq::ScalarOperator(value_0);

    // auto addition = scalar_op + id_op;
    // auto subtraction = scalar_op - id_op;
    // auto multiplication = scalar_op * id_op;
  }

  // // Identity against constant from lambda.
  // {
  //   auto id_op = cudaq::ElementaryOperator::identity(0);
  //   auto scalar_op = cudaq::ScalarOperator(function)
  // }
}