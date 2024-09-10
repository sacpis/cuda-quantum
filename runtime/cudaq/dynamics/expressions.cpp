#include "cudaq/expressions.h"
// #include "cudaq/definition.h"
#include "common/EigenDense.h"

#include <iostream>

namespace cudaq {

// implement everything without `_evaluate` first

// OperatorSum(std::vector<ProductOperator> &terms) {
//   m_terms = terms;
// }

// TODO:
// (1) Elementary Operators
// (2) Scalar Operators

// Definition::Definition(){};
// Definition::~Definition(){};

ElementaryOperator::ElementaryOperator(std::string operator_id,
                                       std::vector<int> degrees)
    : id(operator_id), degrees(degrees) {}

ElementaryOperator ElementaryOperator::identity(int degree) {
  std::string op_id = "identity";
  std::vector<int> degrees = {degree};
  auto op = ElementaryOperator(op_id, degrees);
  // NOTE: I don't actually think I need this if here because this
  // is a static method that creates a new ElementaryOperator (which
  // is what's being checked now) anyways.
  if (op.m_ops.find(op_id) == op.m_ops.end()) {
    auto func = [&](std::vector<int> none, std::vector<VariantArg> _none) {
      // Need to set the degree via the op itself because the
      // argument to the outer function goes out of scope when
      // the user invokes this later on via, e.g, `to_matrix()`.
      auto degree = op.degrees[0];
      auto mat = complex_matrix(degree, degree);
      // Build up the identity matrix.
      for (std::size_t i = 0; i < degree; i++) {
        mat(i, i) = 1.0 + 0.0 * 'j';
      }
      std::cout << "dumping the complex mat: \n";
      mat.dump();
      std::cout << "done\n\n";
      return mat;
    };
    op.define(op_id, degrees, func);
  }
  return op;
}

ElementaryOperator ElementaryOperator::zero(int degree) {
  std::string op_id = "zero";
  std::vector<int> degrees = {degree};
  auto op = ElementaryOperator(op_id, degrees);
  if (op.m_ops.find(op_id) == op.m_ops.end()) {
    auto func = [&](std::vector<int> none, std::vector<VariantArg> _none) {
      // Need to set the degree via the op itself because the
      // argument to the outer function goes out of scope when
      // the user invokes this later on via, e.g, `to_matrix()`.
      auto degree = op.degrees[0];
      auto mat = complex_matrix(degree, degree);
      mat.set_zero();
      std::cout << "dumping the complex mat: \n";
      mat.dump();
      std::cout << "\ndone\n";
      return mat;
    };
    op.define(op_id, degrees, func);
  }
  return op;
}

complex_matrix
ElementaryOperator::to_matrix(std::vector<int> degrees,
                              std::vector<VariantArg> parameters) {
  ReturnType result = m_ops[id].m_generator(degrees, parameters);

  if (std::holds_alternative<complex_matrix>(result)) {
    // Move the complex_matrix from the variant, which avoids copying
    return std::move(std::get<complex_matrix>(result));
  } else {
    // if it's a scalar, convert the scalar to a 1x1 matrix
    std::complex<double> scalar = std::get<std::complex<double>>(result);

    cudaq::complex_matrix scalar_matrix(1, 1);
    scalar_matrix(0, 0) = scalar;

    return scalar_matrix;
  }
}

/// @FIXME: The below function signature can be updated once
/// we support generalized function arguments.
/// @brief Constructor that just takes and returns a complex double value.
ScalarOperator::ScalarOperator(std::complex<double> value) {
  auto func = [&](std::vector<std::complex<double>> _none) { return value; };
  generator = func;
}

std::complex<double>
ScalarOperator::evaluate(std::vector<std::complex<double>> parameters) {
  return generator(parameters);
}

} // namespace cudaq