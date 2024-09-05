#include "cudaq/expressions.h"

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

Definition::Definition() {}
Definition::~Definition() {}

ElementaryOperator::ElementaryOperator(std::string operator_id,
                                       std::vector<int> degrees)
    : id(operator_id), degrees(degrees) {}

ElementaryOperator ElementaryOperator::identity(int degree) {
  std::string op_id = "identity";
  std::vector<int> degrees = {degree};
  auto op = ElementaryOperator(op_id, degrees);
  if (op.m_ops.find(op_id) == op.m_ops.end()) {
    auto func = [&](std::vector<int> none, std::vector<std::complex<double>> _none) {
      auto mat = complex_matrix(degree, degree);
      // Build up the identity matrix.
      for (std::size_t i=0; i<degree; i++) {
        mat(i,i) = 1.0+0.0j;
      }
      std::cout << "dumping the complex mat: \n";
      mat.dump();
      std::cout << "\ndone\n";
      return mat;
    };
      // auto defn = Definition();
      // defn.create_definition(op_id, degrees, func);
      // op.m_ops[op_id] = defn;
  }
  return op;
}

ElementaryOperator ElementaryOperator::zero(int degree) {
  std::string op_id = "zero";
  std::vector<int> degrees = {degree};
  auto op = ElementaryOperator(op_id, degrees);
  if (op.m_ops.find(op_id) == op.m_ops.end()) {
    auto func = [&](std::vector<int> none, std::vector<std::complex<double>> _none) {
      auto mat = complex_matrix(degree, degree);
      mat.set_zero();
      return mat;
    };
    // auto defn = Definition();
    // defn.create_definition(op_id, degrees, func);
    // op.m_ops[op_id] = defn;
  }
  return op;
}

} // namespace cudaq