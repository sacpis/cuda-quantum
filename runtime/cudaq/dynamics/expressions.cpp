#include "cudaq/expressions.h"

#include "common/EigenDense.h"

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
  if (m_ops.find(op_id) == m_ops.end()) {
    auto id_matrix = Eigen::MatrixXcd::Identity(degree, degree);
    auto func = [&]() { return id_matrix; };
    auto defn = Definition();
    defn.create_definition(op_id, degrees, func);
    m_ops[op_id] = defn;
  }
  return ElementaryOperator(op_id, degrees);
}

ElementaryOperator ElementaryOperator::zero(int degree) {
  std::string op_id = "zero";
  std::vector<int> degrees = {degree};
  if (m_ops.find(op_id) == m_ops.end()) {
    auto zero_matrix = Eigen::MatrixXcd::Zero(degree, degree);
    auto func = [&]() { return zero_matrix; };
    auto defn = Definition();
    defn.create_definition(op_id, degrees, func);
    m_ops[op_id] = defn;
  }
  return ElementaryOperator(op_id, degrees);
}

} // namespace cudaq