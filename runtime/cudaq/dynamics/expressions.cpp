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
    const Eigen::MatrixXcd id_matrix =
        Eigen::MatrixXcd::Identity(degree, degree);
    std::cout << op_id << "\n" << id_matrix << "\n\n";
    std::vector<std::complex<double>> std_mat(
        id_matrix.data(), id_matrix.data() + id_matrix.size());
    for (auto elmt : std_mat) {
      std::cout << elmt << " ";
    }
    std::cout << "\n";
    /// FIXME: This could potentially be the most memory unsafe piece of
    /// code I've ever written... I will find out the hard way in testing.
    auto func = [&](std::vector<std::complex<double>> none) {
      auto mat = complex_matrix(std::move(std_mat.data()), degree,
                            degree);
      std::cout << "dumping the complex mat: ";
      mat.dump();
      std::cout << "\ndone\n";
      return mat;
    };
    // For testing purposes, deleting the members of the original std vec
    // to see what happens to the memory.
    std_mat = {6.,7.,8.};
    func({});
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
    const Eigen::MatrixXcd zero_matrix = Eigen::MatrixXcd::Zero(degree, degree);
    /// FIXME: Currently going from eigen -> std::vec -> cudaq::complex_matrix
    /// to avoid eigen bleeding to the user-interface. Ideally, I can get rid of
    /// the std::vec middleman
    std::vector<std::complex<double>> std_mat(
        zero_matrix.data(), zero_matrix.data() + zero_matrix.size());
    std::cout << op_id << "\n" << zero_matrix << "\n\n";
    for (auto elmt : std_mat) {
      std::cout << elmt << " ";
    }
    std::cout << "\n";
    /// FIXME: This could potentially be the most memory unsafe piece of
    /// code I've ever written... I will find out the hard way in testing.
    auto func = [&](std::vector<std::complex<double>> none) {
      auto mat = complex_matrix(std::move(std_mat.data()), degree,
                            degree);
      std::cout << "dumping the complex mat: ";
      mat.dump();
      std::cout << "\ndone\n";
      return mat;
    };
    // For testing purposes, deleting the members of the original std vec
    // to see what happens to the memory.
    std_mat = {1.,2.,3.};
    func({});
    // auto defn = Definition();
    // defn.create_definition(op_id, degrees, func);
    // op.m_ops[op_id] = defn;
  }
  return op;
}

} // namespace cudaq