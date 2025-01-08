/****************************************************************-*- C++ -*-****
 * Copyright (c) 2022 - 2024 NVIDIA Corporation & Affiliates.                  *
 * All rights reserved.                                                        *
 *                                                                             *
 * This source code and the accompanying materials are made available under    *
 * the terms of the Apache License 2.0 which accompanies this distribution.    *
 ******************************************************************************/
#pragma once

#include "details/tensor_impl.h"
#include "extension_point.h"
#include "type_traits.h"

#include <list>

namespace cudaq {

// Forward declarations.
template <typename T>
class tensor;
template <typename T>
tensor<T> operator*(const tensor<T> &, const tensor<T> &);
template <typename T>
tensor<T> operator+(const tensor<T> &, const tensor<T> &);

/// @brief A tensor class implementing the PIMPL idiom.
///
/// The flattened data is stored in row-major layout, where the strides grow
/// from right to left (as in a multi-dimensional C array).
///
/// There are three memory models a client can select from. In all of these, the
/// size of the data must be at least as large as the shape argument.
///
///   - `copy()`: This will make a copy of the data that is passed. The tensor
///     owns the copy of the data.
///   - `take()`: The tensor object is passed a `unique_ptr` to the data and
///     will take ownership of the data. The client's `unique_ptr` will be
///     invalidated with this operation.
///   - `borrow()`: The tensor object is passed a raw pointer to data. The
///     tensor object does \e not have ownership of the data. This means that
///     client code is responsible for ensuring that the pointer remains valid
///     for the entire lifetime of the tensor object.
///
/// Not all of these models will be fully coherent and functional under all
/// scenarios, so wrapping implementation layers may default to `copy()`.
template <typename Scalar = std::complex<double>>
class tensor {
private:
  std::shared_ptr<details::tensor_impl<Scalar>> pimpl;

public:
  /// @brief Type alias for the scalar type used in the tensor
  using scalar_type = typename details::tensor_impl<Scalar>::scalar_type;
  static constexpr auto ScalarAsString = type_to_string<Scalar>();

  /// @brief Construct an empty tensor
  tensor()
      : pimpl(std::shared_ptr<details::tensor_impl<Scalar>>(
            details::tensor_impl<Scalar>::get(
                std::string("xtensor") + std::string(ScalarAsString), {})
                .release())) {}

  /// @brief Construct a tensor with the given shape
  /// @param shape The shape of the tensor
  tensor(const std::vector<std::size_t> &shape)
      : pimpl(std::shared_ptr<details::tensor_impl<Scalar>>(
            details::tensor_impl<Scalar>::get(
                std::string("xtensor") + std::string(ScalarAsString), shape)
                .release())) {}

  /// @brief Construct a tensor with the given data and shape
  /// @param data Pointer to the tensor data
  /// @param shape The shape of the tensor
  tensor(const scalar_type *data, const std::vector<std::size_t> &shape)
      : pimpl(std::shared_ptr<details::tensor_impl<Scalar>>(
            details::tensor_impl<Scalar>::get(std::string("xtensor") +
                                                  std::string(ScalarAsString),
                                              data, shape)
                .release())) {}

  /// Construct a `tensor` from a `tensor_impl` using move semantics.
  tensor(details::tensor_impl<Scalar> *impl) { pimpl.swap(impl); }

  /// @brief Get the rank of the tensor
  /// @return The rank of the tensor
  std::size_t rank() const { return pimpl->rank(); }

  /// @brief Get the total size of the tensor
  /// @return The total number of elements in the tensor
  std::size_t size() const { return pimpl->size(); }

  /// @brief Get the shape of the tensor
  /// @return A vector containing the dimensions of the tensor
  std::vector<std::size_t> shape() const { return pimpl->shape(); }

  /// @brief Access a mutable element of the tensor
  /// @param indices The indices of the element to access
  /// @return A reference to the element at the specified indices
  scalar_type &at(const std::vector<size_t> &indices) {
    if (indices.size() != rank())
      throw std::runtime_error("Invalid indices provided to tensor::at(), size "
                               "must be equal to rank.");
    return pimpl->at(indices);
  }

  /// @brief Access a constant element of the tensor
  /// @param indices The indices of the element to access
  /// @return A constant reference to the element at the specified indices
  const scalar_type &at(const std::vector<size_t> &indices) const {
    return pimpl->at(indices);
  }

  /// @brief Copy data into the tensor
  /// @param data Pointer to the source data
  /// @param shape The shape of the source data
  void copy(const scalar_type *data,
            const std::vector<std::size_t> shape = {}) {
    if (pimpl->shape().empty() && shape.empty())
      throw std::runtime_error(
          "This tensor does not have a shape yet, must provide one to copy()");

    pimpl->copy(data, pimpl->shape().empty() ? shape : pimpl->shape());
  }

  /// @brief Take ownership of the given data
  /// @param data Pointer to the source data
  /// @param shape The shape of the source data
  void take(std::unique_ptr<scalar_type[]> &data,
            const std::vector<std::size_t> shape = {}) {
    if (pimpl->shape().empty() && shape.empty())
      throw std::runtime_error(
          "This tensor does not have a shape yet, must provide one to take()");

    pimpl->take(data, pimpl->shape().empty() ? shape : pimpl->shape());
  }

  /// @brief Borrow the given data without taking ownership
  /// @param data Pointer to the source data
  /// @param shape The shape of the source data
  void borrow(const scalar_type *data,
              const std::vector<std::size_t> shape = {}) {
    if (pimpl->shape().empty() && shape.empty())
      throw std::runtime_error("This tensor does not have a shape yet, must "
                               "provide one to borrow()");

    pimpl->borrow(data, pimpl->shape().empty() ? shape : pimpl->shape());
  }

  /// @brief Get a pointer to the raw data of the tensor.
  ///
  /// This method provides direct access to the underlying data storage of the
  /// tensor. It returns a pointer to the first element of the data array.
  ///
  /// @return `scalar_type*` A pointer to the mutable data of the tensor.
  ///
  /// @note Care should be taken when directly manipulating the raw data to
  /// avoid invalidating the internal state of a tensor or violating its
  /// `invariants`.
  scalar_type *data() { return pimpl->data(); }

  /// @brief Get a `const` pointer to the raw data of the tensor.
  ///
  /// This method provides read-only access to the underlying data storage of
  /// the tensor. It returns a constant pointer to the first element of the data
  /// array.
  ///
  /// @return `const scalar_type *` A constant pointer to the immutable data of
  /// the tensor.
  ///
  /// @note This constant version ensures that the data of a tensor cannot be
  /// modified through the returned pointer, preserving `const` correctness.
  const scalar_type *data() const { return pimpl->data(); }

  void dump() const { pimpl->dump(); }

  template <typename T>
  friend tensor<T> operator*(const tensor<T> &, const tensor<T> &);
  // template <typename T>
  // friend tensor<T> operator*(const tensor<T> &, const T &);
  template <typename T>
  friend tensor<T> operator+(const tensor<T> &, const tensor<T> &);
  template <typename T>
  friend tensor<T> kronecker(const tensor<T> &, const tensor<T> &);
};

/// Multiplication of two tensors.
template <typename T>
tensor<T> operator*(const tensor<T> &left, const tensor<T> &right) {
  if (left.shape() != right.shape()) {
    throw std::runtime_error("expected tensors of equal dimension.");
  }
  // Hard-coding this for 2D tensors only.
  if (left.shape().size() != 2) {
    throw std::runtime_error("expected 2D tensor.");
  }
  tensor<T> returnTensor(left.shape());
  for (std::size_t i = 0; i < left.shape()[0]; i++) {
    for (std::size_t j = 0; j < left.shape()[1]; j++) {
      T sum = 0;
      for (std::size_t k = 0; k < right.shape()[1]; k++) {
        sum += left.at({i, k}) * right.at({k, j});
      }
      returnTensor.at({i, j}) = sum;
    }
  }
  return returnTensor;
  // return tensor((*left.pimpl) * (*right.pimpl));
}

/// Multiplication of scalar against tensor.
template <typename T>
tensor<T> operator*(const T &left, const tensor<T> &right) {
  // Hard-coding this for 2D tensors only.
  if (right.shape().size() != 2) {
    throw std::runtime_error("expected 2D tensor.");
  }
  tensor<T> returnTensor(right.shape());
  for (std::size_t i = 0; i < right.shape()[0]; i++) {
    for (std::size_t j = 0; j < right.shape()[1]; j++) {
      returnTensor.at({i, j}) = (right.at({i, j}) * left);
    }
  }
  return returnTensor;
  // return tensor(left * (*right.pimpl));
}
template <typename T>
tensor<T> operator*(const tensor<T> &left, const T &right) {
  // Hard-coding this for 2D tensors only.
  if (left.shape().size() != 2) {
    throw std::runtime_error("expected 2D tensor.");
  }
  tensor<T> returnTensor(left.shape());
  for (std::size_t i = 0; i < left.shape()[0]; i++) {
    for (std::size_t j = 0; j < left.shape()[1]; j++) {
      returnTensor.at({i, j}) = (left.at({i, j}) * right);
    }
  }
  return returnTensor;
  // return tensor((*left.pimpl) * right);
}

/// Addition of two tensors.
template <typename T>
tensor<T> operator+(const tensor<T> &left, const tensor<T> &right) {
  if (left.shape() != right.shape()) {
    throw std::runtime_error("expected tensors of equal dimension.");
  }
  // Hard-coding this for 2D tensors only.
  if (right.shape().size() != 2) {
    throw std::runtime_error("expected 2D tensor.");
  }
  tensor<T> returnTensor(right.shape());
  for (std::size_t i = 0; i < right.shape()[0]; i++) {
    for (std::size_t j = 0; j < right.shape()[1]; j++) {
      returnTensor.at({i, j}) = right.at({i, j}) + left.at({i, j});
    }
  }
  return returnTensor;
  // return {(*left.pimpl) + (*right.pimpl)};
}

/// Subtraction of two tensors.
template <typename T>
tensor<T> operator-(const tensor<T> &left, const tensor<T> &right) {
  if (left.shape() != right.shape()) {
    throw std::runtime_error("expected tensors of equal dimension.");
  }
  // Hard-coding this for 2D tensors only.
  if (right.shape().size() != 2) {
    throw std::runtime_error("expected 2D tensor.");
  }
  tensor<T> returnTensor(right.shape());
  for (std::size_t i = 0; i < right.shape()[0]; i++) {
    for (std::size_t j = 0; j < right.shape()[1]; j++) {
      returnTensor.at({i, j}) = left.at({i, j}) - right.at({i, j});
    }
  }
  return returnTensor;
  // return {(*left.pimpl) + (*right.pimpl)};
}

/// Kronecker between two tensors.
template <typename T>
tensor<T> kronecker(const tensor<T> &left, const tensor<T> &right) {
  // Hard-coding this for 2D tensors only.
  if (left.shape().size() != 2 || right.shape().size() != 2) {
    throw std::runtime_error("expected 2D tensors.");
  }
  // If either of these is a 1x1 matrix, we will circumvent
  // the full kronecker product loop.
  if (left.size() == 1)
    return left.at({0,0}) * right;
  else if (right.size() == 1)
    return right.at({0,0}) * left;

  auto rowA = left.shape()[0];
  auto colA = left.shape()[1];
  auto rowB = right.shape()[0];
  auto colB = right.shape()[1];

  auto new_shape = {rowA * rowB, colA * colB};
  tensor<T> returnTensor(new_shape);

  std::cout << "rowA = " << left.shape()[0] << " \n";
  std::cout << "colA = " << left.shape()[1] << " \n";
  std::cout << "rowB = " << right.shape()[0] << " \n";
  std::cout << "colB = " << right.shape()[1] << " \n";

  for (std::size_t i = 0; i < rowA; i++) {
    for (std::size_t k = 0; k < colA; k++) {
      for (std::size_t j = 0; j < rowB; j++) {
        for (std::size_t l = 0; l < colB; l++) {
          std::cout << "i,k,j,l = " << i << ", " << k << ", " << j << ", " << l << "\n";
          std::cout << "first index = " << ((i * rowB) + k) << "\n";
          std::cout << "second index = " << ((j * colB) + l) << "\n";
          returnTensor.at(
              {(i * rowB) + k, (j * colB) + l}) =
              left.at({i, j}) * right.at({k, l});
              // left.at({i, k}) * right.at({j, l});
          returnTensor.dump();
        }
      }
    }
  }

  return returnTensor;
  // return {(*left.pimpl).kronecker((*right.pimpl))};
}

/// Kronecker product between a list of tensors.
template <typename T>
tensor<T> kronecker(std::list<tensor<T>> &tensors) {
  // tensor<T> start({1,1});
  // start.at({0,0}) = 1.0+0.0j;
  // auto func = [&](auto &&a, auto &&b) { return kronecker<T>(a,b); };
  // return std::accumulate(begin(tensors), end(tensors), start, func);
  std::vector<tensor<T>> returnTensors = {tensors.front()};
  tensors.pop_front();
  int idx = 0;
  for (auto tensor : tensors) {
    // tenstor<t> returnTensor({
    //   returnTensors[idx].shape()[0] * tensor.shape()[0],
    //   returnTensors[idx].shape()[1] * tensor.shape()[1]
    //  });
    auto returnTensor = kronecker(returnTensors[idx], tensor);
    returnTensors.push_back(returnTensor);
    idx++;
  }
  return returnTensors.back();
}

} // namespace cudaq
