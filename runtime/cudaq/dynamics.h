/****************************************************************-*- C++ -*-****
 * Copyright (c) 2022 - 2024 NVIDIA Corporation & Affiliates.                  *
 * All rights reserved.                                                        *
 *                                                                             *
 * This source code and the accompanying materials are made available under    *
 * the terms of the Apache License 2.0 which accompanies this distribution.    *
 ******************************************************************************/

#pragma once

/// TODO: This file will effectively be deleted once all of its contents
/// are distributed and implemented.

#include "matrix.h"

// #include <cstddef>
// #include <iterator>

namespace cudaq {

class OperatorSum;
class ProductOperator(OperatorSum);
class ScalarOperator(ProductOperator);
class ElementaryOperator(ProductOperator);

// class EvolveResult;
// class AsyncEvolveResult;

// class EvolveResult {
// public:
//   /// @brief Constructor given just the system state.
//   /// @arg state : The final state of the quantum system after time
//   evolution. EvolveResult(state state);

//   /// @brief Constructor given the system state and expectation.
//   /// @arg state : The final state of the quantum system after time
//   evolution.
//   /// @arg expectation : The final expectation values after time evolution.
//   EvolveResult(state state, std::vector<std::complex<double>> expectation);

//   /// @brief Constructor given a sequence of quantum states.
//   /// @arg states : The state of the quantum system at each discrete time
//   step
//   ///               during evolution.
//   EvolveResult(std::vector<state> states);

//   /// @brief Constructor given a sequence of quantum states and expectations.
//   /// @arg states : The state of the quantum system at each discrete time
//   step
//   ///               during evolution.
//   /// @arg expectations : The expectation values of the at each time step.
//   EvolveResult(std::vector<state> &states,
//                std::vector<std::vector<std::complex<double>>> &expectations);

//   // Attributes on the class.

//   /// @brief Stores the final state of the quantum system produced by a call
//   to
//   /// `cudaq::evolve`.
//   state final_state;

//   /// @brief Stores the final expectation value produced by a call to
//   /// `cudaq::evolve`. Each expectation value corresponds to one observable
//   /// provided in the `cudaq::evolve` call. This value will be null if no
//   /// observables were provided.
//   std::vector<std::complex<double>> final_expectation;

//   /// @brief Stores all intermediate states, meaning the state of the system
//   /// after each time step in a defined schedule, produced by a call to
//   /// `cudaq::evolve`. This includes the final state. This value is null
//   unless
//   /// requested in the call to `cudaq::evolve`.
//   std::vector<state> &intermediate_states;

//   /// @brief Stores the expectation values at each step in the schedule
//   produced
//   /// by a call to `cudaq::evolve`, including the final expectation values.
//   Each
//   /// expectation value corresponds to one observable provided in the
//   /// `cudaq.evolve` call. This value is null if no observables were
//   provided,
//   /// or if it wasn't requested in the call to `cudaq::evolve`.
//   std::vector<std::vector<std::complex<double>>> expectation_values;
// };

// class AsyncEvolveResult {

// public:
//   /// @brief Create a class instance that can be used to retrieve the
//   evolution
//   /// result produced from the asynchronous execution, `cudaq::evolve_async`.
//   /// It models a future-like type whose `cudaq::EvolveResult` may be
//   accessed
//   /// via an invocation of the `get` method.
//   AsyncEvolveResult();

//   /// @brief Retrieve the `cudaq::EvolveResult` from the asynchronous evolve
//   /// execution. This causes the current thread to wait until the time
//   evolution
//   /// execution has completed.
//   EvolveResult get();

//   /// @brief Return the result as a string.
//   std::string to_string() const;
// };

// // class OperatorSum;
// // class ProductOperator(OperatorSum);
// // class ScalarOperator(ProductOperator);
// // class ElementaryOperator(ProductOperator);

// // /// @brief Represents an operator expression consisting of a sum of terms,
// where
// // /// each term is a product of elementary and scalar operators. Operator
// // /// expressions cannot be used within quantum kernels, but they provide
// methods
// // /// to convert them to data types that can.
// // class OperatorSum {

// //   /// @brief Construct a `cudaq::OperatorSum` given a sequence of
// //   /// `cudaq::ProductOperator`'s.
// //   /// This operator expression represents a sum of terms, where each term
// //   /// is a product of elementary and scalar operators.
// //   OperatorSum(std::vector<ProductOperator> terms);

// //   // Arithmetic overloads against all other operator types.
// //   OperatorSum operator+(double other);
// //   OperatorSum operator-(double other);
// //   OperatorSum operator*(double other);
// //   OperatorSum operator+=(double other);
// //   OperatorSum operator-=(double other);
// //   OperatorSum operator*=(double other);
// //   OperatorSum operator+(std::complex<double> other);
// //   OperatorSum operator-(std::complex<double> other);
// //   OperatorSum operator*(std::complex<double> other);
// //   OperatorSum operator+=(std::complex<double> other);
// //   OperatorSum operator-=(std::complex<double> other);
// //   OperatorSum operator*=(std::complex<double> other);
// //   OperatorSum operator+(OperatorSum other);
// //   OperatorSum operator-(OperatorSum other);
// //   OperatorSum operator*(OperatorSum other);
// //   OperatorSum operator+=(OperatorSum other);
// //   OperatorSum operator-=(OperatorSum other);
// //   OperatorSum operator*=(OperatorSum other);
// //   OperatorSum operator+(ScalarOperator other);
// //   OperatorSum operator-(ScalarOperator other);
// //   OperatorSum operator*(ScalarOperator other);
// //   OperatorSum operator+=(ScalarOperator other);
// //   OperatorSum operator-=(ScalarOperator other);
// //   OperatorSum operator*=(ScalarOperator other);
// //   OperatorSum operator+(ProductOperator other);
// //   OperatorSum operator-(ProductOperator other);
// //   OperatorSum operator*(ProductOperator other);
// //   OperatorSum operator+=(ProductOperator other);
// //   OperatorSum operator-=(ProductOperator other);
// //   OperatorSum operator*=(ProductOperator other);
// //   OperatorSum operator+(ElementaryOperator other);
// //   OperatorSum operator-(ElementaryOperator other);
// //   OperatorSum operator*(ElementaryOperator other);
// //   OperatorSum operator+=(ElementaryOperator other);
// //   OperatorSum operator-=(ElementaryOperator other);
// //   OperatorSum operator*=(ElementaryOperator other);
// //   /// @brief  True, if the other value is an OperatorSum with equivalent
// terms,
// //   /// and False otherwise. The equality takes into account that operator
// //   /// addition is commutative, as is the product of two operators if they
// //   /// act on different degrees of freedom.
// //   /// The equality comparison does *not* take commutation relations into
// //   /// account, and does not try to reorder terms blockwise; it may hence
// //   /// evaluate to False, even if two operators in reality are the same.
// //   /// If the equality evaluates to True, on the other hand, the operators
// //   /// are guaranteed to represent the same transformation for all
// arguments.
// //   bool operator==(OperatorSum other);

// //   /// @brief Return the OperatorSum as a string.
// //   std::string to_string() const;

// //   /// @brief Return the `OperatorSum` as a matrix.
// //   /// @arg  `dimensions` : A mapping that specifies the number of levels,
// //   ///                      that is, the dimension of each degree of
// freedom
// //   ///                      that the operator acts on. Example for two,
// 2-level
// //   ///                      degrees of freedom: `{0:2, 1:2}`.
// //   /// @arg `parameters` : A map of the paramter names to their concrete,
// complex
// //   /// values.
// //   complex_matrix
// //   to_matrix(std::map<int, int> dimensions,
// //             std::map<std::string, std::complex<double>> parameters);

// //   /// @brief Creates a representation of the operator as a
// `cudaq::pauli_word`
// //   /// that can be passed as an argument to quantum kernels.
// //   pauli_word to_pauli_word();

// //   /// Attributes on the class.

// //   /// @brief The degrees of freedom that the oeprator acts on in canonical
// //   /// order.
// //   std::vector<int> degrees;

// //   /// @brief A map of the paramter names to their concrete, complex
// values.
// //   std::map<std::string, std::complex<double>> parameters;
// // };

// // /// @brief Represents an operator expression consisting of a product of
// // /// elementary and scalar operators. Operator expressions cannot be used
// within
// // /// quantum kernels, but they provide methods to convert them to data
// types that
// // /// can.
// // class ProductOperator {

// //   /// @brief Constructor for an operator expression that represents a
// product
// //   /// of elementary operators.
// //   /// @arg atomic_operators : The operators of which to compute the
// product when
// //   ///                         evaluating the operator expression.
// //   ProductOperator(std::vector<ElementaryOperator> atomic_operators);

// //   /// @brief Constructor for an operator expression that represents a
// product
// //   /// of e scalar operators.
// //   /// @arg atomic_operators : The operators of which to compute the
// product when
// //   ///                         evaluating the operator expression.
// //   ProductOperator(std::vector<ScalarOperator> atomic_operators);

// //   // Arithmetic overloads against all other operator types.
// //   OperatorSum operator+(int other);
// //   OperatorSum operator-(int other);
// //   OperatorSum operator+=(int other);
// //   OperatorSum operator-=(int other);
// //   ProductOperator operator*(int other);
// //   ProductOperator operator*=(int other);
// //   OperatorSum operator+(double other);
// //   OperatorSum operator-(double other);
// //   OperatorSum operator+=(double other);
// //   OperatorSum operator-=(double other);
// //   ProductOperator operator*(double other);
// //   ProductOperator operator*=(double other);
// //   OperatorSum operator+(std::complex<double> other);
// //   OperatorSum operator-(std::complex<double> other);
// //   OperatorSum operator+=(dostd::complex<double> other);
// //   OperatorSum operator-=(std::complex<double> other);
// //   ProductOperator operator*(std::complex<double> other);
// //   ProductOperator operator*=(std::complex<double> other);
// //   OperatorSum operator+(OperatorSum other);
// //   OperatorSum operator-(OperatorSum other);
// //   OperatorSum operator+=(OperatorSum other);
// //   OperatorSum operator-=(OperatorSum other);
// //   ProductOperator operator*(OperatorSum other);
// //   ProductOperator operator*=(OperatorSum other);
// //   OperatorSum operator+(ScalarOperator other);
// //   OperatorSum operator-(ScalarOperator other);
// //   OperatorSum operator+=(ScalarOperator other);
// //   OperatorSum operator-=(ScalarOperator other);
// //   ProductOperator operator*(ScalarOperator other);
// //   ProductOperator operator*=(ScalarOperator other);
// //   OperatorSum operator+(ProductOperator other);
// //   OperatorSum operator-(ProductOperator other);
// //   OperatorSum operator+=(ProductOperator other);
// //   OperatorSum operator-=(ProductOperator other);
// //   ProductOperator operator*(ProductOperator other);
// //   ProductOperator operator*=(ProductOperator other);
// //   OperatorSum operator+(ElementaryOperator other);
// //   OperatorSum operator-(ElementaryOperator other);
// //   OperatorSum operator+=(ElementaryOperator other);
// //   OperatorSum operator-=(ElementaryOperator other);
// //   ProductOperator operator*(ElementaryOperator other);
// //   ProductOperator operator*=(ElementaryOperator other);
// //   /// @brief True, if the other value is an OperatorSum with equivalent
// terms,
// //   ///  and False otherwise. The equality takes into account that operator
// //   ///  addition is commutative, as is the product of two operators if they
// //   ///  act on different degrees of freedom.
// //   ///  The equality comparison does *not* take commutation relations into
// //   ///  account, and does not try to reorder terms blockwise; it may hence
// //   ///  evaluate to False, even if two operators in reality are the same.
// //   ///  If the equality evaluates to True, on the other hand, the operators
// //   ///  are guaranteed to represent the same transformation for all
// arguments.
// //   bool operator==(ProductOperator other);

// //   /// @brief Return the `ProductOperator` as a string.
// //   std::string to_string() const;

// //   /// @brief Return the `OperatorSum` as a matrix.
// //   /// @arg  `dimensions` : A mapping that specifies the number of levels,
// //   ///                      that is, the dimension of each degree of
// freedom
// //   ///                      that the operator acts on. Example for two,
// 2-level
// //   ///                      degrees of freedom: `{0:2, 1:2}`.
// //   /// @arg `parameters` : A map of the paramter names to their concrete,
// complex
// //   /// values.
// //   complex_matrix to_matrix(
// //       std::map<int, int> dimensions,
// //       std::map<std::string, std::complex<double>> parameters);

// //   /// @brief Creates a representation of the operator as a
// `cudaq::pauli_word`
// //   /// that can be passed as an argument to quantum kernels.
// //   pauli_word to_pauli_word();

// //   /// @brief The degrees of freedom that the oeprator acts on in canonical
// //   /// order.
// //   std::vector<int> degrees;

// //   /// @brief A map of the paramter names to their concrete, complex
// values.
// //   std::map<std::string, std::complex<double>> parameters;
// // };

// // class ScalarOperator {

// //   /// @brief Constructor.
// //   /// @arg generator: The value of the scalar operator as a function of
// its
// //   /// parameters. The generator may take any number of complex-valued
// arguments
// //   /// and must return a number.
// //   ScalarOperator(Callable generator,
// //                  std::map<std::string, std::complex<double>> parameters);

// //   // Arithmetic overloads against all other operator types.
// //   ScalarOperator operator+(int other);
// //   ScalarOperator operator-(int other);
// //   ScalarOperator operator+=(int other);
// //   ScalarOperator operator-=(int other);
// //   ScalarOperator operator*(int other);
// //   ScalarOperator operator*=(int other);
// //   ScalarOperator operator/(int other);
// //   ScalarOperator operator/=(int other);
// //   ScalarOperator operator+(double other);
// //   ScalarOperator operator-(double other);
// //   ScalarOperator operator+=(double other);
// //   ScalarOperator operator-=(double other);
// //   ScalarOperator operator*(double other);
// //   ScalarOperator operator*=(double other);
// //   ScalarOperator operator/(double other);
// //   ScalarOperator operator/=(double other);
// //   ScalarOperator operator+(std::complex<double> other);
// //   ScalarOperator operator-(std::complex<double> other);
// //   ScalarOperator operator+=(std::complex<double> other);
// //   ScalarOperator operator-=(std::complex<double> other);
// //   ScalarOperator operator*(std::complex<double> other);
// //   ScalarOperator operator*=(std::complex<double> other);
// //   ScalarOperator operator/(std::complex<double> other);
// //   ScalarOperator operator/=(std::complex<double> other);
// //   ScalarOperator operator+(OperatorSum other);
// //   ScalarOperator operator-(OperatorSum other);
// //   ScalarOperator operator+=(OperatorSum other);
// //   ScalarOperator operator-=(OperatorSum other);
// //   ScalarOperator operator*(OperatorSum other);
// //   ScalarOperator operator*=(OperatorSum other);
// //   ScalarOperator operator/(OperatorSum other);
// //   ScalarOperator operator/=(OperatorSum other);
// //   ScalarOperator operator+(ScalarOperator other);
// //   ScalarOperator operator-(ScalarOperator other);
// //   ScalarOperator operator+=(ScalarOperator other);
// //   ScalarOperator operator-=(ScalarOperator other);
// //   ScalarOperator operator*(ScalarOperator other);
// //   ScalarOperator operator*=(ScalarOperator other);
// //   ScalarOperator operator/(ScalarOperator other);
// //   ScalarOperator operator/=(ScalarOperator other);
// //   ScalarOperator pow(ScalarOperator other);
// //   ScalarOperator operator+(ProductOperator other);
// //   ScalarOperator operator-(ProductOperator other);
// //   ScalarOperator operator+=(ProductOperator other);
// //   ScalarOperator operator-=(ProductOperator other);
// //   ScalarOperator operator*(ProductOperator other);
// //   ScalarOperator operator*=(ProductOperator other);
// //   ScalarOperator operator/(ProductOperator other);
// //   ScalarOperator operator/=(ProductOperator other);
// //   ScalarOperator operator+(ElementaryOperator other);
// //   ScalarOperator operator-(ElementaryOperator other);
// //   ScalarOperator operator+=(ElementaryOperator other);
// //   ScalarOperator operator-=(ElementaryOperator other);
// //   ScalarOperator operator*(ElementaryOperator other);
// //   ScalarOperator operator*=(ElementaryOperator other);
// //   ScalarOperator operator/(ElementaryOperator other);
// //   ScalarOperator operator/=(ElementaryOperator other);

// //   /// @brief Returns true if other is a scalar operator with the same
// //   /// generator.
// //   bool operator==(ScalarOperator other);

// //   /// @brief Return the `OperatorSum` as a matrix.
// //   /// @arg  `dimensions` : A mapping that specifies the number of levels,
// //   ///                      that is, the dimension of each degree of
// freedom
// //   ///                      that the operator acts on. Example for two,
// 2-level
// //   ///                      degrees of freedom: `{0:2, 1:2}`.
// //   /// @arg `parameters` : A map of the paramter names to their concrete,
// complex
// //   /// values.
// //   complex_matrix to_matrix(
// //       std::map<int, int> dimensions,
// //       std::map<std::string, std::complex<double>> parameters);

// //   /// @brief Creates a representation of the operator as `pauli_word` that
// can
// //   /// be passed as an argument to quantum kernels.
// //   pauli_word to_pauli_word ovveride();

// //   /// @brief The number of levels, that is the dimension, for each degree
// of
// //   /// freedom in canonical order that the operator acts on. A value of
// zero or
// //   /// less indicates that the operator is defined for any dimension of
// that
// //   /// degree.
// //   std::vector<int> dimensions;
// //   /// @brief The degrees of freedom that the operator acts on in canonical
// //   /// order.
// //   std::vector<int> degrees;
// //   /// A map of the paramter names to their concrete, complex
// //   /// values.
// //   std::map<std::string, std::complex<double>> parameters;

// //   /// @brief The function that generates the value of the scalar operator.
// //   /// The function can take any number of complex-valued arguments
// //   /// and returns a number.
// //   Callable generator;
// // };

// // class ElementaryOperator {

// //   /// @brief Constructor.
// //   /// @arg operator_id : The ID of the operator as specified when it was
// //   /// defined.
// //   /// @arg degrees : the degrees of freedom that the operator acts upon.
// //   ElementaryOperator(std::string operator_id, std::vector<int> degrees);

// //   // Arithmetic overloads against all other operator types.
// //   OperatorSum operator+(OperatorSum other);
// //   OperatorSum operator-(OperatorSum other);
// //   OperatorSum operator+=(OperatorSum other);
// //   OperatorSum operator-=(OperatorSum other);
// //   ProductOperator operator*(OperatorSum other);
// //   ProductOperator operator*=(OperatorSum other);
// //   OperatorSum operator+(ScalarOperator other);
// //   OperatorSum operator-(ScalarOperator other);
// //   OperatorSum operator+=(ScalarOperator other);
// //   OperatorSum operator-=(ScalarOperator other);
// //   ProductOperator operator*(ScalarOperator other);
// //   ProductOperator operator*=(ScalarOperator other);
// //   OperatorSum operator+(ProductOperator other);
// //   OperatorSum operator-(ProductOperator other);
// //   OperatorSum operator+=(ProductOperator other);
// //   OperatorSum operator-=(ProductOperator other);
// //   ProductOperator operator*(ProductOperator other);
// //   ProductOperator operator*=(ProductOperator other);
// //   OperatorSum operator+(ElementaryOperator other);
// //   OperatorSum operator-(ElementaryOperator other);
// //   OperatorSum operator+=(ElementaryOperator other);
// //   OperatorSum operator-=(ElementaryOperator other);
// //   ProductOperator operator*(ElementaryOperator other);
// //   ProductOperator operator*=(ElementaryOperator other);
// //   /// @brief True, if the other value is an elementary operator with the
// same id
// //   /// acting on the same degrees of freedom, and False otherwise.
// //   bool operator==(ElementaryOperator other);

// //   /// @brief Return the `ElementaryOperator` as a string.
// //   std::string to_string() const;

// //   /// @brief Return the `ElementaryOperator` as a matrix.
// //   /// @arg  `dimensions` : A mapping that specifies the number of levels,
// //   ///                      that is, the dimension of each degree of
// freedom
// //   ///                      that the operator acts on. Example for two,
// 2-level
// //   ///                      degrees of freedom: `{0:2, 1:2}`.
// //   complex_matrix to_matrix(std::map<int, int> dimensions);

// //   ElementaryOperator identity(int degree);
// //   ElementaryOperator zero(int degree);

// //   /// @brief Adds the definition of an elementary operator with the given
// id to
// //   /// the class. After definition, an the defined elementary operator can
// be
// //   /// instantiated by providing the operator id as well as the degree(s)
// of
// //   /// freedom that it acts on. An elementary operator is a parameterized
// object
// //   /// acting on certain degrees of freedom. To evaluate an operator, for
// example
// //   /// to compute its matrix, the level, that is the dimension, for each
// degree
// //   /// of freedom it acts on must be provided, as well as all additional
// //   /// parameters. Additional parameters must be provided in the form of
// keyword
// //   /// arguments. Note: The dimensions passed during operator evaluation
// are
// //   /// automatically validated against the expected dimensions specified
// during
// //   /// definition - the `create` function does not need to do this.
// //   /// @arg operator_id : A string that uniquely identifies the defined
// operator.
// //   /// @arg expected_dimensions : Defines the number of levels, that is the
// //   /// dimension,
// //   ///      for each degree of freedom in canonical (that is sorted) order.
// A
// //   ///      negative or zero value for one (or more) of the expected
// dimensions
// //   ///      indicates that the operator is defined for any dimension of the
// //   ///      corresponding degree of freedom.
// //   /// @arg create : Takes any number of complex-valued arguments and
// returns the
// //   ///      matrix representing the operator in canonical order. If the
// matrix
// //   ///      can be defined for any number of levels for one or more degree
// of
// //   ///      freedom, the `create` function must take an argument called
// //   ///      `dimension` (or `dim` for short), if the operator acts on a
// single
// //   ///      degree of freedom, and an argument called `dimensions` (or
// `dims` for
// //   ///      short), if the operator acts
// //   ///     on multiple degrees of freedom.
// //   void define(std::string operator_id, std::vector<int>
// expected_dimensions,
// //               Callable create);

// //   /// Attributes.

// //   /// @brief The number of levels, that is the dimension, for each degree
// of
// //   /// freedom in canonical order that the operator acts on. A value of
// zero or
// //   /// less indicates that the operator is defined for any dimension of
// that
// //   /// degree.
// //   std::vector<int> expected_dimensions;
// //   /// @brief The degrees of freedom that the operator acts on in canonical
// //   /// order.
// //   std::vector<int> degrees;
// //   /// @brief A map of the paramter names to their concrete, complex
// values.
// //   std::map<std::string, std::complex<double>> parameters;
// //   std::string id;

// //   /// @brief Creates a representation of the operator as `pauli_word` that
// can
// //   /// be passed as an argument to quantum kernels.
// //   pauli_word to_pauli_word ovveride();
// // };

// /// @brief FIXME
// class OperatorArithmetics;
// /// @brief FIXME
// class MatrixArithmetics(OperatorArithmetics);

// /// @brief Create a schedule for evaluating an operator expression at
// different
// /// steps.
// class Schedule {
// public:
//   /// Iterator tags. May be superfluous.
//   using iterator_category = std::forward_iterator_tag;
//   using difference_type = std::ptrdiff_t;
//   using value_type = std::complex<double>;
//   using pointer = std::complex<double> *;
//   using reference = std::complex<double> &;

// private:
//   pointer m_ptr;

// public:
//   Schedule(pointer ptr) : m_ptr(ptr) {}

//   /// @brief Constructor.
//   /// @arg steps: The sequence of steps in the schedule. Restricted to a
//   vector
//   /// of complex values.
//   /// @arg parameters: A sequence of strings representing the parameter names
//   of
//   /// an operator expression.
//   /// @arg value_function: A function that takes the name of a parameter as
//   well
//   /// as an additional value ("step") of arbitrary type as argument and
//   returns
//   /// the complex value for that parameter at the given step.
//   Schedule(std::vector<std::complex<double>> steps,
//            std::vector<std::string> parameters, Callable value_function);

//   // Below, I define what I believe are the minimal necessary methods needed
//   // for this to behave like an iterable. This should be revisited in the
//   // implementation phase.

//   /// @brief Reset the underlying schedule.
//   void reset();

//   /// @brief Pointer
//   reference operator*() const { return *m_ptr; }
//   /// @brief Access
//   pointer operator->() { return m_ptr; }
//   /// @brief Prefix increment.
//   Iterator &operator++() {
//     m_ptr++;
//     return *this;
//   }
//   /// @brief Postfix increment.
//   Iterator operator++(std::complex<double>) {
//     Iterator tmp = *this;
//     ++(*this);
//     return tmp;
//   }
//   /// @brief Comparison.
//   friend bool operator==(const Iterator &a, const Iterator &b) {
//     return a.m_ptr == b.m_ptr;
//   };
//   /// @brief Comparison.
//   friend bool operator!=(const Iterator &a, const Iterator &b) {
//     return a.m_ptr != b.m_ptr;
//   };

//   // TODO: reset option
// };

// /// @brief Quantum Operators to use as building blocks of more
// /// complex Hamiltonians.
// class operators {
//   /// @brief Annihilation Operator.
//   ElementaryOperator annihilate(int degree);
//   /// @brief Creation Operator.
//   ElementaryOperator create(int degree);
//   /// @brief Displacement Operator.
//   ElementaryOperator displace(int degree);
//   /// @brief Momentum Operator.
//   ElementaryOperator momentum(int degree);
//   /// @brief Number Operator.
//   ElementaryOperator number(int degree);
//   /// @brief Parity Operator.
//   ElementaryOperator parity(int degree);
//   /// @brief Position Operator.
//   ElementaryOperator position(int degree);
//   /// @brief Squeeze Operator.
//   ElementaryOperator squeeze(int degree);
//   /// @brief Constant Operator.
//   ScalarOperator constant(NumericType value);
//   /// @brief Identity Operator given a single degree of freedom. May act on
//   any d-level qudit. ElementaryOperator identity(int degree);
//   /// @brief Identity Operator given multiple degrees of freedom.
//   ProductOperator identity(std::vector<int> degrees);
//   /// @brief Zero Operator given a single degree of freedom.
//   ElementaryOperator zero(int degree);
//   /// @brief Zero Operator given multiple degrees of freedom.
//   ProductOperator zero(std::vector<int> degrees);
// };

// /// @brief Pauli Spin Operators to use as building blocks of more
// /// complex Hamiltonians.
// class pauli {
//   /// @brief Pauli Identity. Restricted to qubits.
//   ElementaryOperator i(int degree);
//   /// @brief Pauli X.
//   ElementaryOperator x(int degree);
//   /// @brief Pauli Y.
//   ElementaryOperator y(int degree);
//   /// @brief Pauli Z.
//   ElementaryOperator z(int degree);
//   /// @brief Pauli +
//   ElementaryOperator plus(int degree);
//   /// @brief Pauli -
//   ElementaryOperator minus(int degree);
// };

// // New CUDA-Q API functions.
// /// @brief The minimal call to `cudaq::evolve`.
// /// @arg hamiltonian: The hamiltonian to evolve in time.
// /// @arg dimensions: A mapping that specifies the number of levels, that is,
// the
// /// dimension of each degree of freedom that the operator acts on.
// /// @arg schedule: The discrete time schedule to apply time evolution over.
// /// @arg store_intermediate_states: Store the state of the system at each
// /// discrete time step of evolution. Default is set to false.
// EvolveResult cudaq::evolve(Operator hamiltonian, std::map<int, int>
// dimensions,
//                            Schedule schedule,
//                            bool store_intermediate_states = false);

// /// @brief Evolve the system in time, provided a set of collapse operators
// and
// /// observables to take the expectation value with respect to.
// /// @arg hamiltonian: The hamiltonian to evolve in time.
// /// @arg dimensions: A mapping that specifies the number of levels, that is,
// the
// /// dimension of each degree of freedom that the operator acts on.
// /// @arg schedule: The discrete time schedule to apply time evolution over.
// /// @arg collapse_operators
// /// @arg observables: Operators to take the expectation value of the system
// with
// /// respect to.
// /// @arg store_intermediate_states: Store the state of the system at each
// /// discrete time step of evolution. Default is set to false.
// EvolveResult cudaq::evolve(Operator hamiltonian, std::map<int, int>
// dimensions,
//                            Schedule schedule,
//                            std::vector<Operator> collapse_operators={},
//                            std::vector<Operator> observables={},
//                            bool store_intermediate_states = false);

// /// @brief Evolve the system in time, provided a set of collapse operators
// and
// /// observables to take the expectation value with respect to.
// /// @arg hamiltonian: The hamiltonian to evolve in time.
// /// @arg dimensions: A mapping that specifies the number of levels, that is,
// the
// /// dimension of each degree of freedom that the operator acts on.
// /// @arg schedule: The discrete time schedule to apply time evolution over.
// /// @arg initial_state: The initial quantum state for the system to evolve
// from.
// /// @arg collapse_operators
// /// @arg observables: Operators to take the expectation value of the system
// with
// /// respect to.
// /// @arg store_intermediate_states: Store the state of the system at each
// /// discrete time step of evolution. Default is set to false.
// EvolveResult cudaq::evolve(Operator hamiltonian, std::map<int, int>
// dimensions,
//                            Schedule schedule, state initial_state,
//                            std::vector<Operator> collapse_operators={},
//                            std::vector<Operator> observables={},
//                            bool store_intermediate_states = false);

// /// @brief Evolve the system in time, provided a set of collapse operators
// and
// /// observables to take the expectation value with respect to.
// /// @arg hamiltonian: The hamiltonian to evolve in time.
// /// @arg dimensions: A mapping that specifies the number of levels, that is,
// the
// /// dimension of each degree of freedom that the operator acts on.
// /// @arg schedule: The discrete time schedule to apply time evolution over.
// /// @arg initial_states: A list of initial quantum states for the system to
// /// evolve from. Will result in multiple, independent time evolutions.
// /// @arg collapse_operators
// /// @arg observables: Operators to take the expectation value of the system
// with
// /// respect to.
// /// @arg store_intermediate_states: Store the state of the system at each
// /// discrete time step of evolution. Default is set to false.
// EvolveResult cudaq::evolve(Operator hamiltonian, std::map<int, int>
// dimensions,
//                            Schedule schedule, std::vector<state>
//                            initial_states, std::vector<Operator>
//                            collapse_operators={}, std::vector<Operator>
//                            observables={}, bool store_intermediate_states =
//                            false);

// // New CUDA-Q API functions.
// /// @brief The minimal call to `cudaq::evolve`.
// /// @arg hamiltonian: The hamiltonian to evolve in time.
// /// @arg dimensions: A mapping that specifies the number of levels, that is,
// the
// /// dimension of each degree of freedom that the operator acts on.
// /// @arg schedule: The discrete time schedule to apply time evolution over.
// /// @arg store_intermediate_states: Store the state of the system at each
// /// discrete time step of evolution. Default is set to false.
// AsyncEvolveResult cudaq::evolve_async(Operator hamiltonian,
//                                       std::map<int, int> dimensions,
//                                       Schedule schedule,
//                                       bool store_intermediate_states =
//                                       false);

// /// @brief Evolve the system in time, provided a set of collapse operators
// and
// /// observables to take the expectation value with respect to.
// /// @arg hamiltonian: The hamiltonian to evolve in time.
// /// @arg dimensions: A mapping that specifies the number of levels, that is,
// the
// /// dimension of each degree of freedom that the operator acts on.
// /// @arg schedule: The discrete time schedule to apply time evolution over.
// /// @arg collapse_operators
// /// @arg observables: Operators to take the expectation value of the system
// with
// /// respect to.
// /// @arg store_intermediate_states: Store the state of the system at each
// /// discrete time step of evolution. Default is set to false.
// AsyncEvolveResult cudaq::evolve_async(Operator hamiltonian,
//                                       std::map<int, int> dimensions,
//                                       Schedule schedule,
//                                       std::vector<Operator>
//                                       collapse_operators={},
//                                       std::vector<Operator> observables={},
//                                       bool store_intermediate_states =
//                                       false);

// /// @brief Evolve the system in time, provided a set of collapse operators
// and
// /// observables to take the expectation value with respect to.
// /// @arg hamiltonian: The hamiltonian to evolve in time.
// /// @arg dimensions: A mapping that specifies the number of levels, that is,
// the
// /// dimension of each degree of freedom that the operator acts on.
// /// @arg schedule: The discrete time schedule to apply time evolution over.
// /// @arg initial_state: The initial quantum state for the system to evolve
// from.
// /// @arg collapse_operators
// /// @arg observables: Operators to take the expectation value of the system
// with
// /// respect to.
// /// @arg store_intermediate_states: Store the state of the system at each
// /// discrete time step of evolution. Default is set to false.
// AsyncEvolveResult cudaq::evolve_async(Operator hamiltonian,
//                                       std::map<int, int> dimensions,
//                                       Schedule schedule, state initial_state,
//                                       std::vector<Operator>
//                                       collapse_operators={},
//                                       std::vector<Operator> observables={},
//                                       bool store_intermediate_states =
//                                       false);

// /// @brief Evolve the system in time, provided a set of collapse operators
// and
// /// observables to take the expectation value with respect to.
// /// @arg hamiltonian: The hamiltonian to evolve in time.
// /// @arg dimensions: A mapping that specifies the number of levels, that is,
// the
// /// dimension of each degree of freedom that the operator acts on.
// /// @arg schedule: The discrete time schedule to apply time evolution over.
// /// @arg initial_states: A list of initial quantum states for the system to
// /// evolve from. Will result in multiple, independent time evolutions.
// /// @arg collapse_operators
// /// @arg observables: Operators to take the expectation value of the system
// with
// /// respect to.
// /// @arg store_intermediate_states: Store the state of the system at each
// /// discrete time step of evolution. Default is set to false.
// AsyncEvolveResult cudaq::evolve_async(
//     Operator hamiltonian, std::map<int, int> dimensions, Schedule schedule,
//     std::vector<state> initial_states, std::vector<Operator>
//     collapse_operators={}, std::vector<Operator> observables={}, bool
//     store_intermediate_states = false);

// // Internal facing.
// class PauliWordConversion(OperatorArithmetics);
// // Internal facing.
// class PrettyPrint(OperatorArithmetics);

} // namespace cudaq