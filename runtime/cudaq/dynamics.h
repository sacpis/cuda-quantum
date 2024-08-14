/****************************************************************-*- C++ -*-****
 * Copyright (c) 2022 - 2024 NVIDIA Corporation & Affiliates.                  *
 * All rights reserved.                                                        *
 *                                                                             *
 * This source code and the accompanying materials are made available under    *
 * the terms of the Apache License 2.0 which accompanies this distribution.    *
 ******************************************************************************/

#pragma once

#include "spin_op.h"

namespace cudaq {


// NOTE: Much of this will be dispersed out into more appropriate
// places in the codebase. Just consolidating it all here for now.

class EvolveResult;
class AsyncEvolveResult;

class EvolveResult {
public:

  /// @brief Constructor given just the system state.
  /// @arg state : The final state of the quantum system after time evolution.
  EvolveResult(State state);

  /// @brief Constructor given the system state and expectation.
  /// @arg state : The final state of the quantum system after time evolution.
  /// @arg expectation : The final expectation values after time evolution.
  EvolveResult(State state, std::vector<std::complex<double>> expectation);

  /// @brief Constructor given a sequence of quantum states.
  /// @arg states : The state of the quantum system at each discrete time step
  ///               during evolution.
  EvolveResult(std::vector<State> states);

  /// @brief Constructor given a sequence of quantum states and expectations.
  /// @arg states : The state of the quantum system at each discrete time step
  ///               during evolution.
  /// @arg expectations : The expectation values of the at each time step.
  EvolveResult(std::vector<State> states, std::vector<std::vector<std::complex<double>>> expectations);

  /// Attributes on the class.

  /// @brief Stores the final state of the quantum system produced by a call to `cudaq::evolve`.
  State final_state;

  /// @brief Stores the final expectation value produced by a call to `cudaq::evolve`.
  /// Each expectation value corresponds to one observable provided in the `cudaq::evolve`
  /// call. This value will be null if no observables were provided.
  std::vector<std::complex<double>> final_expectation;

  /// @brief Stores all intermediate states, meaning the state of the system after each time step
  /// in a defined schedule, produced by a call to `cudaq::evolve`. This includes the final state.
  /// This value is null unless requested in the call to `cudaq::evolve`.
  std::vector<State> intermediate_states;

  /// @brief Stores the expectation values at each step in the schedule produced by a call to 
  /// `cudaq::evolve`, including the final expectation values. Each expectation value 
  /// corresponds to one observable provided in the `cudaq.evolve` call. 
  /// This value is null if no observables were provided, or if it wasn't requested in the
  /// call to `cudaq::evolve`.
  std::vector<std::vector<std::complex<double>> expectation_values;

};

class AsyncEvolveResult {
  
public:
  /// @brief Create a class instance that can be used to retrieve the evolution
  /// result produced from the asynchronous execution, `cudaq::evolve_async`.
  /// It models a future-like type whose `cudaq::EvolveResult` may be accessed
  /// via an invocation of the `get` method.
  AsyncEvolveResult();

  /// @brief Retrieve the `cudaq::EvolveResult` from the asynchronous evolve
  /// execution. This causes the current thread to wait until the time evolution
  /// execution has completed.
  EvolveResult get();

  /// @brief Return the result as a string.
  std::string to_string() const;
};


class OperatorSum;
class ProductOperator(OperatorSum);
class ScalarOperator(OperatorSum);
class ElementaryOperator(ProductOperator);

/// @brief Represents an operator expression consisting of a sum of terms, where each
/// term is a product of elementary and scalar operators.
/// Operator expressions cannot be used within quantum kernels, but they provide methods
/// to convert them to data types that can.
class OperatorSum {

  /// @brief Construct a `cudaq::OperatorSum` given a sequence of 
  /// `cudaq::ProductOperator`'s.
  /// This operator expression represents a sum of terms, where each term
  /// is a product of elementary and scalar operators.
  OperatorSum(std::vector<ProductOperator> terms);

  // Arithmetic overloads against all other operator types.
  OperatorSum operator+(OperatorSum other);
  OperatorSum operator-(OperatorSum other);
  OperatorSum operator*(OperatorSum other);
  OperatorSum operator+=(OperatorSum other);
  OperatorSum operator-=(OperatorSum other);
  OperatorSum operator*=(OperatorSum other);
  OperatorSum operator+(ScalarOperator other);
  OperatorSum operator-(ScalarOperator other);
  OperatorSum operator*(ScalarOperator other);
  OperatorSum operator+=(ScalarOperator other);
  OperatorSum operator-=(ScalarOperator other);
  OperatorSum operator*=(ScalarOperator other);
  OperatorSum operator+(ProductOperator other);
  OperatorSum operator-(ProductOperator other);
  OperatorSum operator*(ProductOperator other);
  OperatorSum operator+=(ProductOperator other);
  OperatorSum operator-=(ProductOperator other);
  OperatorSum operator*=(ProductOperator other);
  OperatorSum operator+(ElementaryOperator other);
  OperatorSum operator-(ElementaryOperator other);
  OperatorSum operator*(ElementaryOperator other);
  OperatorSum operator+=(ElementaryOperator other);
  OperatorSum operator-=(ElementaryOperator other);
  OperatorSum operator*=(ElementaryOperator other);
  /// @brief  True, if the other value is an OperatorSum with equivalent terms, 
  /// and False otherwise. The equality takes into account that operator
  /// addition is commutative, as is the product of two operators if they
  /// act on different degrees of freedom.
  /// The equality comparison does *not* take commutation relations into 
  /// account, and does not try to reorder terms blockwise; it may hence 
  /// evaluate to False, even if two operators in reality are the same.
  /// If the equality evaluates to True, on the other hand, the operators 
  /// are guaranteed to represent the same transformation for all arguments.
  bool operator==(OperatorSum other);

  /// @brief Return the OperatorSum as a string.
  std::string to_string() const;

  /// FIXME: Needs to accept extra args...
  /// @brief Return the `OperatorSum` as a matrix.
  /// @arg  `dimensions` : A mapping that specifies the number of levels,
  ///                      that is, the dimension of each degree of freedom
  ///                      that the operator acts on.
  complex_matrix to_matrix(std::tuple<int,int> dimensions);

  /// @brief Creates a representation of the operator as a `cudaq::pauli_word` that
  /// can be passed as an argument to quantum kernels.
  pauli_word to_pauli_word();

  /// Attributes on the class.

  /// @brief The degrees of freedom that the oeprator acts on in canonical order.
  std::vector<int> degrees;

  /// FIXME: Will need an equivalent for `parameters` from python.
};

/// @brief Represents an operator expression consisting of a product of elementary
/// and scalar operators. 
/// Operator expressions cannot be used within quantum kernels, but 
/// they provide methods to convert them to data types that can.
class ProductOperator(OperatorSum) {

  /// @brief Constructor for an operator expression that represents a product
  /// of elementary operators.
  /// @arg atomic_operators : The operators of which to compute the product when
  ///                         evaluating the operator expression.
  ProductOperator(std::vector<ElementaryOperator> atomic_operators);

  /// @brief Constructor for an operator expression that represents a product
  /// of e scalar operators.
  /// @arg atomic_operators : The operators of which to compute the product when
  ///                         evaluating the operator expression.
  ProductOperator(std::vector<ScalarOperator> atomic_operators);

  // Arithmetic overloads against all other operator types.
  OperatorSum operator+(OperatorSum other);
  OperatorSum operator-(OperatorSum other);
  OperatorSum operator+=(OperatorSum other);
  OperatorSum operator-=(OperatorSum other);
  ProductOperator operator*(OperatorSum other);
  ProductOperator operator*=(OperatorSum other);
  OperatorSum operator+(ScalarOperator other);
  OperatorSum operator-(ScalarOperator other);
  OperatorSum operator+=(ScalarOperator other);
  OperatorSum operator-=(ScalarOperator other);
  ProductOperator operator*(ScalarOperator other);
  ProductOperator operator*=(ScalarOperator other);
  OperatorSum operator+(ProductOperator other);
  OperatorSum operator-(ProductOperator other);
  OperatorSum operator+=(ProductOperator other);
  OperatorSum operator-=(ProductOperator other);
  ProductOperator operator*(ProductOperator other);
  ProductOperator operator*=(ProductOperator other);
  OperatorSum operator+(ElementaryOperator other);
  OperatorSum operator-(ElementaryOperator other);
  OperatorSum operator+=(ElementaryOperator other);
  OperatorSum operator-=(ElementaryOperator other);
  ProductOperator operator*(ElementaryOperator other);
  ProductOperator operator*=(ElementaryOperator other);
  /// @brief True, if the other value is an OperatorSum with equivalent terms, 
  ///  and False otherwise. The equality takes into account that operator
  ///  addition is commutative, as is the product of two operators if they
  ///  act on different degrees of freedom.
  ///  The equality comparison does *not* take commutation relations into 
  ///  account, and does not try to reorder terms blockwise; it may hence 
  ///  evaluate to False, even if two operators in reality are the same.
  ///  If the equality evaluates to True, on the other hand, the operators 
  ///  are guaranteed to represent the same transformation for all arguments.
  bool operator==(ProductOperator other);

  /// @brief Return the `ProductOperator` as a string.
  std::string to_string() const;

  /// FIXME: Needs to accept extra args...
  /// @brief Return the `OperatorSum` as a matrix.
  /// @arg  `dimensions` : A mapping that specifies the number of levels,
  ///                      that is, the dimension of each degree of freedom
  ///                      that the operator acts on.
  complex_matrix to_matrix(std::tuple<int,int> dimensions);

  /// @brief Creates a representation of the operator as a `cudaq::pauli_word` that
  /// can be passed as an argument to quantum kernels.
  pauli_word to_pauli_word();

  /// @brief The degrees of freedom that the oeprator acts on in canonical order.
  std::vector<int> degrees;

  /// FIXME: Will need an equivalent for `parameters` from python.
};

class ScalarOperator(OperatorSum) {

};

class ElementaryOperator(ProductOperator) {

  /// @brief Constructor.
  /// @arg operator_id : The ID of the operator as specified when it was defined.
  /// @arg degrees : the degrees of freedom that the operator acts upon.
  ElementaryOperator(std::string operator_id, std::vector<int> degrees);

  // Arithmetic overloads against all other operator types.
  OperatorSum operator+(OperatorSum other);
  OperatorSum operator-(OperatorSum other);
  OperatorSum operator+=(OperatorSum other);
  OperatorSum operator-=(OperatorSum other);
  ProductOperator operator*(OperatorSum other);
  ProductOperator operator*=(OperatorSum other);
  OperatorSum operator+(ScalarOperator other);
  OperatorSum operator-(ScalarOperator other);
  OperatorSum operator+=(ScalarOperator other);
  OperatorSum operator-=(ScalarOperator other);
  ProductOperator operator*(ScalarOperator other);
  ProductOperator operator*=(ScalarOperator other);
  OperatorSum operator+(ProductOperator other);
  OperatorSum operator-(ProductOperator other);
  OperatorSum operator+=(ProductOperator other);
  OperatorSum operator-=(ProductOperator other);
  ProductOperator operator*(ProductOperator other);
  ProductOperator operator*=(ProductOperator other);
  OperatorSum operator+(ElementaryOperator other);
  OperatorSum operator-(ElementaryOperator other);
  OperatorSum operator+=(ElementaryOperator other);
  OperatorSum operator-=(ElementaryOperator other);
  ProductOperator operator*(ElementaryOperator other);
  ProductOperator operator*=(ElementaryOperator other);
  /// @brief True, if the other value is an elementary operator with the same id
  /// acting on the same degrees of freedom, and False otherwise.
  bool operator==(ElementaryOperator other);

  /// @brief Return the `ElementaryOperator` as a string.
  std::string to_string() const;

  /// @brief Return the `ElementaryOperator` as a matrix.
  /// @arg  `dimensions` : A mapping that specifies the number of levels,
  ///                      that is, the dimension of each degree of freedom
  ///                      that the operator acts on.
  complex_matrix to_matrix(std::tuple<int,int> dimensions);

  ElementaryOperator identity(int degree);
  ElementaryOperator zero(int degree);

  /// @brief Adds the definition of an elementary operator with the given id to the class.
  /// After definition, an the defined elementary operator can be instantiated by
  /// providing the operator id as well as the degree(s) of freedom that it acts on.
  /// An elementary operator is a parameterized object acting on certain degrees of 
  /// freedom. To evaluate an operator, for example to compute its matrix, the level, 
  /// that is the dimension, for each degree of freedom it acts on must be provided, 
  /// as well as all additional parameters. Additional parameters must be provided in 
  /// the form of keyword arguments.
  /// Note:
  /// The dimensions passed during operator evaluation are automatically validated 
  /// against the expected dimensions specified during definition - the `create` 
  /// function does not need to do this.
  /// @arg operator_id : A string that uniquely identifies the defined operator.
  /// @arg expected_dimensions : Defines the number of levels, that is the dimension,
  ///      for each degree of freedom in canonical (that is sorted) order. A 
  ///      negative or zero value for one (or more) of the expected dimensions 
  ///      indicates that the operator is defined for any dimension of the 
  ///      corresponding degree of freedom.
  /// @arg create : Takes any number of complex-valued arguments and returns the 
  ///      matrix representing the operator in canonical order. If the matrix can
  ///      be defined for any number of levels for one or more degree of freedom, 
  ///      the `create` function must take an argument called `dimension` (or `dim`
  ///      for short), if the operator acts on a single degree of freedom, and an
  ///      argument called `dimensions` (or `dims` for short), if the operator acts 
   ///     on multiple degrees of freedom.
  void define(std::string operator_id, std::vector<int> expected_dimensions, Callable create);

  /// Attributes.

  /// @brief The number of levels, that is the dimension, for each degree of freedom
  /// in canonical order that the operator acts on. A value of zero or less 
  /// indicates that the operator is defined for any dimension of that degree.
  std::vector<int> expected_dimensions;
  /// @brief The degrees of freedom that the operator acts on in canonical order.
  std::vector<int> degrees;
  std::string id;
  /// FIXME: Need `parameters` equivalent from python.

  pauli_word to_pauli_word ovveride (); 
};


// TODO
class OperatorArithmetics;
// TODO
class MatrixArithmetics(OperatorArithmetics);
// TODO
class PauliWordConversion(OperatorArithmetics);
// TODO
class PrettyPrint(OperatorArithmetics);


// TODO
class Schedule;

// TODO
class operators;

// TODO
class pauli;

// New CUDA-Q API functions.
// TODO
// EvolveResult cudaq::evolve();
// AsyncEvolveResult cudaq::evolve_async();


} // namespace cudaq