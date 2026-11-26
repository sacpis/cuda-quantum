/*******************************************************************************
 * Copyright (c) 2026 NVIDIA Corporation & Affiliates.                         *
 * All rights reserved.                                                        *
 *                                                                             *
 * This source code and the accompanying materials are made available under    *
 * the terms of the Apache License 2.0 which accompanies this distribution.    *
 ******************************************************************************/

// The routines in this file are a qubit-only, in-place replacement for the
// `qpp::apply` / `qpp::applyCTRL` routines of Quantum++
// (https://github.com/softwareqinc/qpp), which is distributed under the
// following license:
//
// Copyright (c) 2013 - 2022 softwareQ Inc. All rights reserved.
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace nvqir {
namespace details {

/// @brief Number of threads to use for a loop over a state vector of the
/// given dimension. Small states run serially, since waking the thread pool
/// costs more than the loop, and mid-sized ones use part of the pool. The
/// count depends only on the dimension so that every loop of one execution
/// uses the same team: resizing the team between loops is expensive.
inline int parallelThreads(std::size_t dimension) {
#if defined(_OPENMP)
  const auto maxThreads =
      static_cast<std::size_t>(std::max(1, omp_get_max_threads()));
  if (dimension < (1ULL << 15))
    return 1;
  return static_cast<int>(
      std::clamp<std::size_t>(dimension >> 13, 1, maxThreads));
#else
  return 1;
#endif
}

/// @brief Insert a 0 bit at each of the given (ascending) bit positions.
inline std::size_t insertZeroBits(std::size_t value,
                                  const std::vector<std::size_t> &sortedBits) {
  for (auto bit : sortedBits) {
    const std::size_t low = value & ((1ULL << bit) - 1);
    value = ((value >> bit) << (bit + 1)) | low;
  }
  return value;
}

/// @brief Complex product without the NaN/Inf handling of `operator*`, which
/// is much slower in the hot loops.
inline std::complex<double> multiply(const std::complex<double> &a,
                                     const std::complex<double> &b) {
  return {a.real() * b.real() - a.imag() * b.imag(),
          a.real() * b.imag() + a.imag() * b.real()};
}

/// @brief A 2x2 complex matrix split into real and imaginary parts.
struct SplitMatrix {
  double m00r, m00i, m01r, m01i, m10r, m10i, m11r, m11i;
  explicit SplitMatrix(const std::complex<double> *matrix)
      : m00r(matrix[0].real()), m00i(matrix[0].imag()), m01r(matrix[1].real()),
        m01i(matrix[1].imag()), m10r(matrix[2].real()), m10i(matrix[2].imag()),
        m11r(matrix[3].real()), m11i(matrix[3].imag()) {}
};

/// @brief Apply the matrix to one amplitude pair, given as (re, im) doubles.
/// The arithmetic is the same as `multiply`.
[[gnu::always_inline]] inline void transformPair(double *a, double *b,
                                                 const SplitMatrix &m) {
  const double ar = a[0], ai = a[1];
  const double br = b[0], bi = b[1];
  a[0] = (m.m00r * ar - m.m00i * ai) + (m.m01r * br - m.m01i * bi);
  a[1] = (m.m00r * ai + m.m00i * ar) + (m.m01r * bi + m.m01i * br);
  b[0] = (m.m10r * ar - m.m10i * ai) + (m.m11r * br - m.m11i * bi);
  b[1] = (m.m10r * ai + m.m10i * ar) + (m.m11r * bi + m.m11i * br);
}

/// @brief Apply a 2x2 matrix to the amplitude pairs (`low[j]`, `high[j]`).
/// Written on real and imaginary parts so the compiler can vectorize it.
[[gnu::always_inline]] inline void
transformPairs(std::complex<double> *low, std::complex<double> *high,
               std::size_t count, const std::complex<double> *matrix) {
  double *__restrict a = reinterpret_cast<double *>(low);
  double *__restrict b = reinterpret_cast<double *>(high);
  const SplitMatrix m(matrix);
  for (std::size_t j = 0; j < count; ++j)
    transformPair(a + 2 * j, b + 2 * j, m);
}

/// @brief Apply a 2x2 matrix to consecutive blocks of 2 * `Run` amplitudes,
/// pairing the first half of each block with its second half.
template <std::size_t Run>
[[gnu::always_inline]] inline void
transformBlocks(std::complex<double> *blocks, std::size_t numBlocks,
                const std::complex<double> *matrix) {
  double *__restrict data = reinterpret_cast<double *>(blocks);
  const SplitMatrix m(matrix);
  for (std::size_t block = 0; block < numBlocks; ++block) {
    double *a = data + 4 * Run * block;
    for (std::size_t j = 0; j < Run; ++j)
      transformPair(a + 2 * j, a + 2 * (Run + j), m);
  }
}

/// @brief Apply a 2x2 matrix to the pairs (`low[j * stride]`,
/// `high[j * stride]`).
inline void transformStrided(std::complex<double> *low,
                             std::complex<double> *high, std::size_t count,
                             std::size_t stride,
                             const std::complex<double> *matrix) {
  double *a = reinterpret_cast<double *>(low);
  double *b = reinterpret_cast<double *>(high);
  const SplitMatrix m(matrix);
  for (std::size_t j = 0; j < count; ++j)
    transformPair(a + 2 * stride * j, b + 2 * stride * j, m);
}

/// @brief X gate on strided pairs: swap the two amplitudes.
inline void swapStrided(std::complex<double> *low, std::complex<double> *high,
                        std::size_t count, std::size_t stride,
                        const std::complex<double> *) {
  for (std::size_t j = 0; j < count; ++j)
    std::swap(low[j * stride], high[j * stride]);
}

/// @brief A function that processes one run of amplitude pairs.
using PairKernel = void (*)(std::complex<double> *, std::complex<double> *,
                            std::size_t, const std::complex<double> *);

inline void transformPairsGeneric(std::complex<double> *low,
                                  std::complex<double> *high, std::size_t count,
                                  const std::complex<double> *matrix) {
  transformPairs(low, high, count, matrix);
}

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#define CUDAQ_QPP_HAS_AVX2_KERNEL
/// @brief True if the CPU supports AVX2.
inline bool hasAvx2() {
  static const bool supported = __builtin_cpu_supports("avx2");
  return supported;
}

/// @brief `transformPairs` compiled for AVX2. No FMA, so results match the
/// generic version exactly.
__attribute__((target("avx2"))) inline void
transformPairsAvx2(std::complex<double> *low, std::complex<double> *high,
                   std::size_t count, const std::complex<double> *matrix) {
  transformPairs(low, high, count, matrix);
}
#endif

/// @brief X gate: swap the two amplitudes.
inline void swapPairs(std::complex<double> *low, std::complex<double> *high,
                      std::size_t count, const std::complex<double> *) {
  std::swap_ranges(low, low + count, high);
}

/// @brief True if the 2x2 matrix is the X gate.
inline bool isXGate(const std::complex<double> *matrix) {
  const std::complex<double> zero = 0.0, one = 1.0;
  return matrix[0] == zero && matrix[1] == one && matrix[2] == one &&
         matrix[3] == zero;
}

/// @brief Pick the kernel for the given 2x2 matrix.
inline PairKernel selectPairKernel(const std::complex<double> *matrix) {
  if (isXGate(matrix))
    return swapPairs;
#if defined(CUDAQ_QPP_HAS_AVX2_KERNEL)
  if (hasAvx2())
    return transformPairsAvx2;
#endif
  return transformPairsGeneric;
}

/// @brief A function that processes consecutive blocks of amplitudes.
using BlockKernel = void (*)(std::complex<double> *, std::size_t,
                             const std::complex<double> *);

template <std::size_t Run>
inline void transformBlocksGeneric(std::complex<double> *blocks,
                                   std::size_t numBlocks,
                                   const std::complex<double> *matrix) {
  transformBlocks<Run>(blocks, numBlocks, matrix);
}

#if defined(CUDAQ_QPP_HAS_AVX2_KERNEL)
template <std::size_t Run>
__attribute__((target("avx2"))) inline void
transformBlocksAvx2(std::complex<double> *blocks, std::size_t numBlocks,
                    const std::complex<double> *matrix) {
  transformBlocks<Run>(blocks, numBlocks, matrix);
}
#endif

/// @brief Pick the block kernel for a target qubit below 4, whose blocks
/// hold 2^(target + 1) amplitudes.
inline BlockKernel selectBlockKernel(std::size_t target) {
#if defined(CUDAQ_QPP_HAS_AVX2_KERNEL)
  static constexpr BlockKernel avx2[] = {
      transformBlocksAvx2<1>, transformBlocksAvx2<2>, transformBlocksAvx2<4>,
      transformBlocksAvx2<8>};
  if (hasAvx2())
    return avx2[target];
#endif
  static constexpr BlockKernel generic[] = {
      transformBlocksGeneric<1>, transformBlocksGeneric<2>,
      transformBlocksGeneric<4>, transformBlocksGeneric<8>};
  return generic[target];
}

/// @brief Longest stretch of consecutive qubits that are not in
/// `sortedInvolved`, as (first qubit, number of qubits).
inline std::pair<std::size_t, std::size_t>
longestFreeStretch(const std::vector<std::size_t> &sortedInvolved,
                   std::size_t numQubits) {
  std::size_t bestBegin = 0, bestBits = 0, begin = 0;
  const auto consider = [&](std::size_t end) {
    if (end - begin > bestBits) {
      bestBegin = begin;
      bestBits = end - begin;
    }
  };
  for (auto qubit : sortedInvolved) {
    consider(qubit);
    begin = qubit + 1;
  }
  consider(numQubits);
  return {bestBegin, bestBits};
}

/// @brief Call `body(outer, first, last)` so that the calls cover every
/// `outer` in [0, numOuter) with [first, last) = [0, 2^innerBits), splitting
/// that index space evenly across the threads.
template <typename Body>
inline void forEachSegment(std::size_t numOuter, std::size_t innerBits,
                           [[maybe_unused]] int numThreads, Body &&body) {
  const std::size_t inner = 1ULL << innerBits;
  const std::size_t total = numOuter << innerBits;
  const auto process = [&](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end;) {
      const std::size_t first = i & (inner - 1);
      const std::size_t last = std::min(inner, first + (end - i));
      body(i >> innerBits, first, last);
      i += last - first;
    }
  };
#if defined(_OPENMP)
#pragma omp parallel num_threads(numThreads) if (numThreads > 1)
  {
    const std::size_t thread = omp_get_thread_num();
    const std::size_t team = omp_get_num_threads();
    process(total * thread / team, total * (thread + 1) / team);
  }
#else
  process(0, total);
#endif
}

/// @brief Apply a (controlled) gate to a state vector, in place.
///
/// Qubit `q` is bit `q` of the amplitude index. `matrix` is row-major, with
/// `targets[0]` as the most significant bit of its index. Only amplitudes
/// whose control bits are all 1 are touched.
inline void applyGateInPlace(std::complex<double> *state, std::size_t numQubits,
                             const std::vector<std::complex<double>> &matrix,
                             const std::vector<std::size_t> &controls,
                             const std::vector<std::size_t> &targets) {
  const std::size_t numTargets = targets.size();
  if (numTargets == 0)
    throw std::invalid_argument("[qpp] gate applied to no target qubits.");
  if (numTargets >= numQubits + 1)
    throw std::invalid_argument("[qpp] gate has more targets than qubits.");
  const std::size_t dim = 1ULL << numTargets;
  if (matrix.size() != dim * dim)
    throw std::invalid_argument(
        "[qpp] gate matrix size does not match the number of targets.");

  // Collect and validate the control and target qubits.
  std::size_t controlMask = 0, involvedMask = 0;
  std::vector<std::size_t> involved;
  involved.reserve(controls.size() + numTargets);
  const auto addQubit = [&](std::size_t qubit) {
    if (qubit >= numQubits)
      throw std::invalid_argument("[qpp] qubit index out of range.");
    const std::size_t bit = 1ULL << qubit;
    if (involvedMask & bit)
      throw std::invalid_argument(
          "[qpp] control and target qubits must be distinct.");
    involvedMask |= bit;
    involved.push_back(qubit);
  };
  for (auto qubit : controls) {
    addQubit(qubit);
    controlMask |= 1ULL << qubit;
  }
  for (auto qubit : targets)
    addQubit(qubit);
  std::sort(involved.begin(), involved.end());

  // One group of amplitudes per assignment of the remaining qubits.
  const std::size_t numGroups = 1ULL << (numQubits - involved.size());
  [[maybe_unused]] const int numThreads = parallelThreads(1ULL << numQubits);

  if (numTargets == 1) {
    const std::size_t targetBit = 1ULL << targets[0];

    const std::size_t target = targets[0];
    const std::size_t lowest = involved.front();

    // Amplitudes below the lowest involved qubit are contiguous: process
    // them one run at a time.
    const std::size_t runLength = 1ULL << lowest;
    if (lowest >= 4 || (lowest >= 2 && target != lowest)) {
      const PairKernel kernel = selectPairKernel(matrix.data());
      const auto processGroups = [&](std::size_t begin, std::size_t end) {
        for (std::size_t group = begin; group < end;) {
          const std::size_t count =
              std::min(runLength - (group & (runLength - 1)), end - group);
          const std::size_t i0 = insertZeroBits(group, involved) | controlMask;
          kernel(state + i0, state + (i0 | targetBit), count, matrix.data());
          group += count;
        }
      };
#if defined(_OPENMP)
#pragma omp parallel num_threads(numThreads) if (numThreads > 1)
      {
        const std::size_t thread = omp_get_thread_num();
        const std::size_t team = omp_get_num_threads();
        processGroups(numGroups * thread / team,
                      numGroups * (thread + 1) / team);
      }
#else
      processGroups(0, numGroups);
#endif
      return;
    }

    // Low target with free qubits right above it: the runs follow each other
    // in memory, so process whole blocks of them at once.
    if (target == lowest) {
      const std::size_t next = involved.size() > 1 ? involved[1] : numQubits;
      const std::size_t blockBits = next - target - 1;
      if (blockBits >= 2) {
        // Everything below `next` is inside the blocks.
        std::vector<std::size_t> outerBits;
        for (std::size_t qubit = 0; qubit < next; ++qubit)
          outerBits.push_back(qubit);
        outerBits.insert(outerBits.end(), involved.begin() + 1, involved.end());
        const BlockKernel kernel = selectBlockKernel(target);
        forEachSegment(
            1ULL << (numQubits - outerBits.size()), blockBits, numThreads,
            [&](std::size_t outer, std::size_t first, std::size_t last) {
              const std::size_t base =
                  insertZeroBits(outer, outerBits) | controlMask;
              kernel(state + base + (first << (target + 1)), last - first,
                     matrix.data());
            });
        return;
      }
    }

    // Otherwise walk the longest stretch of free qubits with a fixed stride.
    if (const auto [stretchBegin, stretchBits] =
            longestFreeStretch(involved, numQubits);
        stretchBits >= 2) {
      std::vector<std::size_t> outerBits = involved;
      for (std::size_t i = 0; i < stretchBits; ++i)
        outerBits.push_back(stretchBegin + i);
      std::sort(outerBits.begin(), outerBits.end());
      const auto kernel =
          isXGate(matrix.data()) ? swapStrided : transformStrided;
      const std::size_t stretchStart = stretchBegin;
      forEachSegment(
          1ULL << (numQubits - outerBits.size()), stretchBits, numThreads,
          [&](std::size_t outer, std::size_t first, std::size_t last) {
            const std::size_t i0 =
                (insertZeroBits(outer, outerBits) | controlMask) +
                (first << stretchStart);
            kernel(state + i0, state + i0 + targetBit, last - first,
                   1ULL << stretchStart, matrix.data());
          });
      return;
    }

    // Too few qubits for either: one pair at a time.
    const auto m00 = matrix[0], m01 = matrix[1], m10 = matrix[2],
               m11 = matrix[3];
#if defined(_OPENMP)
#pragma omp parallel for num_threads(numThreads) if (numThreads > 1)
#endif
    for (std::size_t i = 0; i < numGroups; ++i) {
      const std::size_t i0 = insertZeroBits(i, involved) | controlMask;
      const std::size_t i1 = i0 | targetBit;
      const auto a0 = state[i0], a1 = state[i1];
      state[i0] = multiply(m00, a0) + multiply(m01, a1);
      state[i1] = multiply(m10, a0) + multiply(m11, a1);
    }
    return;
  }

  // Index offset of each target basis state within a group.
  std::vector<std::size_t> offsets(dim, 0);
  for (std::size_t m = 0; m < dim; ++m)
    for (std::size_t j = 0; j < numTargets; ++j)
      if (m & (1ULL << (numTargets - 1 - j)))
        offsets[m] |= 1ULL << targets[j];

#if defined(_OPENMP)
#pragma omp parallel num_threads(numThreads) if (numThreads > 1)
#endif
  {
    std::vector<std::complex<double>> amplitudes(dim);
#if defined(_OPENMP)
#pragma omp for
#endif
    for (std::size_t i = 0; i < numGroups; ++i) {
      const std::size_t base = insertZeroBits(i, involved) | controlMask;
      for (std::size_t m = 0; m < dim; ++m)
        amplitudes[m] = state[base | offsets[m]];
      for (std::size_t row = 0; row < dim; ++row) {
        std::complex<double> sum = 0.0;
        const auto *matrixRow = matrix.data() + row * dim;
        for (std::size_t col = 0; col < dim; ++col)
          sum += multiply(matrixRow[col], amplitudes[col]);
        state[base | offsets[row]] = sum;
      }
    }
  }
}

} // namespace details

} // namespace nvqir
