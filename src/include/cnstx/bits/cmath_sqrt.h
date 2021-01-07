// -*- C++ -*-

// Copyright 2020, 2021 the cnstx authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Square root and friends.

#ifndef CNSTX_BITS_CMATH_SQRT_H_
#define CNSTX_BITS_CMATH_SQRT_H_

#include <cnstx/bits/cmath_float.h>

#include <cstdint>
#include <limits>
#include <tuple>

namespace cnstx {
namespace internal {

template <int max_iter = 8>
constexpr std::uint64_t sqrt_newton(std::uint64_t xx) {
  // Initial approximation of the square root by 0.5*x + 0.5:
  std::uint64_t yy = (xx >> 1) | (1UL << 63);
  for (int i = 0; i < max_iter; ++i) {
    std::uint64_t qq = std::get<0>(div_mod(xx, yy));
    std::uint64_t yy_least = yy & 1;
    std::uint64_t qq_least = qq & 1;
    yy = (yy >> 1) + (qq >> 1);
    if (yy_least & qq_least) {
      yy += 1;
    } else if (yy_least | qq_least) {
      yy += yy & 1;
    }
  }
  return yy;
}

void ignore_this_it_is_only_present_for_clang_format(int *x);

}  // namespace internal

template <typename T,
          typename std::enable_if<std::is_floating_point<T>::value &&
                                  std::numeric_limits<T>::digits <= 64>::type
              * = nullptr>
constexpr T sqrt(T x) {
  if (x < 0) return std::numeric_limits<T>::quiet_NaN();
  if (x == 0 || isinf(x) || isnan(x)) return x;
  auto norm_exp = frexp(x);
  T norm = std::get<0>(norm_exp);
  int exp = std::get<1>(norm_exp);
  std::uint64_t xx = static_cast<std::uint64_t>(ldexp(norm, 64));
  std::uint64_t yy = ~0ULL;
  if (xx < yy - 1) {
    yy = internal::sqrt_newton(xx);
  }
  if (exp & 1) {
    // If the initial strength reduction was by an odd power of 2, we
    // have to divide by the square root of 1/2.
    constexpr std::uint64_t sqrt_half = 0xb504f333f9de6484ULL;  // sqrt(0.5)
    yy = std::get<0>(internal::div_mod(yy, sqrt_half));
    ++exp;
  }
  T y = yy;
  return ldexp(y, (exp >> 1) - 64);
}

}  // namespace cnstx

#endif  // CNSTX_BITS_CMATH_SQRT_H_
