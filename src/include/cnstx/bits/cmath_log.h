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

// Logarithms.

#ifndef CNSTX_BITS_CMATH_LOG_H_
#define CNSTX_BITS_CMATH_LOG_H_

#include <cnstx/bits/cmath_float.h>

namespace cnstx {
namespace internal {

// Computes log(1 + 2**-n) via its Taylor series.
template <typename T>
constexpr T log1p2x(int n) {
  if (n == 0) {
    return 0.6931471805599453094172321214581765681L;
  }
  // Compensated Kahan summation of Taylor series:
  T sum = 0;
  T c = 0;
  for (int i = 1; i < 130; i += 2) {
    T y_i = cnstx::ldexp(T{1} / i, -n * i);
    y_i -= cnstx::ldexp(T{1} / (i + 1), -n * (i + 1));
    y_i -= c;
    T t = sum + y_i;
    c = (t - sum) - y_i;
    if (t == sum) {
      break;
    }
    sum = t;
  }
  return sum;
}

}  // namespace internal
}  // namespace cnstx

#endif  // CNSTX_BITS_CMATH_LOG_H_
