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

// Basic floating-point operations.

#ifndef CNSTX_BITS_CMATH_FLOAT_H_
#define CNSTX_BITS_CMATH_FLOAT_H_

#include <cmath>
#include <tuple>
#include <type_traits>

namespace cnstx {

template <typename T>
struct float_values {};

template <>
struct float_values<double> {
  constexpr static double huge_val = HUGE_VAL;
  constexpr static double e = 2.7182818284590452354;
  constexpr static double log2_e = 1.4426950408889634074;
  constexpr static double log2_10 = 3.3219280948873623479;
  constexpr static double ln_2 = 0.69314718055994530942;
  constexpr static double ln_10 = 2.30258509299404568402;
  constexpr static double log10_2 = 0.30102999566398119521;
  constexpr static double log10_e = 0.43429448190325182765;
  constexpr static double pi = 3.14159265358979323846;
  constexpr static double pi_inv = 0.31830988618379067154;
  constexpr static double sqrt_pi_inv = 1.12837916709551257390;
  constexpr static double sqrt_2 = 1.41421356237309504880;
};

template <>
struct float_values<float> {
  constexpr static float huge_val = HUGE_VALF;
  constexpr static float e = float_values<double>::e;
  constexpr static float pi = float_values<double>::pi;
  constexpr static float sqrt_2 = float_values<double>::sqrt_2;
};

template <>
struct float_values<long double> {
  constexpr static long double huge_val = HUGE_VALL;
  constexpr static long double e = 2.718281828459045235360287471352662498L;
  constexpr static long double pi = 3.141592653589793238462643383279502884L;
  constexpr static long double sqrt_2 = 1.414213562373095048801688724209698079L;
};

template <typename T, typename std::enable_if<
                          std::is_floating_point<T>::value>::type * = nullptr>
constexpr bool isnan(T x) {
  return x != x;
}

template <typename T, typename std::enable_if<
                          std::is_floating_point<T>::value>::type * = nullptr>
constexpr bool isinf(T x) {
  return isnan(x - x);
}

template <typename T, typename std::enable_if<
                          std::is_floating_point<T>::value>::type * = nullptr>
constexpr std::tuple<T, int> frexp(T x) {
  int exp = 0;
  if (x == 0 || isnan(x) || isinf(x)) {
    return std::make_tuple(x, exp);
  }
  while (x >= 1 || x <= -1) {
    x /= 2;
    ++exp;
  }
  while (x < 0.5 && x > -0.5) {
    x *= 2;
    --exp;
  }
  return std::make_tuple(x, exp);
}

constexpr float frexpf(float x, int *exp) {
  auto [y, e] = frexp<float>(x);
  *exp = e;
  return y;
}

constexpr double frexp(double x, int *exp) {
  auto [y, e] = frexp<double>(x);
  *exp = e;
  return y;
}

constexpr long double frexpl(long double x, int *exp) {
  auto [y, e] = frexp<long double>(x);
  *exp = e;
  return y;
}

template <typename T,
          typename std::enable_if<std::is_integral<T>::value>::type * = nullptr>
constexpr double frexp(T x, int *exp) {
  return frexp(static_cast<double>(x), exp);
}

template <typename T, typename std::enable_if<
                          std::is_floating_point<T>::value>::type * = nullptr>
constexpr T ldexp(T x, int exp) {
  if (x == 0 || isnan(x) || isinf(x)) {
    return x;
  }
  if (exp > 0) {
    while (exp--) x *= 2;
    if (isinf(x)) x = float_values<T>::huge_val;
    return x;
  }
  while (exp++) x /= 2;
  return x;
}

constexpr float ldexpf(float x, int exp) { return ldexp<float>(x, exp); }

constexpr long double ldexpl(long double x, int exp) {
  return ldexp<long double>(x, exp);
}

template <typename T,
          typename std::enable_if<std::is_integral<T>::value>::type * = nullptr>
constexpr double ldexp(T x, int exp) {
  return ldexp<double>(static_cast<double>(x), exp);
}

}  // namespace cnstx

#endif  // CNSTX_BITS_CMATH_FLOAT_H_
