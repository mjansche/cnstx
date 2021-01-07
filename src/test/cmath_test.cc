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

#include "cnstx/cmath"

#include <gtest/gtest.h>

#include <cmath>
#include <ios>
#include <limits>
#include <tuple>

namespace {

template <typename T, int value>
void frexp_test() {
  {
    constexpr T x = value;
    constexpr auto y_exp = cnstx::frexp(x);
    constexpr T y = std::get<0>(y_exp);
    constexpr int exp = std::get<1>(y_exp);

    int std_exp;
    T std_y = std::frexp(x, &std_exp);

    EXPECT_EQ(std_y, y);
    EXPECT_EQ(std_exp, exp);
  }
  {
    constexpr T x = static_cast<T>(value) / 16;
    constexpr auto y_exp = cnstx::frexp(x);
    constexpr T y = std::get<0>(y_exp);
    constexpr int exp = std::get<1>(y_exp);

    int std_exp;
    T std_y = std::frexp(x, &std_exp);

    EXPECT_EQ(std_y, y);
    EXPECT_EQ(std_exp, exp);
  }
}

TEST(CMathTest, frexp) {
  frexp_test<float, +3>();
  frexp_test<float, -3>();
  frexp_test<double, +3>();
  frexp_test<double, -3>();
}

template <typename T, int value, int exp>
void ldexp_test() {
  constexpr T x = value;
  constexpr T y = cnstx::ldexp(x, exp);
  T std_y = std::ldexp(x, exp);
  EXPECT_EQ(std_y, y);
}

TEST(CMathTest, ldexp) {
  ldexp_test<float, +3, +20>();
  ldexp_test<float, -3, +20>();
  ldexp_test<float, +3, -20>();
  ldexp_test<float, -3, -20>();
  ldexp_test<float, +3, +50>();
  ldexp_test<float, -3, +50>();
  ldexp_test<double, +3, +20>();
  ldexp_test<double, -3, +20>();
  ldexp_test<double, +3, -20>();
  ldexp_test<double, -3, -20>();
  ldexp_test<double, +3, +333>();
  ldexp_test<double, -3, +333>();
}

template <int num, int den, int exp>
void div_mod_test() {
  constexpr double x = num;
  constexpr double y = den;
  constexpr auto xx = cnstx::internal::fraction64(x);
  constexpr auto yy = cnstx::internal::fraction64(y);
  constexpr auto zz_rem = cnstx::internal::div_mod(xx, yy);
  constexpr auto zz = std::get<0>(zz_rem);
  constexpr double z = cnstx::ldexp(zz, exp);
  EXPECT_EQ(x / y, z);
}

TEST(CMathTest, div_mod) {
  div_mod_test<11, 13, -64>();
  div_mod_test<13, 11, -63>();
  div_mod_test<13, 22, -64>();
}

#define SQRT_TEST(T, value)                                        \
  do {                                                             \
    constexpr T x = value;                                         \
    constexpr T y = cnstx::sqrt(x);                                \
    T z = std::sqrt(x);                                            \
    EXPECT_EQ(z, y) << "Argument: " << std::hexfloat << x          \
                    << "\nExpected: " << z << "\nObserved: " << y; \
  } while (false)

TEST(CMathTest, sqrt) {
  SQRT_TEST(double, 0);
  SQRT_TEST(double, 1);
  SQRT_TEST(double, 2);
  SQRT_TEST(double, 3);
  SQRT_TEST(double, 4);
  SQRT_TEST(double, 5);
  SQRT_TEST(double, 6);
  SQRT_TEST(double, 7);
  SQRT_TEST(double, 8);
  SQRT_TEST(double, 9);
  SQRT_TEST(double, 15);
  SQRT_TEST(double, 16);
  SQRT_TEST(double, 17);
  SQRT_TEST(double, 31);
  SQRT_TEST(double, 32);
  SQRT_TEST(double, 33);
  SQRT_TEST(double, 63);
  SQRT_TEST(double, 64);
  SQRT_TEST(double, 65);
  if constexpr (std::numeric_limits<long double>::digits <= 64) {
    SQRT_TEST(long double, 3);
    SQRT_TEST(long double, 5);
    SQRT_TEST(long double, 15.0L);
    SQRT_TEST(long double, 17.0L);
    SQRT_TEST(long double, 184467440737095516165.0L);
    SQRT_TEST(long double, 184467440737095516167.0L);
    SQRT_TEST(long double, 0xf.ffffffffffffff0p+124L);
    SQRT_TEST(long double, 0xf.ffffffffffffff8p+124L);
    SQRT_TEST(long double, 0xf.ffffffffffffffcp+124L);
    SQRT_TEST(long double, 0xf.ffffffffffffffdp+124L);
    SQRT_TEST(long double, 0xf.ffffffffffffffep+124L);
    SQRT_TEST(long double, 0xf.fffffffffffffffp+124L);
  }
}

}  // namespace
