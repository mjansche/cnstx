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

#include <cfloat>
#include <cmath>
#include <ios>
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
#if LDBL_MANT_DIG <= 64
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
#endif
}

template <typename T, int n>
void log1p2x_test() {
  constexpr T y = cnstx::internal::log1p2x<T>(n);
  T z = std::log1p(std::ldexp(T{1}, -n));
  EXPECT_EQ(z, y) << "Argument: " << n << "\nExpected: " << std::hexfloat << z
                  << "\nObserved: " << y;
}

TEST(CMathTest, log1p2x) {
  log1p2x_test<float, 0>();
  log1p2x_test<float, 1>();
  log1p2x_test<float, 2>();
  log1p2x_test<float, 4>();
  log1p2x_test<float, 8>();
  log1p2x_test<float, 16>();
  log1p2x_test<float, 32>();
  log1p2x_test<float, 64>();
  log1p2x_test<float, 96>();
  log1p2x_test<float, 112>();
  log1p2x_test<float, 120>();
  log1p2x_test<float, 124>();
  log1p2x_test<float, 125>();

  log1p2x_test<double, 0>();
  log1p2x_test<double, 1>();
  log1p2x_test<double, 2>();
  log1p2x_test<double, 4>();
  log1p2x_test<double, 8>();
  log1p2x_test<double, 16>();
  log1p2x_test<double, 32>();
  log1p2x_test<double, 64>();
  log1p2x_test<double, 128>();
  log1p2x_test<double, 256>();
  log1p2x_test<double, 512>();
  log1p2x_test<double, 768>();
  log1p2x_test<double, 896>();
  log1p2x_test<double, 960>();
  log1p2x_test<double, 992>();
  log1p2x_test<double, 1008>();
  log1p2x_test<double, 1016>();
  log1p2x_test<double, 1020>();
  log1p2x_test<double, 1021>();

  log1p2x_test<long double, 0>();
  log1p2x_test<long double, 1>();
  log1p2x_test<long double, 2>();
  log1p2x_test<long double, 4>();
  log1p2x_test<long double, 8>();
  log1p2x_test<long double, 16>();
  log1p2x_test<long double, 32>();
  log1p2x_test<long double, 64>();
  log1p2x_test<long double, 128>();
  log1p2x_test<long double, 256>();
  log1p2x_test<long double, 512>();
#if LDBL_MANT_DIG > 53
  log1p2x_test<long double, 1024>();
  log1p2x_test<long double, 2048>();
  log1p2x_test<long double, 4096>();
  log1p2x_test<long double, 8192>();
  log1p2x_test<long double, 12288>();
  log1p2x_test<long double, 14336>();
  log1p2x_test<long double, 15360>();
  log1p2x_test<long double, 15872>();
  log1p2x_test<long double, 16128>();
  log1p2x_test<long double, 16256>();
  log1p2x_test<long double, 16320>();
  log1p2x_test<long double, 16352>();
  log1p2x_test<long double, 16368>();
  log1p2x_test<long double, 16376>();
  log1p2x_test<long double, 16380>();
  log1p2x_test<long double, 16381>();
#endif
}

}  // namespace
