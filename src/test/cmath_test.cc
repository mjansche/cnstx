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

}  // namespace
