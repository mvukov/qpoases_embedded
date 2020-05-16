/*
 *  This file is part of qpOASES.
 *
 *  qpOASES -- An Implementation of the Online Active Set Strategy.
 *  Copyright (C) 2020 Milan Vukov.
 *
 *  qpOASES is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  qpOASES is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with qpOASES; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301
 *  USA
 *
 */

#ifndef TEST_HELPERS_H_
#define TEST_HELPERS_H_

#include <sstream>
#include <vector>

#include <gtest/gtest.h>

#include "qpoases_embedded/Constants.hpp"

#define EXPECT_STL_VECTORS_EQ(lhs, rhs) \
  EXPECT_TRUE(qpoases_embedded::testing::CompareStlVectors(lhs, rhs));

#define EXPECT_STL_VECTORS_NEAR(lhs, rhs, tolerance) \
  EXPECT_TRUE(                                       \
      qpoases_embedded::testing::CompareStlVectors(lhs, rhs, tolerance));

namespace qpoases_embedded {
namespace testing {

template <typename T>
::testing::AssertionResult CompareStlVectors(const std::vector<T>& lhs,
                                             const std::vector<T>& rhs,
                                             T tolerance) {
  if (lhs.size() != rhs.size()) {
    return ::testing::AssertionFailure() << "Vector sizes must be the same!";
  }

  bool success = true;
  std::stringstream messages;
  for (size_t el = 0; el < lhs.size(); ++el) {
    const T diff = lhs[el] - rhs[el];
    if (std::abs(diff) > tolerance) {
      messages << "Vector elements at " << el << ": lhs = " << lhs[el]
               << ", rhs = " << rhs[el] << ", abs(diff) = " << std::abs(diff)
               << " > " << tolerance << " (tolerance)!" << std::endl;
      success = false;
    }
  }
  if (success) {
    return ::testing::AssertionSuccess();
  }
  return ::testing::AssertionFailure() << std::endl << messages.str();
}

template <typename T>
::testing::AssertionResult CompareStlVectors(const std::vector<T>& lhs,
                                             const std::vector<T>& rhs) {
  return CompareStlVectors(lhs, rhs, T(0));
}

}  // namespace testing

struct QpTestData {
  int num_variables;
  int num_constraints;
  std::vector<real_t> h;
  std::vector<real_t> g;
  std::vector<real_t> a;
  std::vector<real_t> lb;
  std::vector<real_t> ub;
  std::vector<real_t> lba;
  std::vector<real_t> uba;
  std::vector<real_t> x_opt;
  std::vector<real_t> y_opt;
  real_t f_opt;
};

static constexpr real_t inf = INFTY;

}  // namespace qpoases_embedded

#endif  // TEST_HELPERS_H_
