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
#include <vector>

#include <gtest/gtest.h>

#include <qpoases_embedded/MessageHandling.hpp>
#include <qpoases_embedded/QProblem.hpp>
#include <qpoases_embedded/Utils.hpp>

#include "./test_helpers.h"

namespace qpoases_embedded {

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

#include "qpoases_embedded/hanging_chain_test_data.h"

static constexpr real_t kEqualXTolerance = 1.5e-13;  // for primal variables
static constexpr real_t kEqualYTolerance = 1.1e-9;   // for dual variables

class TestHangingChain : public ::testing::TestWithParam<QpTestData> {};

TEST_P(TestHangingChain, ProcessTestData) {
  const auto& test_data = GetParam();
  ASSERT_GT(test_data.num_variables, 0);
  ASSERT_GE(test_data.num_constraints, 0);
  int nWSR = 100;

  QProblem qp(test_data.num_variables, test_data.num_constraints);
  *qp.getMutableH() = test_data.h;
  *qp.getMutableG() = test_data.g;
  *qp.getMutableA() = test_data.a;
  *qp.getMutableLb() = test_data.lb;
  *qp.getMutableUb() = test_data.ub;
  *qp.getMutableLbA() = test_data.lba;
  *qp.getMutableUbA() = test_data.uba;
  const auto success = qp.init(nWSR, printIteration);
  EXPECT_EQ(SUCCESSFUL_RETURN, success) << getErrorString(success);
  EXPECT_STL_VECTORS_NEAR(test_data.x_opt, qp.getX(), kEqualXTolerance);
  EXPECT_STL_VECTORS_NEAR(test_data.y_opt, qp.getY(), kEqualYTolerance);

  if (test_data.num_constraints == 0) {
    nWSR = 100;
    QProblemB qp(test_data.num_variables);
    *qp.getMutableH() = test_data.h;
    *qp.getMutableG() = test_data.g;
    *qp.getMutableLb() = test_data.lb;
    *qp.getMutableUb() = test_data.ub;
    const auto success = qp.init(nWSR, printIterationB);
    EXPECT_EQ(SUCCESSFUL_RETURN, success) << getErrorString(success);
    EXPECT_STL_VECTORS_NEAR(test_data.x_opt, qp.getX(), kEqualXTolerance);
    EXPECT_STL_VECTORS_NEAR(test_data.y_opt, qp.getY(), kEqualYTolerance);
  }
}

#ifdef INSTANTIATE_TEST_SUITE_P
INSTANTIATE_TEST_SUITE_P(HangingChainTests, TestHangingChain,
                         ::testing::ValuesIn(qp_test_data_vectors));
#else
INSTANTIATE_TEST_CASE_P(HangingChainTests, TestHangingChain,
                        ::testing::ValuesIn(qp_test_data_vectors), );
#endif  // INSTANTIATE_TEST_SUITE_P

}  // namespace qpoases_embedded
