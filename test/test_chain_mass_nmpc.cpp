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

#include "qpoases_embedded/chain_mass_nmpc_test_data.h"

namespace qpoases_embedded {

static constexpr real_t kEqualXTolerance = 1e-12;  // for primal variables
static constexpr real_t kEqualYTolerance = 1e-11;  // for dual variables
static constexpr real_t kEqualObjTolerance = 1e-2;

class TestChainMassNmpc : public ::testing::TestWithParam<QpTestData> {};

TEST_P(TestChainMassNmpc, ProcessTestData) {
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
  EXPECT_NEAR(test_data.f_opt, qp.getObjVal(), kEqualObjTolerance);
}

#ifdef INSTANTIATE_TEST_SUITE_P
INSTANTIATE_TEST_SUITE_P(ChainMassNmpcTests, TestChainMassNmpc,
                         ::testing::ValuesIn(qp_test_data_vectors));
#else
INSTANTIATE_TEST_CASE_P(ChainMassNmpcTests, TestChainMassNmpc,
                        ::testing::ValuesIn(qp_test_data_vectors), );
#endif  // INSTANTIATE_TEST_SUITE_P

}  // namespace qpoases_embedded
