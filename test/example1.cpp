/*
 *  This file is part of qpOASES.
 *
 *  qpOASES -- An Implementation of the Online Active Set Strategy.
 *  Copyright (C) 2007-2008 by Hans Joachim Ferreau et al. All rights
 *  reserved.
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

#include <gtest/gtest.h>

#include <qpoases_embedded/QProblem.hpp>
#include <qpoases_embedded/Utils.hpp>

namespace qpoases_embedded {

static constexpr real_t kEqualTolerance = 1e-12;

/** Example for qpOASES using the QProblem class. */
TEST(TestQpOases, Example1) {
  /* Setting up QProblem object. */
  QProblem qp(2, 1);

  /* Setup data of first QP. */
  auto* H = qp.getMutableH();
  *H = {1.0, 0.0, 0.0, 0.5};

  auto* A = qp.getMutableA();
  *A = {1.0, 1.0};

  auto* g = qp.getMutableG();
  *g = {1.5, 1.0};

  auto* lb = qp.getMutableLb();
  *lb = {0.5, -2.0};

  auto* ub = qp.getMutableUb();
  *ub = {5.0, 2.0};

  auto* lbA = qp.getMutableLbA();
  *lbA = {-1.0};

  auto* ubA = qp.getMutableUb();
  *ubA = {2.0};

  /* Solve first QP. */
  int nWSR = 10;
  ASSERT_EQ(SUCCESSFUL_RETURN, qp.init(nWSR, printIteration));
  EXPECT_NEAR(-0.0625, qp.getObjVal(), kEqualTolerance);

  const auto& x = qp.getX();
  const auto& y = qp.getY();

  EXPECT_NEAR(0.5, x[0], kEqualTolerance);
  EXPECT_NEAR(-1.5, x[1], kEqualTolerance);

  EXPECT_NEAR(1.75, y[0], kEqualTolerance);
  EXPECT_NEAR(0, y[1], kEqualTolerance);
  EXPECT_NEAR(0.25, y[2], kEqualTolerance);

  /* Setup data of second QP. */
  *g = {1.0, 1.5};
  *lb = {0.0, -1.0};
  *ub = {5.0, -0.5};
  *lbA = {-2.0};
  *ubA = {1.0};

  /* Solve second QP. */
  nWSR = 10;
  ASSERT_EQ(SUCCESSFUL_RETURN, qp.hotstart(nWSR, printIteration));
  EXPECT_NEAR(-1.25, qp.getObjVal(), kEqualTolerance);

  EXPECT_NEAR(0, x[0], kEqualTolerance);
  EXPECT_NEAR(-1, x[1], kEqualTolerance);

  EXPECT_NEAR(1, y[0], kEqualTolerance);
  EXPECT_NEAR(1, y[1], kEqualTolerance);
  EXPECT_NEAR(0, y[2], kEqualTolerance);
}

}  // namespace qpoases_embedded

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

/*
 *  end of file
 */
