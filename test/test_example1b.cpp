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

#include <qpoases_embedded/QProblemB.hpp>
#include <qpoases_embedded/Utils.hpp>

namespace qpoases_embedded {

static constexpr real_t kEqualTolerance = 1e-12;

/** Example for qpOASES using the QProblemB class. */
TEST(TestQpOases, Example1b) {
  /* Setting up QProblemB object. */
  QProblemB qp(2);

  auto* H = qp.getMutableH();
  *H = {1.0, 0.0, 0.0, 0.5};

  auto* g = qp.getMutableG();
  *g = {1.5, 1.0};

  auto* lb = qp.getMutableLb();
  *lb = {0.5, -2.0};

  auto* ub = qp.getMutableUb();
  *ub = {5.0, 2.0};

  /* Solve first QP. */
  int nWSR = 10;
  ASSERT_EQ(SUCCESSFUL_RETURN, qp.init(nWSR, printIterationB));
  EXPECT_NEAR(-0.125, qp.getObjVal(), kEqualTolerance);

  const auto& x = qp.getX();
  const auto& y = qp.getY();

  EXPECT_NEAR(0.5, x[0], kEqualTolerance);
  EXPECT_NEAR(-2, x[1], kEqualTolerance);

  EXPECT_NEAR(2, y[0], kEqualTolerance);
  EXPECT_NEAR(0, y[1], kEqualTolerance);

  /* Setup data of second QP. */
  *g = {1.0, 1.5};
  *lb = {0.0, -1.0};
  *ub = {5.0, -0.5};

  /* Solve second QP. */
  nWSR = 10;
  ASSERT_EQ(SUCCESSFUL_RETURN, qp.hotstart(nWSR, printIterationB));
  EXPECT_NEAR(-1.25, qp.getObjVal(), kEqualTolerance);

  EXPECT_NEAR(0, x[0], kEqualTolerance);
  EXPECT_NEAR(-1, x[1], kEqualTolerance);

  EXPECT_NEAR(1, y[0], kEqualTolerance);
  EXPECT_NEAR(1, y[1], kEqualTolerance);
}

}  // namespace qpoases_embedded
