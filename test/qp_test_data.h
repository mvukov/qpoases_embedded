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
#include "qpoases_embedded/Constants.hpp"

#ifndef QP_TEST_DATA_H_
#define QP_TEST_DATA_H_

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

#endif  // QP_TEST_DATA_H_

}  // namespace qpoases_embedded
