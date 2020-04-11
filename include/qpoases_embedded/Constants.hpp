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

#ifndef QPOASES_CONSTANTS_HPP
#define QPOASES_CONSTANTS_HPP

#include <cstddef>
#include <limits>

namespace qpoases_embedded {

#ifdef QPOASES_SINGLE_PRECISION
using real_t = float;
#else
using real_t = double;
#endif  // QPOASES_SINGLE_PRECISION

/** Desired KKT tolerance of QP solution; a warning RET_INACCURATE_SOLUTION is
 *  issued if this tolerance is not met.
 *  Note: this value has to be positive! */
static constexpr real_t DESIREDACCURACY = 1.0e-3;

/** Critical KKT tolerance of QP solution; an error is issued if this
 *  tolerance is not met.
 *  Note: this value has to be positive! */
static constexpr real_t CRITICALACCURACY = 1.0e-2;

/** Numerical value of machine precision (min eps, s.t. 1+eps > 1).
 *  Note: this value has to be positive! */
static constexpr real_t EPS = std::numeric_limits<real_t>::epsilon();

/** Numerical value of zero (for situations in which it would be
 *  unreasonable to compare with 0.0).
 *  Note: this value has to be positive! */
static constexpr real_t ZERO = 1.0e-50;

/** Numerical value of infinity (e.g. for non-existing bounds).
 *  Note: this value has to be positive! */
static constexpr real_t INFTY = 1.0e12;

/** Lower/upper (constraints') bound tolerance (an inequality constraint
 *  whose lower and upper bound differ by less than BOUNDTOL is regarded
 *  to be an equality constraint).
 *  Note: this value has to be positive! */
static constexpr real_t BOUNDTOL = 1.0e-10;

/** Offset for relaxing (constraints') bounds at beginning of an initial
 *  homotopy.
 *  Note: this value has to be positive! */
static constexpr real_t BOUNDRELAXATION = 1.0e3;

/** Factor that determines physical lengths of index lists.
 *  Note: this value has to be greater than 1! */
static constexpr size_t INDEXLISTFACTOR = 5;

}  // namespace qpoases_embedded

#endif /* QPOASES_CONSTANTS_HPP */

/*
 *  end of file
 */
