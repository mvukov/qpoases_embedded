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

#include <cmath>

/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

namespace qpoases_embedded {

/*
 *  g e t N F R
 */
inline int QProblemB::getNFR() { return bounds.getNFR(); }

/*
 *  g e t N F X
 */
inline int QProblemB::getNFX() { return bounds.getNFX(); }

/*
 *  g e t N F V
 */
inline int QProblemB::getNFV() const { return bounds.getNFV(); }

/*
 *  g e t S t a t u s
 */
inline QProblemStatus QProblemB::getStatus() const { return status; }

/*
 *  i s I n i t i a l i s e d
 */
inline bool QProblemB::isInitialised() const {
  if (status == QPS_NOTINITIALISED)
    return false;
  else
    return true;
}

/*
 *  i s S o l v e d
 */
inline bool QProblemB::isSolved() const {
  if ((status == QPS_AUXILIARYQPSOLVED) || (status == QPS_HOMOTOPYQPSOLVED) ||
      (status == QPS_SOLVED))
    return true;
  else
    return false;
}

/*
 *  i s I n f e a s i b l e
 */
inline bool QProblemB::isInfeasible() const { return infeasible; }

/*
 *  i s U n b o u n d e d
 */
inline bool QProblemB::isUnbounded() const { return unbounded; }

/*
 *  g e t H e s s i a n T y p e
 */
inline HessianType QProblemB::getHessianType() const { return hessianType; }

/*
 *  s e t H e s s i a n T y p e
 */
inline returnValue QProblemB::setHessianType(HessianType _hessianType) {
  hessianType = _hessianType;
  return SUCCESSFUL_RETURN;
}

/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *  c o m p u t e G i v e n s
 */
inline void QProblemB::computeGivens(real_t xold, real_t yold, real_t& xnew,
                                     real_t& ynew, real_t& c, real_t& s) const {
  if (std::fabs(yold) <= ZERO) {
    c = 1.0;
    s = 0.0;

    xnew = xold;
    ynew = yold;
  } else {
    real_t t, mu;

    mu = std::fabs(xold);
    if (std::fabs(yold) > mu) mu = std::fabs(yold);

    t = mu * sqrt((xold / mu) * (xold / mu) + (yold / mu) * (yold / mu));

    if (xold < 0.0) t = -t;

    c = xold / t;
    s = yold / t;
    xnew = t;
    ynew = 0.0;
  }

  return;
}

/*
 *  a p p l y G i v e n s
 */
inline void QProblemB::applyGivens(real_t c, real_t s, real_t xold, real_t yold,
                                   real_t& xnew, real_t& ynew) const {
  /* Usual Givens plane rotation requiring four multiplications. */
  xnew = c * xold + s * yold;
  ynew = -s * xold + c * yold;
  //   double nu = s/(1.0+c);
  //
  //   xnew = xold*c + yold*s;
  //   ynew = (xnew+xold)*nu - yold;

  return;
}

}  // namespace qpoases_embedded

/*
 *  end of file
 */
