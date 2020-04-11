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

#include <qpoases_embedded/Constraints.hpp>
#include <qpoases_embedded/MessageHandling.hpp>

/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

namespace qpoases_embedded {

/*
 *  C o n s t r a i n t s
 */
Constraints::Constraints(int n, int indexlist_size)
    : SubjectTo(n),
      nC(n),
      nEC(0),
      nIC(0),
      nUC(0),
      active(indexlist_size),
      inactive(indexlist_size) {}

void Constraints::reset() {
  SubjectTo::reset();

  nC = getSize();
  nEC = 0;
  nIC = 0;
  nUC = 0;

  active.reset();
  inactive.reset();
}

/*
 *  s e t u p C o n s t r a i n t
 */
returnValue Constraints::setupConstraint(int _number, SubjectToStatus _status) {
  /* consistency check */
  if ((_number < 0) || (_number >= getNC()))
    return THROWERROR(RET_INDEX_OUT_OF_BOUNDS);

  /* Add constraint index to respective index list. */
  switch (_status) {
    case ST_INACTIVE:
      if (addIndex(getInactive(), _number, _status) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_SETUP_CONSTRAINT_FAILED);
      break;

    case ST_LOWER:
      if (addIndex(getActive(), _number, _status) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_SETUP_CONSTRAINT_FAILED);
      break;

    case ST_UPPER:
      if (addIndex(getActive(), _number, _status) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_SETUP_CONSTRAINT_FAILED);
      break;

    default:
      return THROWERROR(RET_INVALID_ARGUMENTS);
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  s e t u p A l l I n a c t i v e
 */
returnValue Constraints::setupAllInactive() {
  int i;

  /* 1) Place unbounded constraints at the beginning of the index list of
   * inactive constraints. */
  for (i = 0; i < nC; ++i) {
    if (getType(i) == ST_UNBOUNDED) {
      if (setupConstraint(i, ST_INACTIVE) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_SETUP_CONSTRAINT_FAILED);
    }
  }

  /* 2) Add remaining (i.e. "real" inequality) constraints to the index list of
   * inactive constraints. */
  for (i = 0; i < nC; ++i) {
    if (getType(i) == ST_BOUNDED) {
      if (setupConstraint(i, ST_INACTIVE) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_SETUP_CONSTRAINT_FAILED);
    }
  }

  /* 3) Place implicit equality constraints at the end of the index list of
   * inactive constraints. */
  for (i = 0; i < nC; ++i) {
    if (getType(i) == ST_EQUALITY) {
      if (setupConstraint(i, ST_INACTIVE) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_SETUP_CONSTRAINT_FAILED);
    }
  }

  /* 4) Moreover, add all constraints of unknown type. */
  for (i = 0; i < nC; ++i) {
    if (getType(i) == ST_UNKNOWN) {
      if (setupConstraint(i, ST_INACTIVE) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_SETUP_CONSTRAINT_FAILED);
    }
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  m o v e A c t i v e T o I n a c t i v e
 */
returnValue Constraints::moveActiveToInactive(int _number) {
  /* consistency check */
  if ((_number < 0) || (_number >= getNC()))
    return THROWERROR(RET_INDEX_OUT_OF_BOUNDS);

  /* Move index from indexlist of active constraints to that of inactive ones.
   */
  if (removeIndex(getActive(), _number) != SUCCESSFUL_RETURN)
    return THROWERROR(RET_MOVING_BOUND_FAILED);

  if (addIndex(getInactive(), _number, ST_INACTIVE) != SUCCESSFUL_RETURN)
    return THROWERROR(RET_MOVING_BOUND_FAILED);

  return SUCCESSFUL_RETURN;
}

/*
 *  m o v e I n a c t i v e T o A c t i v e
 */
returnValue Constraints::moveInactiveToActive(int _number,
                                              SubjectToStatus _status) {
  /* consistency check */
  if ((_number < 0) || (_number >= getNC()))
    return THROWERROR(RET_INDEX_OUT_OF_BOUNDS);

  /* Move index from indexlist of inactive constraints to that of active ones.
   */
  if (removeIndex(getInactive(), _number) != SUCCESSFUL_RETURN)
    return THROWERROR(RET_MOVING_BOUND_FAILED);

  if (addIndex(getActive(), _number, _status) != SUCCESSFUL_RETURN)
    return THROWERROR(RET_MOVING_BOUND_FAILED);

  return SUCCESSFUL_RETURN;
}

}  // namespace qpoases_embedded

/*
 *  end of file
 */
