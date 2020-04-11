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

#include <qpoases_embedded/QProblemB.hpp>

#include <qpoases_embedded/MessageHandling.hpp>

/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

namespace qpoases_embedded {

QProblemB::QProblemB(size_t _nV) : QProblemB(_nV, 0) {}

QProblemB::QProblemB(size_t _nV, size_t _nC)
    : nV(_nV),
      nC(_nC),
      H(nV * nV, 0),
      g(nV, 0),
      lb(nV, -INFTY),
      ub(nV, INFTY),
      bounds(nV, INDEXLISTFACTOR * nV),
      auxiliaryBounds(nV, INDEXLISTFACTOR * nV),
      R(nV * nV, 0),
      x(nV, 0),
      y(nV + nC, 0),
      g_user(nV, 0),
      lb_user(nV, 0),
      ub_user(nV, 0),
      delta_g(nV, 0),
      delta_lb(nV, 0),
      delta_ub(nV, 0),
      delta_xFR(nV, 0),
      delta_xFX(nV, 0),
      delta_yFX(nV, 0),
      rhs(nV, 0),
      r(nV, 0),
      delta_xFRz_TMP(nV, 0),
      delta_xFRz_RHS(nV, 0),
      HMX_delta_xFX(nV, 0) {
  for (int el = 0; el < nV; ++el) {
    H[el * nV + el] = 1.0;
  }

  reset();
}

void QProblemB::reset() {
  /* 1) Reset bounds. */
  bounds.reset();
  auxiliaryBounds.reset();

  /* 2) Reset Cholesky decomposition. */
  std::fill(R.begin(), R.end(), 0);

  /* 3) Reset steplength and status flags. */
  tau = 0.0;

  hessianType = HST_POSDEF_NULLSPACE; /* Hessian is assumed to be positive
                                         definite by default */
  infeasible = false;
  unbounded = false;

  status = QPS_NOTINITIALISED;
}

/*
 *  i n i t
 */
returnValue QProblemB::init(int& nWSR, const QProblemBCallback& callback) {
  reset();
  return solveInitialQP(nWSR, callback);
}

/*
 *  h o t s t a r t
 */
returnValue QProblemB::hotstart(int& nWSR, const QProblemBCallback& callback) {
  int l;

  /* consistency check */
  if ((getStatus() == QPS_NOTINITIALISED) ||
      (getStatus() == QPS_PREPARINGAUXILIARYQP) ||
      (getStatus() == QPS_PERFORMINGHOMOTOPY)) {
    return THROWERROR(RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED);
  }

  /* I) PREPARATIONS */
  /* 1) Reset status flags. */
  infeasible = false;
  unbounded = false;

  /* 2) Allocate delta vectors of gradient and bounds. */
  returnValue returnvalue;
  bool Delta_bB_isZero;

  int BC_idx;
  SubjectToStatus BC_status;

  /* II) MAIN HOMOTOPY LOOP */
  for (l = 0; l < nWSR; ++l) {
    status = QPS_PERFORMINGHOMOTOPY;

    THROWINFOMSG(RET_ITERATION_STARTED, "%d ...", l);

    /* 1) Setup index arrays. */
    const auto& FR_idx = bounds.getFree()->getNumberArray();
    const auto& FX_idx = bounds.getFixed()->getNumberArray();

    /* 2) Initialize shift direction of the gradient and the bounds. */
    returnvalue = hotstart_determineDataShift(
        FX_idx.data(), g_user.data(), lb_user.data(), ub_user.data(),
        delta_g.data(), delta_lb.data(), delta_ub.data(), Delta_bB_isZero);
    if (returnvalue != SUCCESSFUL_RETURN) {
      nWSR = l;
      THROWERROR(RET_SHIFT_DETERMINATION_FAILED);
      return returnvalue;
    }

    /* 3) Determination of step direction of X and Y. */
    returnvalue = hotstart_determineStepDirection(
        FR_idx.data(), FX_idx.data(), delta_g.data(), delta_lb.data(),
        delta_ub.data(), Delta_bB_isZero, delta_xFX.data(), delta_xFR.data(),
        delta_yFX.data());
    if (returnvalue != SUCCESSFUL_RETURN) {
      nWSR = l;
      THROWERROR(RET_STEPDIRECTION_DETERMINATION_FAILED);
      return returnvalue;
    }

    /* 4) Determination of step length TAU. */
    returnvalue = hotstart_determineStepLength(
        FR_idx.data(), FX_idx.data(), delta_lb.data(), delta_ub.data(),
        delta_xFR.data(), delta_yFX.data(), BC_idx, BC_status);
    if (returnvalue != SUCCESSFUL_RETURN) {
      nWSR = l;
      THROWERROR(RET_STEPLENGTH_DETERMINATION_FAILED);
      return returnvalue;
    }

    /* 5) Realization of the homotopy step. */
    returnvalue = hotstart_performStep(
        FR_idx.data(), FX_idx.data(), delta_g.data(), delta_lb.data(),
        delta_ub.data(), delta_xFX.data(), delta_xFR.data(), delta_yFX.data(),
        BC_idx, BC_status);

    if (returnvalue != SUCCESSFUL_RETURN) {
      nWSR = l;

      /* optimal solution found? */
      if (returnvalue == RET_OPTIMAL_SOLUTION_FOUND) {
        status = QPS_SOLVED;

        THROWINFO(RET_OPTIMAL_SOLUTION_FOUND);

        if (callback) {
          callback(l, tau, getNFX(), BC_idx, BC_status);
        }
        return checkKKTconditions();
      } else {
        /* checks for infeasibility... */
        if (infeasible == true) {
          status = QPS_HOMOTOPYQPSOLVED;
          return THROWERROR(RET_HOTSTART_STOPPED_INFEASIBILITY);
        }

        /* ...unboundedness... */
        if (unbounded ==
            true) /* not necessary since objective function convex! */
          return THROWERROR(RET_HOTSTART_STOPPED_UNBOUNDEDNESS);

        /* ... and throw unspecific error otherwise */
        THROWERROR(RET_HOMOTOPY_STEP_FAILED);
        return returnvalue;
      }
    }

    /* 6) Output information of successful QP iteration. */
    status = QPS_HOMOTOPYQPSOLVED;

    if (callback && !callback(l, tau, getNFX(), BC_idx, BC_status)) {
      THROWINFO(RET_USER_ABORT_REQUESTED);
      return checkKKTconditions();
    }
  }

  THROWERRORMSG(RET_MAX_NWSR_REACHED, "nWSR = %d", nWSR);

  returnValue returnvalueKKTcheck = checkKKTconditions();
  if (returnvalueKKTcheck != SUCCESSFUL_RETURN) return returnvalueKKTcheck;
  return RET_MAX_NWSR_REACHED;
}

/*
 *  g e t N Z
 */
int QProblemB::getNZ() {
  /* if no constraints are present: nZ=nFR */
  return bounds.getFree()->getLength();
}

/*
 *  g e t O b j V a l
 */
real_t QProblemB::getObjVal() const {
  real_t objVal;

  /* calculated optimal objective function value
   * only if current QP has been solved */
  if ((getStatus() == QPS_AUXILIARYQPSOLVED) ||
      (getStatus() == QPS_HOMOTOPYQPSOLVED) || (getStatus() == QPS_SOLVED)) {
    objVal = getObjVal(x.data());
  } else {
    objVal = INFTY;
  }

  return objVal;
}

/*
 *  g e t O b j V a l
 */
real_t QProblemB::getObjVal(const real_t* const _x) const {
  int i, j;

  real_t obj_tmp = 0.0;

  for (i = 0; i < nV; ++i) {
    obj_tmp += _x[i] * g[i];

    for (j = 0; j < nV; ++j) obj_tmp += 0.5 * _x[i] * H[i * nV + j] * _x[j];
  }

  return obj_tmp;
}

/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *  c h e c k F o r I d e n t i t y H e s s i a n
 */
returnValue QProblemB::checkForIdentityHessian() {
  int i, j;

  /* nothing to do as status flag remains unaltered
   * if Hessian differs from identity matrix */
  if (hessianType == HST_IDENTITY) return SUCCESSFUL_RETURN;

  /* 1) If Hessian differs from identity matrix,
   *    return without changing the internal HessianType. */
  for (i = 0; i < nV; ++i)
    if (std::fabs(H[i * nV + i] - 1.0) > EPS) return SUCCESSFUL_RETURN;

  for (i = 0; i < nV; ++i) {
    for (j = 0; j < i; ++j)
      if ((std::fabs(H[i * nV + j]) > EPS) || (std::fabs(H[j * nV + i]) > EPS))
        return SUCCESSFUL_RETURN;
  }

  /* 2) If this point is reached, Hessian equals the idetity matrix. */
  hessianType = HST_IDENTITY;

  return SUCCESSFUL_RETURN;
}

/*
 *  s e t u p S u b j e c t T o T y p e
 */
returnValue QProblemB::setupSubjectToType() {
  int i;

  /* 1) Check if lower bounds are present. */
  bounds.setNoLower(true);
  for (i = 0; i < nV; ++i)
    if (lb[i] > -INFTY) {
      bounds.setNoLower(false);
      break;
    }

  /* 2) Check if upper bounds are present. */
  bounds.setNoUpper(true);
  for (i = 0; i < nV; ++i)
    if (ub[i] < INFTY) {
      bounds.setNoUpper(false);
      break;
    }

  /* 3) Determine implicitly fixed and unbounded variables. */
  int nFV = 0;
  int nUV = 0;

  for (i = 0; i < nV; ++i)
    if ((lb[i] < -INFTY + BOUNDTOL) && (ub[i] > INFTY - BOUNDTOL)) {
      bounds.setType(i, ST_UNBOUNDED);
      ++nUV;
    } else {
      if (lb[i] > ub[i] - BOUNDTOL) {
        bounds.setType(i, ST_EQUALITY);
        ++nFV;
      } else {
        bounds.setType(i, ST_BOUNDED);
      }
    }

  /* 4) Set dimensions of bounds structure. */
  bounds.setNFV(nFV);
  bounds.setNUV(nUV);
  bounds.setNBV(nV - nFV - nUV);

  return SUCCESSFUL_RETURN;
}

/*
 *  c h o l e s k y D e c o m p o s i t i o n
 */
returnValue QProblemB::setupCholeskyDecomposition() {
  int i, j, k, ii, jj;
  int nFR = getNFR();

  /* 1) Initialises R with all zeros. */
  for (i = 0; i < nV; ++i)
    for (j = 0; j < nV; ++j) R[i * nV + j] = 0.0;

  /* 2) Calculate Cholesky decomposition of H (projected to free variables). */
  if (hessianType == HST_IDENTITY) {
    /* if Hessian is identity, so is its Cholesky factor. */
    for (i = 0; i < nFR; ++i) R[i * nV + i] = 1.0;
  } else {
    if (nFR > 0) {
      const auto& FR_idx = bounds.getFree()->getNumberArray();

      /* R'*R = H */
      real_t sum;
      real_t inv;

      for (i = 0; i < nFR; ++i) {
        /* j == i */
        ii = FR_idx[i];
        sum = H[ii * nV + ii];

        for (k = (i - 1); k >= 0; --k) sum -= R[k * nV + i] * R[k * nV + i];

        if (sum > 0.0) {
          R[i * nV + i] = sqrt(sum);
          inv = 1.0 / R[i * nV + i];
        } else {
          hessianType = HST_SEMIDEF;
          return THROWERROR(RET_HESSIAN_NOT_SPD);
        }

        /* j > i */
        for (j = (i + 1); j < nFR; ++j) {
          jj = FR_idx[j];
          sum = H[jj * nV + ii];

          for (k = (i - 1); k >= 0; --k) sum -= R[k * nV + i] * R[k * nV + j];

          R[i * nV + j] = sum * inv;
        }
      }
    }
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  s o l v e I n i t i a l Q P
 */
returnValue QProblemB::solveInitialQP(int& nWSR,
                                      const QProblemBCallback& callback) {
  status = QPS_NOTINITIALISED;

  std::copy(g_user.begin(), g_user.end(), g.begin());
  std::copy(lb_user.begin(), lb_user.end(), lb.begin());
  std::copy(ub_user.begin(), ub_user.end(), ub.begin());

  /* I) ANALYSE QP DATA: */
  /* 1) Check if Hessian happens to be the identity matrix. */
  if (checkForIdentityHessian() != SUCCESSFUL_RETURN)
    return THROWERROR(RET_INIT_FAILED);

  /* 2) Setup type of bounds (i.e. unbounded, implicitly fixed etc.). */
  if (setupSubjectToType() != SUCCESSFUL_RETURN)
    return THROWERROR(RET_INIT_FAILED);

  status = QPS_PREPARINGAUXILIARYQP;

  /* II) SETUP AUXILIARY QP WITH GIVEN OPTIMAL SOLUTION: */
  /* 1) Setup bounds data structure. */
  if (bounds.setupAllFree() != SUCCESSFUL_RETURN)
    return THROWERROR(RET_INIT_FAILED);

  /* 2) Setup optimal primal/dual solution for auxiliary QP. */
  if (setupAuxiliaryQPsolution() != SUCCESSFUL_RETURN)
    return THROWERROR(RET_INIT_FAILED);

  /* 3) Obtain linear independent working set for auxiliary QP. */

  if (obtainAuxiliaryWorkingSet() != SUCCESSFUL_RETURN)
    return THROWERROR(RET_INIT_FAILED);

  /* 4) Setup working set of auxiliary QP and setup cholesky decomposition. */
  if (setupAuxiliaryWorkingSet(true) != SUCCESSFUL_RETURN)
    return THROWERROR(RET_INIT_FAILED);

  /* At the moment we can only provide a Cholesky of the Hessian if
   * the solver is cold-started. */
  if (setupCholeskyDecomposition() != SUCCESSFUL_RETURN)
    return THROWERROR(RET_INIT_FAILED_CHOLESKY);

  /* 5) Setup QP data of an auxiliary QP having an optimal solution
   * as specified by the user (or xOpt = yOpt = 0, by default). */
  if (setupAuxiliaryQPgradient() != SUCCESSFUL_RETURN)
    return THROWERROR(RET_INIT_FAILED);

  if (setupAuxiliaryQPbounds(true) != SUCCESSFUL_RETURN)
    return THROWERROR(RET_INIT_FAILED);

  status = QPS_AUXILIARYQPSOLVED;

  /* III) SOLVE ACTUAL INITIAL QP: */
  /* Use hotstart method to find the solution of the original initial QP,... */
  returnValue returnvalue = hotstart(nWSR, callback);

  /* ... check for infeasibility and unboundedness... */
  if (isInfeasible() == true) return THROWERROR(RET_INIT_FAILED_INFEASIBILITY);

  if (isUnbounded() == true) return THROWERROR(RET_INIT_FAILED_UNBOUNDEDNESS);

  /* ... and internal errors. */
  if ((returnvalue != SUCCESSFUL_RETURN) &&
      (returnvalue != RET_MAX_NWSR_REACHED) &&
      (returnvalue != RET_INACCURATE_SOLUTION) &&
      (returnvalue != RET_NO_SOLUTION))
    return THROWERROR(RET_INIT_FAILED_HOTSTART);

  return returnvalue;
}

/*
 *  o b t a i n A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblemB::obtainAuxiliaryWorkingSet() {
  /* Obtain initial working set in accordance to sign of dual solution
   * vector. */
  for (int i = 0; i < nV; ++i) {
    if (y[i] > ZERO) {
      if (auxiliaryBounds.setupBound(i, ST_LOWER) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_OBTAINING_WORKINGSET_FAILED);
      continue;
    }

    if (y[i] < -ZERO) {
      if (auxiliaryBounds.setupBound(i, ST_UPPER) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_OBTAINING_WORKINGSET_FAILED);
      continue;
    }

    /* Moreover, add all implictly fixed variables if specified. */
    if (bounds.getType(i) == ST_EQUALITY) {
      if (auxiliaryBounds.setupBound(i, ST_LOWER) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_OBTAINING_WORKINGSET_FAILED);
    } else {
      if (auxiliaryBounds.setupBound(i, ST_INACTIVE) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_OBTAINING_WORKINGSET_FAILED);
    }
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  s e t u p A u x i l i a r y W o r k i n g S e t
 */
returnValue QProblemB::setupAuxiliaryWorkingSet(bool setupAfresh) {
  int i;

  /* consistency checks */
  for (i = 0; i < nV; ++i)
    if ((bounds.getStatus(i) == ST_UNDEFINED) ||
        (auxiliaryBounds.getStatus(i) == ST_UNDEFINED))
      return THROWERROR(RET_UNKNOWN_BUG);

  /* I) SETUP CHOLESKY FLAG:
   *    Cholesky decomposition shall only be updated if working set
   *    shall be updated (i.e. NOT setup afresh!) */
  bool updateCholesky;
  if (setupAfresh == true)
    updateCholesky = false;
  else
    updateCholesky = true;

  /* II) REMOVE FORMERLY ACTIVE BOUNDS (IF NECESSARY): */
  if (setupAfresh == false) {
    /* Remove all active bounds that shall be inactive AND
     *  all active bounds that are active at the wrong bound. */
    for (i = 0; i < nV; ++i) {
      if ((bounds.getStatus(i) == ST_LOWER) &&
          (auxiliaryBounds.getStatus(i) != ST_LOWER))
        if (removeBound(i, updateCholesky) != SUCCESSFUL_RETURN)
          return THROWERROR(RET_SETUP_WORKINGSET_FAILED);

      if ((bounds.getStatus(i) == ST_UPPER) &&
          (auxiliaryBounds.getStatus(i) != ST_UPPER))
        if (removeBound(i, updateCholesky) != SUCCESSFUL_RETURN)
          return THROWERROR(RET_SETUP_WORKINGSET_FAILED);
    }
  }

  /* III) ADD NEWLY ACTIVE BOUNDS: */
  /*      Add all inactive bounds that shall be active AND
   *      all formerly active bounds that have been active at the wrong bound.
   */
  for (i = 0; i < nV; ++i) {
    if ((bounds.getStatus(i) == ST_INACTIVE) &&
        (auxiliaryBounds.getStatus(i) != ST_INACTIVE)) {
      if (addBound(i, auxiliaryBounds.getStatus(i), updateCholesky) !=
          SUCCESSFUL_RETURN)
        return THROWERROR(RET_SETUP_WORKINGSET_FAILED);
    }
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  s e t u p A u x i l i a r y Q P s o l u t i o n
 */
returnValue QProblemB::setupAuxiliaryQPsolution() {
  std::fill(x.begin(), x.end(), 0);
  return SUCCESSFUL_RETURN;
}

/*
 *  s e t u p A u x i l i a r y Q P g r a d i e n t
 */
returnValue QProblemB::setupAuxiliaryQPgradient() {
  int i, j;

  /* Setup gradient vector: g = -H*x + y'*Id. */
  for (i = 0; i < nV; ++i) {
    /* y'*Id */
    g[i] = y[i];

    /* -H*x */
    for (j = 0; j < nV; ++j) g[i] -= H[i * nV + j] * x[j];
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  s e t u p A u x i l i a r y Q P b o u n d s
 */
returnValue QProblemB::setupAuxiliaryQPbounds(bool useRelaxation) {
  int i;

  /* Setup bound vectors. */
  for (i = 0; i < nV; ++i) {
    switch (bounds.getStatus(i)) {
      case ST_INACTIVE:
        if (useRelaxation == true) {
          if (bounds.getType(i) == ST_EQUALITY) {
            lb[i] = x[i];
            ub[i] = x[i];
          } else {
            lb[i] = x[i] - BOUNDRELAXATION;
            ub[i] = x[i] + BOUNDRELAXATION;
          }
        }
        break;

      case ST_LOWER:
        lb[i] = x[i];
        if (bounds.getType(i) == ST_EQUALITY) {
          ub[i] = x[i];
        } else {
          if (useRelaxation == true) ub[i] = x[i] + BOUNDRELAXATION;
        }
        break;

      case ST_UPPER:
        ub[i] = x[i];
        if (bounds.getType(i) == ST_EQUALITY) {
          lb[i] = x[i];
        } else {
          if (useRelaxation == true) lb[i] = x[i] - BOUNDRELAXATION;
        }
        break;

      default:
        return THROWERROR(RET_UNKNOWN_BUG);
    }
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  a d d B o u n d
 */
returnValue QProblemB::addBound(int number, SubjectToStatus B_status,
                                bool updateCholesky) {
  int i, j;
  int nFR = getNFR();

  /* consistency check */
  if ((getStatus() == QPS_NOTINITIALISED) ||
      (getStatus() == QPS_AUXILIARYQPSOLVED) ||
      (getStatus() == QPS_HOMOTOPYQPSOLVED) || (getStatus() == QPS_SOLVED)) {
    return THROWERROR(RET_UNKNOWN_BUG);
  }

  /* Perform cholesky updates only if QProblemB has been initialised! */
  if ((getStatus() == QPS_PREPARINGAUXILIARYQP) || (updateCholesky == false)) {
    /* UPDATE INDICES */
    if (bounds.moveFreeToFixed(number, B_status) != SUCCESSFUL_RETURN)
      return THROWERROR(RET_ADDBOUND_FAILED);

    return SUCCESSFUL_RETURN;
  }

  /* I) PERFORM CHOLESKY UPDATE: */
  /* 1) Index of variable to be added within the list of free variables. */
  int number_idx = bounds.getFree()->getIndex(number);

  real_t c, s;

  /* 2) Use row-wise Givens rotations to restore upper triangular form of R. */
  for (i = number_idx + 1; i < nFR; ++i) {
    computeGivens(R[(i - 1) * nV + i], R[i * nV + i], R[(i - 1) * nV + i],
                  R[i * nV + i], c, s);

    for (j = (1 + i); j < nFR; ++j) /* last column of R is thrown away */
      applyGivens(c, s, R[(i - 1) * nV + j], R[i * nV + j], R[(i - 1) * nV + j],
                  R[i * nV + j]);
  }

  /* 3) Delete <number_idx>th column and ... */
  for (i = 0; i < nFR - 1; ++i)
    for (j = number_idx + 1; j < nFR; ++j) R[i * nV + j - 1] = R[i * nV + j];
  /* ... last column of R. */
  for (i = 0; i < nFR; ++i) R[i * nV + nFR - 1] = 0.0;

  /* II) UPDATE INDICES */
  if (bounds.moveFreeToFixed(number, B_status) != SUCCESSFUL_RETURN)
    return THROWERROR(RET_ADDBOUND_FAILED);

  return SUCCESSFUL_RETURN;
}

returnValue QProblemB::removeBound(int number, bool updateCholesky) {
  int i, ii;
  int nFR = getNFR();

  /* consistency check */
  if ((getStatus() == QPS_NOTINITIALISED) ||
      (getStatus() == QPS_AUXILIARYQPSOLVED) ||
      (getStatus() == QPS_HOMOTOPYQPSOLVED) || (getStatus() == QPS_SOLVED)) {
    return THROWERROR(RET_UNKNOWN_BUG);
  }

  /* I) UPDATE INDICES */
  if (bounds.moveFixedToFree(number) != SUCCESSFUL_RETURN)
    return THROWERROR(RET_REMOVEBOUND_FAILED);

  /* Perform cholesky updates only if QProblemB has been initialised! */
  if ((getStatus() == QPS_PREPARINGAUXILIARYQP) || (updateCholesky == false))
    return SUCCESSFUL_RETURN;

  /* II) PERFORM CHOLESKY UPDATE */
  const auto& FR_idx = bounds.getFree()->getNumberArray();

  /* 1) Calculate new column of cholesky decomposition. */
  real_t r0 = H[number * nV + number];

  for (i = 0; i < nFR; ++i) {
    ii = FR_idx[i];
    rhs[i] = H[number * nV + ii];
  }

  if (backsolveR(rhs.data(), true, true, r.data()) != SUCCESSFUL_RETURN)
    return THROWERROR(RET_REMOVEBOUND_FAILED);

  for (i = 0; i < nFR; ++i) r0 -= r[i] * r[i];

  /* 2) Store new column into R. */
  for (i = 0; i < nFR; ++i) R[i * nV + nFR] = r[i];

  if (r0 > 0.0)
    R[nFR * nV + nFR] = sqrt(r0);
  else {
    hessianType = HST_SEMIDEF;
    return THROWERROR(RET_HESSIAN_NOT_SPD);
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  b a c k s o l v e R  (CODE DUPLICATED IN QProblem CLASS!!!)
 */
returnValue QProblemB::backsolveR(const real_t* const b, bool transposed,
                                  real_t* const a) {
  /* Call standard backsolve procedure (i.e. removingBound == false). */
  return backsolveR(b, transposed, false, a);
}

/*
 *  b a c k s o l v e R  (CODE DUPLICATED IN QProblem CLASS!!!)
 */
returnValue QProblemB::backsolveR(const real_t* const b, bool transposed,
                                  bool removingBound, real_t* const a) {
  int i, j;
  int nR = getNZ();

  real_t sum;

  /* if backsolve is called while removing a bound, reduce nZ by one. */
  if (removingBound == true) --nR;

  /* nothing to do */
  if (nR <= 0) return SUCCESSFUL_RETURN;

  /* Solve Ra = b, where R might be transposed. */
  if (transposed == false) {
    /* solve Ra = b */
    for (i = (nR - 1); i >= 0; --i) {
      sum = b[i];
      for (j = (i + 1); j < nR; ++j) sum -= R[i * nV + j] * a[j];

      if (std::fabs(R[i * nV + i]) > ZERO)
        a[i] = sum / R[i * nV + i];
      else
        return THROWERROR(RET_DIV_BY_ZERO);
    }
  } else {
    /* solve R^T*a = b */
    for (i = 0; i < nR; ++i) {
      sum = b[i];

      for (j = 0; j < i; ++j) sum -= R[j * nV + i] * a[j];

      if (std::fabs(R[i * nV + i]) > ZERO)
        a[i] = sum / R[i * nV + i];
      else
        return THROWERROR(RET_DIV_BY_ZERO);
    }
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  h o t s t a r t _ d e t e r m i n e D a t a S h i f t
 */
returnValue QProblemB::hotstart_determineDataShift(
    const int* const FX_idx, const real_t* const g_new,
    const real_t* const lb_new, const real_t* const ub_new,
    real_t* const delta_g, real_t* const delta_lb, real_t* const delta_ub,
    bool& Delta_bB_isZero) {
  int i, ii;
  int nFX = getNFX();

  /* 1) Calculate shift directions. */
  for (i = 0; i < nV; ++i) delta_g[i] = g_new[i] - g[i];

  if (lb_new != 0) {
    for (i = 0; i < nV; ++i) delta_lb[i] = lb_new[i] - lb[i];
  } else {
    /* if no lower bounds exist, assume the new lower bounds to be -infinity */
    for (i = 0; i < nV; ++i) delta_lb[i] = -INFTY - lb[i];
  }

  if (ub_new != 0) {
    for (i = 0; i < nV; ++i) delta_ub[i] = ub_new[i] - ub[i];
  } else {
    /* if no upper bounds exist, assume the new upper bounds to be infinity */
    for (i = 0; i < nV; ++i) delta_ub[i] = INFTY - ub[i];
  }

  /* 2) Determine if active bounds are to be shifted. */
  Delta_bB_isZero = true;

  for (i = 0; i < nFX; ++i) {
    ii = FX_idx[i];

    if ((std::fabs(delta_lb[ii]) > EPS) || (std::fabs(delta_ub[ii]) > EPS)) {
      Delta_bB_isZero = false;
      break;
    }
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  a r e B o u n d s C o n s i s t e n t
 */
bool QProblemB::areBoundsConsistent(const real_t* const delta_lb,
                                    const real_t* const delta_ub) const {
  int i;

  /* Check if delta_lb[i] is greater than delta_ub[i]
   * for a component i whose bounds are already (numerically) equal. */
  for (i = 0; i < nV; ++i)
    if ((lb[i] > ub[i] - BOUNDTOL) && (delta_lb[i] > delta_ub[i] + EPS))
      return false;

  return true;
}

/*****************************************************************************
 *  P R I V A T E                                                            *
 *****************************************************************************/

/*
 *  h o t s t a r t _ d e t e r m i n e S t e p D i r e c t i o n
 */
returnValue QProblemB::hotstart_determineStepDirection(
    const int* const FR_idx, const int* const FX_idx,
    const real_t* const delta_g, const real_t* const delta_lb,
    const real_t* const delta_ub, bool Delta_bB_isZero, real_t* const delta_xFX,
    real_t* const delta_xFR, real_t* const delta_yFX) {
  int i, j, ii, jj;
  int nFR = getNFR();
  int nFX = getNFX();

  /* initialise auxiliary vectors */
  for (i = 0; i < nFR; ++i) HMX_delta_xFX[i] = 0.0;

  /* I) DETERMINE delta_xFX */
  if (nFX > 0) {
    for (i = 0; i < nFX; ++i) {
      ii = FX_idx[i];

      if (bounds.getStatus(ii) == ST_LOWER)
        delta_xFX[i] = delta_lb[ii];
      else
        delta_xFX[i] = delta_ub[ii];
    }
  }

  /* II) DETERMINE delta_xFR */
  if (nFR > 0) {
    /* Determine delta_xFRz. */
    if (Delta_bB_isZero == false) {
      for (i = 0; i < nFR; ++i) {
        ii = FR_idx[i];
        for (j = 0; j < nFX; ++j) {
          jj = FX_idx[j];
          HMX_delta_xFX[i] += H[ii * nV + jj] * delta_xFX[j];
        }
      }
    }

    if (Delta_bB_isZero == true) {
      for (j = 0; j < nFR; ++j) {
        jj = FR_idx[j];
        delta_xFRz_RHS[j] = delta_g[jj];
      }
    } else {
      for (j = 0; j < nFR; ++j) {
        jj = FR_idx[j];
        delta_xFRz_RHS[j] = delta_g[jj] + HMX_delta_xFX[j]; /* *ZFR */
      }
    }

    for (i = 0; i < nFR; ++i) delta_xFR[i] = -delta_xFRz_RHS[i];

    if (backsolveR(delta_xFR, true, delta_xFRz_TMP.data()) != SUCCESSFUL_RETURN)
      return THROWERROR(RET_STEPDIRECTION_FAILED_CHOLESKY);

    if (backsolveR(delta_xFRz_TMP.data(), false, delta_xFR) !=
        SUCCESSFUL_RETURN)
      return THROWERROR(RET_STEPDIRECTION_FAILED_CHOLESKY);
  }

  /* III) DETERMINE delta_yFX */
  if (nFX > 0) {
    for (i = 0; i < nFX; ++i) {
      ii = FX_idx[i];

      delta_yFX[i] = 0.0;
      for (j = 0; j < nFR; ++j) {
        jj = FR_idx[j];
        delta_yFX[i] += H[ii * nV + jj] * delta_xFR[j];
      }

      if (Delta_bB_isZero == false) {
        for (j = 0; j < nFX; ++j) {
          jj = FX_idx[j];
          delta_yFX[i] += H[ii * nV + jj] * delta_xFX[j];
        }
      }

      delta_yFX[i] += delta_g[ii];
    }
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  h o t s t a r t _ d e t e r m i n e S t e p L e n g t h
 */
returnValue QProblemB::hotstart_determineStepLength(
    const int* const FR_idx, const int* const FX_idx,
    const real_t* const delta_lb, const real_t* const delta_ub,
    const real_t* const delta_xFR, const real_t* const delta_yFX, int& BC_idx,
    SubjectToStatus& BC_status) {
  int i, ii;
  int nFR = getNFR();
  int nFX = getNFX();

  real_t tau_tmp;
  real_t tau_new = 1.0;

  BC_idx = 0;
  BC_status = ST_UNDEFINED;

  /* I) DETERMINE MAXIMUM DUAL STEPLENGTH, i.e. ensure that
   *    active dual bounds remain valid (ignoring implicitly fixed variables):
   */
  for (i = 0; i < nFX; ++i) {
    ii = FX_idx[i];

    if (bounds.getType(ii) != ST_EQUALITY) {
      if (bounds.getStatus(ii) == ST_LOWER) {
        /* 1) Active lower bounds. */
        if ((delta_yFX[i] < -ZERO) && (y[ii] >= 0.0)) {
          tau_tmp = y[ii] / (-delta_yFX[i]);
          if (tau_tmp < tau_new) {
            if (tau_tmp >= 0.0) {
              tau_new = tau_tmp;
              BC_idx = ii;
              BC_status = ST_INACTIVE;
            }
          }
        }
      } else {
        /* 2) Active upper bounds. */
        if ((delta_yFX[i] > ZERO) && (y[ii] <= 0.0)) {
          tau_tmp = y[ii] / (-delta_yFX[i]);
          if (tau_tmp < tau_new) {
            if (tau_tmp >= 0.0) {
              tau_new = tau_tmp;
              BC_idx = ii;
              BC_status = ST_INACTIVE;
            }
          }
        }
      }
    }
  }

  /* II) DETERMINE MAXIMUM PRIMAL STEPLENGTH, i.e. ensure that
   *     inactive bounds remain valid (ignoring unbounded variables). */
  /* 1) Inactive lower bounds. */
  if (bounds.isNoLower() == false) {
    for (i = 0; i < nFR; ++i) {
      ii = FR_idx[i];

      if (bounds.getType(ii) != ST_UNBOUNDED) {
        if (delta_lb[ii] > delta_xFR[i]) {
          if (x[ii] > lb[ii])
            tau_tmp = (x[ii] - lb[ii]) / (delta_lb[ii] - delta_xFR[i]);
          else
            tau_tmp = 0.0;

          if (tau_tmp < tau_new) {
            if (tau_tmp >= 0.0) {
              tau_new = tau_tmp;
              BC_idx = ii;
              BC_status = ST_LOWER;
            }
          }
        }
      }
    }
  }

  /* 2) Inactive upper bounds. */
  if (bounds.isNoUpper() == false) {
    for (i = 0; i < nFR; ++i) {
      ii = FR_idx[i];

      if (bounds.getType(ii) != ST_UNBOUNDED) {
        if (delta_ub[ii] < delta_xFR[i]) {
          if (x[ii] < ub[ii])
            tau_tmp = (x[ii] - ub[ii]) / (delta_ub[ii] - delta_xFR[i]);
          else
            tau_tmp = 0.0;

          if (tau_tmp < tau_new) {
            if (tau_tmp >= 0.0) {
              tau_new = tau_tmp;
              BC_idx = ii;
              BC_status = ST_UPPER;
            }
          }
        }
      }
    }
  }

  /* III) SET MAXIMUM HOMOTOPY STEPLENGTH */
  tau = tau_new;

  if (BC_status == ST_UNDEFINED) {
    THROWINFOMSG(RET_STEPSIZE_NONPOSITIVE, "Stepsize is %.6e!", tau);
  } else {
    THROWINFOMSG(RET_STEPSIZE_NONPOSITIVE,
                 "Stepsize is %.6e! (BC_idx = %d, BC_status = %d)", tau, BC_idx,
                 BC_status);
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  h o t s t a r t _ p e r f o r m S t e p
 */
returnValue QProblemB::hotstart_performStep(
    const int* const FR_idx, const int* const FX_idx,
    const real_t* const delta_g, const real_t* const delta_lb,
    const real_t* const delta_ub, const real_t* const delta_xFX,
    const real_t* const delta_xFR, const real_t* const delta_yFX, int BC_idx,
    SubjectToStatus BC_status) {
  int i, ii;
  int nFR = getNFR();
  int nFX = getNFX();

  /* I) CHECK BOUNDS' CONSISTENCY */
  if (areBoundsConsistent(delta_lb, delta_ub) == false) {
    infeasible = true;
    tau = 0.0;

    return THROWERROR(RET_QP_INFEASIBLE);
  }

  /* II) GO TO ACTIVE SET CHANGE */
  if (tau > ZERO) {
    /* 1) Perform step in primal und dual space. */
    for (i = 0; i < nFR; ++i) {
      ii = FR_idx[i];
      x[ii] += tau * delta_xFR[i];
    }

    for (i = 0; i < nFX; ++i) {
      ii = FX_idx[i];
      x[ii] += tau * delta_xFX[i];
      y[ii] += tau * delta_yFX[i];
    }

    /* 2) Shift QP data. */
    for (i = 0; i < nV; ++i) {
      g[i] += tau * delta_g[i];
      lb[i] += tau * delta_lb[i];
      ub[i] += tau * delta_ub[i];
    }
  } else {
    /* print a stepsize error if stepsize is zero */
    THROWERRORMSG(RET_STEPSIZE, "Stepsize is %.6e", tau);
  }

  /* III) UPDATE ACTIVE SET */
  switch (BC_status) {
    /* Optimal solution found as no working set change detected. */
    case ST_UNDEFINED:
      return RET_OPTIMAL_SOLUTION_FOUND;

    /* Remove one variable from active set. */
    case ST_INACTIVE:
      THROWINFOMSG(RET_REMOVE_FROM_ACTIVESET, "bound no. %d.", BC_idx);

      if (removeBound(BC_idx, true) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_REMOVE_FROM_ACTIVESET_FAILED);

      y[BC_idx] = 0.0;
      break;

    /* Add one variable to active set. */
    default:
      if (BC_status == ST_LOWER) {
        THROWINFOMSG(RET_ADD_TO_ACTIVESET, "lower bound no. %d.", BC_idx);
      } else {
        THROWINFOMSG(RET_ADD_TO_ACTIVESET, "upper bound no. %d.", BC_idx);
      }

      if (addBound(BC_idx, BC_status, true) != SUCCESSFUL_RETURN)
        return THROWERROR(RET_ADD_TO_ACTIVESET_FAILED);
      break;
  }

  return SUCCESSFUL_RETURN;
}

/*
 *  c h e c k K K T c o n d i t i o n s
 */
returnValue QProblemB::checkKKTconditions() {
  int i, j;

  real_t tmp;
  real_t maxKKTviolation = 0.0;

  /* 1) Check for Hx + g - y*A' = 0  (here: A = Id). */
  for (i = 0; i < nV; ++i) {
    tmp = g[i];

    for (j = 0; j < nV; ++j) tmp += H[i * nV + j] * x[j];

    tmp -= y[i];

    if (std::fabs(tmp) > maxKKTviolation) maxKKTviolation = std::fabs(tmp);
  }

  /* 2) Check for lb <= x <= ub. */
  for (i = 0; i < nV; ++i) {
    if (lb[i] - x[i] > maxKKTviolation) maxKKTviolation = lb[i] - x[i];

    if (x[i] - ub[i] > maxKKTviolation) maxKKTviolation = x[i] - ub[i];
  }

  /* 3) Check for correct sign of y and for complementary slackness. */
  for (i = 0; i < nV; ++i) {
    switch (bounds.getStatus(i)) {
      case ST_LOWER:
        if (-y[i] > maxKKTviolation) maxKKTviolation = -y[i];
        if (std::fabs((x[i] - lb[i]) * y[i]) > maxKKTviolation)
          maxKKTviolation = std::fabs((x[i] - lb[i]) * y[i]);
        break;

      case ST_UPPER:
        if (y[i] > maxKKTviolation) maxKKTviolation = y[i];
        if (std::fabs((ub[i] - x[i]) * y[i]) > maxKKTviolation)
          maxKKTviolation = std::fabs((ub[i] - x[i]) * y[i]);
        break;

      default: /* inactive */
        if (std::fabs(y[i]) > maxKKTviolation)
          maxKKTviolation = std::fabs(y[i]);
        break;
    }
  }

  if (maxKKTviolation > CRITICALACCURACY) return RET_NO_SOLUTION;

  if (maxKKTviolation > DESIREDACCURACY) return RET_INACCURATE_SOLUTION;

  return SUCCESSFUL_RETURN;
}

}  // namespace qpoases_embedded

/*
 *  end of file
 */
