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

#ifndef QPOASES_QPROBLEMB_HPP
#define QPOASES_QPROBLEMB_HPP

#include <cmath>
#include <functional>

#include <qpoases_embedded/Bounds.hpp>

namespace qpoases_embedded {

/**
 * Class for setting up and solving quadratic programs with (simple) bounds
 * only. The main feature is the possibily to use the newly developed online
 * active set strategy for parametric quadratic programming.
 */
class QProblemB {
  /*
   *  PUBLIC MEMBER FUNCTIONS
   */
 public:
  using QProblemBCallback = std::function<bool(
      int,            /* iteration: Number of current iteration. */
      real_t,         /* tau: last homotopy step length. */
      int,            /* nFX: The number of fixed variables. */
      int,            /* BC_idx: Index of blocking constraint. */
      SubjectToStatus /* BC_status: status of blocking bound/constraint. */
      )>;

  const int nV; /**< Number of variables. */
  const int nC; /**< Number of constraints. */

  QProblemB() = delete;
  /** Constructor which takes the QP dimension only. */
  explicit QProblemB(size_t _nV /**< Number of variables. */);

  /** Clears all data structures of QProblemB except for QP data. */
  void reset();

  /** Initialises a QProblemB with given QP data and solves it
   *  using an initial homotopy with empty working set (at most nWSR
   iterations). *  \return SUCCESSFUL_RETURN \n RET_INIT_FAILED \n
                          RET_INIT_FAILED_CHOLESKY \n
                          RET_INIT_FAILED_HOTSTART \n
                          RET_INIT_FAILED_INFEASIBILITY \n
                          RET_INIT_FAILED_UNBOUNDEDNESS \n
                          RET_MAX_NWSR_REACHED \n
                          RET_INVALID_ARGUMENTS \n
                          RET_INACCURATE_SOLUTION \n
                          RET_NO_SOLUTION */
  returnValue init(
      int& nWSR, /**< Input: Maximum number of working set recalculations when
                    using initial homotopy. \n Output: Number of performed
                    working set recalculations. */
      const QProblemBCallback& callback);

  /** Solves an initialised QProblemB using online active set strategy.
   *  \return SUCCESSFUL_RETURN \n
                          RET_MAX_NWSR_REACHED \n
                          RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED \n
                          RET_HOTSTART_FAILED \n
                          RET_SHIFT_DETERMINATION_FAILED \n
                          RET_STEPDIRECTION_DETERMINATION_FAILED \n
                          RET_STEPLENGTH_DETERMINATION_FAILED \n
                          RET_HOMOTOPY_STEP_FAILED \n
                          RET_HOTSTART_STOPPED_INFEASIBILITY \n
                          RET_HOTSTART_STOPPED_UNBOUNDEDNESS \n
                          RET_INACCURATE_SOLUTION \n
                          RET_NO_SOLUTION */
  returnValue hotstart(int& nWSR, /**< Input: Maximum number of working set
                                    recalculations; \n  Output: Number of
                                    performed working set recalculations. */
                       const QProblemBCallback& callback);

  /** Returns the number of free variables.
   *  \return Number of free variables. */
  inline int getNFR();

  /** Returns the number of fixed variables.
   *  \return Number of fixed variables. */
  inline int getNFX();

  /** Returns the number of implicitly fixed variables.
   *  \return Number of implicitly fixed variables. */
  inline int getNFV() const;

  /** Returns the dimension of null space.
   *  \return Dimension of null space. */
  int getNZ();

  /** Returns the optimal objective function value.
   *  \return finite value: Optimal objective function value (QP was solved)
   \n +infinity:    QP was not yet solved */
  real_t getObjVal() const;

  /** Returns the objective function value at an arbitrary point x.
   *  \return Objective function value at point x */
  real_t getObjVal(const real_t* const _x /**< Point at which the objective
                                             function shall be evaluated. */
                   ) const;

  /** Returns status of the solution process.
   *  \return Status of solution process. */
  inline QProblemStatus getStatus() const;

  /** Returns if the QProblem object is initialised.
   *  \return true:  QProblemB initialised \n
                          false: QProblemB not initialised */
  inline bool isInitialised() const;

  /** Returns if the QP has been solved.
   *  \return true:  QProblemB solved \n
                          false: QProblemB not solved */
  inline bool isSolved() const;

  /** Returns if the QP is infeasible.
   *  \return true:  QP infeasible \n
                          false: QP feasible (or not known to be infeasible!)
 */
  inline bool isInfeasible() const;

  /** Returns if the QP is unbounded.
   *  \return true:  QP unbounded \n
                          false: QP unbounded (or not known to be unbounded!)
 */
  inline bool isUnbounded() const;

  /** Returns Hessian type flag (type is not determined due to this call!).
   *  \return Hessian type. */
  inline HessianType getHessianType() const;

  /** Changes the print level.
   *  \return SUCCESSFUL_RETURN */
  inline returnValue setHessianType(
      HessianType _hessianType /**< New Hessian type. */
  );

  const std::vector<real_t>& getH() const { return H; }
  std::vector<real_t>* getMutableH() { return &H; }

  const std::vector<real_t>& getG() const { return g_user; }
  std::vector<real_t>* getMutableG() { return &g_user; }

  const std::vector<real_t>& getLb() const { return lb_user; }
  std::vector<real_t>* getMutableLb() { return &lb_user; }

  const std::vector<real_t>& getUb() const { return ub_user; }
  std::vector<real_t>* getMutableUb() { return &ub_user; }

  const std::vector<real_t>& getX() const { return x; }
  std::vector<real_t>* getMutableX() { return &x; }

  const std::vector<real_t>& getY() const { return y; }
  std::vector<real_t>* getMutableY() { return &y; }

  /*
   *  PROTECTED MEMBER FUNCTIONS
   */
 protected:
  QProblemB(size_t _nV, /**< Number of variables. */
            size_t _nC /**< Number of constraints. */);

  /** Checks if Hessian happens to be the identity matrix,
   *  and sets corresponding status flag (otherwise the flag remains
   *unaltered!). \return SUCCESSFUL_RETURN */
  returnValue checkForIdentityHessian();

  /** Determines type of constraints and bounds (i.e. implicitly fixed,
   unbounded etc.). *  \return SUCCESSFUL_RETURN \n
                          RET_SETUPSUBJECTTOTYPE_FAILED */
  returnValue setupSubjectToType();

  /** Computes the Cholesky decomposition R of the (simply projected) Hessian
   *(i.e. R^T*R = Z^T*H*Z). It only works in the case where Z is a simple
   *projection matrix! \return SUCCESSFUL_RETURN \n RET_INDEXLIST_CORRUPTED */
  returnValue setupCholeskyDecomposition();

  /** Solves a QProblemB whose QP data is assumed to be stored in the member
   variables.
   *  A guess for its primal/dual optimal solution vectors and the corresponding
   *  optimal working set can be provided.
   *  \return SUCCESSFUL_RETURN \n
                          RET_INIT_FAILED \n
                          RET_INIT_FAILED_CHOLESKY \n
                          RET_INIT_FAILED_HOTSTART \n
                          RET_INIT_FAILED_INFEASIBILITY \n
                          RET_INIT_FAILED_UNBOUNDEDNESS \n
                          RET_MAX_NWSR_REACHED */
  returnValue solveInitialQP(
      int& nWSR, /**< Input: Maximum number of working set recalculations; \n
                  *   Output: Number of performed working set recalculations.
                  */
      const QProblemBCallback& callback);

  /** Obtains the desired working set for the auxiliary initial QP in
   *  accordance with the user specifications
   *  \return SUCCESSFUL_RETURN \n
                          RET_OBTAINING_WORKINGSET_FAILED \n
                          RET_INVALID_ARGUMENTS */
  returnValue obtainAuxiliaryWorkingSet();

  /** Setups bound data structure according to auxiliaryBounds.
   *  (If the working set shall be setup afresh, make sure that
   *  bounds data structure has been resetted!)
   *  \return SUCCESSFUL_RETURN \n
                          RET_SETUP_WORKINGSET_FAILED \n
                          RET_INVALID_ARGUMENTS \n
                          RET_UNKNOWN BUG */
  returnValue setupAuxiliaryWorkingSet(
      bool setupAfresh /**< Flag indicating if given working set shall be
                        *    setup afresh or by updating the current one. */
  );

  /** Setups the optimal primal/dual solution of the auxiliary initial QP.
   *  \return SUCCESSFUL_RETURN */
  returnValue setupAuxiliaryQPsolution();

  /** Setups gradient of the auxiliary initial QP for given
   *  optimal primal/dual solution and given initial working set
   *  (assumes that members X, Y and BOUNDS have already been initialised!).
   *  \return SUCCESSFUL_RETURN */
  returnValue setupAuxiliaryQPgradient();

  /** Setups bounds of the auxiliary initial QP for given
   *  optimal primal/dual solution and given initial working set
   *  (assumes that members X, Y and BOUNDS have already been initialised!).
   *  \return SUCCESSFUL_RETURN \n
                          RET_UNKNOWN BUG */
  returnValue setupAuxiliaryQPbounds(
      bool useRelaxation /**< Flag indicating if inactive bounds shall be
                                   relaxed. */
  );

  /** Adds a bound to active set (specialised version for the case where no
   constraints exist). *  \return SUCCESSFUL_RETURN \n
                          RET_ADDBOUND_FAILED */
  returnValue addBound(
      int number, /**< Number of bound to be added to active set. */
      SubjectToStatus B_status, /**< Status of new active bound. */
      bool updateCholesky       /**< Flag indicating if Cholesky
                                          decomposition shall be updated. */
  );

  /** Removes a bounds from active set (specialised version for the case where
   no constraints exist). *  \return SUCCESSFUL_RETURN \n RET_HESSIAN_NOT_SPD
   \n RET_REMOVEBOUND_FAILED */
  returnValue removeBound(
      int number,         /**< Number of bound to be removed from active set. */
      bool updateCholesky /**< Flag indicating if Cholesky
                                    decomposition shall be updated. */
  );

  /** Solves the system Ra = b or R^Ta = b where R is an upper triangular
   matrix. *  \return SUCCESSFUL_RETURN \n RET_DIV_BY_ZERO */
  returnValue backsolveR(const real_t* const b, /**< Right hand side vector. */
                         bool transposed,       /**< Indicates if the transposed
                                                   system       shall be solved. */
                         real_t* const a        /**< Output: Solution vector */
  );

  /** Solves the system Ra = b or R^Ta = b where R is an upper triangular
   matrix. \n
   *  Special variant for the case that this function is called from within
   "removeBound()". *  \return SUCCESSFUL_RETURN \n RET_DIV_BY_ZERO */
  returnValue backsolveR(
      const real_t* const b, /**< Right hand side vector. */
      bool transposed,       /**< Indicates if the transposed
                                system       shall be solved. */
      bool removingBound,    /**< Indicates if function is
                                called    from "removeBound()". */
      real_t* const a        /**< Output: Solution vector */
  );

  /** Determines step direction of the shift of the QP data.
   *  \return SUCCESSFUL_RETURN */
  returnValue hotstart_determineDataShift(
      const int* const FX_idx,    /**< Index array of fixed variables. */
      const real_t* const g_new,  /**< New gradient vector. */
      const real_t* const lb_new, /**< New lower bounds. */
      const real_t* const ub_new, /**< New upper bounds. */
      real_t* const delta_g,  /**< Output: Step direction of gradient vector. */
      real_t* const delta_lb, /**< Output: Step direction of lower bounds. */
      real_t* const delta_ub, /**< Output: Step direction of upper bounds. */
      bool& Delta_bB_isZero   /**< Output: Indicates if active bounds are
                                        to be shifted. */
  );

  /** Checks if lower/upper bounds remain consistent
   *  (i.e. if lb <= ub) during the current step.
   *  \return true iff bounds remain consistent
   */
  bool areBoundsConsistent(
      const real_t* const delta_lb, /**< Step direction of lower bounds. */
      const real_t* const delta_ub  /**< Step direction of upper bounds. */
      ) const;

  /** Computes parameters for the Givens matrix G for which [x,y]*G = [z,0]
   *  \return SUCCESSFUL_RETURN */
  inline void computeGivens(
      real_t xold,  /**< Matrix entry to be normalised. */
      real_t yold,  /**< Matrix entry to be annihilated. */
      real_t& xnew, /**< Output: Normalised matrix entry. */
      real_t& ynew, /**< Output: Annihilated matrix entry. */
      real_t& c,    /**< Output: Cosine entry of Givens matrix. */
      real_t& s     /**< Output: Sine entry of Givens matrix. */
      ) const;

  /** Applies Givens matrix determined by c and s (cf. computeGivens).
   *  \return SUCCESSFUL_RETURN */
  inline void applyGivens(
      real_t c,     /**< Cosine entry of Givens matrix. */
      real_t s,     /**< Sine entry of Givens matrix. */
      real_t xold,  /**< Matrix entry to be transformed corresponding to
                     *   the normalised entry of the original matrix. */
      real_t yold,  /**< Matrix entry to be transformed corresponding to
                     *   the annihilated entry of the original matrix. */
      real_t& xnew, /**< Output: Transformed matrix entry corresponding to
                     *   the normalised entry of the original matrix. */
      real_t& ynew  /**< Output: Transformed matrix entry corresponding to
                     *   the annihilated entry of the original matrix. */
      ) const;

  /*
   *  PRIVATE MEMBER FUNCTIONS
   */
 private:
  /** Determines step direction of the homotopy path.
   *  \return SUCCESSFUL_RETURN \n
                          RET_STEPDIRECTION_FAILED_CHOLESKY */
  returnValue hotstart_determineStepDirection(
      const int* const FR_idx,      /**< Index array of free variables. */
      const int* const FX_idx,      /**< Index array of fixed variables. */
      const real_t* const delta_g,  /**< Step direction of gradient vector. */
      const real_t* const delta_lb, /**< Step direction of lower bounds. */
      const real_t* const delta_ub, /**< Step direction of upper bounds. */
      bool
          Delta_bB_isZero, /**< Indicates if active bounds are to be shifted. */
      real_t* const delta_xFX, /**< Output: Primal homotopy step direction of
                                  fixed variables. */
      real_t* const delta_xFR, /**< Output: Primal homotopy step direction of
                                  free variables. */
      real_t* const delta_yFX /**< Output: Dual homotopy step direction of fixed
                                 variables' multiplier. */
  );

  /** Determines the maximum possible step length along the homotopy path.
   *  \return SUCCESSFUL_RETURN */
  returnValue hotstart_determineStepLength(
      const int* const FR_idx,      /**< Index array of free variables. */
      const int* const FX_idx,      /**< Index array of fixed variables. */
      const real_t* const delta_lb, /**< Step direction of lower bounds. */
      const real_t* const delta_ub, /**< Step direction of upper bounds. */
      const real_t* const
          delta_xFR, /**< Primal homotopy step direction of free variables. */
      const real_t* const delta_yFX, /**< Dual homotopy step direction of fixed
                                        variables' multiplier. */
      int& BC_idx,               /**< Output: Index of blocking constraint. */
      SubjectToStatus& BC_status /**< Output: Status of blocking constraint. */
  );

  /** Performs a step along the homotopy path (and updates active set).
   *  \return  SUCCESSFUL_RETURN \n
                           RET_OPTIMAL_SOLUTION_FOUND \n
                           RET_REMOVE_FROM_ACTIVESET_FAILED \n
                           RET_ADD_TO_ACTIVESET_FAILED \n
                           RET_QP_INFEASIBLE */
  returnValue hotstart_performStep(
      const int* const FR_idx,      /**< Index array of free variables. */
      const int* const FX_idx,      /**< Index array of fixed variables. */
      const real_t* const delta_g,  /**< Step direction of gradient vector. */
      const real_t* const delta_lb, /**< Step direction of lower bounds. */
      const real_t* const delta_ub, /**< Step direction of upper bounds. */
      const real_t* const
          delta_xFX, /**< Primal homotopy step direction of fixed variables. */
      const real_t* const
          delta_xFR, /**< Primal homotopy step direction of free variables. */
      const real_t* const delta_yFX, /**< Dual homotopy step direction of fixed
                                        variables' multiplier. */
      int BC_idx,                    /**< Index of blocking constraint. */
      SubjectToStatus BC_status      /**< Status of blocking constraint. */
  );

  /** Determines the maximum violation of the KKT optimality conditions
   *  of the current iterate within the QProblemB object.
   *  \return SUCCESSFUL_RETURN \n
   *       RET_INACCURATE_SOLUTION \n
   *       RET_NO_SOLUTION */
  returnValue checkKKTconditions();

  /*
   *  PROTECTED MEMBER VARIABLES
   */
 protected:
  std::vector<real_t> H;  /**< Hessian matrix. */
  std::vector<real_t> g;  /**< Gradient. */
  std::vector<real_t> lb; /**< Lower bound vector (on variables). */
  std::vector<real_t> ub; /**< Upper bound vector (on variables). */

  Bounds bounds; /**< Data structure for problem's bounds. */
  Bounds auxiliaryBounds;

  std::vector<real_t> R; /**< Cholesky decomposition of H (i.e. H = R^T*R). */

  std::vector<real_t> x; /**< Primal solution vector. */
  std::vector<real_t> y; /**< Dual solution vector. */

  real_t tau; /**< Last homotopy step length. */

  QProblemStatus status; /**< Current status of the solution process. */

  bool infeasible; /**< QP infeasible? */
  bool unbounded;  /**< QP unbounded? */

  HessianType hessianType; /**< Type of Hessian matrix. */

  std::vector<real_t> g_user;
  std::vector<real_t> lb_user;
  std::vector<real_t> ub_user;

  std::vector<real_t> delta_g;
  std::vector<real_t> delta_lb;
  std::vector<real_t> delta_ub;

  std::vector<real_t> delta_xFR;
  std::vector<real_t> delta_xFX;
  std::vector<real_t> delta_yFX;

  std::vector<real_t> rhs;
  std::vector<real_t> r;

  std::vector<real_t> delta_xFRz_TMP;
  std::vector<real_t> delta_xFRz_RHS;

 private:
  std::vector<real_t> HMX_delta_xFX;
};

}  // namespace qpoases_embedded

#include <qpoases_embedded/QProblemB.ipp>

#endif /* QPOASES_QPROBLEMB_HPP */

/*
 *  end of file
 */
