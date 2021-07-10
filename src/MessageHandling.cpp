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

#include <qpoases_embedded/MessageHandling.hpp>

#include <map>

#ifdef QPOASES_DEBUG

#include <cstdarg>
#include <cstdio>

#endif  // QPOASES_DEBUG

namespace qpoases_embedded {

static const std::map<returnValue, std::string> returnValueErrorToString = {
    /* miscellaneous */
    {SUCCESSFUL_RETURN, "Successful return"},
    {RET_DIV_BY_ZERO, "Division by zero"},
    {RET_INDEX_OUT_OF_BOUNDS, "Index out of bounds"},
    {RET_INVALID_ARGUMENTS, "At least one of the arguments is invalid"},
    {RET_ERROR_UNDEFINED, "Error number undefined"},
    {RET_WARNING_UNDEFINED, "Warning number undefined"},
    {RET_INFO_UNDEFINED, "Info number undefined"},
    {RET_EWI_UNDEFINED, "Error/warning/info number undefined"},
    {RET_AVAILABLE_WITH_LINUX_ONLY,
     "This function is available under Linux only"},
    {RET_UNKNOWN_BUG, "The error occured is not yet known"},
    {RET_PRINTLEVEL_CHANGED, "Print level changed"},
    {
        RET_NOT_YET_IMPLEMENTED,
        "Requested function is not yet implemented.",
    },
    /* Indexlist */
    {
        RET_INDEXLIST_MUST_BE_REORDERD,
        "Index list has to be reordered",
    },
    {RET_INDEXLIST_EXCEEDS_MAX_LENGTH,
     "Index list exceeds its maximal physical length"},
    {RET_INDEXLIST_CORRUPTED, "Index list corrupted"},
    {RET_INDEXLIST_OUTOFBOUNDS, "Physical index is out of bounds"},
    {RET_INDEXLIST_ADD_FAILED, "Adding indices from another index set failed"},
    {RET_INDEXLIST_INTERSECT_FAILED,
     "Intersection with another index set failed"},
    /* SubjectTo / Bounds / Constraints */
    {RET_INDEX_ALREADY_OF_DESIRED_STATUS, "Index is already of desired status"},
    {RET_SWAPINDEX_FAILED, "Cannot swap between different indexsets"},
    {RET_ADDINDEX_FAILED, "Adding index to index set failed"},
    {RET_NOTHING_TO_DO, "Nothing to do"},
    {RET_SETUP_BOUND_FAILED, "Setting up bound index failed"},
    {RET_SETUP_CONSTRAINT_FAILED, "Setting up constraint index failed"},
    {RET_MOVING_BOUND_FAILED, "Moving bound between index sets failed"},
    {RET_MOVING_CONSTRAINT_FAILED,
     "Moving constraint between index sets failed"},
    /* QProblem */
    {RET_QP_ALREADY_INITIALISED, "QProblem has already been initialised"},
    {RET_NO_INIT_WITH_STANDARD_SOLVER,
     "Initialisation via extern QP solver is not yet implemented"},
    {RET_RESET_FAILED, "Reset failed"},
    {RET_INIT_FAILED, "Initialisation failed"},
    {RET_INIT_FAILED_TQ, "Initialisation failed due to TQ factorisation"},
    {RET_INIT_FAILED_CHOLESKY,
     "Initialisation failed due to Cholesky decomposition"},
    {RET_INIT_FAILED_HOTSTART,
     "Initialisation failed! QP could not be solved!"},
    {RET_INIT_FAILED_INFEASIBILITY,
     "Initial QP could not be solved due to infeasibility!"},
    {RET_INIT_FAILED_UNBOUNDEDNESS,
     "Initial QP could not be solved due to unboundedness!"},
    {RET_INIT_SUCCESSFUL, "Initialisation done"},
    {RET_OBTAINING_WORKINGSET_FAILED,
     "Failed to obtain working set for auxiliary QP"},
    {RET_SETUP_WORKINGSET_FAILED,
     "Failed to setup working set for auxiliary QP"},
    {RET_SETUP_AUXILIARYQP_FAILED,
     "Failed to setup auxiliary QP for initialised homotopy"},
    {RET_NO_EXTERN_SOLVER, "No extern QP solver available"},
    {RET_QP_UNBOUNDED, "QP is unbounded"},
    {RET_QP_INFEASIBLE, "QP is infeasible"},
    {RET_QP_NOT_SOLVED,
     "Problems occured while solving QP with standard solver"},
    {RET_QP_SOLVED, "QP successfully solved"},
    {RET_UNABLE_TO_SOLVE_QP, "Problems occured while solving QP"},
    {RET_INITIALISATION_STARTED, "Starting problem initialisation..."},
    {RET_HOTSTART_FAILED, "Unable to perform homotopy due to internal error"},
    {RET_HOTSTART_FAILED_TO_INIT, "Unable to initialise problem"},
    {RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED,
     "Unable to perform homotopy as previous QP is not solved"},
    {RET_ITERATION_STARTED, "Iteration"},
    {RET_SHIFT_DETERMINATION_FAILED,
     "Determination of shift of the QP data failed"},
    {RET_STEPDIRECTION_DETERMINATION_FAILED,
     "Determination of step direction failed"},
    {RET_STEPLENGTH_DETERMINATION_FAILED,
     "Determination of step direction failed"},
    {RET_OPTIMAL_SOLUTION_FOUND, "Optimal solution of neighbouring QP found"},
    {RET_HOMOTOPY_STEP_FAILED, "Unable to perform homotopy step"},
    {RET_HOTSTART_STOPPED_INFEASIBILITY,
     "Premature homotopy termination because QP is infeasible"},
    {RET_HOTSTART_STOPPED_UNBOUNDEDNESS,
     "Premature homotopy termination because QP is unbounded"},
    {RET_WORKINGSET_UPDATE_FAILED,
     "Unable to update working sets according to initial guesses"},
    {RET_MAX_NWSR_REACHED,
     "Maximum number of working set recalculations performed"},
    {RET_CONSTRAINTS_NOT_SPECIFIED,
     "Problem does comprise constraints! You have to specify new constraints' "
     "bounds"},
    {RET_INVALID_FACTORISATION_FLAG, "Invalid factorisation flag"},
    {RET_UNABLE_TO_SAVE_QPDATA, "Unable to save QP data"},
    {RET_STEPDIRECTION_FAILED_TQ,
     "Abnormal termination due to TQ factorisation"},
    {RET_STEPDIRECTION_FAILED_CHOLESKY,
     "Abnormal termination due to Cholesky factorisation"},
    {RET_STEPSIZE, ""},
    {RET_STEPSIZE_NONPOSITIVE, ""},
    {RET_SETUPSUBJECTTOTYPE_FAILED, "Setup of SubjectToTypes failed"},
    {RET_ADDCONSTRAINT_FAILED, "Addition of constraint to working set failed"},
    {RET_ADDCONSTRAINT_FAILED_INFEASIBILITY,
     "Addition of constraint to working set failed"},
    {RET_ADDBOUND_FAILED, "Addition of bound to working set failed"},
    {RET_ADDBOUND_FAILED_INFEASIBILITY,
     "Addition of bound to working set failed"},
    {RET_REMOVECONSTRAINT_FAILED,
     "Removal of constraint from working set failed"},
    {RET_REMOVEBOUND_FAILED, "Removal of bound from working set failed"},
    {RET_REMOVE_FROM_ACTIVESET, "Removing from active set:"},
    {RET_ADD_TO_ACTIVESET, "Adding to active set:"},
    {RET_REMOVE_FROM_ACTIVESET_FAILED, "Removing from active set failed"},
    {RET_ADD_TO_ACTIVESET_FAILED, "Adding to active set failed"},
    {RET_CONSTRAINT_ALREADY_ACTIVE, "Constraint is already active"},
    {RET_ALL_CONSTRAINTS_ACTIVE,
     "All constraints are active, no further constraint can be added"},
    {RET_LINEARLY_DEPENDENT, "New bound/constraint is linearly dependent"},
    {RET_LINEARLY_INDEPENDENT, "New bound/constraint is linearly independent"},
    {RET_LI_RESOLVED,
     "Linear independence of active contraint matrix successfully resolved"},
    {RET_ENSURELI_FAILED,
     "Failed to ensure linear indepence of active contraint matrix"},
    {RET_ENSURELI_FAILED_TQ, "Abnormal termination due to TQ factorisation"},
    {RET_ENSURELI_FAILED_NOINDEX, "No index found, QP is probably infeasible"},
    {RET_BOUND_ALREADY_ACTIVE, "Bound is already active"},
    {RET_ALL_BOUNDS_ACTIVE,
     "All bounds are active, no further bound can be added"},
    {RET_CONSTRAINT_NOT_ACTIVE, "Constraint is not active"},
    {RET_BOUND_NOT_ACTIVE, "Bound is not active"},
    {RET_HESSIAN_NOT_SPD, "Projected Hessian matrix not positive definite"},
    {RET_MATRIX_SHIFT_FAILED,
     "Unable to update matrices or to transform vectors"},
    {RET_MATRIX_FACTORISATION_FAILED,
     "Unable to calculate new matrix factorisations"},
    {RET_USER_ABORT_REQUESTED, "User requested early termination"},
    /* SolutionAnalysis */
    {RET_NO_SOLUTION, "QP solution does not satisfy KKT optimality conditions"},
    {RET_INACCURATE_SOLUTION,
     "KKT optimality conditions not satisfied to sufficient accuracy"}};

const std::string& getErrorString(returnValue error) {
  return returnValueErrorToString.at(error);
}

#ifdef QPOASES_DEBUG

namespace internal {

static const std::map<MsgType, std::string> kMsgTypeEnumToString = {
    {MsgType::kInfo, "INFO"}, {MsgType::kError, "ERROR"}};

returnValue throwMessage(internal::MsgType msgType, returnValue retval,
                         const char* filename, const unsigned long linenumber,
                         const char* functionname) {
  const std::string& returnValueString = returnValueErrorToString.at(retval);
  printf("%s %s:%4d:%s: %s.\n", kMsgTypeEnumToString.at(msgType).c_str(),
         filename, (int)linenumber, functionname, returnValueString.c_str());
  return retval;
}

returnValue throwMessage(internal::MsgType msgType, returnValue retval,
                         const char* filename, const unsigned long linenumber,
                         const char* functionname, const char* fmt, ...) {
  const std::string& returnValueString = returnValueErrorToString.at(retval);
  printf("%s %s:%4d:%s: %s ", kMsgTypeEnumToString.at(msgType).c_str(),
         filename, (int)linenumber, functionname, returnValueString.c_str());
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  printf("\n");
  return retval;
}

}  // namespace internal

#endif /* QPOASES_DEBUG */

}  // namespace qpoases_embedded

/*
 *  end of file
 */
