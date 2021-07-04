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

#ifndef QPOASES_UTILS_HPP
#define QPOASES_UTILS_HPP

#include <qpoases_embedded/Constants.hpp>
#include <qpoases_embedded/Types.hpp>

namespace qpoases_embedded {

/** Prints concise information on the current iteration. Use with QProblemB.
 *  \returns Always returns true. */
bool printIterationB(int iteration, /**< Number of current iteration. */
                     real_t tau,    /**< The last homotopy step length. */
                     int nFX,       /**< The number of fixed variables. */
                     int BC_idx,    /**< Index of blocking bound. */
                     SubjectToStatus BC_status /**< Status of blocking bound. */
);

/** Prints concise information on the current iteration. Use with QProblem.
 *  \returns Always returns true. */
bool printIteration(
    int iteration,             /**< Number of current iteration. */
    real_t tau,                /**< The last homotopy step length. */
    int nFX,                   /**< The number of fixed variables. */
    int nAC,                   /**< The number of active constraints. */
    int BC_idx,                /**< Index of blocking constraint. */
    SubjectToStatus BC_status, /**< Status of blocking constraint. */
    bool BC_isBound /**< Indicates if blocking constraint is a bound. */
);

/** Prints a vector.
 * \return SUCCESSFUL_RETURN */
returnValue print(const real_t* const v, /**< Vector to be printed. */
                  int n                  /**< Length of vector. */
);

/** Prints a permuted vector.
 * \return SUCCESSFUL_RETURN */
returnValue print(const real_t* const v, /**< Vector to be printed. */
                  int n,                 /**< Length of vector. */
                  const int* const V_idx /**< Pemutation vector. */
);

/** Prints a named vector.
 * \return SUCCESSFUL_RETURN */
returnValue print(const real_t* const v, /**< Vector to be printed. */
                  int n,                 /**< Length of vector. */
                  const char* name       /** Name of vector. */
);

/** Prints a matrix.
 * \return SUCCESSFUL_RETURN */
returnValue print(const real_t* const M, /**< Matrix to be printed. */
                  int nrow,              /**< Row number of matrix. */
                  int ncol               /**< Column number of matrix. */
);

/** Prints a permuted matrix.
 * \return SUCCESSFUL_RETURN */
returnValue print(const real_t* const M,    /**< Matrix to be printed. */
                  int nrow,                 /**< Row number of matrix. */
                  int ncol,                 /**< Column number of matrix. */
                  const int* const ROW_idx, /**< Row pemutation vector. */
                  const int* const COL_idx  /**< Column pemutation vector. */
);

/** Prints a named matrix.
 * \return SUCCESSFUL_RETURN */
returnValue print(const real_t* const M, /**< Matrix to be printed. */
                  int nrow,              /**< Row number of matrix. */
                  int ncol,              /**< Column number of matrix. */
                  const char* name       /** Name of matrix. */
);

/** Prints an index array.
 * \return SUCCESSFUL_RETURN */
returnValue print(const int* const index, /**< Index array to be printed. */
                  int n                   /**< Length of index array. */
);

/** Prints a named index array.
 * \return SUCCESSFUL_RETURN */
returnValue print(const int* const index, /**< Index array to be printed. */
                  int n,                  /**< Length of index array. */
                  const char* name        /**< Name of index array. */
);

}  // namespace qpoases_embedded

#endif /* QPOASES_UTILS_HPP */

/*
 *  end of file
 */
