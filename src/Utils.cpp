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

#include <qpoases_embedded/Utils.hpp>

#ifdef QPOASES_DEBUG

#include <cassert>
#include <cstdio>

namespace qpoases_embedded {

/*
 *  p r i n t I t e r a t i o n
 */
bool printIterationB(int iteration, real_t tau, int nFX, int BC_idx,
                     SubjectToStatus BC_status) {
  assert(iteration >= 0);
  assert(tau >= 0);
  assert(nFX >= 0);
  assert(BC_idx >= 0);

  /* 1) Print header. */
  if (iteration % 10 == 0) {
    printf("  Iter  |  StepLength   |     Info      |  nFX   \n");
  }

  /* 2) Print iteration line. */
  if (BC_status == ST_UNDEFINED) {
    printf("  %4.1d  |  %1.5e  |   QP SOLVED   | %4.1d   \n", iteration, tau,
           nFX);

  } else {
    char info[8];

    if (BC_status == ST_INACTIVE)
      sprintf(info, "REM BND");
    else
      sprintf(info, "ADD BND");

    printf("  %4.1d  |  %1.5e  |  %s%4.1d  | %4.1d   \n", iteration, tau, info,
           BC_idx, nFX);
  }

  return true;
}

/*
 *  p r i n t I t e r a t i o n
 */
bool printIteration(int iteration, real_t tau, int nFX, int nAC, int BC_idx,
                    SubjectToStatus BC_status, bool BC_isBound) {
  assert(iteration >= 0);
  assert(tau >= 0);
  assert(nFX >= 0);
  assert(nAC >= 0);
  assert(BC_idx >= 0);

  /* 1) Print header. */
  if (iteration % 10 == 0) {
    printf("  Iter  |  StepLength   |     Info      |  nFX  |  nAC  \n");
  }

  /* 2) Print iteration line. */
  if (BC_status == ST_UNDEFINED) {
    printf("  %4.1d  |  %1.5e  |   QP SOLVED   | %4.1d  | %4.1d  \n", iteration,
           tau, nFX, nAC);

  } else {
    char info[8];

    if (BC_status == ST_INACTIVE)
      sprintf(info, "REM ");
    else
      sprintf(info, "ADD ");

    if (BC_isBound == true)
      sprintf(&(info[4]), "BND");
    else
      sprintf(&(info[4]), "CON");

    printf("  %4.1d  |  %1.5e  |  %s%4.1d  | %4.1d  | %4.1d  \n", iteration,
           tau, info, BC_idx, nFX, nAC);
  }

  return true;
}

/*
 *  p r i n t
 */
returnValue print(const real_t* const v, int n) {
  int i;

  /* Print a vector. */

  printf("%s", "[\t");
  for (i = 0; i < n; ++i) {
    printf(" %.16e\t", v[i]);
  }
  printf("%s", "]\n");

  return SUCCESSFUL_RETURN;
}

/*
 *  p r i n t
 */
returnValue print(const real_t* const v, int n, const int* const V_idx) {
  int i;

  /* Print a permuted vector. */
  printf("%s", "[\t");
  for (i = 0; i < n; ++i) {
    printf(" %.16e\t", v[V_idx[i]]);
  }
  printf("%s", "]\n");

  return SUCCESSFUL_RETURN;
}

/*
 *  p r i n t
 */
returnValue print(const real_t* const v, int n, const char* name) {
  /* Print vector name ... */
  printf("%s = ", name);

  /* ... and the vector itself. */
  return print(v, n);
}

/*
 *  p r i n t
 */
returnValue print(const real_t* const M, int nrow, int ncol) {
  int i;

  /* Print a matrix as a collection of row vectors. */
  for (i = 0; i < nrow; ++i) print(&(M[i * ncol]), ncol);
  printf("%s", "\n");

  return SUCCESSFUL_RETURN;
}

/*
 *  p r i n t
 */
returnValue print(const real_t* const M, int nrow, int ncol,
                  const int* const ROW_idx, const int* const COL_idx) {
  int i;

  /* Print a permuted matrix as a collection of permuted row vectors. */
  for (i = 0; i < nrow; ++i) print(&(M[ROW_idx[i] * ncol]), ncol, COL_idx);
  printf("%s", "\n");

  return SUCCESSFUL_RETURN;
}

/*
 *  p r i n t
 */
returnValue print(const real_t* const M, int nrow, int ncol, const char* name) {
  /* Print matrix name ... */
  printf("%s = ", name);

  /* ... and the matrix itself. */
  return print(M, nrow, ncol);
}

/*
 *  p r i n t
 */
returnValue print(const int* const index, int n) {
  int i;

  /* Print a indexlist. */
  printf("%s", "[\t");
  for (i = 0; i < n; ++i) {
    printf(" %d\t", index[i]);
  }
  printf("%s", "]\n");

  return SUCCESSFUL_RETURN;
}

/*
 *  p r i n t
 */
returnValue print(const int* const index, int n, const char* name) {
  /* Print indexlist name ... */
  printf("%s = ", name);

  /* ... and the indexlist itself. */
  return print(index, n);
}

}  // namespace qpoases_embedded

#endif /* QPOASES_DEBUG */

/*
 *  end of file
 */
