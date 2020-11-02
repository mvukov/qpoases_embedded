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

#ifndef QPOASES_SUBJECTTO_HPP
#define QPOASES_SUBJECTTO_HPP

#include <qpoases_embedded/Indexlist.hpp>

namespace qpoases_embedded {

/** This class manages working sets of constraints and bounds by storing
 *  index sets and other status information.
 */
class SubjectTo {
  /*
   *  PUBLIC MEMBER FUNCTIONS
   */
 public:
  SubjectTo() = delete;
  explicit SubjectTo(int size);

  void reset();

  int getSize() const { return type.size(); }

  /** Returns type of (constraints') bound.
   *  \return Type of (constraints') bound \n
                          RET_INDEX_OUT_OF_BOUNDS */
  inline SubjectToType getType(int i /**< Number of (constraints') bound. */
                               ) const;

  /** Returns status of (constraints') bound.
   *  \return Status of (constraints') bound \n
                          ST_UNDEFINED */
  inline SubjectToStatus getStatus(int i /**< Number of (constraints') bound. */
                                   ) const;

  /** Sets type of (constraints') bound.
   *  \return SUCCESSFUL_RETURN \n
                          RET_INDEX_OUT_OF_BOUNDS */
  inline returnValue setType(
      int i,              /**< Number of (constraints') bound. */
      SubjectToType value /**< Type of (constraints') bound. */
  );

  /** Sets status of (constraints') bound.
   *  \return SUCCESSFUL_RETURN \n
                          RET_INDEX_OUT_OF_BOUNDS */
  inline returnValue setStatus(
      int i,                /**< Number of (constraints') bound. */
      SubjectToStatus value /**< Status of (constraints') bound. */
  );

  /** Sets status of lower (constraints') bounds. */
  inline void setNoLower(
      bool _status /**< Status of lower (constraints') bounds. */
  );

  /** Sets status of upper (constraints') bounds. */
  inline void setNoUpper(
      bool _status /**< Status of upper (constraints') bounds. */
  );

  /** Returns status of lower (constraints') bounds.
   *  \return true if there is no lower (constraints') bound on any
   *variable. */
  inline bool isNoLower() const;

  /** Returns status of upper bounds.
   *  \return true if there is no upper (constraints') bound on any
   *variable. */
  inline bool isNoUpper() const;

  /*
   *  PROTECTED MEMBER FUNCTIONS
   */
 protected:
  /** Adds the index of a new constraint or bound to index set.
   *  \return SUCCESSFUL_RETURN \n
                          RET_ADDINDEX_FAILED */
  returnValue addIndex(
      Indexlist* const
          indexlist, /**< Index list to which the new index shall be added. */
      int newnumber, /**< Number of new constraint or bound. */
      SubjectToStatus newstatus /**< Status of new constraint or bound. */
  );

  /** Removes the index of a constraint or bound from index set.
   *  \return SUCCESSFUL_RETURN \n
                          RET_UNKNOWN_BUG */
  returnValue removeIndex(
      Indexlist* const indexlist, /**< Index list from which the new index shall
                                     be removed. */
      int removenumber /**< Number of constraint or bound to be removed. */
  );

  /** Swaps the indices of two constraints or bounds within the index set.
   *  \return SUCCESSFUL_RETURN \n
                          RET_SWAPINDEX_FAILED */
  returnValue swapIndex(
      Indexlist* const
          indexlist, /**< Index list in which the indices shold be swapped. */
      int number1,   /**< Number of first constraint or bound. */
      int number2    /**< Number of second constraint or bound. */
  );

  /*
   *  PROTECTED MEMBER VARIABLES
   */
 protected:
  std::vector<SubjectToType> type;     /**< Type of constraints/bounds. */
  std::vector<SubjectToStatus> status; /**< Status of constraints/bounds. */

  bool noLower; /**< This flag indicates if there is no lower bound on
                          any variable. */
  bool noUpper; /**< This flag indicates if there is no upper bound on
                          any variable. */
};

}  // namespace qpoases_embedded

#include <qpoases_embedded/SubjectTo.ipp>

#endif /* QPOASES_SUBJECTTO_HPP */

/*
 *  end of file
 */
