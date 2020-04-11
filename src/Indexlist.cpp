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

#include <qpoases_embedded/Indexlist.hpp>
#include <qpoases_embedded/MessageHandling.hpp>

/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

namespace qpoases_embedded {

/*
 *  I n d e x l i s t
 */
Indexlist::Indexlist(int physicallength)
    : number(physicallength),
      next(physicallength),
      previous(physicallength),
      numberArray(physicallength) {
  reset();
}

void Indexlist::reset() {
  std::fill(number.begin(), number.end(), -1);
  std::fill(next.begin(), next.end(), -1);
  std::fill(previous.begin(), previous.end(), -1);
  std::fill(numberArray.begin(), numberArray.end(), -1);

  length = 0;
  first = -1;
  last = -1;
  lastusedindex = -1;
}

/*
 *  g e t I n d e x
 */
int Indexlist::getIndex(int givennumber) const {
  int i;
  int n = first;
  int index = -1; /* return -1 by default */

  /* Run trough indexlist until number is found, if so return it index. */
  for (i = 0; i < length; ++i) {
    if (number[n] == givennumber) {
      index = i;
      break;
    }

    n = next[n];
  }

  return index;
}

/*
 *  g e t P h y s i c a l I n d e x
 */
int Indexlist::getPhysicalIndex(int givennumber) const {
  int i;
  int n = first;
  int index = -1; /* return -1 by default */

  /* Run trough indexlist until number is found, if so return it physicalindex.
   */
  for (i = 0; i < length; ++i) {
    if (number[n] == givennumber) {
      index = n;
      break;
    }

    n = next[n];
  }

  return index;
}

/*
 *  a d d N u m b e r
 */
returnValue Indexlist::addNumber(int addnumber) {
  int i;

  const int physicallength = number.size();
  if (lastusedindex + 1 < physicallength) {
    /* If there is enough storage, add number to indexlist. */
    ++lastusedindex;
    number[lastusedindex] = addnumber;
    next[lastusedindex] = 0;

    if (length == 0) {
      first = lastusedindex;
      previous[lastusedindex] = 0;
    } else {
      next[last] = lastusedindex;
      previous[lastusedindex] = last;
    }

    last = lastusedindex;
    ++length;

    return updateNumberArray();
  }
  /* Rearrangement of index list necessary! */
  if (length == physicallength) {
    return THROWERROR(RET_INDEXLIST_EXCEEDS_MAX_LENGTH);
  }

  /* copy existing elements */
  for (i = 0; i < length; ++i) {
    number[i] = numberArray[i];
    next[i] = i + 1;
    previous[i] = i - 1;
  }

  /* add new number at end of list */
  number[length] = addnumber;
  next[length] = -1;
  previous[length] = length - 1;

  /* and set remaining entries to empty */
  for (i = length + 1; i < physicallength; ++i) {
    number[i] = -1;
    next[i] = -1;
    previous[i] = -1;
  }

  first = 0;
  last = length;
  lastusedindex = length;
  ++length;

  const auto updateSuccess = updateNumberArray();
  if (updateSuccess != SUCCESSFUL_RETURN) {
    return updateSuccess;
  }
  return THROWERROR(RET_INDEXLIST_MUST_BE_REORDERD);
}

/*
 *  r e m o v e N u m b e r
 */
returnValue Indexlist::removeNumber(int removenumber) {
  int i = getPhysicalIndex(removenumber);

  /* nothing to be done if number is not contained in index set */
  if (i < 0) return SUCCESSFUL_RETURN;

  int p = previous[i];
  int n = next[i];

  if (i == last)
    last = p;
  else
    previous[n] = p;

  if (i == first)
    first = n;
  else
    next[p] = n;

  number[i] = -1;
  next[i] = -1;
  previous[i] = -1;
  --length;

  return updateNumberArray();
}

/*
 *  s w a p N u m b e r s
 */
returnValue Indexlist::swapNumbers(int number1, int number2) {
  int index1 = getPhysicalIndex(number1);
  int index2 = getPhysicalIndex(number2);

  /* consistency check */
  if ((index1 < 0) || (index2 < 0)) return THROWERROR(RET_INDEXLIST_CORRUPTED);

  int tmp = number[index1];
  number[index1] = number[index2];
  number[index2] = tmp;

  return updateNumberArray();
}

returnValue Indexlist::updateNumberArray() {
  int i;
  int n = first;

  /* Run trough indexlist and store numbers in numberArray. */
  for (i = 0; i < length; ++i) {
    if ((n >= 0) && (number[n] >= 0))
      numberArray[i] = number[n];
    else
      return THROWERROR(RET_INDEXLIST_CORRUPTED);

    n = next[n];
  }

  return SUCCESSFUL_RETURN;
}

}  // namespace qpoases_embedded

/*
 *  end of file
 */
