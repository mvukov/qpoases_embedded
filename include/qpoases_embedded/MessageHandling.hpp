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

#ifndef QPOASES_MESSAGEHANDLING_HPP
#define QPOASES_MESSAGEHANDLING_HPP

#include <string>

#include <qpoases_embedded/Types.hpp>

namespace qpoases_embedded {

const std::string& getErrorString(returnValue error);

#ifdef QPOASES_DEBUG

namespace internal {

enum class MsgType { kInfo, kError };

returnValue throwMessage(MsgType msgType, returnValue retval,
                         const char* filename, const unsigned long linenumber,
                         const char* functionname);

returnValue throwMessage(MsgType msgType, returnValue retval,
                         const char* filename, const unsigned long linenumber,
                         const char* functionname, const char* fmt, ...);

}  // namespace internal

#ifndef __FUNCTION__
/** Ensures that __FUNCTION__ macro is defined. */
#elif defined(__func__)
#define __FUNCTION__ __func__
#else
#define __FUNCTION__ 0
#endif

#ifndef __FILE__
/** Ensures that __FILE__ macro is defined. */
#define __FILENAME__ 0
#else
#define __FILENAME__ strrchr("/" __FILE__, '/') + 1
#endif

#ifndef __LINE__
/** Ensures that __LINE__ macro is defined. */
#define __LINE__ 0
#endif

#define THROWINFO(retval)                                                \
  internal::throwMessage(internal::MsgType::kInfo, retval, __FILENAME__, \
                         __LINE__, __FUNCTION__)

#define THROWINFOMSG(retval, ...)                                        \
  internal::throwMessage(internal::MsgType::kInfo, retval, __FILENAME__, \
                         __LINE__, __FUNCTION__, __VA_ARGS__)

#define THROWERROR(retval)                                                \
  internal::throwMessage(internal::MsgType::kError, retval, __FILENAME__, \
                         __LINE__, __FUNCTION__)

#define THROWERRORMSG(retval, ...)                                        \
  internal::throwMessage(internal::MsgType::kError, retval, __FILENAME__, \
                         __LINE__, __FUNCTION__, __VA_ARGS__)

#else /* QPOASES_DEBUG */

namespace internal {
inline returnValue throwMessage(returnValue retval) { return retval; }
}  // namespace internal

#define THROWINFO(retval) internal::throwMessage(retval);

#define THROWINFOMSG(retval, ...)

#define THROWERROR(retval) internal::throwMessage(retval);

#define THROWERRORMSG(retval, ...)

#endif /* QPOASES_DEBUG */

}  // namespace qpoases_embedded

#endif /* QPOASES_MESSAGEHANDLING_HPP */

/*
 *  end of file
 */
