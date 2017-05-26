//
// CDDL HEADER START
//
// The contents of this file are subject to the terms of the Common
// Development and Distribution License Version 1.0 (the "License").
//
// You can obtain a copy of the license at
// http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
// specific language governing permissions and limitations under the License.
//
// When distributing Covered Code, include this CDDL HEADER in each file and
// include the License file in a prominent location with the name
// LICENSE.CDDL.
// If applicable, add the following below this CDDL HEADER, with the fields
// enclosed by brackets "[]" replaced with your own identifying information:
//
// Portions Copyright (c) [yyyy] [name of copyright owner].
// All rights reserved.
//
// CDDL HEADER END
//

//
// Copyright (c) 2016--2017, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Ryan S. Elliott
//

//
// Release: This file is part of the kim-api.git repository.
//

#ifndef KIM_LOG_LEVEL_HPP_
#include "KIM_LogLevel.hpp"
#endif
extern "C"
{
#ifndef KIM_LOG_LEVEL_H_
#include "KIM_LogLevel.h"
#endif
}  // extern "C"

#ifndef KIM_LOG_HPP_
#include "KIM_Log.hpp"
#endif
extern "C"
{
#ifndef KIM_LOG_H_
#include "KIM_Log.h"
#endif
}  // extern "C"


namespace
{
KIM::LogLevel const makeLogLevelCpp(KIM_LogLevel const logLevel)
{
  return KIM::LogLevel(logLevel.logLevelID);
}
}  // namespace

extern "C"
{
/* @@@@@ to be removed */
void KIM_report_error(int const line, char const * const file,
                      char const * const userMessage, int const statusCode)
{
  KIM_Log(KIM_LOG_LEVEL_Error, userMessage, line, file);
}

void KIM_Log(KIM_LogLevel const logLevel, char const * const message,
             int const lineNumber, char const * const fileName)
{
  KIM::Log(makeLogLevelCpp(logLevel), message, lineNumber, fileName);
}
}  // extern "C"
