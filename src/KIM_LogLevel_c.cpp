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


namespace
{
KIM::LogLevel const makeLogLevelCpp(KIM_LogLevel const logLevel)
{
  return KIM::LogLevel(logLevel.logLevelID);
}
}  // namespace

extern "C"
{
int KIM_LogLevelEqual(KIM_LogLevel const left, KIM_LogLevel const right)
{
  return (left.logLevelID == right.logLevelID);
}

int KIM_LogLevelNotEqual(KIM_LogLevel const left, KIM_LogLevel const right)
{
  return (!KIM_LogLevelEqual(left, right));
}

char const * const KIM_LogLevelString(KIM_LogLevel const logLevel)
{
  return (makeLogLevelCpp(logLevel)).string().c_str();
}

KIM_LogLevel const KIM_LOG_LEVEL_Silent = {KIM::LOG_LEVEL::Silent.logLevelID};
KIM_LogLevel const KIM_LOG_LEVEL_Fatal = {KIM::LOG_LEVEL::Fatal.logLevelID};
KIM_LogLevel const KIM_LOG_LEVEL_Error = {KIM::LOG_LEVEL::Error.logLevelID};
KIM_LogLevel const KIM_LOG_LEVEL_Warning = {KIM::LOG_LEVEL::Warning.logLevelID};
KIM_LogLevel const KIM_LOG_LEVEL_Information
= {KIM::LOG_LEVEL::Information.logLevelID};
KIM_LogLevel const KIM_LOG_LEVEL_Debug = {KIM::LOG_LEVEL::Debug.logLevelID};

}  // extern "C"
