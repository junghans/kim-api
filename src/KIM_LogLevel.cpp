//
// CDDL HEADER START
//
// The contents of this file are subject to the terms of the Common Development
// and Distribution License Version 1.0 (the "License").
//
// You can obtain a copy of the license at
// http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
// specific language governing permissions and limitations under the License.
//
// When distributing Covered Code, include this CDDL HEADER in each file and
// include the License file in a prominent location with the name LICENSE.CDDL.
// If applicable, add the following below this CDDL HEADER, with the fields
// enclosed by brackets "[]" replaced with your own identifying information:
//
// Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
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

namespace KIM
{

LogLevel::LogLevel() : logLevelID(0){}
LogLevel::LogLevel(int const id) : logLevelID(id){}

bool LogLevel::operator<(LogLevel const & rhs) const
{return logLevelID < rhs.logLevelID;}
bool LogLevel::operator>(LogLevel const & rhs) const
{return logLevelID > rhs.logLevelID;}
bool LogLevel::operator<=(LogLevel const & rhs) const
{return logLevelID <= rhs.logLevelID;}
bool LogLevel::operator>=(LogLevel const & rhs) const
{return logLevelID >= rhs.logLevelID;}
bool LogLevel::operator==(LogLevel const & rhs) const
{return logLevelID == rhs.logLevelID;}
bool LogLevel::operator!=(LogLevel const & rhs) const
{return logLevelID != rhs.logLevelID;}

std::string LogLevel::string() const
{
  if (*this == LOG_LEVEL::Silent)
    return "Silent";
  else if (*this == LOG_LEVEL::Fatal)
    return "Fatal";
  else if (*this == LOG_LEVEL::Error)
    return "Error";
  else if (*this == LOG_LEVEL::Warning)
    return "Warning";
  else if (*this == LOG_LEVEL::Information)
    return "Information";
  else if (*this == LOG_LEVEL::Debug)
    return "Debug";
  else
    return "unknown";
}

// Order is important
namespace LOG_LEVEL
{
LogLevel const Silent(0);
LogLevel const Fatal(1);
LogLevel const Error(2);
LogLevel const Warning(3);
LogLevel const Information(4);
LogLevel const Debug(5);
}  // namespace LOG_LEVEL

}  // namespace KIM
