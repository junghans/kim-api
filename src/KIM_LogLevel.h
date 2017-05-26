/*                                                                            */
/* CDDL HEADER START                                                          */
/*                                                                            */
/* The contents of this file are subject to the terms of the Common           */
/* Development and Distribution License Version 1.0 (the "License").          */
/*                                                                            */
/* You can obtain a copy of the license at                                    */
/* http://www.opensource.org/licenses/CDDL-1.0.  See the License for the      */
/* specific language governing permissions and limitations under the License. */
/*                                                                            */
/* When distributing Covered Code, include this CDDL HEADER in each file and  */
/* include the License file in a prominent location with the name             */
/* LICENSE.CDDL.                                                              */
/* If applicable, add the following below this CDDL HEADER, with the fields   */
/* enclosed by brackets "[]" replaced with your own identifying information:  */
/*                                                                            */
/* Portions Copyright (c) [yyyy] [name of copyright owner].                   */
/* All rights reserved.                                                       */
/*                                                                            */
/* CDDL HEADER END                                                            */
/*                                                                            */

/*                                                                            */
/* Copyright (c) 2016--2017, Regents of the University of Minnesota.          */
/* All rights reserved.                                                       */
/*                                                                            */
/* Contributors:                                                              */
/*    Ryan S. Elliott                                                         */
/*                                                                            */

/*                                                                            */
/* Release: This file is part of the kim-api.git repository.                  */
/*                                                                            */


#ifndef KIM_LOG_LEVEL_H_
#define KIM_LOG_LEVEL_H_

#include "KIM_LOG_DEFINES.inc"


struct KIM_LogLevel
{
  int logLevelID;
};
#ifndef KIM_LOG_LEVEL_DEFINED_
#define KIM_LOG_LEVEL_DEFINED_
typedef struct KIM_LogLevel KIM_LogLevel;
#endif

int KIM_LogLevelLessThan(KIM_LogLevel const left, KIM_LogLevel const right);
int KIM_LogLevelGreaterThan(KIM_LogLevel const left, KIM_LogLevel const right);
int KIM_LogLevelLessThanEqual(KIM_LogLevel const left,
                              KIM_LogLevel const right);
int KIM_LogLevelGreaterThanEqual(KIM_LogLevel const left,
                                 KIM_LogLevel const right);
int KIM_LogLevelEqual(KIM_LogLevel const left, KIM_LogLevel const right);
int KIM_LogLevelNotEqual(KIM_LogLevel const left, KIM_LogLevel const right);
char const * const KIM_LogLevelString(KIM_LogLevel const logLevel);

extern KIM_LogLevel const KIM_LOG_LEVEL_Silent;
extern KIM_LogLevel const KIM_LOG_LEVEL_Fatal;
extern KIM_LogLevel const KIM_LOG_LEVEL_Error;
extern KIM_LogLevel const KIM_LOG_LEVEL_Warning;
extern KIM_LogLevel const KIM_LOG_LEVEL_Information;
extern KIM_LogLevel const KIM_LOG_LEVEL_Debug;

#endif  /* KIM_LOG_LEVEL_H_ */
