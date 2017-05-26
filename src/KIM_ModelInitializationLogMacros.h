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


#ifdef LOG_SILENT
#undef LOG_SILENT
#endif
#if (KIM_LOG_MAXIMUM_LEVEL >= KIM_LOG_LEVEL_SILENT_)
#define LOG_SILENT(message)                                             \
  KIM_ModelInitialization_Log(modelInitialization,                      \
                              KIM_LOG_LEVEL_Silent, message,            \
                              __LINE__, __FILE__)
#else
#define LOG_SILENT(message)
#endif

#ifdef LOG_FATAL
#undef LOG_FATAL
#endif
#if (KIM_LOG_MAXIMUM_LEVEL >= KIM_LOG_LEVEL_FATAL_)
#define LOG_FATAL(message)                                            \
  KIM_ModelInitialization_Log(modelInitialization,                    \
                              KIM_LOG_LEVEL_Fatal, message,           \
                              __LINE__, __FILE__)
#else
#define LOG_FATAL(message)
#endif

#ifdef LOG_ERROR
#undef LOG_ERROR
#endif
#if (KIM_LOG_MAXIMUM_LEVEL >= KIM_LOG_LEVEL_ERROR_)
#define LOG_ERROR(message)                                           \
  KIM_ModelInitialization_Log(modelInitialization,                   \
                              KIM_LOG_LEVEL_Error, message,          \
                              __LINE__, __FILE__)
#else
#define LOG_ERROR(message)
#endif

#ifdef LOG_WARNING
#undef LOG_WARNING
#endif
#if (KIM_LOG_MAXIMUM_LEVEL >= KIM_LOG_LEVEL_WARNING_)
#define LOG_WARNING(message)                                          \
  KIM_ModelInitialization_Log(modelInitialization,                    \
                              KIM_LOG_LEVEL_Warning, message,         \
                              __LINE__, __FILE__)
#else
#define LOG_WARNING(message)
#endif

#ifdef LOG_INFORMATION
#undef LOG_INFORMATION
#endif
#if (KIM_LOG_MAXIMUM_LEVEL >= KIM_LOG_LEVEL_INFORMATION_)
#define LOG_INFORMATION(message)                                        \
  KIM_ModelInitialization_Log(modelInitialization,                      \
                              KIM_LOG_LEVEL_Information, message,       \
                              __LINE__, __FILE__)
#else
#define LOG_INFORMATION(message)
#endif

#ifdef LOG_DEBUG
#undef LOG_DEBUG
#endif
#if (KIM_LOG_MAXIMUM_LEVEL >= KIM_LOG_LEVEL_DEBUG_)
#define LOG_DEBUG(message)                                      \
  KIM_ModelInitialization_Log(modelInitialization,              \
                              KIM_LOG_LEVEL_Debug, message,     \
                              __LINE__, __FILE__)
#else
#define LOG_DEBUG(message)
#endif
