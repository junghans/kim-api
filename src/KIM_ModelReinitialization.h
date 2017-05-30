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


#ifndef KIM_MODEL_REINITIALIZATION_H_
#define KIM_MODEL_REINITIALIZATION_H_

struct KIM_ModelReinitialization {
  void * p;
};

#ifndef KIM_MODEL_REINITIALIZATION_DEFINED_
#define KIM_MODEL_REINITIALIZATION_DEFINED_
typedef struct KIM_ModelReinitialization KIM_ModelReinitialization;
#endif

void KIM_ModelReinitialization_set_influence_distance(
    KIM_ModelReinitialization * const modelReinitialization,
    double * const influenceDistance);

void KIM_ModelReinitialization_set_cutoffs(
    KIM_ModelReinitialization * const modelReinitialization,
    int const numberOfCutoffs, double const * const cutoffs);

void KIM_ModelReinitialization_Log(
    KIM_ModelReinitialization const * const modelReinitialization,
    KIM_LogLevel const logLevel, char const * const message,
    int const lineNumber, char const * const fileName);

char const * const KIM_ModelReinitialization_string(
    KIM_ModelReinitialization const * const modelReinitialization);

void KIM_ModelReinitialization_get_model_buffer(
    KIM_ModelReinitialization const * const modelReinitialization,
    void ** const ptr);

#endif  /* KIM_MODEL_REINITIALIZATION_H_ */
