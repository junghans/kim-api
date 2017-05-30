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


#ifndef KIM_MODEL_COMPUTE_H_
#define KIM_MODEL_COMPUTE_H_

/* Forward declarations */
#ifndef KIM_LOG_LEVEL_DEFINED_
#define KIM_LOG_LEVEL_DEFINED_
typedef struct KIM_LogLevel KIM_LogLevel;
#endif

#ifndef KIM_ARGUMENT_NAME_DEFINED_
#define KIM_ARGUMENT_NAME_DEFINED_
typedef struct KIM_ArgumentName KIM_ArgumentName;
#endif

#ifndef KIM_CALL_BACK_NAME_DEFINED_
#define KIM_CALL_BACK_NAME_DEFINED_
typedef struct KIM_CallBackName KIM_CallBackName;
#endif


struct KIM_ModelCompute {
  void * p;
};

#ifndef KIM_MODEL_COMPUTE_DEFINED_
#define KIM_MODEL_COMPUTE_DEFINED_
typedef struct KIM_ModelCompute KIM_ModelCompute;
#endif

int KIM_ModelCompute_get_neigh(KIM_ModelCompute const * const modelCompute,
                               int const neighborListIndex,
                               int const particleNumber,
                               int * const numberOfNeighbors,
                               int const ** const neighborsOfParticle);

int KIM_ModelCompute_process_dEdr(KIM_ModelCompute const * const modelCompute,
                                  double const de, double const r,
                                  double const * const dx, int const i,
                                  int const j);

int KIM_ModelCompute_process_d2Edr2(KIM_ModelCompute const * const modelCompute,
                                    double const de, double const * const r,
                                    double const * const dx,
                                    int const * const i,
                                    int const * const j);

/* *data functions */
int KIM_ModelCompute_get_data_int(KIM_ModelCompute const * const modelCompute,
                                  KIM_ArgumentName const argumentName,
                                  int ** const ptr);

int KIM_ModelCompute_get_data_double(
    KIM_ModelCompute const * const modelCompute,
    KIM_ArgumentName const argumentName, double ** const ptr);

int KIM_ModelCompute_is_call_back_present(
    KIM_ModelCompute const * const modelCompute,
    KIM_CallBackName const callBackName, int * const present);

void KIM_ModelCompute_get_model_buffer(
    KIM_ModelCompute const * const modelCompute, void ** const ptr);

void KIM_ModelCompute_Log(
    KIM_ModelCompute const * const modelCompute,
    KIM_LogLevel const logLevel, char const * const message,
    int const lineNumber, char const * const fileName);

char const * const KIM_ModelCompute_string(
    KIM_ModelCompute const * const modelCompute);

#endif  /* KIM_MODEL_COMPUTE_H_ */
