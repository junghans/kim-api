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


#ifndef KIM_CHARGE_UNIT_H_
#define KIM_CHARGE_UNIT_H_

struct KIM_ChargeUnit
{
  int chargeUnitID;
};
#ifndef KIM_CHARGE_UNIT_DEFINED_
#define KIM_CHARGE_UNIT_DEFINED_
typedef struct KIM_ChargeUnit KIM_ChargeUnit;
#endif

int KIM_ChargeUnitEqual(KIM_ChargeUnit const left, KIM_ChargeUnit right);
int KIM_ChargeUnitNotEqual(KIM_ChargeUnit const left, KIM_ChargeUnit right);
char const * const KIM_ChargeUnitString(KIM_ChargeUnit const chargeUnit);

extern KIM_ChargeUnit const KIM_CHARGE_UNIT_any;
extern KIM_ChargeUnit const KIM_CHARGE_UNIT_C;
extern KIM_ChargeUnit const KIM_CHARGE_UNIT_e;
extern KIM_ChargeUnit const KIM_CHARGE_UNIT_statC;

#endif  /* KIM_CHARGE_UNIT_H_ */
