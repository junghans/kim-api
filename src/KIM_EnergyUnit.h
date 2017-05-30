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


#ifndef KIM_ENERGY_UNIT_H_
#define KIM_ENERGY_UNIT_H_

struct KIM_EnergyUnit
{
  int energyUnitID;
};
#ifndef KIM_ENERGY_UNIT_DEFINED_
#define KIM_ENERGY_UNIT_DEFINED_
typedef struct KIM_EnergyUnit KIM_EnergyUnit;
#endif

int KIM_EnergyUnitEqual(KIM_EnergyUnit const left, KIM_EnergyUnit const right);
int KIM_EnergyUnitNotEqual(KIM_EnergyUnit const left,
                           KIM_EnergyUnit const right);
char const* const KIM_EnergyUnitString(KIM_EnergyUnit const energyUnit);

extern KIM_EnergyUnit const KIM_ENERGY_UNIT_any;
extern KIM_EnergyUnit const KIM_ENERGY_UNIT_amu_A2_per_ps2;
extern KIM_EnergyUnit const KIM_ENERGY_UNIT_erg;
extern KIM_EnergyUnit const KIM_ENERGY_UNIT_eV;
extern KIM_EnergyUnit const KIM_ENERGY_UNIT_Hartree;
extern KIM_EnergyUnit const KIM_ENERGY_UNIT_J;
extern KIM_EnergyUnit const KIM_ENERGY_UNIT_kcal_mol;

#endif  /* KIM_ENERGY_UNIT_H_ */
