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

#ifndef KIM_ENERGY_UNIT_HPP_
#include "KIM_EnergyUnit.hpp"
#endif

namespace KIM
{

EnergyUnit::EnergyUnit() : energyUnitID(0){}
EnergyUnit::EnergyUnit(int const id) : energyUnitID(id){}
bool EnergyUnit::operator==(EnergyUnit const & rhs) const
{return energyUnitID==rhs.energyUnitID;}
bool EnergyUnit::operator!=(EnergyUnit const & rhs) const
{return energyUnitID!=rhs.energyUnitID;}

std::string EnergyUnit::string() const
{
  if (*this == ENERGY_UNIT::any) return "any";
  else if (*this == ENERGY_UNIT::amu_A2_per_ps2) return "amu_A2_per_ps2";
  else if (*this == ENERGY_UNIT::erg) return "erg";
  else if (*this == ENERGY_UNIT::eV) return "eV";
  else if (*this == ENERGY_UNIT::Hartree) return "Hartree";
  else if (*this == ENERGY_UNIT::J) return "J";
  else if (*this == ENERGY_UNIT::kcal_mol) return "kcal_mol";
  else return "unknown";
}

namespace ENERGY_UNIT
{
EnergyUnit const any(0);
EnergyUnit const amu_A2_per_ps2(1);
EnergyUnit const erg(2);
EnergyUnit const eV(3);
EnergyUnit const Hartree(4);
EnergyUnit const J(5);
EnergyUnit const kcal_mol(6);
}  // namespace ENERGY_UNIT

}  // namespace KIM
