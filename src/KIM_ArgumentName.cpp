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

#include <vector>

#ifndef KIM_DATA_TYPE_HPP_
#include "KIM_DataType.hpp"
#endif

#ifndef KIM_ARGUMENT_NAME_HPP_
#include "KIM_ArgumentName.hpp"
#endif

namespace KIM
{

ArgumentName::ArgumentName() : argumentNameID(0){}
ArgumentName::ArgumentName(int const id) : argumentNameID(id){}

bool ArgumentName::operator==(ArgumentName const & rhs) const
{return argumentNameID == rhs.argumentNameID;}
bool ArgumentName::operator!=(ArgumentName const & rhs) const
{return argumentNameID != rhs.argumentNameID;}

std::string ArgumentName::string() const
{
  if (*this == ARGUMENT_NAME::numberOfParticles) return "numberOfParticles";
  else if (*this == ARGUMENT_NAME::particleSpecies) return "particleSpecies";
  else if (*this == ARGUMENT_NAME::particleContributing)
    return "particleContributing";
  else if (*this == ARGUMENT_NAME::coordinates) return "coordinates";
  else if (*this == ARGUMENT_NAME::energy) return "energy";
  else if (*this == ARGUMENT_NAME::forces) return "forces";
  else if (*this == ARGUMENT_NAME::particleEnergy) return "particleEnergy";
  else if (*this == ARGUMENT_NAME::virial) return "virial";
  else if (*this == ARGUMENT_NAME::particleVirial) return "particleVirial";
  else if (*this == ARGUMENT_NAME::hessian) return "hessian";

  return "unknown";
}

// Order doesn't matter as long as all values are unique
namespace ARGUMENT_NAME
{
ArgumentName const numberOfParticles(0);
ArgumentName const particleSpecies(1);
ArgumentName const particleContributing(2);
ArgumentName const coordinates(3);
ArgumentName const energy(4);
ArgumentName const forces(5);
ArgumentName const particleEnergy(6);
ArgumentName const virial(7);
ArgumentName const particleVirial(8);
ArgumentName const hessian(9);

extern std::vector<ArgumentName> const mandatoryArguments = {
  numberOfParticles,
  particleSpecies,
  particleContributing,
  coordinates};

void get_number_of_arguments(int * const numberOfArguments)
{
  *numberOfArguments = 10;
}

int get_argument_name(int const index, ArgumentName * const argumentName)
{
  switch (index)
  {
    case 0:
      *argumentName = ARGUMENT_NAME::numberOfParticles;
      break;
    case 1:
      *argumentName = ARGUMENT_NAME::particleSpecies;
      break;
    case 2:
      *argumentName = ARGUMENT_NAME::particleContributing;
      break;
    case 3:
      *argumentName = ARGUMENT_NAME::coordinates;
      break;
    case 4:
      *argumentName = ARGUMENT_NAME::energy;
      break;
    case 5:
      *argumentName = ARGUMENT_NAME::forces;
      break;
    case 6:
      *argumentName = ARGUMENT_NAME::particleEnergy;
      break;
    case 7:
      *argumentName = ARGUMENT_NAME::virial;
      break;
    case 8:
      *argumentName = ARGUMENT_NAME::particleVirial;
      break;
    case 9:
      *argumentName = ARGUMENT_NAME::hessian;
      break;
    default:
      return true;  // invalid index
      break;
  }

  return false;  // no error
}

int get_argument_data_type(ArgumentName const argumentName,
                           DataType * const dataType)
{
  if (argumentName == ARGUMENT_NAME::numberOfParticles)
    *dataType = DATA_TYPE::Integer;
  if (argumentName == ARGUMENT_NAME::particleSpecies)
    *dataType = DATA_TYPE::Integer;
  if (argumentName == ARGUMENT_NAME::particleContributing)
    *dataType = DATA_TYPE::Integer;
  if (argumentName == ARGUMENT_NAME::coordinates)
    *dataType = DATA_TYPE::Double;
  if (argumentName == ARGUMENT_NAME::energy)
    *dataType = DATA_TYPE::Double;
  if (argumentName == ARGUMENT_NAME::forces)
    *dataType = DATA_TYPE::Double;
  if (argumentName == ARGUMENT_NAME::particleEnergy)
    *dataType = DATA_TYPE::Double;
  if (argumentName == ARGUMENT_NAME::virial)
    *dataType = DATA_TYPE::Double;
  if (argumentName == ARGUMENT_NAME::particleVirial)
    *dataType = DATA_TYPE::Double;
  if (argumentName == ARGUMENT_NAME::hessian)
    *dataType = DATA_TYPE::Double;
  else
    return true;  // unknown argument name

  return false;  // no error
}

}  // namespace ARGUMENT_NAME

}  // namespace KIM
