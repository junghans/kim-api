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

#ifndef KIM_NUMBERING_HPP_
#include "KIM_Numbering.hpp"
#endif

namespace KIM
{

Numbering::Numbering() : numberingID(0) {}
Numbering::Numbering(int const id) : numberingID(id) {}
bool Numbering::operator==(Numbering const & rhs) const
{return numberingID==rhs.numberingID;}
bool Numbering::operator!=(Numbering const & rhs) const
{return numberingID!=rhs.numberingID;}

std::string Numbering::string() const
{
  if (*this == NUMBERING::zeroBased) return "zeroBased";
  else if (*this == NUMBERING::oneBased) return "oneBased";
  else return "unknown";
}

namespace NUMBERING
{
Numbering const zeroBased(0);
Numbering const oneBased(1);
}  // namespace NUMBERING
}  // namespace KIM
