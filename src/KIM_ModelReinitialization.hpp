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


#ifndef KIM_MODEL_REINITIALIZATION_HPP_
#define KIM_MODEL_REINITIALIZATION_HPP_

#include <string>

namespace KIM
{
// Forward declarations
class LengthUnit;
class EnergyUnit;
class ChargeUnit;
class TemperatureUnit;
class TimeUnit;


class ModelReinitialization{
 public:
  // stores pointer value
  void set_influence_distance(double const * const influenceDistance);

  // stores pointer value
  void set_cutoffs(int const numberOfCutoffs, double const * const cutoffs);

  void get_model_buffer(void ** const ptr) const;

  void Log(LogLevel const logLevel, std::string const & message,
           int const lineNumber, std::string const & fileName) const;
  std::string string() const;

 private:
  // do not allow copy constructor or operator=
  ModelReinitialization(ModelReinitialization const &);
  void operator=(ModelReinitialization const &);

  ModelReinitialization();
  ~ModelReinitialization();

  class ModelReinitializationImplementation;
  ModelReinitializationImplementation * pimpl;
};  // class ModelReinitialization
}  // namespace KIM
#endif  // KIM_MODEL_REINITIALIZATION_HPP_
