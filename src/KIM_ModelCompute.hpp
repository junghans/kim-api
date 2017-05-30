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


#ifndef KIM_MODEL_COMPUTE_HPP_
#define KIM_MODEL_COMPUTE_HPP_

namespace KIM
{
// Forward declarations
class LogLevel;
class ArgumentName;
class CallBackName;


class ModelCompute{
 public:
  // call back functions
  int get_neigh(int const neighborListIndex, int const particleNumber,
                int * const numberOfNeighbors,
                int const ** const neighborsOfParticle) const;

  int process_dEdr(double const de, double const r, double const * const dx,
                   int const i, int const j) const;

  int process_d2Edr2(double const de, double const * const r,
                     double const * const dx, int const * const i,
                     int const * const j) const;

  // data functions
  int get_data(ArgumentName const argumentName, int const ** const ptr) const;
  int get_data(ArgumentName const argumentName, int ** const ptr) const;
  int get_data(ArgumentName const argumentName, double const ** const ptr)
      const;
  int get_data(ArgumentName const argumentName, double ** const ptr) const;

  int is_call_back_present(CallBackName const callBackName, int * const present)
      const;

  void get_model_buffer(void ** const ptr) const;

  void Log(LogLevel const logLevel, std::string const & message,
           int const lineNumber, std::string const & fileName) const;
  std::string string() const;

 private:
  // do not allow copy constructor or operator=
  ModelCompute(ModelCompute const &);
  void operator=(ModelCompute const &);

  ModelCompute();
  ~ModelCompute();

  class ModelComputeImplementation;
  ModelComputeImplementation * pimpl;
};  // class ModelCompute
}  // namespace KIM
#endif  // KIM_MODEL_COMPUTE_HPP_
