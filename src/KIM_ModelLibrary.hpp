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


#ifndef KIM_MODEL_LIBRARY_HPP_
#define KIM_MODEL_LIBRARY_HPP_

#include <string>

#ifndef KIM_FUNC_HPP_
#include "KIM_func.hpp"
#endif

namespace KIM
{
// Forward declarations
class LanguageName;


class ModelLibrary
{
 public:
  ModelLibrary();
  ~ModelLibrary();

  int open(std::string const & modelName);
  int getModelInitializationFunctionPointer(LanguageName * const languageName,
                                            func ** const functionPointer);
  int getModelCompiledWithVersion(std::string * const versionString);

 private:
  std::string modelName_;
  std::string libraryPath_;
  void * libraryHandle_;
};
}  // namespace KIM
#endif  // KIM_MODEL_LIBRARY_HPP_
