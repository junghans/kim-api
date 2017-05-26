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

#include <iostream>
#include <iomanip>
#include <vector>

#include <dlfcn.h>

#include "old_KIM_API_DIRS.h"

#ifndef KIM_MODEL_LIBRARY_HPP_
#include "KIM_ModelLibrary.hpp"
#endif

#ifndef KIM_LANGUAGE_NAME_HPP_
#include "KIM_LanguageName.hpp"
#endif

namespace KIM
{

ModelLibrary::ModelLibrary() :
    libraryHandle_(0)
{}

ModelLibrary::~ModelLibrary()
{
  if (libraryHandle_ != 0)
  {
    int error = dlclose(libraryHandle_);
    if (error)
    {
      // @@@@ report error in log....
    }
  }
}

int ModelLibrary::open(std::string const & modelName)
{
  if (libraryHandle_ != 0) return true;  // already open


  modelName_ = modelName;

  std::vector<std::string> item;
  bool accessible = findItem(OLD_KIM::KIM_MODELS_DIR, modelName_, &item);
  if (!accessible) return true;  // cannot find modelName
  libraryPath_ = item[1] + "/" + item[0] + "/" + MODELLIBFILE + ".so";

  libraryHandle_ = dlopen(libraryPath_.c_str(), RTLD_NOW);
  if (libraryHandle_ == 0)
  {
    std::cout << dlerror() << std::endl;
    return true;
  }

  return false;
}

int ModelLibrary::getModelInitializationFunctionPointer(
    LanguageName * const languageName, func ** const functionPointer)
{
  if (libraryHandle_ == 0) return true;  // not open

  std::string languageSymbol(modelName_ + "_language");
  LanguageName * pLanguageName
      = reinterpret_cast<LanguageName *>(dlsym(libraryHandle_,
                                               languageSymbol.c_str()));
  if (pLanguageName == 0)
  {
    std::cout << dlerror() << std::endl;
    return true;
  }
  else
  {
    *languageName = *pLanguageName;
  }

  std::string initFunctionSymbol(modelName_ + "_init_pointer");
  func ** pointerToFunctionPointer
      = reinterpret_cast<func **>(dlsym(libraryHandle_,
                                        initFunctionSymbol.c_str()));

  if (pointerToFunctionPointer == 0)
  {
    std::cout << dlerror() << std::endl;
    return true;
  }

  *functionPointer = *(pointerToFunctionPointer);
  return false;
}

int ModelLibrary::getModelCompiledWithVersion(
    std::string * const versionString)
{
  if (libraryHandle_ == 0) return true;  // not open

  std::string versionSymbol(modelName_ + "_compiled_with_version");
  char const * versionCharString
      = static_cast<char const *>(dlsym(libraryHandle_, versionSymbol.c_str()));
  if (versionCharString == 0)
  {
    std::cout << dlerror() << std::endl;
    return true;
  }

  *versionString = versionCharString;
  return false;
}

}  // namespace KIM
