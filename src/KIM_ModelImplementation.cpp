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

#ifndef KIM_MODEL_IMPLEMENTATION_HPP_
#include "KIM_ModelImplementation.hpp"
#endif

#ifndef KIM_MODEL_LIBRARY_HPP_
#include "KIM_ModelLibrary.hpp"
#endif

#ifndef KIM_LANGUAGE_NAME_HPP_
#include "KIM_LanguageName.hpp"
#endif

#ifndef KIM_UNIT_SYSTEM_H_
extern "C"
{
#include "KIM_UnitSystem.h"
}  // extern "C"
#endif

#ifndef KIM_MODEL_INITIALIZATION_H_
extern "C"
{
#include "KIM_ModelInitialization.h"
}  // extern "C"
#endif

#ifndef KIM_MODEL_REINITIALIZATION_H_
extern "C"
{
#include "KIM_ModelReinitialization.h"
}  // extern "C"
#endif

#ifndef KIM_MODEL_COMPUTE_H_
extern "C"
{
#include "KIM_ModelCompute.h"
}  // extern "C"
#endif

#ifndef KIM_MODEL_DESTROY_H_
extern "C"
{
#include "KIM_ModelDestroy.h"
}  // extern "C"
#endif


namespace KIM
{
namespace ARGUMENT_NAME
{
extern std::vector<ArgumentName> const mandatoryArguments;
}  // namespace ARGUMENT_NAME

namespace CALL_BACK_NAME
{
extern std::vector<CallBackName> const mandatoryCallBacks;
}  // namespace CALL_BACK_NAME
}  // namespace KIM

namespace std
{
size_t hash<KIM::SpeciesName>::operator()(KIM::SpeciesName const & speciesName)
    const
{
  return speciesName.speciesNameID;
}

size_t hash<KIM::ArgumentName>::operator()(
    KIM::ArgumentName const & argumentName) const
{
  return argumentName.argumentNameID;
}


size_t hash<KIM::CallBackName>::operator()(
    KIM::CallBackName const & callBackName) const
{
  return callBackName.callBackNameID;
}
}  // namespace std


namespace
{
KIM_LengthUnit makeLengthUnitC(KIM::LengthUnit const lengthUnit)
{
  KIM_LengthUnit lengthUnitC = {lengthUnit.lengthUnitID};
  return lengthUnitC;
}

KIM_EnergyUnit makeEnergyUnitC(KIM::EnergyUnit const energyUnit)
{
  KIM_EnergyUnit energyUnitC = {energyUnit.energyUnitID};
  return energyUnitC;
}

KIM_ChargeUnit makeChargeUnitC(KIM::ChargeUnit const chargeUnit)
{
  KIM_ChargeUnit chargeUnitC = {chargeUnit.chargeUnitID};
  return chargeUnitC;
}

KIM_TemperatureUnit makeTemperatureUnitC(
    KIM::TemperatureUnit const temperatureUnit)
{
  KIM_TemperatureUnit temperatureUnitC = {temperatureUnit.temperatureUnitID};
  return temperatureUnitC;
}

KIM_TimeUnit makeTimeUnitC(KIM::TimeUnit const timeUnit)
{
  KIM_TimeUnit timeUnitC = {timeUnit.timeUnitID};
  return timeUnitC;
}
}  // namespace

namespace KIM
{
// Forward declarations
class ModelInitialization;
class ModelReinitialization;
class ModelCompute;
class ModelDestroy;

int ModelImplementation::create(
    Numbering const numbering,
    LengthUnit const requestedLengthUnit,
    EnergyUnit const requestedEnergyUnit,
    ChargeUnit const requestedChargeUnit,
    TemperatureUnit const requestedTemperatureUnit,
    TimeUnit const requestedTimeUnit,
    std::string const & modelName,
    int * const requestedUnitsAccepted,
    ModelImplementation ** const modelImplementation)
{
  ModelImplementation * pModelImplementation;
  pModelImplementation = new ModelImplementation(new ModelLibrary());

  int error = pModelImplementation->ModelInitialization(
      numbering, requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
      requestedTemperatureUnit, requestedTimeUnit, modelName);
  if (error)
  {
    delete pModelImplementation;
    return true;
  }

  LengthUnit finalLengthUnit;
  EnergyUnit finalEnergyUnit;
  ChargeUnit finalChargeUnit;
  TemperatureUnit finalTemperatureUnit;
  TimeUnit finalTimeUnit;
  pModelImplementation->get_units(&finalLengthUnit, &finalEnergyUnit,
                                  &finalChargeUnit, &finalTemperatureUnit,
                                  &finalTimeUnit);

  if (((finalLengthUnit == LENGTH_UNIT::any) ||
       (finalLengthUnit == requestedLengthUnit))
      &&
      ((finalEnergyUnit == ENERGY_UNIT::any) ||
       (finalEnergyUnit == requestedEnergyUnit))
      &&
      ((finalChargeUnit == CHARGE_UNIT::any) ||
       (finalChargeUnit == requestedChargeUnit))
      &&
      ((finalTemperatureUnit == TEMPERATURE_UNIT::any) ||
       (finalTemperatureUnit == requestedTemperatureUnit))
      &&
      ((finalTimeUnit == TIME_UNIT::any) ||
       (finalTimeUnit == requestedTimeUnit)))
  {
    *requestedUnitsAccepted = true;
  }
  else
  {
    *requestedUnitsAccepted = false;
  }

  *modelImplementation = pModelImplementation;

  return false;
}


void ModelImplementation::destroy(
    ModelImplementation ** const modelImplementation)
{

  int error = (*modelImplementation)->ModelDestroy();
  if (error)
    ;// @@@@@ log message

  delete *modelImplementation;
  *modelImplementation = 0;
}

void ModelImplementation::set_influence_distance(
    double const * const influenceDistance)
{
  influenceDistance_ = influenceDistance;
}
void ModelImplementation::get_influence_distance(
    double * const influenceDistance) const
{
  *influenceDistance = *influenceDistance_;
}


void ModelImplementation::set_cutoffs(
    int const numberOfCutoffs, double const * const cutoffs)
{
  numberOfCutoffs_ = numberOfCutoffs;
  cutoffs_ = cutoffs;
}

// allows NULL as value of cutoffs (to get just numberOfCutoffs)
void ModelImplementation::get_cutoffs(int * const numberOfCutoffs,
                                      double const ** const cutoffs) const
{
  *numberOfCutoffs = numberOfCutoffs_;
  *cutoffs = cutoffs_;
}


int ModelImplementation::set_reinit(LanguageName const languageName,
                                    func * const fptr)
{
  reinitializationLanguage_ = languageName;
  reinitializationFunction_ = fptr;

  return false;
}

int ModelImplementation::set_destroy(LanguageName const languageName,
                                     func * const fptr)
{
  destroyLanguage_ = languageName;
  destroyFunction_ = fptr;

  return false;
}

int ModelImplementation::set_compute_func(LanguageName const languageName,
                                          func * const fptr)
{
  computeLanguage_ = languageName;
  computeFunction_ = fptr;

  return false;
}


int ModelImplementation::set_species_code(SpeciesName const speciesName,
                                          int const code)
{
  supportedSpecies_[speciesName] = code;

  return false;
}

int ModelImplementation::get_species_support_and_code(
    KIM::SpeciesName const speciesName,
    int * const speciesIsSupported,
    int * const code) const
{
  auto result = supportedSpecies_.find(speciesName);

  if (result == supportedSpecies_.end())
  {
    *speciesIsSupported = false;
  }
  else
  {
    *speciesIsSupported = true;
    *code = result->second;
  }

  return false;
}


int ModelImplementation::set_argument_attribute(ArgumentName const argumentName,
                                                Attribute const attribute)
{
  argumentAttribute_[argumentName] = attribute;

  return false;
}

int ModelImplementation::get_argument_attribute(ArgumentName const argumentName,
                                                Attribute * const attribute)
    const
{
  auto result = argumentAttribute_.find(argumentName);

  if (result == argumentAttribute_.end())
  {
    return true;
  }
  else
  {
    *attribute = result->second;
    return false;
  }
}


int ModelImplementation::set_call_back_attribute(
    CallBackName const callBackName,
    Attribute const attribute)
{
  callBackAttribute_[callBackName] = attribute;

  return false;
}

int ModelImplementation::get_call_back_attribute(
    CallBackName const callBackName,
    Attribute * const attribute) const
{
  auto result = callBackAttribute_.find(callBackName);

  if (result == callBackAttribute_.end())
  {
    return true;
  }
  else
  {
    *attribute = result->second;
    return false;
  }
}


int ModelImplementation::set_model_numbering(Numbering const numbering)
{
  modelNumbering_ = numbering;

  return false;
}

int ModelImplementation::set_simulator_numbering(Numbering const numbering)
{
  simulatorNumbering_ = numbering;
  return false;
}


int ModelImplementation::set_units(LengthUnit const lengthUnit,
                                   EnergyUnit const energyUnit,
                                   ChargeUnit const chargeUnit,
                                   TemperatureUnit const temperatureUnit,
                                   TimeUnit const timeUnit)
{
  lengthUnit_ = lengthUnit;
  energyUnit_ = energyUnit;
  chargeUnit_ = chargeUnit;
  temperatureUnit_ = temperatureUnit;
  timeUnit_ = timeUnit;

  return false;
}

void ModelImplementation::get_units(LengthUnit * const lengthUnit,
                                    EnergyUnit * const energyUnit,
                                    ChargeUnit * const chargeUnit,
                                    TemperatureUnit * const temperatureUnit,
                                    TimeUnit * const timeUnit) const
{
  *lengthUnit = lengthUnit_;
  *energyUnit = energyUnit_;
  *chargeUnit = chargeUnit_;
  *temperatureUnit = temperatureUnit_;
  *timeUnit = timeUnit_;
}


int ModelImplementation::set_parameter(int const extent, int * const ptr,
                                       std::string const & description)
{
  parameterDescription_.push_back(description);
  parameterDataType_.push_back(DATA_TYPE::Integer);
  parameterExtent_.push_back(extent);
  parameterPointer_.push_back(ptr);

  return false;
}

int ModelImplementation::set_parameter(int const extent, double * const ptr,
                                       std::string const & description)
{
  parameterDescription_.push_back(description);
  parameterDataType_.push_back(DATA_TYPE::Double);
  parameterExtent_.push_back(extent);
  parameterPointer_.push_back(ptr);

  return false;
}

void ModelImplementation::get_num_params(int * const numberOfParameters) const
{
  *numberOfParameters = parameterPointer_.size();
}

int ModelImplementation::get_parameter_data_type_and_description(
    int const index, DataType * const dataType,
    std::string * const description) const
{
  *dataType = parameterDataType_[index];
  *description = parameterDescription_[index];

  return false;
}

int ModelImplementation::get_parameter_extent_and_pointer(
    int const index, int * extent, int ** const ptr)
{
  *extent = parameterExtent_[index];
  *ptr = reinterpret_cast<int *>(parameterPointer_[index]);

  return false;
}

int ModelImplementation::get_parameter_extent_and_pointer(
    int const index, int * extent, int const ** const ptr) const
{
  *extent = parameterExtent_[index];
  *ptr = reinterpret_cast<int const *>(parameterPointer_[index]);

  return false;
}

int ModelImplementation::get_parameter_extent_and_pointer(
    int const index, int * extent, double ** const ptr)
{
  *extent = parameterExtent_[index];
  *ptr = reinterpret_cast<double *>(parameterPointer_[index]);

  return false;
}

int ModelImplementation::get_parameter_extent_and_pointer(
    int const index, int * extent, double const ** const ptr) const
{
  *extent = parameterExtent_[index];
  *ptr = reinterpret_cast<double const *>(parameterPointer_[index]);

  return false;
}


int ModelImplementation::set_data(ArgumentName const argumentName,
                                  int const * const ptr)
{
  argumentPointer_[argumentName]
      = reinterpret_cast<void *>(const_cast<int *>(ptr));

  return false;
}

int ModelImplementation::set_data(ArgumentName const argumentName,
                                  double const * const ptr)
{
  argumentPointer_[argumentName]
      = reinterpret_cast<void *>(const_cast<double *>(ptr));

  return false;
}

int ModelImplementation::get_data(ArgumentName const argumentName,
                                  int const ** const ptr) const
{
  auto result = argumentPointer_.find(argumentName);

  if (result == argumentPointer_.end())
  {
    *ptr = 0;
    return false;
  }
  else
  {
    *ptr = reinterpret_cast<int const *>(result->second);
    return false;
  }
}

int ModelImplementation::get_data(ArgumentName const argumentName,
                                  int ** const ptr) const
{
  auto result = argumentPointer_.find(argumentName);

  if (result == argumentPointer_.end())
  {
    *ptr = 0;
    return false;
  }
  else
  {
    *ptr = reinterpret_cast<int *>(result->second);
    return false;
  }
}

int ModelImplementation::get_data(ArgumentName const argumentName,
                                  double const ** const ptr) const
{
  auto result = argumentPointer_.find(argumentName);

  if (result == argumentPointer_.end())
  {
    *ptr = 0;
    return false;
  }
  else
  {
    *ptr = reinterpret_cast<double const *>(result->second);
    return false;
  }
}

int ModelImplementation::get_data(ArgumentName const argumentName,
                                  double ** const ptr) const
{
  auto result = argumentPointer_.find(argumentName);

  if (result == argumentPointer_.end())
  {
    *ptr = 0;
    return false;
  }
  else
  {
    *ptr = reinterpret_cast<double *>(result->second);
    return false;
  }
}



int ModelImplementation::set_call_back(CallBackName const callBackName,
                                       LanguageName const languageName,
                                       func * const fptr,
                                       void const * const dataObject)
{
  auto result = callBackAttribute_.find(callBackName);

  if ((result == callBackAttribute_.end())
      ||
      (result->second == ATTRIBUTE::notSupported))
  {
    return true;
  }
  else
  {
    callBackLanguage_[callBackName] = languageName;
    callBackFunctionPointer_[callBackName] = fptr;
    callBackDataObjectPointer_[callBackName] = dataObject;
    return false;
  }
}

int ModelImplementation::is_call_back_present(
    CallBackName const callBackName, int * const present) const
{
  auto result = callBackFunctionPointer_.find(callBackName);

  if ((result == callBackFunctionPointer_.end())
      ||
      (result->second == 0))
  {
    *present = false;
    return false;
  }
  else
  {
    *present = true;
    return false;
  }
}


int ModelImplementation::compute() const
{
  typedef int ModelComputeCpp(KIM::ModelCompute * const);
  ModelComputeCpp * CppCompute
      = reinterpret_cast<ModelComputeCpp *>(computeFunction_);
  typedef int ModelComputeC(KIM_ModelCompute * const);
  ModelComputeC * CCompute
      = reinterpret_cast<ModelComputeC *>(computeFunction_);
  typedef void ModelComputeF(KIM_ModelCompute * const, int * const);
  ModelComputeF * FCompute
      = reinterpret_cast<ModelComputeF *>(computeFunction_);

  int error;
  struct Mdl {void const * p;};
  Mdl M;
  M.p = this;
  if (computeLanguage_ == LANGUAGE_NAME::Cpp)
  {
    error = CppCompute(reinterpret_cast<KIM::ModelCompute *>(&M));
  }
  else if (computeLanguage_ == LANGUAGE_NAME::C)
  {
    KIM_ModelCompute cM;
    cM.p = &M;
    error = CCompute(&cM);
  }
  else if (computeLanguage_ == LANGUAGE_NAME::Fortran)
  {
    KIM_ModelCompute cM;
    cM.p = &M;
    FCompute(&cM, &error);
  }
  else
  {
    return true;
  }

  if (error)
    return true;
  else
    return false;
}

int ModelImplementation::clear_pointers_and_reinitialize_model()
{
  influenceDistance_ = 0;
  numberOfCutoffs_ = 0;
  cutoffs_ = 0;
  argumentPointer_.clear();
  callBackDataObjectPointer_.clear();
  callBackFunctionPointer_.clear();

  typedef int ModelReinitializationCpp(KIM::ModelReinitialization * const);
  ModelReinitializationCpp * CppReinitialization
      = reinterpret_cast<ModelReinitializationCpp *>(reinitializationFunction_);
  typedef int ModelReinitializationC(KIM_ModelReinitialization * const);
  ModelReinitializationC * CReinitialization
      = reinterpret_cast<ModelReinitializationC *>(reinitializationFunction_);
  typedef void ModelReinitializationF(KIM_ModelReinitialization * const,
                                      int * const);
  ModelReinitializationF * FReinitialization
      = reinterpret_cast<ModelReinitializationF *>(reinitializationFunction_);

  int error;
  struct Mdl {void * p;};
  Mdl M;
  M.p = this;
  if (reinitializationLanguage_ == LANGUAGE_NAME::Cpp)
  {
    error = CppReinitialization(
        reinterpret_cast<KIM::ModelReinitialization *>(&M));
  }
  else if (reinitializationLanguage_ == LANGUAGE_NAME::C)
  {
    KIM_ModelReinitialization cM;
    cM.p = &M;
    error = CReinitialization(&cM);
  }
  else if (reinitializationLanguage_ == LANGUAGE_NAME::Fortran)
  {
    KIM_ModelReinitialization cM;
    cM.p = &M;
    FReinitialization(&cM, &error);
  }
  else
  {
    return true;
  }

  if (error)
    return true;
  else
    return false;
}

int ModelImplementation::get_neigh(int const neighborListIndex,
                                   int const particleNumber,
                                   int * const numberOfNeighbors,
                                   int const ** const neighborsOfParticle)
    const
{
  auto languageResult = callBackLanguage_.find(CALL_BACK_NAME::get_neigh);
  if (languageResult == callBackLanguage_.end())
  {
    // @@@@ log message
    return true;
  }
  LanguageName languageName = languageResult->second;
  void const * dataObject
      = (callBackDataObjectPointer_.find(CALL_BACK_NAME::get_neigh))->second;

  func * functionPointer
      = (callBackFunctionPointer_.find(CALL_BACK_NAME::get_neigh))->second;
  typedef int get_NeighCpp(void const * const dataObject,
                           int const neighborListIndex,
                           int const particleNumber,
                           int * const numberOfNeighbors,
                           int const ** const neighborsOfParticle);
  get_NeighCpp * CppGet_Neigh
      = reinterpret_cast<get_NeighCpp *>(functionPointer);
  typedef int get_NeighC(void const * const dataObject,
                         int const neighborListIndex,
                         int const particleNumber,
                         int * const numberOfNeighbors,
                         int const ** const neighborsOfParticle);
  get_NeighC * CGet_Neigh = reinterpret_cast<get_NeighC *>(functionPointer);
  typedef void get_NeighF(void const * const dataObject,
                          int const neighborListIndex,
                          int const particleNumber,
                          int * const numberOfNeighbors,
                          int const ** const neighborsOfParticle,
                          int * const ierr);
  get_NeighF * FGet_Neigh = reinterpret_cast<get_NeighF *>(functionPointer);


  int simulatorParticleNumber = particleNumber +
      ((simulatorNumbering_ == modelNumbering_) ? 0 : -numberingOffset_);
  int const * simulatorNeighborsOfParticle;
  int error;
  if (languageName == LANGUAGE_NAME::Cpp)
  {
    error = CppGet_Neigh(dataObject, neighborListIndex,
                         simulatorParticleNumber, numberOfNeighbors,
                         &simulatorNeighborsOfParticle);
  }
  else if (languageName == LANGUAGE_NAME::C)
  {
    error = CGet_Neigh(dataObject, neighborListIndex, simulatorParticleNumber,
                       numberOfNeighbors, &simulatorNeighborsOfParticle);
  }
  else if (languageName == LANGUAGE_NAME::Fortran)
  {
    FGet_Neigh(dataObject, neighborListIndex+1, simulatorParticleNumber,
               numberOfNeighbors, &simulatorNeighborsOfParticle, &error);
  }
  else
  {
    return true;
  }

  if (error) return true;

  // account for numbering differences if needed
  if (simulatorNumbering_ != modelNumbering_)
  {
    std::vector<int> & list = getNeighborListStorage_[neighborListIndex];
    list.resize(*numberOfNeighbors);
    for (int i=0; i<*numberOfNeighbors; ++i)
      list[i] = simulatorNeighborsOfParticle[i] + numberingOffset_;

    *neighborsOfParticle = list.data();
  }
  else
  {
    *neighborsOfParticle = simulatorNeighborsOfParticle;
  }

  return false;
}

int ModelImplementation::process_dEdr(double const de, double const r,
                                      double const * const dx,
                                      int const i, int const j) const
{
  auto languageResult = callBackLanguage_.find(CALL_BACK_NAME::process_dEdr);
  if (languageResult == callBackLanguage_.end())
  {
    // @@@@ log message
    return true;
  }
  LanguageName languageName = languageResult->second;
  void const * dataObject
      = (callBackDataObjectPointer_.find(CALL_BACK_NAME::process_dEdr))->second;

  func * functionPointer
      = (callBackFunctionPointer_.find(CALL_BACK_NAME::process_dEdr))->second;
  typedef int process_dEdrCpp(void const * const dataObject, double const de,
                              double const r, double const * const dx,
                              int const i, int const j);
  process_dEdrCpp * CppProcess_dEdr
      = reinterpret_cast<process_dEdrCpp *>(functionPointer);
  typedef int process_dEdrC(void const * const dataObject, double const de,
                            double const r, double const * const dx,
                            int const i, int const j);
  process_dEdrC * CProcess_dEdr
      = reinterpret_cast<process_dEdrC *>(functionPointer);
  typedef void process_dEdrF(void const * const dataObject, double const de,
                             double const r, double const * const dx,
                             int const i, int const j, int * const ierr);
  process_dEdrF * FProcess_dEdr
      = reinterpret_cast<process_dEdrF *>(functionPointer);

  int error;
  if (languageName == LANGUAGE_NAME::Cpp)
  {
    error = CppProcess_dEdr(dataObject, de, r, dx, i, j);
  }
  else if (languageName == LANGUAGE_NAME::C)
  {
    error = CProcess_dEdr(dataObject, de, r, dx, i, j);
  }
  else if (languageName == LANGUAGE_NAME::Fortran)
  {
    FProcess_dEdr(dataObject, de, r, dx, i, j, &error);
  }
  else
  {
    return true;
  }

  if (error)
    return true;
  else
    return false;
}

int ModelImplementation::process_d2Edr2(double const de, double const * const r,
                                        double const * const dx,
                                        int const * const i,
                                        int const * const j) const
{
  auto languageResult = callBackLanguage_.find(CALL_BACK_NAME::process_d2Edr2);
  if (languageResult == callBackLanguage_.end())
  {
    // @@@@ log message
    return true;
  }
  LanguageName languageName = languageResult->second;
  void const * dataObject = (callBackDataObjectPointer_
                             .find(CALL_BACK_NAME::process_d2Edr2))->second;

  func * functionPointer
      = (callBackFunctionPointer_.find(CALL_BACK_NAME::process_d2Edr2))->second;
  typedef int process_d2Edr2Cpp(void const * const dataObject, double const de,
                                double const * const r, double const * const dx,
                                int const * const i, int const * const j);
  process_d2Edr2Cpp * CppProcess_d2Edr2
      = reinterpret_cast<process_d2Edr2Cpp *>(functionPointer);
  typedef int process_d2Edr2C(void const * const dataObject, double const de,
                              double const * const r, double const * const dx,
                              int const * const i, int const * const j);
  process_d2Edr2C * CProcess_d2Edr2
      = reinterpret_cast<process_d2Edr2C *>(functionPointer);
  typedef void process_d2Edr2F(void const * const dataObject, double const de,
                               double const * const r, double const * const dx,
                               int const * const i, int const * const j,
                               int * const ierr);
  process_d2Edr2F * FProcess_d2Edr2
      = reinterpret_cast<process_d2Edr2F *>(functionPointer);

  int error;
  if (languageName == LANGUAGE_NAME::Cpp)
  {
    error = CppProcess_d2Edr2(dataObject, de, r, dx, i, j);
  }
  else if (languageName == LANGUAGE_NAME::C)
  {
    error = CProcess_d2Edr2(dataObject, de, r, dx, i, j);
  }
  else if (languageName == LANGUAGE_NAME::Fortran)
  {
    FProcess_d2Edr2(dataObject, de, r, dx, i, j, &error);
  }
  else
  {
    return true;
  }

  if (error)
    return true;
  else
    return false;
}

void ModelImplementation::set_model_buffer(void * const ptr)
{
  modelBuffer_ = ptr;
}

void ModelImplementation::get_model_buffer(void ** const ptr) const
{
  *ptr = modelBuffer_;
}


void ModelImplementation::set_sim_buffer(void * const ptr)
{
  simulatorBuffer_ = ptr;
}

void ModelImplementation::get_sim_buffer(void ** const ptr) const
{
  *ptr = simulatorBuffer_;
}


int ModelImplementation::convert_unit(
    LengthUnit const fromLengthUnit,
    EnergyUnit const fromEnergyUnit,
    ChargeUnit const fromChargeUnit,
    TemperatureUnit const fromTemperatureUnit,
    TimeUnit const fromTimeUnit,
    LengthUnit const toLengthUnit,
    EnergyUnit const toEnergyUnit,
    ChargeUnit const toChargeUnit,
    TemperatureUnit const toTemperatureUnit,
    TimeUnit const toTimeUnit,
    double const lengthExponent,
    double const energyExponent,
    double const chargeExponent,
    double const temperatureExponent,
    double const timeExponent,
    double * const conversionFactor) const
{
  // @@@@@ currently does nothing...
  return false;
}

void ModelImplementation::Log(LogLevel const logLevel,
                              std::string const & message,
                              int const lineNumber,
                              std::string const & fileName) const
{
  KIM::Log(logLevel, message, lineNumber, fileName);
}

std::string ModelImplementation::string() const
{
  // @@@@@ currently does nothing...
  return "empty!!!!";
}


ModelImplementation::ModelImplementation(ModelLibrary * const modelLibrary) :
    modelLibrary_(modelLibrary),
    influenceDistance_(0),
    numberOfCutoffs_(0),
    cutoffs_(0),
    reinitializationFunction_(0),
    destroyFunction_(0),
    computeFunction_(0),
    modelBuffer_(0),
    simulatorBuffer_(0)
{
  // populate Arguments
  int numberOfArguments;
  ARGUMENT_NAME::get_number_of_arguments(&numberOfArguments);
  for (int i=0; i<numberOfArguments; ++i)
  {
    ArgumentName argumentName;
    ARGUMENT_NAME::get_argument_name(i, &argumentName);
    argumentAttribute_[argumentName] = ATTRIBUTE::notSupported;
  }
  // populate mandatory Arguments
  for (auto mandatoryArgument = ARGUMENT_NAME::mandatoryArguments.begin();
       mandatoryArgument != ARGUMENT_NAME::mandatoryArguments.end();
       ++mandatoryArgument)
  {
    argumentAttribute_[*mandatoryArgument] = ATTRIBUTE::mandatory;
  }

  // populate CallBacks
  int numberOfCallBacks;
  CALL_BACK_NAME::get_number_of_call_backs(&numberOfCallBacks);
  for (int i=0; i<numberOfCallBacks; ++i)
  {
    CallBackName callBackName;
    CALL_BACK_NAME::get_call_back_name(i, &callBackName);
    callBackAttribute_[callBackName] = ATTRIBUTE::notSupported;
  }
  // populate mandatory CallBacks
  for (auto mandatoryCallBack = CALL_BACK_NAME::mandatoryCallBacks.begin();
       mandatoryCallBack != CALL_BACK_NAME::mandatoryCallBacks.end();
       ++mandatoryCallBack)
  {
    callBackAttribute_[*mandatoryCallBack] = ATTRIBUTE::mandatory;
  }
}

ModelImplementation::~ModelImplementation()
{
  delete modelLibrary_;
}

int ModelImplementation::ModelInitialization(
    Numbering const numbering,
    LengthUnit const requestedLengthUnit,
    EnergyUnit const requestedEnergyUnit,
    ChargeUnit const requestedChargeUnit,
    TemperatureUnit const requestedTemperatureUnit,
    TimeUnit const requestedTimeUnit,
    std::string const & modelName)
{
  int error = set_simulator_numbering(numbering);
  if (error) return true;

  error = modelLibrary_->open(modelName);
  if (error) return true;

  LanguageName languageName;
  func * functionPointer = 0;
  error = modelLibrary_->getModelInitializationFunctionPointer(
      &languageName, &functionPointer);
  if (error) return true;

  typedef int ModelInitializationCpp(
      KIM::ModelInitialization * const modelInitialization,
      LengthUnit const requestedLengthUnit,
      EnergyUnit const requestedEnergyUnit,
      ChargeUnit const requestedChargeUnit,
      TemperatureUnit const requestedTemperatureUnit,
      TimeUnit const requestedTimeUnit);
  ModelInitializationCpp * CppInitialization
      = reinterpret_cast<ModelInitializationCpp *>(functionPointer);
  typedef int ModelInitializationC(
      KIM_ModelInitialization * const modelInitialization,
      KIM_LengthUnit const requestedLengthUnit,
      KIM_EnergyUnit const requestedEnergyUnit,
      KIM_ChargeUnit const requestedChargeUnit,
      KIM_TemperatureUnit const requestedTemperatureUnit,
      KIM_TimeUnit const requestedTimeUnit);
  ModelInitializationC * CInitialization
      = reinterpret_cast<ModelInitializationC *>(functionPointer);
  typedef void ModelInitializationF(
      KIM_ModelInitialization * const modelInitialization,
      KIM_LengthUnit const requestedLengthUnit,
      KIM_EnergyUnit const requestedEnergyUnit,
      KIM_ChargeUnit const requestedChargeUnit,
      KIM_TemperatureUnit const requestedTemperatureUnit,
      KIM_TimeUnit const requestedTimeUnit,
      int * const);
  ModelInitializationF * FInitialization
      = reinterpret_cast<ModelInitializationF *>(functionPointer);

  struct Mdl {void * p;};
  Mdl M;
  M.p = this;
  KIM_LengthUnit requestedLengthUnitC = makeLengthUnitC(requestedLengthUnit);
  KIM_EnergyUnit requestedEnergyUnitC = makeEnergyUnitC(requestedEnergyUnit);
  KIM_ChargeUnit requestedChargeUnitC = makeChargeUnitC(requestedChargeUnit);
  KIM_TemperatureUnit requestedTemperatureUnitC
      = makeTemperatureUnitC(requestedTemperatureUnit);
  KIM_TimeUnit requestedTimeUnitC = makeTimeUnitC(requestedTimeUnit);
  if (languageName == LANGUAGE_NAME::Cpp)
  {
    error = CppInitialization(reinterpret_cast<KIM::ModelInitialization *>(&M),
                              requestedLengthUnit, requestedEnergyUnit,
                              requestedChargeUnit, requestedTemperatureUnit,
                              requestedTimeUnit);
  }
  else if (languageName == LANGUAGE_NAME::C)
  {
    KIM_ModelInitialization cM;
    cM.p = &M;
    error = CInitialization(&cM, requestedLengthUnitC, requestedEnergyUnitC,
                            requestedChargeUnitC, requestedTemperatureUnitC,
                            requestedTimeUnitC);

  }
  else if (languageName == LANGUAGE_NAME::Fortran)
  {
    KIM_ModelInitialization cM;
    cM.p = &M;
    FInitialization(&cM, requestedLengthUnitC, requestedEnergyUnitC,
                    requestedChargeUnitC, requestedTemperatureUnitC,
                    requestedTimeUnitC, &error);
  }
  else
  {
    return true;
  }
  if (error) return true;

  // set numberingOffset_
  if (simulatorNumbering_ == modelNumbering_)
    numberingOffset_ = 0;
  else if (simulatorNumbering_ == NUMBERING::zeroBased)
    numberingOffset_ = 1;
  else
    numberingOffset_ = -1;

  // resize getNeighborListStorage_
  if (simulatorNumbering_ != modelNumbering_)
    getNeighborListStorage_.resize(numberOfCutoffs_);

  return false;
}

int ModelImplementation::ModelDestroy()
{
  typedef int ModelDestroyCpp(KIM::ModelDestroy * const);
  ModelDestroyCpp * CppDestroy
      = reinterpret_cast<ModelDestroyCpp *>(destroyFunction_);
  typedef int ModelDestroyC(KIM_ModelDestroy * const);
  ModelDestroyC * CDestroy
      = reinterpret_cast<ModelDestroyC *>(destroyFunction_);
  typedef void ModelDestroyF(KIM_ModelDestroy * const, int * const);
  ModelDestroyF * FDestroy
      = reinterpret_cast<ModelDestroyF *>(destroyFunction_);

  int error;
  struct Mdl {void * p;};
  Mdl M;
  M.p = this;
  if (destroyLanguage_ == LANGUAGE_NAME::Cpp)
  {
    error = CppDestroy(reinterpret_cast<KIM::ModelDestroy *>(&M));
  }
  else if (destroyLanguage_ == LANGUAGE_NAME::C)
  {
    KIM_ModelDestroy cM;
    cM.p = &M;
    error = CDestroy(&cM);
  }
  else if (destroyLanguage_ == LANGUAGE_NAME::Fortran)
  {
    KIM_ModelDestroy cM;
    cM.p = &M;
    FDestroy(&cM, &error);
  }
  else
  {
    return true;
  }

  if (error)
    return true;
  else
    return false;
}

}  // namespace KIM
