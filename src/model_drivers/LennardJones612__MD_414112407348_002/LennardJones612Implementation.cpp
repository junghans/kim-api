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
// Copyright (c) 2013--2015, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Ryan S. Elliott
//    Stephen M. Whalen
//    Andrew Akerson
//


#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include "LennardJones612.hpp"
#include "LennardJones612Implementation.hpp"
#include "KIM_LanguageName.hpp"
#include "KIM_SpeciesName.hpp"
#include "KIM_UnitSystem.hpp"
#include "KIM_ArgumentName.hpp"
#include "KIM_CallBackName.hpp"

#define MAXLINE 1024
#define IGNORE_RESULT(fn) if(fn){}


//==============================================================================
//
// Implementation of LennardJones612Implementation public member functions
//
//==============================================================================

//******************************************************************************
LennardJones612Implementation::LennardJones612Implementation(
    KIM::ModelInitialization * const modelInitialization,
    char const * const  parameterFileNames,
    int const parameterFileNameLength,
    int const numberParameterFiles,
    int * const ier)
    : numberModelSpecies_(0),
      numberUniqueSpeciesPairs_(0),
      shift_(0),
      cutoffs_(0),
      epsilons_(0),
      sigmas_(0),
      influenceDistance_(0.0),
      cutoffsSq2D_(0),
      fourEpsilonSigma6_2D_(0),
      fourEpsilonSigma12_2D_(0),
      twentyFourEpsilonSigma6_2D_(0),
      fortyEightEpsilonSigma12_2D_(0),
      oneSixtyEightEpsilonSigma6_2D_(0),
      sixTwentyFourEpsilonSigma12_2D_(0),
      shifts2D_(0),
      cachedNumberOfParticles_(0)
{
  *ier = SetConstantValues(modelInitialization);
  if (*ier) return;

  AllocateFreeParameterMemory();

  FILE* parameterFilePointers[MAX_PARAMETER_FILES];
  *ier = OpenParameterFiles(parameterFileNames, parameterFileNameLength,
                            numberParameterFiles, parameterFilePointers);
  if (*ier) return;

  *ier = ProcessParameterFiles(modelInitialization, parameterFilePointers,
                               numberParameterFiles);
  CloseParameterFiles(parameterFilePointers, numberParameterFiles);
  if (*ier) return;

  *ier = ConvertUnits(modelInitialization);
  if (*ier) return;

  *ier = SetReinitMutableValues(modelInitialization);
  if (*ier) return;

  *ier = RegisterKIMParameters(modelInitialization);
  if (*ier) return;

  *ier = RegisterKIMFunctions(modelInitialization);
  if (*ier) return;

  // everything is good
  *ier = false;
  return;
}

//******************************************************************************
LennardJones612Implementation::~LennardJones612Implementation()
{ // note: it is ok to delete a null pointer and we have ensured that
  // everything is initialized to null

  delete [] cutoffs_;
  Deallocate2DArray(cutoffsSq2D_);
  delete [] epsilons_;
  delete [] sigmas_;
  Deallocate2DArray(fourEpsilonSigma6_2D_);
  Deallocate2DArray(fourEpsilonSigma12_2D_);
  Deallocate2DArray(twentyFourEpsilonSigma6_2D_);
  Deallocate2DArray(fortyEightEpsilonSigma12_2D_);
  Deallocate2DArray(oneSixtyEightEpsilonSigma6_2D_);
  Deallocate2DArray(sixTwentyFourEpsilonSigma12_2D_);
  Deallocate2DArray(shifts2D_);
}

//******************************************************************************
int LennardJones612Implementation::Reinit(
    KIM::ModelReinitialization * const modelReinitialization)
{
  int ier;

  ier = SetReinitMutableValues(modelReinitialization);
  if (ier) return ier;

  // nothing else to do for this case

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int LennardJones612Implementation::Compute(
    KIM::ModelCompute const * const modelCompute)
{
  int ier;

  // KIM API Model Input compute flags
  bool isComputeProcess_dEdr = false;
  bool isComputeProcess_d2Edr2 = false;
  //
  // KIM API Model Output compute flags
  bool isComputeEnergy = false;
  bool isComputeForces = false;
  bool isComputeParticleEnergy = false;
  //
  // KIM API Model Input
  int const* particleSpecies = 0;
  int const* particleContributing = 0;
  VectorOfSizeDIM const* coordinates = 0;
  //
  // KIM API Model Output
  double* energy = 0;
  double* particleEnergy = 0;
  VectorOfSizeDIM* forces = 0;
  ier = SetComputeMutableValues(modelCompute, isComputeProcess_dEdr,
                                isComputeProcess_d2Edr2, isComputeEnergy,
                                isComputeForces, isComputeParticleEnergy,
                                particleSpecies, particleContributing,
                                coordinates, energy, particleEnergy, forces);
  if (ier) return ier;

  // Skip this check for efficiency
  //
  // ier = CheckParticleSpecies(pkim, particleSpecies);
  // if (ier) return ier;

  bool const isShift = (1 == shift_);

#include "LennardJones612ImplementationComputeDispatch.cpp"
  return ier;
}

//==============================================================================
//
// Implementation of LennardJones612Implementation private member functions
//
//==============================================================================

//******************************************************************************
int LennardJones612Implementation::
SetConstantValues(KIM::ModelInitialization const * const modelInitialization)
{
  // set numberModelSpecies & numberUniqueSpeciesPairs
  modelInitialization->get_num_model_species(&numberModelSpecies_);
  numberUniqueSpeciesPairs_ = ((numberModelSpecies_+1)*numberModelSpecies_)/2;

  // everything is good
  return false;
}

//******************************************************************************
int LennardJones612Implementation::OpenParameterFiles(
    char const* const parameterFileNames,
    int const parameterFileNameLength,
    int const numberParameterFiles,
    FILE* parameterFilePointers[MAX_PARAMETER_FILES])
{
  int ier;

  if (numberParameterFiles > MAX_PARAMETER_FILES)
  {
    ier = true;
    KIM::report_error(__LINE__, __FILE__, "LennardJones612 given too many"
                       " parameter files", ier);
    return ier;
  }

  for (int i = 0; i < numberParameterFiles; ++i)
  {
    parameterFilePointers[i]
        = fopen(&parameterFileNames[i * parameterFileNameLength], "r");
    if (parameterFilePointers[i] == 0)
    {
      char message[MAXLINE];
      sprintf(message,
              "LennardJones612 parameter file number %d cannot be opened",
              i);
      ier = true;
      KIM::report_error(__LINE__, __FILE__, message, ier);
      for (int j = i - 1; i <= 0; --i)
      {
        fclose(parameterFilePointers[j]);
      }
      return ier;
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int LennardJones612Implementation::ProcessParameterFiles(
    KIM::ModelInitialization const * const modelInitialization,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles)
{
  int N, ier;
  int endOfFileFlag = 0;
  char spec1[MAXLINE], spec2[MAXLINE], nextLine[MAXLINE];
  char *nextLinePtr;
  int iIndex, jIndex , indx, iiIndex, jjIndex;
  double nextCutoff, nextEpsilon, nextSigma;
  nextLinePtr = nextLine;

  getNextDataLine(parameterFilePointers[0], nextLinePtr,
                  MAXLINE, &endOfFileFlag);
  ier = sscanf(nextLine, "%d %d", &N, &shift_);
  if (ier != 2)
  {
    sprintf(nextLine, "unable to read first line of the parameter file");
    ier = true;
    KIM::report_error(__LINE__, __FILE__, nextLine, ier);
    fclose(parameterFilePointers[0]);
    return ier;
  }
  if (N != numberModelSpecies_)
  {
    sprintf(nextLine, "The value for N from the parameter file is inconsistent "
            "with numberModelSpecies_");
    ier = true;
    KIM::report_error(__LINE__, __FILE__, nextLine, ier);
    fclose(parameterFilePointers[0]);
    return ier;
  }

  // get and correctly order the particle names
  std::string* const particleNames = new std::string[numberModelSpecies_];
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    KIM::SpeciesName kimModelParticleSpecies;
    ier = modelInitialization->get_model_species(i, &kimModelParticleSpecies);
    modelSpeciesNameList_.push_back(kimModelParticleSpecies);
    if (ier)
    {
      KIM::report_error(__LINE__, __FILE__, "get_model_species", ier);
      delete [] particleNames;
      return ier;
    }
    int index;
    ier = modelInitialization->get_species_code(kimModelParticleSpecies,
                                                &index);
    modelSpeciesCodeList_.push_back(index);
    if (index >= numberModelSpecies_)
    {
      KIM::report_error(__LINE__, __FILE__, "get_species_code",
                        true);
      delete [] particleNames;
      return true;
    }
    particleNames[index] = kimModelParticleSpecies.string();
  }

  // set all values in the arrays to -1 for mixing later
  for (int i = 0; i < ((N+1)*N/2); i++)
  {
    cutoffs_[i]  = -1;
    epsilons_[i] = -1;
    sigmas_[i] = -1;
  }

  // Read and process data lines
  getNextDataLine(parameterFilePointers[0], nextLinePtr,
                  MAXLINE, &endOfFileFlag);
  while (endOfFileFlag == 0)
  {
    ier = sscanf(nextLine, "%s  %s %lg %lg %lg",
	         spec1, spec2, &nextCutoff, &nextEpsilon, &nextSigma);
    if (ier != 5)
    {
      sprintf(nextLine, "error reading lines of the parameter file");
      KIM::report_error(__LINE__, __FILE__, nextLine, true);
      delete [] particleNames;
      return true;
    }
    iIndex = jIndex = -1;
    for (int i = 0; i <  N; i++)
    {
      if (strcmp(spec1, particleNames[i].c_str()) == 0)
      {
        iIndex = i;
      }
      if (strcmp(spec2, particleNames[i].c_str()) == 0)
      {
        jIndex = i;
      }
    }
    if ((iIndex == -1) || (jIndex == -1))
    {
      sprintf(nextLine, "Unsupported Species name found in parameter file");
      KIM::report_error(__LINE__, __FILE__, nextLine, true);
      delete [] particleNames;
      return true;
    }
    if (iIndex >= jIndex)
    {
      indx = jIndex*N + iIndex - (jIndex*jIndex + jIndex)/2;
    }
    else
    {
      indx = iIndex*N + jIndex - (iIndex*iIndex + iIndex)/2;
    }
    cutoffs_[indx] = nextCutoff;
    epsilons_[indx] = nextEpsilon;
    sigmas_[indx] = nextSigma;

    getNextDataLine(parameterFilePointers[0], nextLinePtr,
                    MAXLINE, &endOfFileFlag);
  }

  // check that we got all like - like pairs
  sprintf(nextLine, "There are not values for like-like pairs of:");
  for (int i = 0; i < N; i++)
  {
    if (cutoffs_[(i*N + i - (i*i + i)/2)] == -1)
    {
      strcat(nextLine, "  ");
      strcat(nextLine, particleNames[i].c_str());
      ier = -1;
    }
  }
  if (ier == -1)
  {
    KIM::report_error(__LINE__, __FILE__, nextLine, true);
    delete [] particleNames;
    return true;
  }

  // Perform Mixing if nessisary
  for (int jIndex = 0; jIndex < N; jIndex++)
  {
    jjIndex = (jIndex*N + jIndex - (jIndex*jIndex + jIndex)/2);
    for (int iIndex = (jIndex+1) ; iIndex < N; iIndex++)
    {
      indx = jIndex*N + iIndex - (jIndex*jIndex + jIndex)/2;
      if (cutoffs_[indx] == -1)
      {
        iiIndex = (iIndex*N + iIndex - (iIndex*iIndex + iIndex)/2);
        epsilons_[indx] = sqrt(epsilons_[iiIndex]*epsilons_[jjIndex]);
        sigmas_[indx] = (sigmas_[iiIndex] + sigmas_[jjIndex])/2.0;
        cutoffs_[indx] = (cutoffs_[iiIndex] + cutoffs_[jjIndex])/2.0;
      }
    }
  }
  delete [] particleNames;

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
void LennardJones612Implementation::getNextDataLine(
    FILE* const filePtr, char* nextLinePtr, int const maxSize,
    int *endOfFileFlag)
{
  do
  {
    if(fgets(nextLinePtr, maxSize, filePtr) == NULL)
    {
       *endOfFileFlag = 1;
       break;
    }
    while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t') ||
           (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r' ))
    {
      nextLinePtr = (nextLinePtr + 1);
    }
  }
  while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
}

//******************************************************************************
void LennardJones612Implementation::CloseParameterFiles(
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles)
{
  for (int i = 0; i < numberParameterFiles; ++i)
    fclose(parameterFilePointers[i]);
}

//******************************************************************************
void LennardJones612Implementation::AllocateFreeParameterMemory()
{ // allocate memory for data
  cutoffs_ = new double[numberUniqueSpeciesPairs_];
  AllocateAndInitialize2DArray(cutoffsSq2D_, numberModelSpecies_,
                               numberModelSpecies_);

  epsilons_ = new double[numberUniqueSpeciesPairs_];
  sigmas_ = new double[numberUniqueSpeciesPairs_];
  AllocateAndInitialize2DArray(fourEpsilonSigma6_2D_, numberModelSpecies_,
                               numberModelSpecies_);
  AllocateAndInitialize2DArray(fourEpsilonSigma12_2D_, numberModelSpecies_,
                               numberModelSpecies_);
  AllocateAndInitialize2DArray(twentyFourEpsilonSigma6_2D_, numberModelSpecies_,
                               numberModelSpecies_);
  AllocateAndInitialize2DArray(fortyEightEpsilonSigma12_2D_,
                               numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray(oneSixtyEightEpsilonSigma6_2D_,
                               numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray(sixTwentyFourEpsilonSigma12_2D_,
                               numberModelSpecies_, numberModelSpecies_);

  AllocateAndInitialize2DArray(shifts2D_, numberModelSpecies_,
                               numberModelSpecies_);
}

//******************************************************************************
int LennardJones612Implementation::ConvertUnits(
    KIM::ModelInitialization const * const modelInitialization)
{
  int ier;

  // define default base units
  KIM::LengthUnit length = KIM::LENGTH_UNIT::A;
  KIM::EnergyUnit energy = KIM::ENERGY_UNIT::eV;
  KIM::ChargeUnit charge = KIM::CHARGE_UNIT::e;
  KIM::TemperatureUnit temperature = KIM::TEMPERATURE_UNIT::K;
  KIM::TimeUnit time = KIM::TIME_UNIT::ps;

  // changing units of cutoffs and sigmas
  double convertLength;
  ier = modelInitialization->convert_to_act_unit(length, energy, charge,
                                                 temperature, time,
                                                 1.0, 0.0, 0.0, 0.0, 0.0,
                                                 &convertLength);
  if (ier)
  {
    KIM::report_error(__LINE__, __FILE__, "convert_to_act_unit", ier);
    return ier;
  }
  if (convertLength != ONE)
  {
    for (int i = 0; i < numberUniqueSpeciesPairs_; ++i)
    {
      cutoffs_[i] *= convertLength;  // convert to active units
      sigmas_[i] *= convertLength;  // convert to active units
    }
  }
  // changing units of epsilons
  double convertEnergy;
  ier = modelInitialization->convert_to_act_unit(length, energy, charge,
                                                 temperature, time,
                                                 0.0, 1.0, 0.0, 0.0, 0.0,
                                                 &convertEnergy);
  if (ier)
  {
    KIM::report_error(__LINE__, __FILE__, "convert_to_act_unit", ier);
    return ier;
  }
  if (convertEnergy != ONE)
  {
    for (int i = 0; i < numberUniqueSpeciesPairs_; ++i)
    {
      epsilons_[i] *= convertEnergy;  // convert to active units
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int LennardJones612Implementation::RegisterKIMParameters(
    KIM::ModelInitialization * const modelInitialization) const
{
  int ier = false;

  // publish parameters
  ier = modelInitialization->set_parameter(PARAM_SHIFT_INDEX, 1,
                                           (void *) &shift_);
  if (ier)
  {
    KIM::report_error(__LINE__, __FILE__, "set_parameter", ier);
    return ier;
  }
  ier = modelInitialization->set_parameter(PARAM_CUTOFFS_INDEX,
                                           numberUniqueSpeciesPairs_,
                                           (void *) cutoffs_);
  if (ier)
  {
    KIM::report_error(__LINE__, __FILE__, "set_parameter", ier);
    return ier;
  }
  ier = modelInitialization->set_parameter(PARAM_EPSILONS_INDEX,
                                           numberUniqueSpeciesPairs_,
                                           (void *) epsilons_);
  if (ier)
  {
    KIM::report_error(__LINE__, __FILE__, "set_parameter", ier);
    return ier;
  }
  ier = modelInitialization->set_parameter(PARAM_SIGMAS_INDEX,
                                           numberUniqueSpeciesPairs_,
                                           (void *) sigmas_);
  if (ier)
  {
    KIM::report_error(__LINE__, __FILE__, "set_parameter", ier);
    return ier;
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int LennardJones612Implementation::RegisterKIMFunctions(
    KIM::ModelInitialization * const modelInitialization)
    const
{
  int ier;

  // register the destroy() and reinit() functions
  modelInitialization->set_destroy(
      KIM::LANGUAGE_NAME::Cpp, (KIM::func*) &(LennardJones612::Destroy));
  modelInitialization->set_reinit(
      KIM::LANGUAGE_NAME::Cpp, (KIM::func*) &(LennardJones612::Reinit));
  modelInitialization->set_compute_func(
      KIM::LANGUAGE_NAME::Cpp, (KIM::func*) &(LennardJones612::Compute));

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
template<class ModelObj>
int LennardJones612Implementation::SetReinitMutableValues(
    ModelObj * const modelObj)
{ // use (possibly) new values of free parameters to compute other quantities
  int ier;

  // update cutoffsSq, epsilons, and sigmas
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    for (int j = 0; j <= i ; ++j)
    {
      int const index = j*numberModelSpecies_ + i - (j*j + j)/2;
      cutoffsSq2D_[i][j] = cutoffsSq2D_[j][i]
          = (cutoffs_[index]*cutoffs_[index]);
      fourEpsilonSigma6_2D_[i][j] = fourEpsilonSigma6_2D_[j][i]
          = 4.0*epsilons_[index]*pow(sigmas_[index],6.0);
      fourEpsilonSigma12_2D_[i][j] = fourEpsilonSigma12_2D_[j][i]
          = 4.0*epsilons_[index]*pow(sigmas_[index],12.0);
      twentyFourEpsilonSigma6_2D_[i][j] = twentyFourEpsilonSigma6_2D_[j][i]
          = 6.0*fourEpsilonSigma6_2D_[i][j];
      fortyEightEpsilonSigma12_2D_[i][j] = fortyEightEpsilonSigma12_2D_[j][i]
          = 12.0*fourEpsilonSigma12_2D_[i][j];
      oneSixtyEightEpsilonSigma6_2D_[i][j]
          = oneSixtyEightEpsilonSigma6_2D_[j][i]
          = 7.0*twentyFourEpsilonSigma6_2D_[i][j];
      sixTwentyFourEpsilonSigma12_2D_[i][j]
          = sixTwentyFourEpsilonSigma12_2D_[j][i]
          = 13.0*fortyEightEpsilonSigma12_2D_[i][j];
    }
  }

  // update cutoff value in KIM API object
  influenceDistance_ = 0.0;

  for (int i = 0; i < numberModelSpecies_; i++)
  {
    int indexI = modelSpeciesCodeList_[i];

    for (int j = 0; j < numberModelSpecies_; j++)
    {
      int indexJ = modelSpeciesCodeList_[j];

      if (influenceDistance_ < cutoffsSq2D_[indexI][indexJ])
      {
        influenceDistance_ = cutoffsSq2D_[indexI][indexJ];
      }
    }
  }
  influenceDistance_ = sqrt(influenceDistance_);
  modelObj->set_influence_distance(&influenceDistance_);
  modelObj->set_cutoffs(1, &influenceDistance_);

  // update shifts
  // compute and set shifts2D_ check if minus sign
  double const* const* const  constFourEpsSig6_2D = fourEpsilonSigma6_2D_;
  double const* const* const  constFourEpsSig12_2D = fourEpsilonSigma12_2D_;
  if (1 == shift_)
  {
    double phi;
    for (int iSpecies = 0; iSpecies < numberModelSpecies_; iSpecies++)
    {
      for(int jSpecies = 0; jSpecies <= iSpecies; jSpecies++)
      {
        int const index = jSpecies*numberModelSpecies_ + iSpecies
                           - (jSpecies*jSpecies + jSpecies)/2;
        double const rij2 = cutoffs_[index]*cutoffs_[index];
        double const r2iv = 1.0/rij2;
        double const r6iv = r2iv*r2iv*r2iv;
        LENNARD_JONES_PHI(;);
        shifts2D_[iSpecies][jSpecies] = shifts2D_[jSpecies][iSpecies] = phi;
      }
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int LennardJones612Implementation::SetComputeMutableValues(
    KIM::ModelCompute const * const modelCompute,
    bool& isComputeProcess_dEdr,
    bool& isComputeProcess_d2Edr2,
    bool& isComputeEnergy,
    bool& isComputeForces,
    bool& isComputeParticleEnergy,
    int const*& particleSpecies,
    int const*& particleContributing,
    VectorOfSizeDIM const*& coordinates,
    double*& energy,
    double*& particleEnergy,
    VectorOfSizeDIM*& forces)
{
  int ier = true;

  // get compute flags
  int compProcess_dEdr;
  int compProcess_d2Edr2;

  modelCompute->is_call_back_present(KIM::CALL_BACK_NAME::process_dEdr,
                                     &compProcess_dEdr);
  modelCompute->is_call_back_present(KIM::CALL_BACK_NAME::process_d2Edr2,
                                     &compProcess_d2Edr2);

  isComputeProcess_dEdr = compProcess_dEdr;
  isComputeProcess_d2Edr2 = compProcess_d2Edr2;

  // double const* cutoff;            // currently unused
  // int const* numberOfSpecies;  // currently unused
  int const* numberOfParticles;
  ier =
      modelCompute->get_data(KIM::ARGUMENT_NAME::numberOfParticles,
                             &numberOfParticles)
      || modelCompute->get_data(KIM::ARGUMENT_NAME::particleSpecies,
                                &particleSpecies)
      || modelCompute->get_data(KIM::ARGUMENT_NAME::particleContributing,
                                &particleContributing)
      || modelCompute->get_data(KIM::ARGUMENT_NAME::coordinates,
                                (double const ** const) &coordinates)
      || modelCompute->get_data(KIM::ARGUMENT_NAME::energy, &energy)
      || modelCompute->get_data(KIM::ARGUMENT_NAME::particleEnergy,
                                &particleEnergy)
      || modelCompute->get_data(KIM::ARGUMENT_NAME::forces,
                                (double const ** const) &forces);
  if (ier)
  {
    KIM::report_error(__LINE__, __FILE__, "get_data", ier);
    return ier;
  }

  isComputeEnergy = (energy != 0);
  isComputeParticleEnergy = (particleEnergy != 0);
  isComputeForces = (forces != 0);

  // update values
  cachedNumberOfParticles_ = *numberOfParticles;

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int LennardJones612Implementation::CheckParticleSpecies(
    int const* const particleSpecies)
    const
{
  int ier;
  for (int i = 0; i < cachedNumberOfParticles_; ++i)
  {
    if ((particleSpecies[i] < 0) || (particleSpecies[i] >= numberModelSpecies_))
    {
      ier = true;
      KIM::report_error(__LINE__, __FILE__,
                        "unsupported particle species detected", ier);
      return ier;
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int LennardJones612Implementation::GetComputeIndex(
    const bool& isComputeProcess_dEdr,
    const bool& isComputeProcess_d2Edr2,
    const bool& isComputeEnergy,
    const bool& isComputeForces,
    const bool& isComputeParticleEnergy,
    const bool& isShift) const
{
  //const int processdE = 2;
  const int processd2E = 2;
  const int energy = 2;
  const int force = 2;
  const int particleEnergy = 2;
  const int shift = 2;


  int index = 0;

  // processdE
  index += (int(isComputeProcess_dEdr))
      * processd2E * energy * force * particleEnergy * shift;

  // processd2E
  index += (int(isComputeProcess_d2Edr2))
      * energy * force * particleEnergy * shift;

  // energy
  index += (int(isComputeEnergy))
      * force * particleEnergy * shift;

  // force
  index += (int(isComputeForces))
      * particleEnergy * shift;

  // particleEnergy
  index += (int(isComputeParticleEnergy))
      * shift;

  // shift
  index += (int(isShift));

  return index;
}

//==============================================================================
//
// Implementation of helper functions
//
//==============================================================================

//******************************************************************************
void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
                                  int const extentOne)
{ // allocate memory and set pointers
  arrayPtr = new double*[extentZero];
  arrayPtr[0] = new double[extentZero * extentOne];
  for (int i = 1; i < extentZero; ++i)
  {
    arrayPtr[i] = arrayPtr[i-1] + extentOne;
  }

  // initialize
  for (int i = 0; i < extentZero; ++i)
  {
    for (int j = 0; j < extentOne; ++j)
    {
      arrayPtr[i][j] = 0.0;
    }
  }
}

//******************************************************************************
void Deallocate2DArray(double**& arrayPtr)
{ // deallocate memory
  if (arrayPtr != 0) delete [] arrayPtr[0];
  delete [] arrayPtr;

  // nullify pointer
  arrayPtr = 0;
}
