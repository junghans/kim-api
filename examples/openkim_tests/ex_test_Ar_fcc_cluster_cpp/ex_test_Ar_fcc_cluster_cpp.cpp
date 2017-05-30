//
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
//
// Copyright (c) 2013--2017, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Ryan S. Elliott
//    Stephen M. Whalen
//
//

//
// Release: This file is part of the kim-api.git repository.
//

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include "KIM_LogLevel.hpp"
#include "KIM_Log.hpp"
#include "KIM_LanguageName.hpp"
#include "KIM_DataType.hpp"
#include "KIM_SpeciesName.hpp"
#include "KIM_Numbering.hpp"
#include "KIM_Model.hpp"
#include "KIM_ArgumentName.hpp"
#include "KIM_CallBackName.hpp"
#include "KIM_Attribute.hpp"
#include "KIM_UnitSystem.hpp"
#include "KIM_Log.hpp"

#define NAMESTRLEN    128

#define FCCSPACING    5.260
#define DIM           3
#define NCELLSPERSIDE 2
#define NCLUSTERPARTS (4*(NCELLSPERSIDE*NCELLSPERSIDE*NCELLSPERSIDE) + \
                       6*(NCELLSPERSIDE*NCELLSPERSIDE)                 \
                       + 3*(NCELLSPERSIDE) + 1)

#define MY_ERROR(message)                                       \
  {                                                             \
    std::cout << "* Error : \"" << message << "\" : "           \
              << __LINE__ << ":" << __FILE__ << std::endl;      \
    exit(1);                                                    \
  }

#define MY_WARNING(message)                                             \
  {                                                                     \
    std::cout << "* Error : \"" << message << "\" : "                   \
              << __LINE__ << ":" << __FILE__ << std::endl;              \
  }


/* Define neighborlist structure */
typedef struct
{
  int numberOfParticles;
  int iteratorId;
  int* NNeighbors;
  int* neighborList;
} NeighList;

/* Define prototypes */
void fcc_cluster_neighborlist(int half, int numberOfParticles, double* coords,
                              double cutoff, NeighList* nl);

int get_cluster_neigh(void const * const dataObject,
                      int const neighborListIndex, int const particleNumber,
                      int * const numberOfNeighbors,
                      int const ** const neighborsOfParticle);

void create_FCC_cluster(double FCCspacing, int nCellsPerSide, double *coords);


/* Main program */
int main()
{
  /* Local variable declarations */
  double const MinSpacing = 0.8*FCCSPACING;
  double const MaxSpacing = 1.2*FCCSPACING;
  double const SpacingIncr = 0.025*FCCSPACING;
  double CurrentSpacing;
  double cutpad = 0.75; /* Angstroms */
  int i;
  int error;


  /* model inputs */
  int numberOfParticles_cluster = NCLUSTERPARTS;
  int particleSpecies_cluster_model[NCLUSTERPARTS];
  int particleContributing_cluster_model[NCLUSTERPARTS];
  double coords_cluster[NCLUSTERPARTS][DIM];
  NeighList nl_cluster_model;
  /* model outputs */
  double influence_distance_cluster_model;
  int number_of_cutoffs;
  double const * cutoff_cluster_model;
  double energy_cluster_model;

  std::string modelname;

  /* Get KIM Model names */
  std::cout << "Please enter valid KIM Model name: \n";
  std::cin >> modelname;


  /* initialize the model */
  KIM::Model * kim_cluster_model;
  int requestedUnitsAccepted;
  error = KIM::Model::create(
      KIM::NUMBERING::zeroBased,
      KIM::LENGTH_UNIT::A,
      KIM::ENERGY_UNIT::eV,
      KIM::CHARGE_UNIT::e,
      KIM::TEMPERATURE_UNIT::K,
      KIM::TIME_UNIT::ps,
      modelname,
      &requestedUnitsAccepted,
      &kim_cluster_model);
  if (error)
  {
    MY_ERROR("KIM_create_model_interface()");
  }

  // Check for compatibility with the model
  if (!requestedUnitsAccepted)
  {
    MY_ERROR("Must Adapt to model units");
  }

  // print model units
  KIM::LengthUnit lengthUnit;
  KIM::EnergyUnit energyUnit;
  KIM::ChargeUnit chargeUnit;
  KIM::TemperatureUnit temperatureUnit;
  KIM::TimeUnit timeUnit;

  kim_cluster_model->get_units(&lengthUnit, &energyUnit, &chargeUnit,
                               &temperatureUnit, &timeUnit);

  std::cout << "LengthUnit is \"" << lengthUnit.string() << "\"" << std::endl
            << "EnergyUnit is \"" << energyUnit.string() << "\"" << std::endl
            << "ChargeUnit is \"" << chargeUnit.string() << "\"" << std::endl
            << "TemperatureUnit is \"" << temperatureUnit.string()
            << "\"" << std::endl
            << "TimeUnit is \"" << timeUnit.string() << "\"" << std::endl;

  // check species
  int speciesIsSupported;
  int modelArCode;
  error = kim_cluster_model->get_species_support_and_code(
      KIM::SPECIES_NAME::Ar, &speciesIsSupported, &modelArCode);
  if ((error) || (!speciesIsSupported))
  {
    MY_ERROR("Species Ar not supported");
  }

  // check arguments
  int numberOfArguments;
  KIM::ARGUMENT_NAME::get_number_of_arguments(&numberOfArguments);
  for (int i=0; i<numberOfArguments; ++i)
  {
    KIM::ArgumentName argumentName;
    KIM::Attribute attribute;
    KIM::ARGUMENT_NAME::get_argument_name(i, &argumentName);
    KIM::DataType dataType;
    KIM::ARGUMENT_NAME::get_argument_data_type(argumentName, &dataType);
    error = kim_cluster_model->get_argument_attribute(argumentName, &attribute);
    if (error)
      MY_ERROR("unable to get argument attribute");

    std::cout << "Argument Name \""
              << argumentName.string() << "\""
              << " is of type \""
              << dataType.string() << "\""
              << " and has attribute \""
              << attribute.string() << "\""
              << std::endl;

    // can only handle energy as a required arg
    if (attribute == KIM::ATTRIBUTE::required)
    {
      if (argumentName != KIM::ARGUMENT_NAME::energy)
      {
        MY_ERROR("unsupported required argument");
      }
    }

    // must have energy
    if (argumentName == KIM::ARGUMENT_NAME::energy)
    {
      if (! ((attribute == KIM::ATTRIBUTE::required)
             ||
             (attribute == KIM::ATTRIBUTE::optional)))
      {
        MY_ERROR("energy not available");
      }
    }
  }

  // check call backs
  int numberOfCallBacks;
  KIM::CALL_BACK_NAME::get_number_of_call_backs(&numberOfCallBacks);
  for (int i=0; i<numberOfCallBacks; ++i)
  {
    KIM::CallBackName callBackName;
    KIM::CALL_BACK_NAME::get_call_back_name(i, &callBackName);
    KIM::Attribute attribute;
    kim_cluster_model->get_call_back_attribute(callBackName, &attribute);

    std::cout << "CallBack Name \""
              << callBackName.string() << "\""
              << " has attribute \""
              << attribute.string() << "\""
              << std::endl;

    // cannot handle any "required" call backs
    if (attribute == KIM::ATTRIBUTE::required)
    {
      MY_ERROR("unsupported required call back");
    }
  }

  // We're compatible with the model. Let's do it.

  int numberOfParameters;
  kim_cluster_model->get_num_params(&numberOfParameters);
  for (int i=0; i<numberOfParameters; ++i)
  {
    KIM::DataType dataType;
    std::string str;
    int extent;
    kim_cluster_model->get_parameter_data_type_and_description(
        i, &dataType, &str);
    if (dataType == KIM::DATA_TYPE::Integer)
    {
      int const * intPtr;
      kim_cluster_model->get_parameter_extent_and_pointer(i, &extent, &intPtr);
      std::cout << "Parameter No. " << i
                << " has data type \"" << dataType.string() << "\""
                << " with has extent " << extent
                << " and description : " << str << std::endl;
    }
    else
    {
      double const * doublePtr;
      kim_cluster_model->get_parameter_extent_and_pointer(i, &extent,
                                                          &doublePtr);
      std::cout << "Parameter No. " << i
                << " has data type \"" << dataType.string() << "\""
                << " with has extent " << extent
                << " and description : " << str << std::endl;
    }
  }

  error = kim_cluster_model->set_data(
      KIM::ARGUMENT_NAME::numberOfParticles, (int *) &numberOfParticles_cluster)
      || kim_cluster_model->set_data(
          KIM::ARGUMENT_NAME::particleSpecies, particleSpecies_cluster_model)
      || kim_cluster_model->set_data(
          KIM::ARGUMENT_NAME::particleContributing, particleContributing_cluster_model)
      || kim_cluster_model->set_data(
          KIM::ARGUMENT_NAME::coordinates, (double*) coords_cluster)
      || kim_cluster_model->set_data(
          KIM::ARGUMENT_NAME::energy, &energy_cluster_model);
  if (error) MY_ERROR("KIM_API_set_data");
  error = kim_cluster_model->set_call_back(KIM::CALL_BACK_NAME::get_neigh,
                                           KIM::LANGUAGE_NAME::Cpp,
                                           (KIM::func *) &get_cluster_neigh,
                                           &nl_cluster_model);
  if (error) MY_ERROR("set_call_back");

  kim_cluster_model->get_influence_distance(&influence_distance_cluster_model);
  kim_cluster_model->get_cutoffs(&number_of_cutoffs, &cutoff_cluster_model);
  if (number_of_cutoffs != 1) MY_ERROR("too many cutoffs");

  /* setup particleSpecies */
  int isSpeciesSupported;
  error = kim_cluster_model->get_species_support_and_code(
      KIM::SPECIES_NAME::Ar,
      &isSpeciesSupported,
      &(particleSpecies_cluster_model[0]));
  if (error) MY_ERROR("get_species_code");
  for (i = 1; i < NCLUSTERPARTS; ++i)
    particleSpecies_cluster_model[i] = particleSpecies_cluster_model[0];
  /* setup particleContributing */
  for (i = 0; i < NCLUSTERPARTS; ++i)
    particleContributing_cluster_model[i] = 1;  /* every particle contributes */

  /* setup neighbor lists */
  /* allocate memory for list */
  nl_cluster_model.numberOfParticles = NCLUSTERPARTS;
  nl_cluster_model.NNeighbors = new int[NCLUSTERPARTS];
  if (NULL==nl_cluster_model.NNeighbors) MY_ERROR("new unsuccessful");

  nl_cluster_model.neighborList = new int[NCLUSTERPARTS*NCLUSTERPARTS];
  if (NULL==nl_cluster_model.neighborList) MY_ERROR("new unsuccessful");

  /* ready to compute */
  std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(10);
  std::cout << "--------------------------------------------------------------------------------\n";
  std::cout << "This is Test  : ex_test_Ar_fcc_cluster\n";
  std::cout << "MODEL is : " << modelname << std::endl;

  for (CurrentSpacing = MinSpacing; CurrentSpacing < MaxSpacing; CurrentSpacing += SpacingIncr)
  {
    /* update coordinates for cluster */
    create_FCC_cluster(CurrentSpacing, NCELLSPERSIDE, &(coords_cluster[0][0]));
    /* compute neighbor lists */
    fcc_cluster_neighborlist(0, NCLUSTERPARTS, &(coords_cluster[0][0]),
                             (*cutoff_cluster_model + cutpad), &nl_cluster_model);

    /* call compute functions */
    error = kim_cluster_model->compute();
    if (error) MY_ERROR("compute");

    /* print the results */
    std::cout << "Energy for " << NCLUSTERPARTS << " parts = "
              << std::setw(20) << energy_cluster_model
              << std::setw(20) << CurrentSpacing
              << std::endl;
  }


  /* call model destroy */
  KIM::Model::destroy(&kim_cluster_model);

  /* free memory of neighbor lists */
  delete [] nl_cluster_model.NNeighbors;
  delete [] nl_cluster_model.neighborList;

  /* everything is great */
  return 0;
}

void create_FCC_cluster(double FCCspacing, int nCellsPerSide, double *coords)
{
  /* local variables */
  double FCCshifts[4][DIM];
  double latVec[DIM];
  int a;
  int i;
  int j;
  int k;
  int m;
  int n;

  /* create a cubic FCC cluster of parts */
  FCCshifts[0][0] = 0.0;            FCCshifts[0][1] = 0.0;            FCCshifts[0][2] = 0.0;
  FCCshifts[1][0] = 0.5*FCCspacing; FCCshifts[1][1] = 0.5*FCCspacing; FCCshifts[1][2] = 0.0;
  FCCshifts[2][0] = 0.5*FCCspacing; FCCshifts[2][1] = 0.0;            FCCshifts[2][2] = 0.5*FCCspacing;
  FCCshifts[3][0] = 0.0;            FCCshifts[3][1] = 0.5*FCCspacing; FCCshifts[3][2] = 0.5*FCCspacing;

  a = 0;
  for (i = 0; i < nCellsPerSide; ++i)
  {
    latVec[0] = ((double) i)*FCCspacing;
    for (j = 0; j < nCellsPerSide; ++j)
    {
      latVec[1] = ((double) j)*FCCspacing;
      for (k = 0; k < nCellsPerSide; ++k)
      {
        latVec[2] = ((double) k)*FCCspacing;
        for (m = 0; m < 4; ++m)
        {
          for (n = 0; n < DIM; ++n)
          {
            coords[a*DIM + n] = latVec[n] + FCCshifts[m][n];
          }
          a++;
        }
      }
      /* add in the remaining three faces */
      /* pos-x face */
      latVec[0] = NCELLSPERSIDE*FCCspacing;
      latVec[1] = ((double) i)*FCCspacing;
      latVec[2] = ((double) j)*FCCspacing;
      for (n = 0; n < DIM; ++n)
      {
        coords[a*DIM + n] = latVec[n];
      }
      a++;
      for (n = 0; n < DIM; ++n)
      {
        coords[a*DIM + n] = latVec[n] + FCCshifts[3][n];
      }
      a++;
      /* pos-y face */
      latVec[0] = ((double) i)*FCCspacing;
      latVec[1] = NCELLSPERSIDE*FCCspacing;
      latVec[2] = ((double) j)*FCCspacing;
      for (n = 0; n < DIM; ++n)
      {
        coords[a*DIM + n] = latVec[n];
      }
      a++;
      for (n = 0; n < DIM; ++n)
      {
        coords[a*DIM + n] = latVec[n] + FCCshifts[2][n];
      }
      a++;
      /* pos-z face */
      latVec[0] = ((double) i)*FCCspacing;
      latVec[1] = ((double) j)*FCCspacing;
      latVec[2] = NCELLSPERSIDE*FCCspacing;
      for (n = 0; n < DIM; ++n)
      {
        coords[a*DIM + n] = latVec[n];
      }
      a++;
      for (n = 0; n < DIM; ++n)
      {
        coords[a*DIM + n] = latVec[n] + FCCshifts[1][n];
      }
      a++;
    }
    /* add in the remaining three edges */
    latVec[0] = ((double) i)*FCCspacing;
    latVec[1] = NCELLSPERSIDE*FCCspacing;
    latVec[2] = NCELLSPERSIDE*FCCspacing;
    for (n = 0; n < DIM; ++n)
    {
      coords[a*DIM + n] = latVec[n];
    }
    a++;
    latVec[0] = NCELLSPERSIDE*FCCspacing;
    latVec[1] = ((double) i)*FCCspacing;
    latVec[2] = NCELLSPERSIDE*FCCspacing;
    for (n = 0; n < DIM; ++n)
    {
      coords[a*DIM + n] = latVec[n];
    }
    a++;
    latVec[0] = NCELLSPERSIDE*FCCspacing;
    latVec[1] = NCELLSPERSIDE*FCCspacing;
    latVec[2] = ((double) i)*FCCspacing;
    for (n = 0; n < DIM; ++n)
    {
      coords[a*DIM + n] = latVec[n];
    }
    a++;
  }
  /* add in the remaining corner */
  for (n = 0; n < DIM; ++n)
  {
    coords[a*DIM + n] = NCELLSPERSIDE*FCCspacing;
  }
  a++;

  return;
}


void fcc_cluster_neighborlist(int half, int numberOfParticles, double* coords,
                              double cutoff, NeighList* nl)
{
  /* local variables */
  int i;
  int j;
  int k;
  int a;

  double dx[DIM];
  double r2;
  double cutoff2;

  cutoff2 = cutoff*cutoff;

  for (i = 0; i < numberOfParticles; ++i)
  {
    a = 0;
    for (j = 0; j < numberOfParticles; ++j)
    {
      r2 = 0.0;
      for (k = 0; k < DIM; ++k)
      {
        dx[k] = coords[j*DIM + k] - coords[i*DIM + k];
        r2 += dx[k]*dx[k];
      }

      if (r2 < cutoff2)
      {
        if ((half && i < j) || (!half && i != j))
        {
          /* part j is a neighbor of part i */
          (*nl).neighborList[i*NCLUSTERPARTS + a] = j;
          a++;
        }
      }
    }
    /* part i has `a' neighbors */
    (*nl).NNeighbors[i] = a;
  }

  return;
}

int get_cluster_neigh(void const * const dataObject,
                      int const neighborListIndex, int const particleNumber,
                      int * const numberOfNeighbors,
                      int const ** const neighborsOfParticle)
{
  /* local variables */
  int error = true;
  NeighList* nl = (NeighList*) dataObject;
  int numberOfParticles = nl->numberOfParticles;

  if (neighborListIndex != 0) return error;

  /* initialize numNeigh */
  *numberOfNeighbors = 0;

  if ((particleNumber >= numberOfParticles) || (particleNumber < 0)) /* invalid id */
  {
    MY_WARNING("Invalid part ID in get_cluster_neigh");
    return true;
  }

  /* set the returned number of neighbors for the returned part */
  *numberOfNeighbors = (*nl).NNeighbors[particleNumber];

  /* set the location for the returned neighbor list */
  *neighborsOfParticle = &((*nl).neighborList[(particleNumber)*numberOfParticles]);

  return false;
}
