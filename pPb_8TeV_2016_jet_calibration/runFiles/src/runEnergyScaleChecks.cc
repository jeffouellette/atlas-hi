#include <EnergyScaleChecks.h>
#include <iostream>
#include <string>
#include <TString.h>

using namespace pPb8TeV2016JetCalibration;

/**
 * Runs an instance of the EnergyScaleChecks analysis routine.
 * Author: Jeff Ouellette
 * Dated: 8/20/2018
 */
int main (int argc, char** argv) {

  if (argc == 2) { // if just the dataset is specified
   const int dataSet = atoi (argv[1]);
   EnergyScaleChecks (dataSet);
  }
  else if (argc == 3 && argv[1] && argv[2]) { // if period is also specified
   const int dataSet = atoi (argv[1]);
   const bool isPeriodA = (string(argv[2]) == "true" ? true : false);
   EnergyScaleChecks (dataSet, isPeriodA);
  }
  else if (argc == 4 && argv[1] && argv[2] && argv[3]) { // if filename is also specified
   const int dataSet = atoi (argv[1]);
   const bool isPeriodA = (string(argv[2]) == "true" ? true : false);
   const TString fileName = argv[3];
   EnergyScaleChecks (dataSet, isPeriodA, fileName);
  }

  return 0;
}
