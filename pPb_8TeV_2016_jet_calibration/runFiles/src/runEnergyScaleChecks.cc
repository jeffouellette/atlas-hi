#include <EnergyScaleChecks.h>
#include <iostream>
#include <string>
#include <TString.h>

using namespace std;
using namespace pPb8TeV2016JetCalibration;

/**
 * Runs an instance of the EnergyScaleChecks analysis routine.
 * Author: Jeff Ouellette
 * Dated: 8/20/2018
 */
int main (int argc, char** argv) {

  int dataSet, nevt;
  bool isPeriodA;
  TString fileName;
  double xs, feff;

  switch (argc) {
   case 7:
    if (argv[6]) nevt = atoi (argv[6]);
   case 6:
    if (argv[5]) feff = atof (argv[5]);
   case 5:
    if (argv[4]) xs = atof (argv[4]);
   case 4:
    if (argv[3]) fileName = TString (argv[3]);
   case 3:
    if (argv[2]) isPeriodA = (string(argv[2]) == "true" ? true : false);
   case 2:
    if (argv[1]) dataSet = atoi (argv[1]);
  }

  switch (argc) {
   case 2:
    EnergyScaleChecks (dataSet);
    break;
   case 3:
    EnergyScaleChecks (dataSet, isPeriodA);
    break;
   case 4:
    EnergyScaleChecks (dataSet, isPeriodA, fileName);
    break;
   case 5:
    EnergyScaleChecks (dataSet, isPeriodA, fileName, xs);
    break;
   case 6:
    EnergyScaleChecks (dataSet, isPeriodA, fileName, xs, feff);
    break;
   case 7:
    EnergyScaleChecks (dataSet, isPeriodA, fileName, xs, feff, nevt);
    break;
  }

  return 0;
}
