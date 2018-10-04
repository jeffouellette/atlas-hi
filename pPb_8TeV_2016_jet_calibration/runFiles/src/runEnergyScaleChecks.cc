#include <EnergyScaleChecks.h>
#include <iostream>
#include <string>

using namespace std;
using namespace pPb8TeV2016JetCalibration;

/**
 * Runs an instance of the EnergyScaleChecks analysis routine.
 * Author: Jeff Ouellette
 * Dated: 8/20/2018
 */
int main (int argc, char** argv) {
  if (argc < 3 || !argv[1] || !argv[2]) { // needs at least 2 argument (+ executable)
   cout << "Not enough arguments specified. Quitting." << endl;
   return 0;
  }

  int dataSet, nevt;
  bool isPeriodA;
  char* directory,* fileName;
  double xs, feff;

  int argn = argc-1;
  switch (argc) {
   case 8:
    if (argv[argn]) nevt = atoi (argv[argn--]);
   case 7:
    if (argv[argn]) feff = atof (argv[argn--]);
   case 6:
    if (argv[argn]) xs = atof (argv[argn--]);
   case 5:
    if (argv[argn]) fileName = argv[argn--];
   case 4:
    if (argv[argn]) isPeriodA = (string(argv[argn--]) == "true" ? true : false);
   case 3:
    if (argv[argn]) dataSet = atoi (argv[argn--]);
   case 2:
    if (argv[argn]) directory = argv[argn--];
  }

  switch (argc) {
   case 3:
    EnergyScaleChecks (directory, dataSet);
    break;
   case 4:
    EnergyScaleChecks (directory, dataSet, isPeriodA);
    break;
   case 5:
    EnergyScaleChecks (directory, dataSet, isPeriodA, fileName);
    break;
   case 6:
    EnergyScaleChecks (directory, dataSet, isPeriodA, fileName, xs);
    break;
   case 7:
    EnergyScaleChecks (directory, dataSet, isPeriodA, fileName, xs, feff);
    break;
   case 8:
    EnergyScaleChecks (directory, dataSet, isPeriodA, fileName, xs, feff, nevt);
    break;
  }

  return 0;
}
