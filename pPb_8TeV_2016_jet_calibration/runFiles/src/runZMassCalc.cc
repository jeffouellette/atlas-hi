#include <ZMassCalc.h>
#include <iostream>
#include <string>
#include <TString.h>

using namespace std;
using namespace pPb8TeV2016JetCalibration;

/**
 * Runs an instance of the ZMassCalc analysis routine.
 * Author: Jeff Ouellette
 * Dated: 9/10/2018
 */
int main (int argc, char** argv) {
  if (argc < 2 || !argv[1]) { // needs at least 1 argument (+ executable)
   cout << "Not enough arguments specified. Quitting." << endl;
   return 0;
  }

  if (argc == 2 && argv[1]) {
   const int dataSet = atoi (argv[1]);
   ZMassCalc (dataSet);
  }
  else if (argc == 3 && argv[1] && argv[2]) {
   const int dataSet = atoi (argv[1]);
   const double luminosity = atof (argv[2]);
   ZMassCalc (dataSet, luminosity);
  }
  else if (argc == 4 && argv[1] && argv[2] && argv[3]) {
   const int dataSet = atoi (argv[1]);
   const double luminosity = atof (argv[2]);
   const bool isMC = (string (argv[3]) == "true" ? true : false);
   ZMassCalc (dataSet, luminosity, isMC);
  }
  else if (argc == 5 && argv[1] && argv[2] && argv[3] && argv[4]) {
   const int dataSet = atoi (argv[1]);
   const double luminosity = atof (argv[2]);
   const bool isMC = (string (argv[3]) == "true" ? true : false);
   const bool isPeriodA = (string (argv[4]) == "true" ? true : false);
   ZMassCalc (dataSet, luminosity, isMC, isPeriodA);
  }
  else if (argc == 6 && argv[1] && argv[2] && argv[3] && argv[4] && argv[5]) {
   const int dataSet = atoi (argv[1]);
   const double luminosity = atof (argv[2]);
   const bool isMC = (string (argv[3]) == "true" ? true : false);
   const bool isPeriodA = (string (argv[4]) == "true" ? true : false);
   const TString fileName = argv[5];
   ZMassCalc (dataSet, luminosity, isMC, isPeriodA, fileName);
  }

  return 0;
}
