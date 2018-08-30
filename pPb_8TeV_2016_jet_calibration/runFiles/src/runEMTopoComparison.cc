#include <EMTopoComparison.h>
#include <iostream>
#include <string>
#include <TString.h>

using namespace std;
using namespace pPb8TeV2016JetCalibration;

/**
 * Runs an instance of the EMTopoComparison analysis routine.
 * Author: Jeff Ouellette
 * Dated: 8/20/2018
 */
int main (int argc, char** argv) {
  if (!argv[0] || !argv[1] || !argv[2] || !argv[3] || !argv[4]) {
   cout << "Not enough arguments specified. Quitting." << endl;
   return 0;
  }

  const int dataSet = atoi(argv[1]);
  const double luminosity = atof(argv[2]);
  const bool isMC = (string(argv[3]) == "true" ? true : false);
  const bool isPeriodA = (string(argv[4]) == "true" ? true : false);

  if (argc == 5) { // if file name is not given
   EMTopoComparison (dataSet, luminosity, isMC, isPeriodA);
  }
  else if (argc == 6 && argv[5]) { // if file name is given
   const TString fileName = argv[5];
   EMTopoComparison (dataSet, luminosity, isMC, isPeriodA, fileName);
  }
  return 0;
}
