#include <ElectronContaminationStudy.h>
#include <iostream>
#include <string>
#include <TString.h>

using namespace std;
using namespace pPb8TeV2016JetCalibration;

/**
 * Runs an instance of the ElectronContaminationStudy analysis routine.
 * Author: Jeff Ouellette
 * Dated: 8/20/2018
 */
int main (int argc, char** argv) {
  std::cout << argc << " arguments given." << endl;

  if (argc == 1) {
   ElectronContaminationStudy ();
  }
  else if (argc == 2 && argv[1]) {
   const bool isPeriodA = (string(argv[1]) == "true" ? true : false);
   ElectronContaminationStudy (isPeriodA);
  }
  else if (argc == 3 && argv[1] && argv[2]) {
   const bool isPeriodA = (string(argv[1]) == "true" ? true : false);
   const TString fileName = argv[2];
   ElectronContaminationStudy (isPeriodA, fileName);
  }

  return 0;
}
