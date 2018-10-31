#include <JetInsituCorrectionCheck.h>
#include <iostream>
#include <string>

using namespace std;
using namespace JetCalibration;

/**
 * Runs an instance of the JetInsituCorrectionCheck analysis routine.
 * Author: Jeff Ouellette
 * Dated: 8/20/2018
 */
int main (int argc, char** argv) {
  if (argc < 3 || !argv[1] || !argv[2]) { // needs at least 2 argument (+ executable)
   cout << "Not enough arguments specified. Quitting." << endl;
   return 0;
  }

  int dataSet;
  bool isPeriodA;
  char* directory,* inFileName;

  int argn = argc-1;
  switch (argc) {
   case 5:
    if (argv[argn]) inFileName = argv[argn--];
   case 4:
    if (argv[argn]) isPeriodA = (string (argv[argn--]) == "true");
   case 3:
    if (argv[argn]) dataSet = atoi (argv[argn--]);
    if (argv[argn]) directory = argv[argn--];
  }

  switch (argc) {
   case 3:
    JetInsituCorrectionCheck (directory, dataSet);
    break;
   case 4:
    JetInsituCorrectionCheck (directory, dataSet, isPeriodA);
    break;
   case 5:
    JetInsituCorrectionCheck (directory, dataSet, isPeriodA, inFileName);
    break;
  }

  return 0;
}
