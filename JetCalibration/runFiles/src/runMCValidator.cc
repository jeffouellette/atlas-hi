#include <MCValidator.h>
#include <iostream>
#include <string>

using namespace std;
using namespace JetCalibration;

/**
 * Runs an instance of the MCValidator analysis routine.
 * Author: Jeff Ouellette
 * Dated: 10/15/2018
 */
int main (int argc, char** argv) {
  if (argc < 7 || !argv[1] || !argv[2] || !argv[3] || !argv[4] || !argv[5] || !argv[6]) { // needs at least 6 arguments (+ executable)
   cout << "Not enough arguments specified. Quitting." << endl;
   return 0;
  }

  int dataSet, nevt;
  bool isPeriodA;
  char *directory, *fileName;
  double pt_lo, pt_hi, xs, feff;

  int argn = argc-1;
  switch (argc) {
   case 10:
    if (argv[argn]) nevt = atoi (argv[argn--]);
   case 9:
    if (argv[argn]) feff = atof (argv[argn--]);
   case 8:
    if (argv[argn]) xs = atof (argv[argn--]);
   case 7:
    if (argv[argn]) fileName = argv[argn--];
    if (argv[argn]) pt_hi = atof (argv[argn--]);
    if (argv[argn]) pt_lo = atof (argv[argn--]);
    if (argv[argn]) isPeriodA = (string (argv[argn--]) == "true" ? true : false);
    if (argv[argn]) dataSet = atoi (argv[argn--]);
    if (argv[argn]) directory = argv[argn--];
  }

  switch (argc) {
   case 7:
    MCValidator (directory, dataSet, isPeriodA, pt_lo, pt_hi, fileName);
    break;
   case 8:
    MCValidator (directory, dataSet, isPeriodA, pt_lo, pt_hi, fileName, xs);
    break;
   case 9:
    MCValidator (directory, dataSet, isPeriodA, pt_lo, pt_hi, fileName, xs, feff);
    break;
   case 10:
    MCValidator (directory, dataSet, isPeriodA, pt_lo, pt_hi, fileName, xs, feff, nevt);
    break;
  }

  return 0;
}
