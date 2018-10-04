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
  if (argc < 3 || !argv[1] || !argv[2]) { // needs at least 2 argument (+ executable)
   cout << "Not enough arguments specified. Quitting." << endl;
   return 0;
  }

  int dataSet, nevt;
  bool isMC, isPeriodA;
  char* directory,* fileName;
  double lumi, xs, feff;

  int argn = argc-1;
  switch (argc) {
   case 10:
    if (argv[argn]) nevt = atoi (argv[argn]);
    argn--;
   case 9:
    if (argv[argn]) feff = atof (argv[argn]);
    argn--;
   case 8:
    if (argv[argn]) xs = atof (argv[argn]);
    argn--;
   case 7:
    if (argv[argn]) fileName = argv[argn];
    argn--;
   case 6:
    if (argv[argn]) isPeriodA = (string(argv[argn]) == "true" ? true : false);
    argn--;
   case 5:
    if (argv[argn]) isMC = (string(argv[argn]) == "true" ? true : false);
    argn--;
   case 4:
    if (argv[argn]) lumi = atof (argv[argn]);
    argn--;
   case 3:
    if (argv[argn]) dataSet = atoi (argv[argn]);
    argn--;
    if (argv[argn]) directory = argv[argn];
    argn--;
  }

  switch (argc) {
   case 3:
    EMTopoComparison (directory, dataSet);
    break;
   case 4:
    EMTopoComparison (directory, dataSet, lumi);
    break;
   case 5:
    EMTopoComparison (directory, dataSet, lumi, isMC);
    break;
   case 6:
    EMTopoComparison (directory, dataSet, lumi, isMC, isPeriodA);
    break;
   case 7:
    EMTopoComparison (directory, dataSet, lumi, isMC, isPeriodA, fileName);
    break;
   case 8:
    EMTopoComparison (directory, dataSet, lumi, isMC, isPeriodA, fileName, xs);
    break;
   case 9:
    EMTopoComparison (directory, dataSet, lumi, isMC, isPeriodA, fileName, xs, feff);
    break;
   case 10:
    EMTopoComparison (directory, dataSet, lumi, isMC, isPeriodA, fileName, xs, feff, nevt);
    break;
  }

  return 0;
}
