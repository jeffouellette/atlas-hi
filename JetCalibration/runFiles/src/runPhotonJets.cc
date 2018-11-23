#include <PhotonJets.h>
#include <iostream>
#include <string>

using namespace std;
using namespace JetCalibration;

/**
 * Runs an instance of the PhotonJets analysis routine.
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
    if (argv[argn]) nevt = atoi (argv[argn--]);
   case 9:
    if (argv[argn]) feff = atof (argv[argn--]);
   case 8:
    if (argv[argn]) xs = atof (argv[argn--]);
   case 7:
    if (argv[argn]) fileName = argv[argn--];
   case 6:
    if (argv[argn]) isPeriodA = string(argv[argn--]) == "true";
   case 5:
    if (argv[argn]) isMC = string(argv[argn--]) == "true";
   case 4:
    if (argv[argn]) lumi = atof (argv[argn--]);
   case 3:
    if (argv[argn]) dataSet = atoi (argv[argn--]);
    if (argv[argn]) directory = argv[argn--];
  }

  switch (argc) {
   case 3:
    PhotonJets (directory, dataSet);
    break;
   case 4:
    PhotonJets (directory, dataSet, lumi);
    break;
   case 5:
    PhotonJets (directory, dataSet, lumi, isMC);
    break;
   case 6:
    PhotonJets (directory, dataSet, lumi, isMC, isPeriodA);
    break;
   case 7:
    PhotonJets (directory, dataSet, lumi, isMC, isPeriodA, fileName);
    break;
   case 8:
    PhotonJets (directory, dataSet, lumi, isMC, isPeriodA, fileName, xs);
    break;
   case 9:
    PhotonJets (directory, dataSet, lumi, isMC, isPeriodA, fileName, xs, feff);
    break;
   case 10:
    PhotonJets (directory, dataSet, lumi, isMC, isPeriodA, fileName, xs, feff, nevt);
    break;
  }

  return 0;
}
