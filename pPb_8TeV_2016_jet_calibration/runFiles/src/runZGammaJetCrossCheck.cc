#include <ZGammaJetCrossCheck.h>
#include <iostream>
#include <string>
#include <TString.h>

using namespace std;
using namespace pPb8TeV2016JetCalibration;

/**
 * Runs an instance of the ZGammaJetCrossCheck analysis routine.
 * Author: Jeff Ouellette
 * Dated: 8/20/2018
 */
int main (int argc, char** argv) {
  if (argc < 2 || !argv[1]) { // needs at least 1 argument (+ executable)
   cout << "Not enough arguments specified. Quitting." << endl;
   return 0;
  }

  int dataSet, nevt;
  bool isMC, isPeriodA;
  TString fileName;
  double lumi, xs, feff;

  switch (argc) {
   case 9:
    if (argv[8]) nevt = atoi (argv[8]);
   case 8:
    if (argv[7]) feff = atof (argv[7]);
   case 7:
    if (argv[6]) xs = atof (argv[6]);
   case 6:
    if (argv[5]) fileName = TString (argv[5]);
   case 5:
    if (argv[4]) isPeriodA = (string(argv[4]) == "true" ? true : false);
   case 4:
    if (argv[3]) isMC = (string(argv[3]) == "true" ? true : false);
   case 3:
    if (argv[2]) lumi = atof (argv[2]);
   case 2:
    if (argv[1]) dataSet = atoi (argv[1]);
  }

  switch (argc) {
   case 2:
    ZGammaJetCrossCheck (dataSet);
    break;
   case 3:
    ZGammaJetCrossCheck (dataSet, lumi);
    break;
   case 4:
    ZGammaJetCrossCheck (dataSet, lumi, isMC);
    break;
   case 5:
    ZGammaJetCrossCheck (dataSet, lumi, isMC, isPeriodA);
    break;
   case 6:
    ZGammaJetCrossCheck (dataSet, lumi, isMC, isPeriodA, fileName);
    break;
   case 7:
    ZGammaJetCrossCheck (dataSet, lumi, isMC, isPeriodA, fileName, xs);
    break;
   case 8:
    ZGammaJetCrossCheck (dataSet, lumi, isMC, isPeriodA, fileName, xs, feff);
    break;
   case 9:
    ZGammaJetCrossCheck (dataSet, lumi, isMC, isPeriodA, fileName, xs, feff, nevt);
    break;
  }

  return 0;
}
