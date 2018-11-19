#include <PhotonAnalysis.h>
#include <iostream>
#include <string>

using namespace std;
using namespace offlineAnalyses;

int main (int argc, char** argv) {
  int dataSet = 0;
  bool isMC = false;
  float lumi = 0;

  if (argc >= 4 && argv[1] && argv[2] && argv[3]) {
   dataSet = atoi (argv[1]);
   isMC = (string (argv[2]) == "true");
   lumi = atof (argv[3]);
  }

  if (argc == 5 && argv[4])
   PhotonAnalysis (dataSet, isMC, lumi, argv[4]);
  else if (argc == 4)
   PhotonAnalysis (dataSet, isMC, lumi);

  return 0;
}
