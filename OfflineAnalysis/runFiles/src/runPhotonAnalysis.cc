#include <PhotonAnalysis.h>
#include <iostream>
#include <string>

using namespace std;
using namespace offlineAnalyses;

int main (int argc, char** argv) {
  int dataSet = 0;
  bool isMC = false;

  if (argc >= 3 && argv[1] && argv[2]) {
   dataSet = atoi (argv[1]);
   isMC = (string (argv[2]) == "true");
  }

  if (argc == 4 && argv[3])
   PhotonAnalysis (dataSet, isMC, argv[3]);
  else if (argc == 3)
   PhotonAnalysis (dataSet, isMC);

  return 0;
}
