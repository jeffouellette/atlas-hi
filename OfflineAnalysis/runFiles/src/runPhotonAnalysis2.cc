#include <PhotonAnalysis2.h>
#include <iostream>
#include <string>

using namespace std;
using namespace offlineAnalyses;

int main (int argc, char** argv) {
  int dataSet;
  bool isMC;

  if (argc == 3 && argv[1] && argv[2]) {
   dataSet = atoi (argv[1]);
   isMC = (string (argv[2]) == "true");
   PhotonAnalysis2 (dataSet, isMC);
  }

  return 0;
}
