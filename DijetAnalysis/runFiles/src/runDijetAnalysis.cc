#include <DijetAnalysis.h>

using namespace JetAnalysis;

int main (int argc, char** argv) {
  //if (argc < 3 || !argv[1] || !argv[2]) // needs 3 arguments (run ID & luminosity, plus executable)
 
  int runNumber = atoi (argv[1]);
  double lumi = (double)atof (argv[2]);
  bool isMC = (string (argv[3]) == "true");
  bool isPeriodA = (string (argv[4]) == "true");

  DijetAnalysis (runNumber, lumi, isMC, isPeriodA);

  return 0;
}
