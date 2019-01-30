#include <RtrkTreeMaker.h>

#include <string>

using namespace std;

using namespace JetTrackAnalysis;

int main (int argc, char** argv) {

  int argn = 1;
  const char* subdir = argv[argn++];
  const int dataSet = atoi (argv[argn++]);
  const bool isMC = string (argv[argn++]) == "true";
  const bool isPeriodA = string (argv[argn++]) == "true";
  const char* inFileName = (argc > argn && argv[argn] ? argv[argn++] : "");
  const double xs = (argc > argn && argv[argn] ? atof (argv[argn++]) : 0.);
  const double feff = (argc > argn && argv[argn] ? atof (argv[argn++]) : 0.);
  const int nevt = (argc > argn && argv[argn] ? atoi (argv[argn++]) : 0);

  RtrkTreeMaker (subdir, dataSet, isMC, isPeriodA, inFileName, xs, feff, nevt);

  return 0;
}
