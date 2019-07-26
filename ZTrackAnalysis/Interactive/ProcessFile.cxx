#include "MCAnalysis.h"
#include "MinbiasAnalysis.h"

int main (int argc, char** argv) {

  if (argc < 4) {
    cout << "Insufficient arguments! Exiting." << endl;
    return 1;
  }

  string algo = argv[1];
  string inFileName = argv[2];
  string outFileName = argv[3];

  if (algo == "mc") {
    MCAnalysis* mc = new MCAnalysis ();
    mc->Execute (inFileName.c_str (), outFileName.c_str ());
    delete mc;
  }
  else if (algo == "minbias") {
    MinbiasAnalysis* bkg = new MinbiasAnalysis (); 
    bkg->Execute (inFileName.c_str (), outFileName.c_str ());
    delete bkg;
  }

  return 0;
}
