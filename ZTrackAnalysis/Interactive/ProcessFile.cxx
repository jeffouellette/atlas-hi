#include "MCAnalysis.h"
#include "TruthAnalysis.h"
#include "MinbiasAnalysis.h"

int main (int argc, char** argv) {

  if (argc < 5) {
    cout << "Insufficient arguments! Exiting." << endl;
    return 1;
  }

  string algo = argv[1];
  string inFileName = argv[2];
  string outFileName = argv[3];
  bool doHITightVar = (string (argv[4]) == "true");

  if (algo == "mc") {
    MCAnalysis* mc = nullptr;
    if (!doHITightVar)
      mc = new MCAnalysis ();
    else
      mc = new MCAnalysis ("mc_trackHITightVar", "Variations/TrackHITightWPVariation", true);
    mc->Execute (inFileName.c_str (), outFileName.c_str ());
    delete mc;
  }
  else if (algo == "truth") {
    TruthAnalysis* truth = nullptr;
    truth = new TruthAnalysis ();
    truth->Execute (inFileName.c_str (), outFileName.c_str ());
    delete truth;
  }
  else if (algo == "minbias") {
    MinbiasAnalysis* bkg = nullptr;
    if (!doHITightVar)
      bkg = new MinbiasAnalysis (); 
    else
      bkg = new MinbiasAnalysis ("bkg_trackHITightVar", "Variations/TrackHITightWPVariation", true);
    bkg->Execute (inFileName.c_str (), outFileName.c_str ());
    delete bkg;
  }

  return 0;
}
