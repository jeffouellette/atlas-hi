#include "MCAnalysis.h"
#include "TruthAnalysis.h"
#include "MinbiasAnalysis.h"

int main (int argc, char** argv) {

  if (argc < 6) {
    cout << "Insufficient arguments! Exiting." << endl;
    return 1;
  }

  string algo = argv[1];
  string inFileName = argv[2];
  string outFileName = argv[3];
  bool doHITightVar = (string (argv[4]) == "true");
  bool doTrkPurVar = (string (argv[5]) == "true");

  if (!doHITightVar && !doTrkPurVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Nominal/" + outFileName;
  }
  else if (doHITightVar) {
    inFileName = "Variations/TrackHITightWPVariation/" + inFileName;
    outFileName = "Variations/TrackHITightWPVariation/" + outFileName;
  }
  else if (doTrkPurVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Variations/TrackPurityVariation/" + outFileName;
  }

  if (algo == "mc") {
    MCAnalysis* mc = nullptr;
    if (doHITightVar)
      mc = new MCAnalysis ("mc_trackHITightVar", "");
    else if (doTrkPurVar)
      mc = new MCAnalysis ("mc_trackPurityVar", "");
    else
      mc = new MCAnalysis ("mc", "");
    mc->useHITight = doHITightVar;
    mc->doTrackPurVar = doTrkPurVar;
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
    if (doHITightVar)
      bkg = new MinbiasAnalysis ("bkg_trackHITightVar", "");
    else if (doTrkPurVar)
      bkg = new MinbiasAnalysis ("bkg_trackPurityVar", "");
    else
      bkg = new MinbiasAnalysis ("minbias", ""); 
    bkg->useHITight = doHITightVar;
    bkg->doTrackPurVar = doTrkPurVar;
    bkg->Execute (inFileName.c_str (), outFileName.c_str ());
    delete bkg;
  }

  return 0;
}
