#include "MCAnalysis.h"
#include "TruthAnalysis.h"
#include "MinbiasAnalysis.h"

int main (int argc, char** argv) {

  if (argc < 7) {
    cout << "Insufficient arguments! Exiting." << endl;
    return 1;
  }

  string algo = argv[1];
  string inFileName = argv[2];
  string outFileName = argv[3];
  bool doHITightVar = (string (argv[4]) == "true");
  bool doTrkPurVar = (string (argv[5]) == "true");
  bool doRunVar = (string (argv[6]) == "true");

  string subDir = "";
  if (!doHITightVar && !doTrkPurVar && !doRunVar) {
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
  else if (doRunVar) {
    subDir = "Variations/RunVariation";
  }

  if (algo == "mc") {
    MCAnalysis* mc = nullptr;
    if (doHITightVar)
      mc = new MCAnalysis ("mc_trackHITightVar", subDir.c_str ());
    else if (doTrkPurVar)
      mc = new MCAnalysis ("mc_trackPurityVar", subDir.c_str ());
    else
      mc = new MCAnalysis ("mc", subDir.c_str ());
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
      bkg = new MinbiasAnalysis ("bkg_trackHITightVar", subDir.c_str ());
    else if (doTrkPurVar)
      bkg = new MinbiasAnalysis ("bkg_trackPurityVar", subDir.c_str ());
    else if (doRunVar)
      bkg = new MinbiasAnalysis ("bkg_runVar", subDir.c_str ()); 
    else
      bkg = new MinbiasAnalysis ("minbias", subDir.c_str ()); 
    bkg->useHITight = doHITightVar;
    bkg->doTrackPurVar = doTrkPurVar;
    bkg->Execute (inFileName.c_str (), outFileName.c_str ());
    delete bkg;
  }

  return 0;
}
