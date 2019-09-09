#include "MCAnalysis.h"
#include "TruthAnalysis.h"
#include "MinbiasAnalysis.h"

int main (int argc, char** argv) {

  if (argc < 12) {
    cout << "Insufficient arguments! Exiting." << endl;
    return 1;
  }

  string algo = argv[1];
  string inFileName = argv[2];
  string outFileName = argv[3];

  bool doElectronPtUpVar = (string (argv[4]) == "true");
  bool doElectronPtDownVar = (string (argv[5]) == "true");
  bool doMuonPtUpVar = (string (argv[6]) == "true");
  bool doMuonPtDownVar = (string (argv[7]) == "true");
  bool doElectronLHMediumVar = (string (argv[8]) == "true");
  bool doMuonTightVar = (string (argv[9]) == "true");
  bool doHITightVar = (string (argv[10]) == "true");
  bool doTrkPurVar = (string (argv[11]) == "true");


  if (!doElectronPtUpVar && !doElectronPtDownVar && !doMuonPtUpVar && !doMuonPtDownVar && !doElectronLHMediumVar && !doMuonTightVar && !doHITightVar && !doTrkPurVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Nominal/" + outFileName;
  }
  else if (doElectronPtUpVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Variations/ElectronPtUpVariation/" + outFileName;
  }
  else if (doElectronPtDownVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Variations/ElectronPtDownVariation/" + outFileName;
  }
  else if (doMuonPtUpVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Variations/MuonPtUpVariation/" + outFileName;
  }
  else if (doMuonPtDownVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Variations/MuonPtDownVariation/" + outFileName;
  }
  else if (doElectronLHMediumVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Variations/ElectronLHMediumVariation/" + outFileName;
  }
  else if (doMuonTightVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Variations/MuonTightVariation/" + outFileName;
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
    inFileName = "MCAnalysis/" + inFileName;
    outFileName = "MCAnalysis/" + outFileName;
    if (doHITightVar)
      mc = new MCAnalysis ("mc_trackHITightVar");
    else if (doTrkPurVar)
      mc = new MCAnalysis ("mc_trackPurityVar");
    else
      mc = new MCAnalysis ("mc");
    mc->useHITight = doHITightVar;
    mc->doTrackPurVar = doTrkPurVar;
    mc->Execute (inFileName.c_str (), outFileName.c_str ());
    delete mc;
  }

  else if (algo == "truth") {
    TruthAnalysis* truth = nullptr;
    inFileName = "TruthAnalysis/" + inFileName;
    outFileName = "TruthAnalysis/" + outFileName;
    truth = new TruthAnalysis ("truth");
    truth->Execute (inFileName.c_str (), outFileName.c_str ());
    delete truth;
  }

  else if (algo == "minbias") {
    MinbiasAnalysis* bkg = nullptr;
    inFileName = "MinbiasAnalysis/" + inFileName;
    outFileName = "MinbiasAnalysis/" + outFileName;
    if (doHITightVar) {
      bkg = new MinbiasAnalysis ("bkg_trackHITightVar");
      bkg->eventWeightsExt = "minbias";
    }
    else if (doTrkPurVar) {
      bkg = new MinbiasAnalysis ("bkg_trackPurityVar");
      bkg->eventWeightsExt = "minbias";
    }
    else {
      bkg = new MinbiasAnalysis ("bkg");
      bkg->eventWeightsExt = "minbias";
    }
    bkg->useHITight = doHITightVar;
    bkg->doTrackPurVar = doTrkPurVar;
    bkg->Execute (inFileName.c_str (), outFileName.c_str ());
    delete bkg;
  }


  return 0;
}
