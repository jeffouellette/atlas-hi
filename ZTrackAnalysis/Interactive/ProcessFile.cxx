#include "MCAnalysis.h"
#include "TruthAnalysis.h"
#include "MinbiasAnalysis.h"

int main (int argc, char** argv) {

  if (argc < 13) {
    cout << "Insufficient arguments! Exiting." << endl;
    return 1;
  }

  string algo = argv[1];
  string inFileName = argv[2];
  string outFileName = argv[3];

  bool use2015conds = (string(argv[4]) == "true");
  bool doElectronPtUpVar = (string (argv[5]) == "true");
  bool doElectronPtDownVar = (string (argv[6]) == "true");
  bool doMuonPtUpVar = (string (argv[7]) == "true");
  bool doMuonPtDownVar = (string (argv[8]) == "true");
  bool doHITightVar = (string (argv[9]) == "true");
  bool doTrkEffVar = (string (argv[10]) == "true");
  bool doTrkPurUpVar = (string (argv[11]) == "true");
  bool doTrkPurDownVar = (string (argv[12]) == "true");


  if (!doElectronPtUpVar && !doElectronPtDownVar && !doMuonPtUpVar && !doMuonPtDownVar && !doHITightVar && !doTrkEffVar && !doTrkPurUpVar && !doTrkPurDownVar) {
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
  else if (doHITightVar) {
    inFileName = "Variations/TrackHITightWPVariation/" + inFileName;
    outFileName = "Variations/TrackHITightWPVariation/" + outFileName;
  }
  else if (doTrkEffVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Variations/TrackEffPionsVariation/" + outFileName;
  }
  else if (doTrkPurUpVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Variations/TrackPurityUpVariation/" + outFileName;
  }
  else if (doTrkPurDownVar) {
    inFileName = "Nominal/" + inFileName;
    outFileName = "Variations/TrackPurityDownVariation/" + outFileName;
  }


  if (algo == "mc") {
    MCAnalysis* mc = nullptr;
    inFileName = "MCAnalysis/" + inFileName;
    outFileName = "MCAnalysis/" + outFileName;
    if (doHITightVar)
      mc = new MCAnalysis ("mc_trackHITightVar");
    else if (doTrkEffVar)
      mc = new MCAnalysis ("mc_trackEffVar");
    else if (doTrkPurUpVar)
      mc = new MCAnalysis ("mc_trkPurUpVar");
    else if (doTrkPurDownVar)
      mc = new MCAnalysis ("mc_trkPurDownVar");
    else
      mc = new MCAnalysis ("mc");
    mc->useHITight = doHITightVar;
    mc->doTrackEffVar = doTrkEffVar;
    mc->doTrackPurVar = (doTrkPurUpVar || doTrkPurDownVar);
    mc->trkPurNSigma = (doTrkPurUpVar ? 1. : (doTrkPurDownVar ? -1 : 0));
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
    if (doHITightVar)
      bkg = new MinbiasAnalysis ("bkg_trackHITightVar");
    else if (doTrkEffVar)
      bkg = new MinbiasAnalysis ("bkg_trackEffVar");
    else if (doTrkPurUpVar)
      bkg = new MinbiasAnalysis ("bkg_trkPurUpVar");
    else if (doTrkPurDownVar)
      bkg = new MinbiasAnalysis ("bkg_trkPurDownVar");
    else
      bkg = new MinbiasAnalysis ("bkg");
    bkg->is2015Conds = use2015conds;
    bkg->useHijingEffs = use2015conds;
    bkg->useHITight = doHITightVar;
    bkg->doTrackEffVar = doTrkEffVar;
    bkg->doTrackPurVar = (doTrkPurUpVar || doTrkPurDownVar);
    bkg->trkPurNSigma = (doTrkPurUpVar ? 1. : (doTrkPurDownVar ? -1 : 0));
    bkg->Execute (inFileName.c_str (), outFileName.c_str ());
    delete bkg;
  }


  return 0;
}
