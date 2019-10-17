#include "MCAnalysis.h"
#include "MCClosureAnalysis.h"
#include "MixedMCAnalysis.h"
#include "TruthAnalysis.h"
#include "MinbiasAnalysis.h"

int main (int argc, char** argv) {

  if (argc < 7) {
    cout << "Insufficient arguments! Exiting." << endl;
    return 1;
  }

  string algo = argv[1];

  string inFileName = argv[2];
  string mbInFileName = (string (argv[3]) == "0" ? inFileName : argv[3]); // defaults to 1st file name if "0"
  string outFileName = argv[4];

  string mixdir = "DataAnalysis/";
  if (algo == "mcminbias") {
    mixdir = "MCAnalysis/"; 
    //outFileName = "Background/" + outFileName;
  }

  bool use2015conds = (string(argv[5]) == "true");
  bool doHITightVar = (string (argv[6]) == "true");


  if (!doHITightVar) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Nominal/" + outFileName;
  }
  else if (doHITightVar) {
    inFileName = "Variations/TrackHITightWPVariation/" + inFileName;
    mbInFileName = "Variations/TrackHITightWPVariation/" + mbInFileName;
    outFileName = "Variations/TrackHITightWPVariation/" + outFileName;
  }


  //if (algo == "mc" || algo == "mcminbias") {
  if (algo == "mc") {
    MCAnalysis* mc = nullptr;
    inFileName = "MCAnalysis/" + inFileName;
    outFileName = "MCAnalysis/" + outFileName;
    if (doHITightVar)
      mc = new MCAnalysis ("mc_trackHITightVar");
    else if (algo == "mcminbias")
      mc = new MCAnalysis ("mc_bkg");
    else
      mc = new MCAnalysis ("mc");

    if (algo == "mcminbias")
      mc->takeNonTruthTracks = true;

    mc->is2015Conds = use2015conds;
    mc->useHijingEffs = use2015conds;
    mc->useHITight = doHITightVar;

    mc->Execute (inFileName.c_str (), outFileName.c_str ());
    delete mc;
  }

  else if (algo == "mcclosure") {
    MCClosureAnalysis* mc = nullptr;
    inFileName = "MCAnalysis/" + inFileName;
    mbInFileName = "MinbiasAnalysis/" + mbInFileName;
    outFileName = "MCAnalysis/" + outFileName;
    mc = new MCClosureAnalysis ("mc_closure");

    mc->is2015Conds = use2015conds;
    mc->useHijingEffs = use2015conds;
    mc->useHITight = doHITightVar;

    mc->Execute (inFileName.c_str (), mbInFileName.c_str (), outFileName.c_str ());
    delete mc;
  }

  else if (algo == "mixmc") {
    MixedMCAnalysis* mc = nullptr;
    inFileName = "MCAnalysis/" + inFileName;
    outFileName = "MixedMCAnalysis/" + outFileName;
    mc = new MixedMCAnalysis ("mixed_mc");

    mc->is2015Conds = use2015conds;
    mc->useHijingEffs = use2015conds;
    mc->useHITight = doHITightVar;

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

  //else if (algo == "minbias") {
  else if (algo == "minbias" || algo == "mcminbias") {
    MinbiasAnalysis* bkg = nullptr;

    inFileName = mixdir + inFileName;
    mbInFileName = "MinbiasAnalysis/" + mbInFileName;
    outFileName = "MinbiasAnalysis/" + outFileName;
    if (doHITightVar)
      bkg = new MinbiasAnalysis ("bkg_trackHITightVar");
    else
      bkg = new MinbiasAnalysis ("bkg");

    bkg->is2015Conds = use2015conds;
    bkg->useHijingEffs = use2015conds;
    bkg->useHITight = doHITightVar;

    //bkg->Execute (mbInFileName.c_str (), outFileName.c_str ()); // old code
    bkg->Execute (inFileName.c_str (), mbInFileName.c_str (), outFileName.c_str ()); // new code
    delete bkg;
  }


  return 0;
}
