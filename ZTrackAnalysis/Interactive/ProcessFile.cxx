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

  bool use2015conds = (string (argv[5]) == "true");
  bool doHITightVar = (string (argv[6]) == "true");
  //bool doQ2Mixing   = (string (argv[7]) == "true");
  //bool doPsi2Mixing = (string (argv[8]) == "true");
  bool doMixVarA    = (string (argv[7]) == "true");
  bool doMixVarB    = (string (argv[8]) == "true");
  bool doMixVarC    = (string (argv[9]) == "true");
  bool doMixVarD    = (string (argv[10]) == "true");
  bool doMixVarE    = (string (argv[11]) == "true");
  bool doMixVarF    = (string (argv[12]) == "true");


  if (!doHITightVar && /*!doQ2Mixing && !doPsi2Mixing*/ !doMixVarA && !doMixVarB && !doMixVarC && !doMixVarD && !doMixVarE && !doMixVarF) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Nominal/" + outFileName;
  }
  else if (doHITightVar) {
    inFileName = "Variations/TrackHITightWPVariation/" + inFileName;
    mbInFileName = "Variations/TrackHITightWPVariation/" + mbInFileName;
    outFileName = "Variations/TrackHITightWPVariation/" + outFileName;
  }
  //else if (doQ2Mixing) {
  //  inFileName = "Nominal/" + inFileName;
  //  mbInFileName = "Nominal/" + mbInFileName;
  //  outFileName = "Variations/Q2MixingVariation/" + outFileName;
  //}
  //else if (doPsi2Mixing) {
  //  inFileName = "Nominal/" + inFileName;
  //  mbInFileName = "Nominal/" + mbInFileName;
  //  outFileName = "Variations/Psi2MixingVariation/" + outFileName;
  //}
  else if (doMixVarA) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationA/" + outFileName;
  }
  else if (doMixVarB) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationB/" + outFileName;
  }
  else if (doMixVarC) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationC/" + outFileName;
  }
  else if (doMixVarD) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationD/" + outFileName;
  }
  else if (doMixVarE) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationE/" + outFileName;
  }
  else if (doMixVarF) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationF/" + outFileName;
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

    mc->takeNonTruthTracks = true;

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
    else if (doMixVarA)
      bkg = new MinbiasAnalysis ("bkg_mixVarA");
    else if (doMixVarB)
      bkg = new MinbiasAnalysis ("bkg_mixVarB");
    else if (doMixVarC)
      bkg = new MinbiasAnalysis ("bkg_mixVarC");
    else if (doMixVarD)
      bkg = new MinbiasAnalysis ("bkg_mixVarD");
    else if (doMixVarE)
      bkg = new MinbiasAnalysis ("bkg_mixVarE");
    else if (doMixVarF)
      bkg = new MinbiasAnalysis ("bkg_mixVarF");
    //else if (doQ2Mixing)
    //  bkg = new MinbiasAnalysis ("bkg_q2mixed");
    //else if (doPsi2Mixing)
    //  bkg = new MinbiasAnalysis ("bkg_psi2mixed");
    else
      bkg = new MinbiasAnalysis ("bkg");

    bkg->is2015Conds = use2015conds;
    bkg->useHijingEffs = use2015conds;
    bkg->useHITight = doHITightVar;
    if (doMixVarA) {
      bkg->doPsi2Mixing = false;
      bkg->doQ2Mixing = false;
    }
    if (doMixVarB) {
      bkg->doPsi2Mixing = true;
      bkg->numPsi2MixBins = 2;
      bkg->doQ2Mixing = false;
    }
    if (doMixVarC) {
      bkg->doPsi2Mixing = true;
      bkg->numPsi2MixBins = 4;
      bkg->doQ2Mixing = false;
    }
    if (doMixVarD) {
      bkg->doPsi2Mixing = true;
      bkg->numPsi2MixBins = 8;
      bkg->doQ2Mixing = false;
    }
    if (doMixVarE) {
      bkg->doPsi2Mixing = true;
      bkg->numPsi2MixBins = 16;
      bkg->doQ2Mixing = false;
    }
    if (doMixVarF) {
      bkg->doPsi2Mixing = true;
      bkg->numPsi2MixBins = 32;
      bkg->doQ2Mixing = false;
    }
    if (!doMixVarA && !doMixVarB && !doMixVarC && !doMixVarD && !doMixVarE && !doMixVarF) {
      bkg->doPsi2Mixing = true;
      bkg->numPsi2MixBins = 8;
      bkg->doQ2Mixing = false;
    }
    //bkg->doQ2Mixing = doQ2Mixing;
    //bkg->doPsi2Mixing = doPsi2Mixing;

    //bkg->Execute (mbInFileName.c_str (), outFileName.c_str ()); // old code
    bkg->Execute (inFileName.c_str (), mbInFileName.c_str (), outFileName.c_str ()); // new code
    delete bkg;
  }


  return 0;
}
