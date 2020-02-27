#include "MCAnalysis.h"
#include "MCClosureAnalysis.h"
#include "MixedMCAnalysis.h"
#include "TruthAnalysis.h"
#include "MinbiasAnalysis.h"

int mixingFraction = 40;

int main (int argc, char** argv) {

  if (argc < 9) {
    cout << "Insufficient arguments! Exiting." << endl;
    return 1;
  }

  string algo = argv[1];

  string inFileName = argv[2];
  string mbInFileName = (string (argv[3]) == "0" ? inFileName : argv[3]); // defaults to 1st file name if "0"
  string outFileName = argv[4];
  string mixedFileName = (string (argv[5]) == "0" ? "" : string ("MixedEvents/") + string (argv[5]));
  

  string inDir = "DataAnalysis/";
  string histDir = "DataAnalysis/";
  if (algo == "mcminbias" || algo == "mc") {
    inDir = "MCAnalysis/"; 
    histDir = "MCAnalysis/";
    //outFileName = "Background/" + outFileName;
  }
  if (algo == "minbias") {
    inDir = "DataAnalysis/";
    histDir = "MinbiasAnalysis/";
  }

  const bool isPbPb       = (string (argv[6]) == "true");
  const bool use2015conds = (string (argv[7]) == "true");
  const bool doHITightVar = (string (argv[8]) == "true");

  const bool doMixVarA    = (string (argv[9]) == "doMixVarA");
  const bool doMixVarB    = (string (argv[9]) == "doMixVarB");
  const bool doMixVarC    = (string (argv[9]) == "doMixVarC");
  const bool doMixVarD    = (string (argv[9]) == "doMixVarD");
  const bool doMixVarE    = (string (argv[9]) == "doMixVarE");
  const bool doMixVarF    = (string (argv[9]) == "doMixVarF");
  const bool doMixVarG    = (string (argv[9]) == "doMixVarG");
  const bool doMixVarH    = (string (argv[9]) == "doMixVarH");
  const bool doPPMixVar   = (string (argv[9]) == "doPPMixVar");
  const bool doPbPbMixVar  = (doMixVarA || doMixVarB || doMixVarC || doMixVarD || doMixVarE || doMixVarF || doMixVarG || doMixVarH);


  if (!doHITightVar && !doPbPbMixVar && !doPPMixVar) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Nominal/" + outFileName;
  } else if (doHITightVar) {
    inFileName = "Variations/TrackHITightWPVariation/" + inFileName;
    mbInFileName = "Variations/TrackHITightWPVariation/" + mbInFileName;
    outFileName = "Variations/TrackHITightWPVariation/" + outFileName;
  } else if (doMixVarA) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationA/" + outFileName;
  } else if (doMixVarB) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationB/" + outFileName;
  } else if (doMixVarC) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationC/" + outFileName;
  } else if (doMixVarD) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationD/" + outFileName;
  } else if (doMixVarE) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationE/" + outFileName;
  } else if (doMixVarF) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationF/" + outFileName;
  } else if (doMixVarG) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationG/" + outFileName;
  } else if (doMixVarH) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/MixingVariationH/" + outFileName;
  } else if (doPPMixVar) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/PPMixingVariation/" + outFileName;
  }


  //if (algo == "mc" || algo == "mcminbias") {
  if (algo == "mc") {
    MCAnalysis* mc = nullptr;
    inFileName = inDir + inFileName;
    outFileName = histDir + outFileName;
    if (doHITightVar)
      mc = new MCAnalysis ("mc_trackHITightVar");
    //else if (algo == "mcminbias")
    //  mc = new MCAnalysis ("mc_bkg");
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

  //else if (algo == "mcclosure") {
  //  MCClosureAnalysis* mc = nullptr;
  //  inFileName = "MCAnalysis/" + inFileName;
  //  mbInFileName = "MinbiasAnalysis/" + mbInFileName;
  //  outFileName = "MCAnalysis/" + outFileName;
  //  mc = new MCClosureAnalysis ("mc_closure");

  //  mc->is2015Conds = use2015conds;
  //  mc->useHijingEffs = use2015conds;
  //  mc->useHITight = doHITightVar;

  //  mc->Execute (inFileName.c_str (), mbInFileName.c_str (), outFileName.c_str ());
  //  delete mc;
  //}

  //else if (algo == "mixmc") {
  //  MixedMCAnalysis* mc = nullptr;
  //  inFileName = "MCAnalysis/" + inFileName;
  //  outFileName = "MixedMCAnalysis/" + outFileName;
  //  mc = new MixedMCAnalysis ("mixed_mc");

  //  mc->is2015Conds = use2015conds;
  //  mc->useHijingEffs = use2015conds;
  //  mc->useHITight = doHITightVar;

  //  mc->takeNonTruthTracks = true;

  //  mc->doPsi2Mixing = true;
  //  mc->nPsi2MixBins = 16;
  //  mc->doQ2Mixing = false;

  //  mc->Execute (inFileName.c_str (), outFileName.c_str ());
  //  delete mc;
  //}

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

    inFileName = inDir + inFileName;

    if (!doPPMixVar && !isPbPb && algo == "minbias")
      mbInFileName = "DataAnalysis/" + mbInFileName;
    else if (!doPPMixVar && !isPbPb && algo == "mcminbias")
      mbInFileName = inFileName;
    else
      mbInFileName = "MinbiasAnalysis/" + mbInFileName;

    outFileName = histDir + outFileName;

    if (algo == "mcminbias")  bkg = new MinbiasAnalysis ("mc_bkg");
    else if (doHITightVar)    bkg = new MinbiasAnalysis ("bkg_trackHITightVar");
    else if (doMixVarA)       bkg = new MinbiasAnalysis ("bkg_mixVarA");
    else if (doMixVarB)       bkg = new MinbiasAnalysis ("bkg_mixVarB");
    else if (doMixVarC)       bkg = new MinbiasAnalysis ("bkg_mixVarC");
    else if (doMixVarD)       bkg = new MinbiasAnalysis ("bkg_mixVarD");
    else if (doMixVarE)       bkg = new MinbiasAnalysis ("bkg_mixVarE");
    else if (doMixVarF)       bkg = new MinbiasAnalysis ("bkg_mixVarF");
    else if (doMixVarG)       bkg = new MinbiasAnalysis ("bkg_mixVarG");
    else if (doMixVarH)       bkg = new MinbiasAnalysis ("bkg_mixVarH");
    else if (doPPMixVar)      bkg = new MinbiasAnalysis ("bkg_ppMixVar");
    else                      bkg = new MinbiasAnalysis ("bkg");

    if ((!isPbPb && !doPPMixVar))
      mixingFraction = 1;
    else if (algo == "mcminbias")
      mixingFraction = 10;

    bkg->is2015Conds = use2015conds;
    bkg->useHijingEffs = use2015conds; // for now
    bkg->useHITight = doHITightVar;

    if (doMixVarA) {
      bkg->doPsi2Mixing = false;
      bkg->nPsi2MixBins = 1;
      bkg->doPsi3Mixing = false;
      bkg->doQ2Mixing = false;
    } else if (doMixVarB) {
      bkg->nPsi2MixBins = 2;
    } else if (doMixVarC) {
      bkg->nPsi2MixBins = 4;
    } else if (doMixVarD) {
      bkg->nPsi2MixBins = 8;
    } else if (doMixVarE) {
      bkg->nPsi2MixBins = 16;
    } else if (doMixVarF) {
      bkg->nPsi2MixBins = 32;
    } else if (doMixVarG) {
      bkg->nPsi2MixBins = 64;
    } else if (doMixVarH) {
      bkg->nPsi2MixBins = 16;
      bkg->doPsi3Mixing = true;
      bkg->nPsi3MixBins = 3;
    } else if (doPPMixVar) {
      bkg->doPPTransMinMixing = false;
    }

    bkg->Execute (isPbPb, inFileName.c_str (), mbInFileName.c_str (), outFileName.c_str (), mixedFileName.c_str ()); // new code
    delete bkg;
  }


  return 0;
}
