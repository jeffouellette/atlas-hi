#include "MCAnalysis.h"
//#include "MCClosureAnalysis.h"
//#include "MixedMCAnalysis.h"
#include "TruthAnalysis.h"
#include "MixingAnalysis.h"
#include "HijingAnalysis.h"
#include "HijingMixedAnalysis.h"
#include "HijingTruthAnalysis.h"

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

  const bool isPbPb               = (string (argv[6]) == "true");
  const bool use2015conds         = (string (argv[7]) == "true");
  

  string inDir = "DataAnalysis/";
  string histDir = "DataAnalysis/";
  if (algo == "mcminbias" || algo == "mc"/* || algo == "hijing_mixing"*/) {
    inDir = "MCAnalysis/"; 
    histDir = "MCAnalysis/";
    //outFileName = "Background/" + outFileName;
  }
  if (algo == "hijing_truth") {
    inDir = isPbPb ? "MixedEvents/" : "TruthAnalysis/";
    histDir = "TruthAnalysis/";
  }
  if (algo == "hijing") {
    inDir = isPbPb ? "MixedEvents/" : "MCAnalysis/";
    histDir = "MCAnalysis/";
  }
  if (algo == "minbias") {
    inDir = "DataAnalysis/";
    histDir = "MixingAnalysis/";
  }
  if (algo == "hijing_mixing") {
    inDir = isPbPb ? "MixedEvents/" : "MCAnalysis/";
    histDir = "MCAnalysis/";
  }

  const bool doHITightVar         = (string (argv[8]) == "doHITightVar");
  const bool doMuonPtUpVar        = (string (argv[8]) == "doMuonPtUpVar");
  const bool doMuonPtDownVar      = (string (argv[8]) == "doMuonPtDownVar");
  const bool doElectronPtUpVar    = (string (argv[8]) == "doElectronPtUpVar");
  const bool doElectronPtDownVar  = (string (argv[8]) == "doElectronPtDownVar");
  const bool doMixVarA            = (string (argv[8]) == "doMixVarA");
  const bool doMixVarB            = (string (argv[8]) == "doMixVarB");
  const bool doMixVarC            = (string (argv[8]) == "doMixVarC");
  const bool doMixVarD            = (string (argv[8]) == "doMixVarD");
  const bool doMixVarE            = (string (argv[8]) == "doMixVarE");
  const bool doMixVarF            = (string (argv[8]) == "doMixVarF");
  const bool doMixVarG            = (string (argv[8]) == "doMixVarG");
  const bool doMixVarH            = (string (argv[8]) == "doMixVarH");
  const bool doPPMBMixVar         = (string (argv[8]) == "doPPMBMixVar");
  const bool doPPTransMaxVar      = (string (argv[8]) == "doPPTransMaxVar");
  const bool doPPMixVar           = (doPPMBMixVar || doPPTransMaxVar);
  const bool doPbPbMixVar         = (doMixVarA || doMixVarB || doMixVarC || doMixVarD || doMixVarE || doMixVarF || doMixVarG || doMixVarH);


  if (!doHITightVar && !doElectronPtUpVar && !doElectronPtDownVar && !doMuonPtUpVar && !doMuonPtDownVar && !doPbPbMixVar && !doPPMixVar) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Nominal/" + outFileName;
  } else if (doHITightVar) {
    inFileName = "Variations/TrackHITightWPVariation/" + inFileName;
    mbInFileName = "Variations/TrackHITightWPVariation/" + mbInFileName;
    outFileName = "Variations/TrackHITightWPVariation/" + outFileName;
  } else if (doElectronPtUpVar) {
    inFileName = "Variations/ElectronPtUpVariation/" + inFileName;
    mbInFileName = "Variations/ElectronPtUpVariation/" + mbInFileName;
    outFileName = "Variations/ElectronPtUpVariation/" + outFileName;
  } else if (doElectronPtDownVar) {
    inFileName = "Variations/ElectronPtDownVariation/" + inFileName;
    mbInFileName = "Variations/ElectronPtDownVariation/" + mbInFileName;
    outFileName = "Variations/ElectronPtDownVariation/" + outFileName;
  } else if (doMuonPtUpVar) {
    inFileName = "Variations/MuonPtUpVariation/" + inFileName;
    mbInFileName = "Variations/MuonPtUpVariation/" + mbInFileName;
    outFileName = "Variations/MuonPtUpVariation/" + outFileName;
  } else if (doMuonPtDownVar) {
    inFileName = "Variations/MuonPtDownVariation/" + inFileName;
    mbInFileName = "Variations/MuonPtDownVariation/" + mbInFileName;
    outFileName = "Variations/MuonPtDownVariation/" + outFileName;
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
  } else if (doPPMBMixVar) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/PPMinBiasVariation/" + outFileName;
  } else if (doPPTransMaxVar) {
    inFileName = "Nominal/" + inFileName;
    mbInFileName = "Nominal/" + mbInFileName;
    outFileName = "Variations/PPTransMaxVariation/" + outFileName;
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
    else if (doElectronPtUpVar) {
      mc = new MCAnalysis ("mc_electronPtUpVar");
      mc->eventWeightsFileName = "MCAnalysis/Variations/ElectronPtUpVariation/eventWeightsFile.root";
    }
    else if (doElectronPtDownVar) {
      mc = new MCAnalysis ("mc_electronPtDownVar");
      mc->eventWeightsFileName = "MCAnalysis/Variations/ElectronPtDownVariation/eventWeightsFile.root";
    }
    else if (doMuonPtUpVar) {
      mc = new MCAnalysis ("mc_muonPtUpVar");
      mc->eventWeightsFileName = "MCAnalysis/Variations/MuonPtUpVariation/eventWeightsFile.root";
    }
    else if (doMuonPtDownVar) {
      mc = new MCAnalysis ("mc_muonPtDownVar");
      mc->eventWeightsFileName = "MCAnalysis/Variations/MuonPtDownVariation/eventWeightsFile.root";
    }
    else
      mc = new MCAnalysis ("mc");

    if (algo == "mcminbias")
      mc->takeNonTruthTracks = true;

    mc->is2015Conds = use2015conds;
    mc->useHITight = doHITightVar;
    mc->useCentWgts = true;

    mc->Execute (inFileName.c_str (), outFileName.c_str ());
    delete mc;
  }
  else if (algo == "hijing") {

    HijingAnalysis* mc = nullptr;
    inFileName = inDir + inFileName;
    outFileName = histDir + outFileName;

    mc = new HijingAnalysis ("hijing");

    mc->is2015Conds = use2015conds;
    mc->useHITight = doHITightVar;
    //mc->useCentWgts = true;
    mc->useHijingEffs = true;

    //mc->ewExt = "hijing";
    mixingFraction = 10;

    mc->Execute (inFileName.c_str (), outFileName.c_str ());
    delete mc;
  }
  else if (algo == "hijing_truth") {

    HijingTruthAnalysis* truth = nullptr;
    inFileName = inDir + inFileName;
    outFileName = histDir + outFileName;

    truth = new HijingTruthAnalysis ("hijing_truth");

    truth->is2015Conds = use2015conds;
    truth->useHITight = doHITightVar;
    //truth->useCentWgts = true;
    truth->useHijingEffs = true;

    mixingFraction = 10;

    truth->Execute (inFileName.c_str (), outFileName.c_str ());
    delete truth;
  }

  //else if (algo == "mcclosure") {
  //  MCClosureAnalysis* mc = nullptr;
  //  inFileName = "MCAnalysis/" + inFileName;
  //  mbInFileName = "MixingAnalysis/" + mbInFileName;
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

  else if (algo == "hijing_mixing") {
    HijingMixedAnalysis* bkg = nullptr;

    inFileName = inDir + inFileName;
    outFileName = histDir + outFileName;

    bkg = new HijingMixedAnalysis ("hijing_bkg");

    bkg->is2015Conds = use2015conds;
    bkg->useHITight = doHITightVar;

    bkg->Execute (inFileName.c_str (), outFileName.c_str ());
    delete bkg;
  }

  else if (algo == "minbias" || algo == "mcminbias"/* || algo == "hijing_mixing"*/) {
    MixingAnalysis* bkg = nullptr;

    inFileName = inDir + inFileName;

    if (!doPPMBMixVar && !isPbPb && algo == "minbias")
      mbInFileName = "DataAnalysis/" + mbInFileName;
    else if (!doPPMBMixVar && !isPbPb && (algo == "mcminbias"/* || algo == "hijing_mixing"*/))
      mbInFileName = inFileName;
    else
      mbInFileName = "MixingAnalysis/" + mbInFileName;

    outFileName = histDir + outFileName;

    if (algo == "mcminbias")  bkg = new MixingAnalysis ("mc_bkg");
    else if (doHITightVar)    bkg = new MixingAnalysis ("bkg_trackHITightVar");
    else if (doMixVarA)       bkg = new MixingAnalysis ("bkg_mixVarA");
    else if (doMixVarB)       bkg = new MixingAnalysis ("bkg_mixVarB");
    else if (doMixVarC)       bkg = new MixingAnalysis ("bkg_mixVarC");
    else if (doMixVarD)       bkg = new MixingAnalysis ("bkg_mixVarD");
    else if (doMixVarE)       bkg = new MixingAnalysis ("bkg_mixVarE");
    else if (doMixVarF)       bkg = new MixingAnalysis ("bkg_mixVarF");
    else if (doMixVarG)       bkg = new MixingAnalysis ("bkg_mixVarG");
    else if (doMixVarH)       bkg = new MixingAnalysis ("bkg_mixVarH");
    else if (doPPMBMixVar)    bkg = new MixingAnalysis ("bkg_ppMBMixVar");
    else if (doPPTransMaxVar) bkg = new MixingAnalysis ("bkg_ppTransMaxVar");
    else                      bkg = new MixingAnalysis ("bkg");

    //if (algo == "hijing_mixing")
    //  mixingFraction = 5;
    if ((!isPbPb && !doPPMBMixVar))
      mixingFraction = 1;
    else if (!isPbPb && !doPPMBMixVar) {
      bkg->doPPMBMixing = true;
      mixingFraction = 40;
    }
    else if (algo == "mcminbias")
      mixingFraction = 10;
    else if (isPbPb) {
      const int _first = inFileName.find ("iCent") + 5;
      const int _last = inFileName.find (".root");
      const int iCent = atoi (inFileName.substr (_first, _last-_first).c_str ());
      if (1 <= iCent && iCent <= 8)         mixingFraction = 160;
      else if (9 <= iCent && iCent <= 18)   mixingFraction = 80;
      else if (19 <= iCent && iCent <= 28)  mixingFraction = 40;
      else {
        cout << "Invalid centrality bin = " << iCent << "??? Check input & string parsing!" << endl;
        return -1;
      }
    }
    else
      mixingFraction = 40;

    bkg->is2015Conds = use2015conds;
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
    } else if (doPPMBMixVar) {
      bkg->doPPTransMinMixing = false;
      bkg->doPPMBMixing = true;
    } else if (doPPTransMaxVar) {
      bkg->doPPTransMinMixing = false;
      bkg->doPPTransMaxMixing = true;
    }

    bkg->Execute (isPbPb, inFileName.c_str (), mbInFileName.c_str (), outFileName.c_str (), mixedFileName.c_str ()); // new code
    delete bkg;
  }


  return 0;
}
