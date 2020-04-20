#ifndef __Run_C__
#define __Run_C__

#include "PhysicsAnalysis.h"
#include "FullAnalysis.h"
#include "MCAnalysis.h"
#include "MixingAnalysis.h"
#include "TruthAnalysis.h"
#include "Auxiliary.h"

#include "Systematic.h"
#include "TrackIDSystematic.h"
#include "ReweightingVariation.h"
#include "SystematicFits.h"
#include "LowPtVariation.h"
#include "BinWidthSystematic.h"

#include "PlotHybridModel.h"

const bool doSys = false;
int mixingFraction = 0;

// nominal analyses
FullAnalysis* data18 = nullptr;
//FullAnalysis* data15 = nullptr;
MixingAnalysis* bkg18 = nullptr;
//MixingAnalysis* bkg15 = nullptr;
MCAnalysis* mc = nullptr;
MixingAnalysis* mc_bkg = nullptr;
TruthAnalysis* truth = nullptr;

// master systematics objects
Systematic* combSys = nullptr;
Systematic* bkgStatSys = nullptr;
Systematic* bkgMixSys = nullptr;
Systematic* ppMixSys = nullptr;
TrackIDSystematic* trkSys = nullptr;
ReweightingVariation* trkEffSysWeights = nullptr;
Systematic* trkEffSys = nullptr;
ReweightingVariation* trkPurSysWeights = nullptr;
Systematic* trkPurSys = nullptr;
Systematic* leptonRejSys = nullptr;
SystematicFits* data_electronPtUpSysFits = nullptr, *data_electronPtDownSysFits = nullptr;
Systematic* electronPtSys = nullptr;
SystematicFits* data_muonPtUpSysFits = nullptr, *data_muonPtDownSysFits = nullptr;
Systematic* muonPtSys = nullptr;
Systematic* lowPtSys = nullptr;
LowPtVariation* lowPtSysWeights = nullptr;
Systematic* binCenterSys = nullptr;

//Systematic* mcBkgStatSys = nullptr;

// variations for systematics
MixingAnalysis* bkg_statUpVar = nullptr, *bkg_statDownVar = nullptr;
PhysicsAnalysis* data_bkgStatUpVar = nullptr, *data_bkgStatDownVar = nullptr;

MixingAnalysis* bkg_mixVarA = nullptr, *bkg_mixVarB = nullptr, *bkg_mixVarC = nullptr, *bkg_mixVarD = nullptr, *bkg_mixVarE = nullptr, *bkg_mixVarF = nullptr, *bkg_mixVarG = nullptr, *bkg_mixVarH = nullptr;
//MixingAnalysis* mcBkg_statUpVar = nullptr, *mcBkg_statDownVar = nullptr;
//PhysicsAnalysis* mc_bkgStatUpVar = nullptr, *mc_bkgStatDownVar = nullptr;
PhysicsAnalysis* data_mixVarA = nullptr, *data_mixVarB = nullptr, *data_mixVarC = nullptr, *data_mixVarD = nullptr, *data_mixVarE = nullptr, *data_mixVarF = nullptr, *data_mixVarG = nullptr, *data_mixVarH = nullptr;

MixingAnalysis* bkg_ppMBMixVar = nullptr, *bkg_ppTransMaxVar = nullptr;
PhysicsAnalysis* data_ppMBMixVar = nullptr, *data_ppTransMaxVar = nullptr;

PhysicsAnalysis* data_trackHItight = nullptr, *data_trkIDUpVar = nullptr, *data_trkIDDownVar = nullptr;
MixingAnalysis* bkg_trkIDUpVar = nullptr, *bkg_trkIDDownVar = nullptr;

PhysicsAnalysis* data_partComp = nullptr, *data_partCompUpVar = nullptr, *data_partCompDownVar = nullptr;
MixingAnalysis* bkg_partCompUpVar = nullptr, *bkg_partCompDownVar = nullptr;

PhysicsAnalysis* data_trkPurUpVar = nullptr, *data_trkPurDownVar = nullptr;
MixingAnalysis* bkg_trkPurUpVar = nullptr, *bkg_trkPurDownVar = nullptr;

PhysicsAnalysis* data_leptonRejVar = nullptr;

PhysicsAnalysis* data_electronPtUp = nullptr, *data_electronPtDown = nullptr;
MCAnalysis* mc_electronPtUp = nullptr, *mc_electronPtDown = nullptr;

PhysicsAnalysis* data_muonPtUp = nullptr, *data_muonPtDown = nullptr;
MCAnalysis* mc_muonPtUp = nullptr, *mc_muonPtDown = nullptr;

PhysicsAnalysis* data_lowPtUpVar = nullptr, *data_lowPtDownVar = nullptr;

PhysicsAnalysis* data_binCenterUpVar = nullptr, *data_binCenterDownVar = nullptr;



void Run () {
  data18    = new FullAnalysis ("data18");
  //data15  = new FullAnalysis ("data15");
  //data15->is2015Conds = true;
  //data15->useHijingEffs = true;
  bkg18     = new MixingAnalysis ("bkg");
  //bkg15->is2015Conds = true;
  //bkg15->useHijingEffs = true;

  mc      = new MCAnalysis ();
  mc_bkg  = new MixingAnalysis ("mc_bkg");
  //truth   = new TruthAnalysis ();

  //if (doSys) {
  //  mc_bkgStatUpVar         = new PhysicsAnalysis ("mc_bkgStatUpVar");
  //  mc_bkgStatDownVar       = new PhysicsAnalysis ("mc_bkgStatDownVar");
  //  mcBkg_statUpVar         = new MixingAnalysis ("mcBkg_statUpVar");
  //  mcBkg_statDownVar       = new MixingAnalysis ("mcBkg_statDownVar");
  //}

  //data18->Execute   ("DataAnalysis/Nominal/data18hi.root",   "DataAnalysis/Nominal/data18hi_hists.root");
  //data15->Execute ("DataAnalysis/Nominal/data15hi.root",  "DataAnalysis/Nominal/data15hi_hists.root");

  //if (doSys) {
  //  data_leptonRejVar       = new PhysicsAnalysis ("data_leptonRejVar");
  //  data_leptonRejVar->doLeptonRejVar = true;
  //  data_leptonRejVar->Execute      ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/LeptonRejVariation/data18hi_hists.root");

  //  data_trackHItight       = new PhysicsAnalysis ("data_trackHITightVar");
  //  data_trackHItight->useHITight = true;
  //  data_trackHItight->Execute      ("DataAnalysis/Variations/TrackHITightWPVariation/data18hi.root",      "DataAnalysis/Variations/TrackHITightWPVariation/data18hi_hists.root");

  //  data_partComp        = new PhysicsAnalysis ("data_trackEffVar");
  //  data_partComp->doTrackEffVar = true;
  //  data_partComp->Execute          ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackEffPionsVariation/data18hi_hists.root");

  //  data_trkPurUpVar      = new PhysicsAnalysis ("data_trkPurUpVar");
  //  data_trkPurUpVar->doTrackPurVar = true; data_trkPurUpVar->trkPurNSigma = 1;
  //  data_trkPurUpVar->Execute       ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackPurityUpVariation/data18hi_hists.root");

  //  data_trkPurDownVar    = new PhysicsAnalysis ("data_trkPurDownVar");
  //  data_trkPurDownVar->doTrackPurVar = true; data_trkPurDownVar->trkPurNSigma = -1;
  //  data_trkPurDownVar->Execute     ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackPurityDownVariation/data18hi_hists.root");
  //}


  //mc->LoadHists ("MCAnalysis/Nominal/savedHists.root");
  //mc->CombineHists ();

  mc->LoadHists ("MCAnalysis/Nominal/PbPb18_pp_Zee_Hijing/savedHists.root");
  //mc_bkg->LoadHists ("MCAnalysis/Variations/PPMixingVariation/mixedHists.root");
  mc_bkg->LoadHists ("MCAnalysis/Nominal/PbPb18_pp_Zee_Hijing/mixedHists.root");
  mc->SubtractBackground (mc_bkg);
  //mc->CalculateIAA ();
  //SaferDelete (&mc_bkg);


  data18->LoadHists ("DataAnalysis/Nominal/data18hi_hists.root");
  //data15->LoadHists ("DataAnalysis/Nominal/data15hi_hists.root");
  bkg18->LoadHists ("MixingAnalysis/Nominal/data18hi_hists.root");
  //bkg15->LoadHists ("MixingAnalysis/Nominal/data15hi_hists.root");
  //truth->LoadHists ("TruthAnalysis/Nominal/savedHists.root");

  data18->SubtractBackground (bkg18);
  //data18->CalculateTrackMeans (data18, data18->h_z_pt);
  //data18->CalculateIAA ();

  //data15->SubtractBackground (bkg15);
  //data15->CalculateIAA ();

  //truth->SubtractBackground ();
  //truth->CalculateIAA ();

  if (doSys) {
    cout << "Initializing systematic objects. " << endl;


    cout << "Calculating bkg. stat. systematic errors." << endl;
    bkgStatSys = new Systematic (data18, "bkgStatSys", "Bkg. Modeling Unc.");
    bkgStatSys->cancelIAA = false;
    bkgStatSys->meanTrackUncStoredAtCentralValues = false;
    data_bkgStatUpVar = new PhysicsAnalysis ("data_bkgStatUpVar");
    data_bkgStatUpVar->CopyAnalysis (data18, false);
    data_bkgStatDownVar = new PhysicsAnalysis ("data_bkgStatDownVar");
    data_bkgStatDownVar->CopyAnalysis (data18, false);
    bkg_statUpVar = new MixingAnalysis ("bkg_statUpVar");
    bkg_statUpVar->CopyAnalysis (bkg18, false);
    bkg_statUpVar->ConvertToStatVariation (true, 1);
    bkg_statDownVar = new MixingAnalysis ("bkg_statDownVar");
    bkg_statDownVar->CopyAnalysis (bkg18, false);
    bkg_statDownVar->ConvertToStatVariation (false, 1);
    data_bkgStatUpVar->SubtractBackground (bkg_statUpVar);
    data_bkgStatDownVar->SubtractBackground (bkg_statDownVar);
    data_bkgStatUpVar->CalculateTrackMeans (data18, data18->h_z_pt, bkg_statUpVar);
    data_bkgStatDownVar->CalculateTrackMeans (data18, data18->h_z_pt, bkg_statDownVar);
    data_bkgStatUpVar->CalculateIAA ();
    data_bkgStatDownVar->CalculateIAA ();
    bkgStatSys->AddVariation (data_bkgStatUpVar, true);
    bkgStatSys->AddVariation (data_bkgStatDownVar, true);
    bkgStatSys->AddVariations ();
    SaferDelete (&data_bkgStatUpVar);
    SaferDelete (&data_bkgStatDownVar);
    SaferDelete (&bkg_statUpVar);
    SaferDelete (&bkg_statDownVar);



    //cout << "Calculating mixed event systematic errors." << endl;
    //bkgMixSys = new Systematic (data18, "bkgMixSys", "Mixed event");
    //data_mixVarA    = new PhysicsAnalysis ("data_mixVarA");
    //data_mixVarB    = new PhysicsAnalysis ("data_mixVarB");
    //data_mixVarC    = new PhysicsAnalysis ("data_mixVarC");
    //data_mixVarD    = new PhysicsAnalysis ("data_mixVarD");
    //data_mixVarE    = new PhysicsAnalysis ("data_mixVarE");
    //data_mixVarF    = new PhysicsAnalysis ("data_mixVarF");
    //data_mixVarG    = new PhysicsAnalysis ("data_mixVarG");
    //data_mixVarH    = new PhysicsAnalysis ("data_mixVarH");
    //data_mixVarA->CopyAnalysis (data18, false);
    //data_mixVarB->CopyAnalysis (data18, false);
    //data_mixVarC->CopyAnalysis (data18, false);
    //data_mixVarD->CopyAnalysis (data18, false);
    //data_mixVarE->CopyAnalysis (data18, false);
    //data_mixVarF->CopyAnalysis (data18, false);
    //data_mixVarG->CopyAnalysis (data18, false);
    //data_mixVarH->CopyAnalysis (data18, false);
    //bkg_mixVarA     = new MixingAnalysis ("bkg_mixVarA");
    //bkg_mixVarB     = new MixingAnalysis ("bkg_mixVarB");
    //bkg_mixVarC     = new MixingAnalysis ("bkg_mixVarC");
    //bkg_mixVarD     = new MixingAnalysis ("bkg_mixVarD");
    //bkg_mixVarE     = new MixingAnalysis ("bkg_mixVarE");
    //bkg_mixVarF     = new MixingAnalysis ("bkg_mixVarF");
    //bkg_mixVarG     = new MixingAnalysis ("bkg_mixVarG");
    //bkg_mixVarH     = new MixingAnalysis ("bkg_mixVarH");
    //bkg_mixVarA->LoadHists ("MixingAnalysis/Variations/MixingVariationA/data18hi_hists.root");
    //bkg_mixVarB->LoadHists ("MixingAnalysis/Variations/MixingVariationB/data18hi_hists.root");
    //bkg_mixVarC->LoadHists ("MixingAnalysis/Variations/MixingVariationC/data18hi_hists.root");
    //bkg_mixVarD->LoadHists ("MixingAnalysis/Variations/MixingVariationD/data18hi_hists.root");
    //bkg_mixVarE->LoadHists ("MixingAnalysis/Variations/MixingVariationE/data18hi_hists.root");
    //bkg_mixVarF->LoadHists ("MixingAnalysis/Variations/MixingVariationF/data18hi_hists.root");
    //bkg_mixVarG->LoadHists ("MixingAnalysis/Variations/MixingVariationG/data18hi_hists.root");
    //bkg_mixVarH->LoadHists ("MixingAnalysis/Variations/MixingVariationH/data18hi_hists.root");
    //data_mixVarA->SubtractBackground (bkg_mixVarA);
    //data_mixVarA->CalculateIAA ();
    //data_mixVarB->SubtractBackground (bkg_mixVarB);
    //data_mixVarB->CalculateIAA ();
    //data_mixVarC->SubtractBackground (bkg_mixVarC);
    //data_mixVarC->CalculateIAA ();
    //data_mixVarD->SubtractBackground (bkg_mixVarD);
    //data_mixVarD->CalculateIAA ();
    //data_mixVarE->SubtractBackground (bkg_mixVarE);
    //data_mixVarE->CalculateIAA ();
    //data_mixVarF->SubtractBackground (bkg_mixVarF);
    //data_mixVarF->CalculateIAA ();
    //data_mixVarG->SubtractBackground (bkg_mixVarG);
    //data_mixVarG->CalculateIAA ();
    //data_mixVarH->SubtractBackground (bkg_mixVarH);
    //data_mixVarH->CalculateIAA ();
    //bkgMixSys->AddVariation (data_mixVarA);
    //bkgMixSys->AddVariation (data_mixVarB);
    //bkgMixSys->AddVariation (data_mixVarC);
    //bkgMixSys->AddVariation (data_mixVarD);
    //bkgMixSys->AddVariation (data_mixVarE);
    //bkgMixSys->AddVariation (data_mixVarF);
    //bkgMixSys->AddVariation (data_mixVarG);
    //bkgMixSys->AddVariation (data_mixVarH);
    //bkgMixSys->AddVarDesc (data_mixVarA, "1 #Psi_{2} bin");
    //bkgMixSys->AddVarDesc (data_mixVarB, "2 #Psi_{2} bins");
    //bkgMixSys->AddVarDesc (data_mixVarC, "4 #Psi_{2} bins");
    //bkgMixSys->AddVarDesc (data_mixVarD, "8 #Psi_{2} bins");
    //bkgMixSys->AddVarDesc (data_mixVarE, "16 #Psi_{2} bins");
    //bkgMixSys->AddVarDesc (data_mixVarF, "32 #Psi_{2} bins");
    //bkgMixSys->AddVarDesc (data_mixVarG, "64 #Psi_{2} bins");
    //bkgMixSys->AddVarDesc (data_mixVarH, "16 #Psi_{2} bins + 3 #Psi_{3} bins");
    //bkgMixSys->AddVariations ();
    ////bkgMixSys->AddVariationsUsingStdDev ();
    ////SaferDelete (&data_mixVarA);
    ////SaferDelete (&data_mixVarB);
    ////SaferDelete (&data_mixVarC);
    ////SaferDelete (&data_mixVarD);
    ////SaferDelete (&data_mixVarE);
    ////SaferDelete (&data_mixVarF);
    ////SaferDelete (&data_mixVarF);
    ////SaferDelete (&data_mixVarH);
    //SaferDelete (&bkg_mixVarA);
    //SaferDelete (&bkg_mixVarB);
    //SaferDelete (&bkg_mixVarC);
    //SaferDelete (&bkg_mixVarD);
    //SaferDelete (&bkg_mixVarE);
    //SaferDelete (&bkg_mixVarF);
    //SaferDelete (&bkg_mixVarG);
    //SaferDelete (&bkg_mixVarH);



    //cout << "Calculating pp mixing systematic errors" << endl;
    //ppMixSys = new Systematic (data18, "ppMixSys", "#it{pp} mixing");
    //ppMixSys->AddVariation (data18);
    //ppMixSys->AddVarDesc (data18, "Z-Z #it{pp} mixing (trans-min)");

    //data_ppTransMaxVar = new PhysicsAnalysis ("data_ppTransMaxVar");
    //data_ppTransMaxVar->CopyAnalysis (data18, false);
    //bkg_ppTransMaxVar = new MixingAnalysis ("bkg_ppTransMaxVar");
    //bkg_ppTransMaxVar->doPPTransMinMixing = false;
    //bkg_ppTransMaxVar->doPPTransMaxMixing = true;
    //bkg_ppTransMaxVar->doPPMBMixing = false;
    //bkg_ppTransMaxVar->LoadHists ("MixingAnalysis/Variations/PPTransMaxVariation/data18hi_hists.root");
    //data_ppTransMaxVar->SubtractBackground (bkg_ppTransMaxVar);
    //data_ppTransMaxVar->CalculateIAA ();
    //ppMixSys->AddVariation (data_ppTransMaxVar);
    //ppMixSys->AddVarDesc (data_ppTransMaxVar, "Z-Z #it{pp} mixing (trans-max)");
    //SaferDelete (&bkg_ppTransMaxVar);

    //data_ppMBMixVar = new PhysicsAnalysis ("data_ppMBMixVar");
    //data_ppMBMixVar->CopyAnalysis (data18, false);
    //bkg_ppMBMixVar = new MixingAnalysis ("bkg_ppMBMixVar");
    //bkg_ppMBMixVar->doPPTransMinMixing = false;
    //bkg_ppMBMixVar->doPPMBMixing = true;
    //bkg_ppMBMixVar->LoadHists ("MixingAnalysis/Variations/PPMinBiasVariation/data18hi_hists.root");
    //data_ppMBMixVar->SubtractBackground (bkg_ppMBMixVar);
    //data_ppMBMixVar->CalculateIAA ();
    //ppMixSys->AddVariation (data_ppMBMixVar);
    //ppMixSys->AddVarDesc (data_ppMBMixVar, "Z-MinBias #it{pp} mixing");
    //SaferDelete (&bkg_ppMBMixVar);

    //ppMixSys->AddVariations ();
    ////SaferDelete (&data_ppMBMixVar);
    ////SaferDelete (&data_ppTransMaxVar);


    cout << "Calculating track ID relative systematic errors." << endl;
    trkSys = new TrackIDSystematic (data18, "trkSys", "Track ID Cuts");
    data_trackHItight = new PhysicsAnalysis ("data_trackHITightVar");
    data_trackHItight->useHITight = true;
    data_trackHItight->LoadHists ("DataAnalysis/Variations/TrackHITightWPVariation/data18hi_hists.root");
    trkSys->GetRelativeVariation (data18, data_trackHItight);
    SaferDelete (&data_trackHItight);
    cout << "Applying track ID systematic errors." << endl;
    data_trkIDUpVar = new PhysicsAnalysis ("data_trkIDUpVar");
    data_trkIDUpVar->useHITight = true;
    data_trkIDUpVar->CopyAnalysis (data18, false);
    data_trkIDDownVar = new PhysicsAnalysis ("data_trkIDDownVar");
    data_trkIDDownVar->useHITight = true;
    data_trkIDDownVar->CopyAnalysis (data18, false);
    bkg_trkIDUpVar = new MixingAnalysis ("bkg_trkIDUpVar");
    bkg_trkIDUpVar->useHITight = true;
    bkg_trkIDUpVar->CopyAnalysis (bkg18, false);
    bkg_trkIDDownVar = new MixingAnalysis ("bkg_trkIDDownVar");
    bkg_trkIDDownVar->useHITight = true;
    bkg_trkIDDownVar->CopyAnalysis (bkg18, false);

    data_trkIDUpVar->ApplyRelativeVariation (trkSys->relVar, true); // track yields go up
    data_trkIDDownVar->ApplyRelativeVariation (trkSys->relVar, false); // track yields go down
    bkg_trkIDUpVar->ApplyRelativeVariation (trkSys->relVar, true); // track yields go up
    bkg_trkIDDownVar->ApplyRelativeVariation (trkSys->relVar, false); // track yields go down

    data_trkIDUpVar->SubtractBackground (bkg_trkIDUpVar);
    data_trkIDDownVar->SubtractBackground (bkg_trkIDDownVar);
    data_trkIDUpVar->CalculateTrackMeans (data_trkIDUpVar, data18->h_z_pt);
    data_trkIDDownVar->CalculateTrackMeans (data_trkIDDownVar, data18->h_z_pt);
    data_trkIDUpVar->CalculateIAA ();
    data_trkIDDownVar->CalculateIAA ();
    trkSys->AddVariation (data_trkIDUpVar, true);
    trkSys->AddVariation (data_trkIDDownVar, true);
    trkSys->AddVariations ();
    SaferDelete (&data_trkIDUpVar);
    SaferDelete (&data_trkIDDownVar);
    SaferDelete (&bkg_trkIDUpVar);
    SaferDelete (&bkg_trkIDDownVar);



    cout << "Calculating particle composition systematic errors." << endl;
    trkEffSysWeights = new ReweightingVariation ("trkEffSysWeights");
    trkEffSys = new Systematic (data18, "trkEffSys", "Particle Composition");
    trkEffSys->cancelIAA = false;
    data_partComp        = new PhysicsAnalysis ("data_trackEffVar");
    data_partComp->doTrackEffVar = true;
    data_partComp->LoadHists ("DataAnalysis/Variations/TrackEffPionsVariation/data18hi_hists.root");
    trkEffSysWeights->GetRelativeVariations (data18, data_partComp);
    data_partCompUpVar = new PhysicsAnalysis ("data_partCompUpVar");
    data_partCompUpVar->doTrackEffVar = true;
    data_partCompUpVar->CopyAnalysis (data18, false);
    data_partCompDownVar = new PhysicsAnalysis ("data_partCompDownVar");
    data_partCompDownVar->doTrackEffVar = true;
    data_partCompDownVar->CopyAnalysis (data18, false);
    bkg_partCompUpVar = new MixingAnalysis ("bkg_partCompUpVar");
    bkg_partCompUpVar->doTrackEffVar = true;
    bkg_partCompUpVar->CopyAnalysis (bkg18, false);
    bkg_partCompDownVar = new MixingAnalysis ("bkg_partCompDownVar");
    bkg_partCompDownVar->doTrackEffVar = true;
    bkg_partCompDownVar->CopyAnalysis (bkg18, false);
    trkEffSysWeights->ApplyRelativeVariations (data_partCompUpVar, true);
    trkEffSysWeights->ApplyRelativeVariations (bkg_partCompUpVar, true);
    trkEffSysWeights->ApplyRelativeVariations (data_partCompDownVar, false);
    trkEffSysWeights->ApplyRelativeVariations (bkg_partCompDownVar, false);
    SaferDelete (&trkEffSysWeights);

    data_partCompUpVar->SubtractBackground (bkg_partCompUpVar);
    data_partCompDownVar->SubtractBackground (bkg_partCompDownVar);
    data_partCompUpVar->CalculateTrackMeans (data_partCompUpVar, data18->h_z_pt);
    data_partCompDownVar->CalculateTrackMeans (data_partCompDownVar, data18->h_z_pt);
    data_partCompUpVar->CalculateIAA ();
    data_partCompDownVar->CalculateIAA ();
    trkEffSys->AddVariation (data_partCompUpVar);
    trkEffSys->AddVariation (data_partCompDownVar);
    trkEffSys->AddVariations ();
    SaferDelete (&data_partCompUpVar);
    SaferDelete (&data_partCompDownVar);
    SaferDelete (&bkg_partCompUpVar);
    SaferDelete (&bkg_partCompDownVar);



    cout << "Calculating track purity systematic errors." << endl;
    trkPurSys = new Systematic (data18, "trkPurSys", "Purity Correction");
    data_trkPurUpVar = new PhysicsAnalysis ("data_trkPurUpVar");
    data_trkPurUpVar->doTrackPurVar = true; data_trkPurUpVar->trkPurNSigma = 1;
    data_trkPurUpVar->LoadHists ("DataAnalysis/Variations/TrackPurityUpVariation/data18hi_hists.root");
    trkPurSysWeights = new ReweightingVariation ("trkPurSysWeights");
    trkPurSysWeights->GetRelativeVariations (data18, data_trkPurUpVar);
    data_trkPurUpVar->CopyAnalysis (data18, false);
    bkg_trkPurUpVar = new MixingAnalysis ("bkg_trkPurUpVar");
    bkg_trkPurUpVar->doTrackPurVar = true; bkg_trkPurUpVar->trkPurNSigma = 1;
    bkg_trkPurUpVar->CopyAnalysis (bkg18, false);
    trkPurSysWeights->ApplyRelativeVariations (data_trkPurUpVar, true);
    trkPurSysWeights->ApplyRelativeVariations (bkg_trkPurUpVar, true);
    SaferDelete (&trkPurSysWeights);

    data_trkPurDownVar = new PhysicsAnalysis ("data_trkPurDownVar");
    data_trkPurDownVar->doTrackPurVar = true; data_trkPurDownVar->trkPurNSigma = -1;
    data_trkPurDownVar->LoadHists ("DataAnalysis/Variations/TrackPurityDownVariation/data18hi_hists.root");
    trkPurSysWeights = new ReweightingVariation ("trkPurSysWeights");
    trkPurSysWeights->GetRelativeVariations (data18, data_trkPurDownVar);
    data_trkPurDownVar->CopyAnalysis (data18, false);
    bkg_trkPurDownVar = new MixingAnalysis ("bkg_trkPurDownVar");
    bkg_trkPurDownVar->doTrackPurVar = true; bkg_trkPurDownVar->trkPurNSigma = -1;
    bkg_trkPurDownVar->CopyAnalysis (bkg18, false);
    trkPurSysWeights->ApplyRelativeVariations (data_trkPurDownVar, false);
    trkPurSysWeights->ApplyRelativeVariations (bkg_trkPurDownVar, false);
    SaferDelete (&trkPurSysWeights);

    data_trkPurUpVar->SubtractBackground (bkg_trkPurUpVar);
    data_trkPurDownVar->SubtractBackground (bkg_trkPurDownVar);
    data_trkPurUpVar->CalculateTrackMeans (data_trkPurUpVar, data18->h_z_pt);
    data_trkPurDownVar->CalculateTrackMeans (data_trkPurDownVar, data18->h_z_pt);
    data_trkPurUpVar->CalculateIAA ();
    data_trkPurDownVar->CalculateIAA ();
    trkPurSys->AddVariation (data_trkPurUpVar);
    trkPurSys->AddVariation (data_trkPurDownVar);
    trkPurSys->AddVariations ();
    SaferDelete (&data_trkPurUpVar);
    SaferDelete (&data_trkPurDownVar);
    SaferDelete (&bkg_trkPurUpVar);
    SaferDelete (&bkg_trkPurDownVar);



    cout << "Calculating lepton rejection systematic errors." << endl;
    leptonRejSys = new Systematic (data18, "leptonRejSys", "Lepton Rejection");
    data_leptonRejVar = new PhysicsAnalysis ("data_leptonRejVar");
    data_leptonRejVar->doLeptonRejVar = true;
    data_leptonRejVar->LoadHists ("DataAnalysis/Variations/LeptonRejVariation/data18hi_hists.root");
    data_leptonRejVar->SubtractBackground (bkg18);
    data_leptonRejVar->CalculateTrackMeans (data_leptonRejVar, data18->h_z_pt);
    data_leptonRejVar->CalculateIAA ();
    leptonRejSys->AddVariation (data_leptonRejVar, -1);
    leptonRejSys->AddVariations ();
    SaferDelete (&data_leptonRejVar);



    cout << "Calculating electron ES systematic errors." << endl;
    data_electronPtUpSysFits = new SystematicFits ("data_electronPtUpSysFit");
    data_electronPtDownSysFits = new SystematicFits ("data_electronPtDownSysFit");
    electronPtSys = new Systematic (data18, "electronPtSys", "Electron Energy Scale");
    electronPtSys->cancelIAA = false;

    mc_electronPtUp = new MCAnalysis ("mc_electronPtUpVar");
    mc_electronPtUp->LoadHists ("MCAnalysis/Variations/ElectronPtUpVariation/savedHists.root");
    mc_electronPtUp->CombineHists ();
    data_electronPtUpSysFits->GetRelativeVariations (mc, mc_electronPtUp);
    data_electronPtUp = new PhysicsAnalysis ("data_electronPtUpVar");
    data_electronPtUp->CopyAnalysis (data18, true, false);
    data_electronPtUpSysFits->ApplyRelativeVariations (data_electronPtUp, true);
    data_electronPtUp->CalculateTrackMeans (data_electronPtUp, mc_electronPtUp->h_z_pt);
    data_electronPtUp->CalculateIAA ();
    SaferDelete (&mc_electronPtUp);
    SaferDelete (&data_electronPtUpSysFits);

    mc_electronPtDown = new MCAnalysis ("mc_electronPtDownVar");
    mc_electronPtDown->LoadHists ("MCAnalysis/Variations/ElectronPtDownVariation/savedHists.root");
    mc_electronPtDown->CombineHists ();
    data_electronPtDownSysFits->GetRelativeVariations (mc, mc_electronPtDown);
    data_electronPtDown = new PhysicsAnalysis ("data_electronPtDownVar");
    data_electronPtDown->CopyAnalysis (data18, true, false);
    data_electronPtDownSysFits->ApplyRelativeVariations (data_electronPtDown, true);
    data_electronPtDown->CalculateTrackMeans (data_electronPtDown, mc_electronPtDown->h_z_pt);
    data_electronPtDown->CalculateIAA ();
    SaferDelete (&mc_electronPtDown);
    SaferDelete (&data_electronPtDownSysFits);

    electronPtSys->AddVariation (data_electronPtUp);
    electronPtSys->AddVariation (data_electronPtDown);
    electronPtSys->AddVariations ();
    SaferDelete (&data_electronPtUp);
    SaferDelete (&data_electronPtDown);



    cout << "Calculating muon ES systematic errors." << endl;
    data_muonPtUpSysFits = new SystematicFits ("data_muonPtUpSysFit");
    data_muonPtDownSysFits = new SystematicFits ("data_muonPtDownSysFit");
    muonPtSys = new Systematic (data18, "muonPtSys", "Muon Energy Scale");
    muonPtSys->cancelIAA = false;

    mc_muonPtUp = new MCAnalysis ("mc_muonPtUpVar");
    mc_muonPtUp->LoadHists ("MCAnalysis/Variations/MuonPtUpVariation/savedHists.root");
    mc_muonPtUp->CombineHists ();
    data_muonPtUpSysFits->GetRelativeVariations (mc, mc_muonPtUp);
    data_muonPtUp = new PhysicsAnalysis ("data_muonPtUpVar");
    data_muonPtUp->CopyAnalysis (data18, true, false);
    data_muonPtUpSysFits->ApplyRelativeVariations (data_muonPtUp, true);
    data_muonPtUp->CalculateTrackMeans (data_muonPtUp, mc_muonPtUp->h_z_pt);
    data_muonPtUp->CalculateIAA ();
    SaferDelete (&mc_muonPtUp);
    SaferDelete (&data_muonPtUpSysFits);

    mc_muonPtDown = new MCAnalysis ("mc_muonPtDownVar");
    mc_muonPtDown->LoadHists ("MCAnalysis/Variations/MuonPtDownVariation/savedHists.root");
    mc_muonPtDown->CombineHists ();
    data_muonPtDownSysFits->GetRelativeVariations (mc, mc_muonPtDown);
    data_muonPtDown = new PhysicsAnalysis ("data_muonPtDownVar");
    data_muonPtDown->CopyAnalysis (data18, true, false);
    data_muonPtDownSysFits->ApplyRelativeVariations (data_muonPtDown, true);
    data_muonPtDown->CalculateTrackMeans (data_muonPtDown, mc_muonPtDown->h_z_pt);
    data_muonPtDown->CalculateIAA ();
    SaferDelete (&mc_muonPtDown);
    SaferDelete (&data_muonPtDownSysFits);

    muonPtSys->AddVariation (data_muonPtUp);
    muonPtSys->AddVariation (data_muonPtDown);
    muonPtSys->AddVariations ();
    SaferDelete (&data_muonPtUp);
    SaferDelete (&data_muonPtDown);



    cout << "Calculating low pT additional systematic errors." << endl;
    lowPtSys = new Systematic (data18, "lowPtSys", "ee/#mu#mu channel difference");
    lowPtSysWeights = new LowPtVariation ("lowPtSysWeights");
    lowPtSysWeights->GetRelativeVariations (data18, bkg18);
    data_lowPtUpVar = new PhysicsAnalysis ("data_lowPtUpVar");
    data_lowPtDownVar = new PhysicsAnalysis ("data_lowPtDownVar");
    data_lowPtUpVar->CopyAnalysis (data18, true, false);
    data_lowPtDownVar->CopyAnalysis (data18, true, false);
    lowPtSysWeights->ApplyRelativeVariations (data_lowPtUpVar, true);
    lowPtSysWeights->ApplyRelativeVariations (data_lowPtDownVar, false);
    SaferDelete (&lowPtSysWeights);

    data_lowPtUpVar->CalculateTrackMeans (data_lowPtUpVar, data18->h_z_pt);
    data_lowPtDownVar->CalculateTrackMeans (data_lowPtDownVar, data18->h_z_pt);
    data_lowPtUpVar->CalculateIAA ();
    data_lowPtDownVar->CalculateIAA ();
    lowPtSys->AddVariation (data_lowPtUpVar);
    lowPtSys->AddVariation (data_lowPtDownVar);
    lowPtSys->AddVariations ();
    SaferDelete (&data_lowPtUpVar);
    SaferDelete (&data_lowPtDownVar);



    cout << "Calculating bin width systematic error in <pTch>, <xhZ> measurements." << endl;
    binCenterSys = new BinWidthSystematic (data18, "binCenterSys", "Bin center variation");
    data_binCenterUpVar = new PhysicsAnalysis ("data_binCenterUpVar");
    data_binCenterDownVar = new PhysicsAnalysis ("data_binCenterDownVar");
    data_binCenterUpVar->CalculateTrackMeans (data18, data18->h_z_pt, nullptr, 1);
    data_binCenterDownVar->CalculateTrackMeans (data18, data18->h_z_pt, nullptr, -1);
    binCenterSys->AddVariation (data_binCenterUpVar);
    binCenterSys->AddVariation (data_binCenterDownVar);
    binCenterSys->AddVariations ();
    SaferDelete (&data_binCenterUpVar);
    SaferDelete (&data_binCenterDownVar);


    cout << "Adding errors in quadrature." << endl;
    combSys = new Systematic (data18, "combSys", "Total");
    combSys->AddSystematic (trkSys);
    combSys->AddSystematic (trkEffSys);
    combSys->AddSystematic (trkPurSys);
    combSys->AddSystematic (leptonRejSys);
    combSys->AddSystematic (electronPtSys);
    combSys->AddSystematic (muonPtSys);
    combSys->AddSystematic (lowPtSys);
    combSys->AddSystematic (bkgStatSys);
    combSys->AddSystematic (binCenterSys);
    combSys->AddSystematics ();

    // the following lines are necessary to avoid messing up xhZ systematic plots!
    //cout << "Truncating excess bins." << endl;
    //trkSys->TruncatePhysicsPlots ();
    //trkEffSys->TruncatePhysicsPlots ();
    //trkPurSys->TruncatePhysicsPlots ();
    //leptonRejSys->TruncatePhysicsPlots ();
    //electronPtSys->TruncatePhysicsPlots ();
    //muonPtSys->TruncatePhysicsPlots ();
    //lowPtSys->TruncatePhysicsPlots ();
    //bkgStatSys->TruncatePhysicsPlots ();
    //binCenterSys->TruncatePhysicsPlots ();
    //combSys->TruncatePhysicsPlots ();


    /*
    trkSys->SaveGraphs ("Systematics/TrackIDSys.root");
    trkEffSys->SaveGraphs ("Systematics/PartCompSys.root");
    trkPurSys->SaveGraphs ("Systematics/TrackPurSys.root");
    leptonRejSys->SaveGraphs ("Systematics/LeptonRejSys.root");
    electronPtSys->SaveGraphs ("Systematics/ElectronPtSys.root");
    muonPtSys->SaveGraphs ("Systematics/MuonPtSys.root");
    bkgStatSys->SaveGraphs ("Systematics/BkgStatSys.root");
    */
    combSys->SaveGraphs ("Systematics/CombinedSys.root"); 
  }

  data18->CalculateTrackMeans (data18, data18->h_z_pt);
  data18->CalculateIAA ();



  SetupDirectories ("", "ZTrackAnalysis/");

}




void MakePhysicsPlots () {

  bkg18->PlotAllYields_dPhi (1, 2, 4);
  //combSys->PlotAllYields_dPhi (1, 2, 4);
  data18->PlotAllYields_dPhi (1, 2, 4);

  max_iaa=4.4;
  //combSys->PlotIAA_dPhi (1, 2, 4);
  data18->PlotIAA_dPhi (1, 2, 4);
 
  bkg18->PlotAllYields_dPtZ (1, 2); 
  if (doSys) combSys->PlotAllYields_dPtZ (1, 2);
  data18->PlotAllYields_dPtZ (1, 2);

  bkg18->PlotAllYields_dPtZ (0, 2); 
  if (doSys) combSys->PlotAllYields_dPtZ (0, 2);
  data18->PlotAllYields_dPtZ (0, 2);

  max_iaa = 10;
  if (doSys) combSys->PlotSingleIAA_dPtZ (1);
  data18->PlotSingleIAA_dPtZ (1);

  if (doSys) combSys->PlotSingleIAA_dPtZ (0);
  data18->PlotSingleIAA_dPtZ (0);

  max_iaa = 3.2;
  if (doSys) combSys->PlotIAA_dPtZ (1);
  data18->PlotIAA_dPtZ (1);

  max_iaa = 3.2;
  if (doSys) combSys->PlotIAA_dPtZ (0);
  data18->PlotIAA_dPtZ (0);

  combSys->PlotSubYields_dPtZ_Fits (0);
  data18->PlotSubYields_dPtZ_Fits (0);

  combSys->PlotSubYields_dPtZ_Fits (1);
  data18->PlotSubYields_dPtZ_Fits (1);

}




void CompareTrkYieldSpcComp (const short iCent = 1, const short iPtZ = nPtZBins-1) {
  TCanvas* c = new TCanvas ("c", "", 900, 600);
  const double dPadY = 0.5;
  const double uPadY = 1. - dPadY;
  const int axisTextSize = 23;

  const double llMargin = 0.17;
  const double lrMargin = 0.024;
  const double rlMargin = 0.100;
  const double rrMargin = 0.02;

  const double a = (double) 1./(2. + (llMargin+lrMargin)/(1.-llMargin-lrMargin) + (rlMargin+rrMargin)/(1.-rlMargin-rrMargin));
  const double xPadMiddle = a * (1 + (llMargin+lrMargin)/(1.-llMargin-lrMargin));

  const double yPadMiddle = 0.5;

  c->cd ();

  const char* luPadName = Form ("luPad_%i", iCent);
  const char* ldPadName = Form ("ldPad_%i", iCent);
  const char* ruPadName = Form ("ruPad_%i", iCent);
  const char* rdPadName = Form ("rdPad_%i", iCent);

  TPad* luPad = new TPad (luPadName, "", 0, dPadY, xPadMiddle, 1);
  TPad* ldPad = new TPad (ldPadName, "", 0, 0, xPadMiddle, dPadY);
  TPad* ruPad = new TPad (ruPadName, "", xPadMiddle, dPadY, 1, 1);
  TPad* rdPad = new TPad (rdPadName, "", xPadMiddle, 0, 1, dPadY);

  luPad->SetLeftMargin (llMargin);
  luPad->SetRightMargin (lrMargin);
  ruPad->SetLeftMargin (rlMargin);
  ruPad->SetRightMargin (rrMargin);
  luPad->SetTopMargin (0.02);
  ruPad->SetTopMargin (0.02);
  luPad->SetBottomMargin (0);
  ruPad->SetBottomMargin (0);

  ldPad->SetLeftMargin (llMargin);
  ldPad->SetRightMargin (lrMargin);
  rdPad->SetLeftMargin (rlMargin);
  rdPad->SetRightMargin (rrMargin);
  ldPad->SetTopMargin (0.02);
  rdPad->SetTopMargin (0.02);
  ldPad->SetBottomMargin (0.25);
  rdPad->SetBottomMargin (0.25);

  luPad->Draw ();
  ruPad->Draw ();
  ldPad->Draw ();
  rdPad->Draw ();


  PhysicsAnalysis* a1 = data18;
  PhysicsAnalysis* a2 = bkg18;
  TH1D* h_ee = nullptr, *h_ee_bkg = nullptr, *h_ee_sub = nullptr;
  TH1D* h_mumu = nullptr, *h_mumu_bkg = nullptr, *h_mumu_sub = nullptr;

  h_ee = a1->h_trk_pt_ptz[0][iPtZ][iCent];
  h_mumu = a1->h_trk_pt_ptz[1][iPtZ][iCent];

  h_ee_bkg = a2->h_trk_pt_ptz[0][iPtZ][iCent];
  h_mumu_bkg = a2->h_trk_pt_ptz[1][iPtZ][iCent];
  
  h_ee_sub = a1->h_trk_pt_ptz_sub[0][iPtZ][iCent];
  h_mumu_sub = a1->h_trk_pt_ptz_sub[1][iPtZ][iCent];

  h_ee->SetMarkerColor (kRed);
  h_ee->SetLineColor (kRed);
  h_ee_bkg->SetMarkerColor (kRed);
  h_ee_bkg->SetLineColor (kRed);
  h_ee_sub->SetMarkerColor (kRed);
  h_ee_sub->SetLineColor (kRed);
  h_mumu->SetMarkerColor (kBlue);
  h_mumu->SetLineColor (kBlue);
  h_mumu_bkg->SetMarkerColor (kBlue);
  h_mumu_bkg->SetLineColor (kBlue);
  h_mumu_sub->SetMarkerColor (kBlue);
  h_mumu_sub->SetLineColor (kBlue);
  h_ee_bkg->SetLineStyle (2);
  h_mumu_bkg->SetLineStyle (2);
  h_ee->SetMarkerStyle (kOpenSquare);
  h_ee_bkg->SetMarkerStyle (kOpenSquare);
  h_ee_sub->SetMarkerStyle (kOpenSquare);
  h_mumu->SetMarkerStyle (kOpenCircle);
  h_mumu_bkg->SetMarkerStyle (kOpenCircle);
  h_mumu_sub->SetMarkerStyle (kOpenCircle);

  luPad->cd ();
  luPad->SetLogx ();
  luPad->SetLogy ();
  TH1D* h = new TH1D ("h", "", nPtchBins[iPtZ], pTchBins[iPtZ]);
  h->GetYaxis ()->SetRangeUser (3e-5, 200);
  h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
  h->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T} d#Delta#phi");
  h->GetXaxis ()->SetMoreLogLabels ();
  h->GetXaxis ()->SetTitleFont (43);
  h->GetYaxis ()->SetTitleFont (43);
  h->GetXaxis ()->SetLabelFont (43);
  h->GetYaxis ()->SetLabelFont (43);
  h->GetXaxis ()->SetTitleSize (26);
  h->GetYaxis ()->SetTitleSize (26);
  h->GetXaxis ()->SetLabelSize (26);
  h->GetYaxis ()->SetLabelSize (26);
  h->GetXaxis ()->SetTitleOffset (2.4);
  h->GetYaxis ()->SetTitleOffset (1.8);
  //h->GetXaxis ()->SetLabelOffset (1);
  //h->GetYaxis ()->SetLabelOffset (1);

  h->Draw ("");
  h_ee->Draw ("E1 same");
  h_mumu->Draw ("E1 same");
  h_ee_bkg->Draw ("hist ][ same");
  h_mumu_bkg->Draw ("hist ][ same");

  myText             (-0.090+0.30,   -0.55+0.88-0.01, kBlack, "ee", 0.056);
  myText             (-0.040+0.30,   -0.55+0.88-0.01, kBlack, "#mu#mu", 0.056);
  myText             (-0.050+0.38,   -0.55+0.78-0.01, kBlack, "Z-tagged", 0.054);
  myText             (-0.050+0.38,   -0.55+0.68-0.01, kBlack, "Bkg.", 0.054);
  myMarkerTextNoLine (-0.050+0.34,   -0.55+0.78+0.0001, kBlue, markerStyles[0], "", 1.3 * 1.5, 0.036);
  myMarkerTextNoLine (-0.100+0.34,   -0.55+0.78+0.0001, kRed, markerStyles[1], "", 1.3 * 1.5, 0.036);
  myLineText         (-0.050+0.3466, -0.55+0.68, kBlue, 2, "", 1.0, 0.024);
  myLineText         (-0.100+0.3466, -0.55+0.68, kRed, 2, "", 1.0, 0.024);

  myText             (0.65, 0.86, kBlack, "Pb+Pb, 30-80%", 0.054);
  myText             (0.65, 0.76, kBlack, "#it{p}_{T}^{Z} > 60 GeV", 0.054);


  ldPad->cd ();
  ldPad->SetLogx ();
  h = new TH1D ("h", "", nPtchBins[iPtZ], pTchBins[iPtZ]);
  h->GetYaxis ()->SetRangeUser (0, 2);
  h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
  h->GetYaxis ()->SetTitle ("Electrons / muons");
  h->GetXaxis ()->SetMoreLogLabels ();
  h->GetXaxis ()->SetTitleFont (43);
  h->GetYaxis ()->SetTitleFont (43);
  h->GetXaxis ()->SetLabelFont (43);
  h->GetYaxis ()->SetLabelFont (43);
  h->GetXaxis ()->SetTitleSize (26);
  h->GetYaxis ()->SetTitleSize (26);
  h->GetXaxis ()->SetLabelSize (26);
  h->GetYaxis ()->SetLabelSize (26);
  h->GetXaxis ()->SetTitleOffset (2.4);
  h->GetYaxis ()->SetTitleOffset (1.8);
  //h->GetXaxis ()->SetLabelOffset (1);
  //h->GetYaxis ()->SetLabelOffset (1);
  h->Draw ();

  TH1D* h_rat = (TH1D*) h_ee->Clone ("h_rat");
  h_rat->Divide (h_mumu);
  h_rat->SetMarkerStyle (kOpenCircle);
  h_rat->SetMarkerColor (kBlack);
  h_rat->SetLineColor (kBlack);
  TH1D* h_rat_bkg = (TH1D*) h_ee_bkg->Clone ("h_rat_bkg");
  h_rat_bkg->Divide (h_mumu_bkg);
  h_rat_bkg->SetLineStyle (2);
  h_rat_bkg->SetLineColor (kBlack);
  h_rat->Draw ("E1 same");
  h_rat_bkg->Draw ("hist ][ same");

  TLine* l = new TLine (pTchBins[iPtZ][0], 1, pTchBins[iPtZ][nPtchBins[iPtZ]], 1);
  l->SetLineColor (kPink-8);
  l->SetLineStyle (2);
  l->Draw ("same");


  ruPad->cd ();
  ruPad->SetLogx ();
  ruPad->SetLogy ();
  h = new TH1D ("h", "", nPtchBins[iPtZ], pTchBins[iPtZ]);
  h->GetYaxis ()->SetRangeUser (3e-5, 20);
  h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
  h->GetYaxis ()->SetTitle ("d^{2}Y_{sub} / d#it{p}_{T} d#Delta#phi");
  h->GetXaxis ()->SetMoreLogLabels ();
  h->GetXaxis ()->SetTitleFont (43);
  h->GetYaxis ()->SetTitleFont (43);
  h->GetXaxis ()->SetLabelFont (43);
  h->GetYaxis ()->SetLabelFont (43);
  h->GetXaxis ()->SetTitleSize (26);
  h->GetYaxis ()->SetTitleSize (26);
  h->GetXaxis ()->SetLabelSize (26);
  h->GetYaxis ()->SetLabelSize (26);
  h->GetXaxis ()->SetTitleOffset (2.4);
  h->GetYaxis ()->SetTitleOffset (1.8);
  //h->GetXaxis ()->SetLabelOffset (1);
  //h->GetYaxis ()->SetLabelOffset (1);
  h->Draw ("");

  h_ee_sub->Draw ("E1 same");
  h_mumu_sub->Draw ("E1 same");

  rdPad->cd ();
  rdPad->SetLogx ();
  h = new TH1D ("h", "", nPtchBins[iPtZ], pTchBins[iPtZ]);
  h->GetYaxis ()->SetRangeUser (0, 2);
  h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
  h->GetYaxis ()->SetTitle ("d^{2}Y_{sub} / d#it{p}_{T} d#Delta#phi");
  h->GetXaxis ()->SetMoreLogLabels ();
  h->GetXaxis ()->SetTitleFont (43);
  h->GetYaxis ()->SetTitleFont (43);
  h->GetXaxis ()->SetLabelFont (43);
  h->GetYaxis ()->SetLabelFont (43);
  h->GetXaxis ()->SetTitleSize (26);
  h->GetYaxis ()->SetTitleSize (26);
  h->GetXaxis ()->SetLabelSize (26);
  h->GetYaxis ()->SetLabelSize (26);
  h->GetXaxis ()->SetTitleOffset (2.4);
  h->GetYaxis ()->SetTitleOffset (1.8);
  //h->GetXaxis ()->SetLabelOffset (1);
  //h->GetYaxis ()->SetLabelOffset (1);
  h->Draw ("");

  TH1D* h_rat_sub = (TH1D*) h_ee_sub->Clone ("h_rat_sub");
  h_rat_sub->Divide (h_mumu_sub);
  h_rat_sub->SetMarkerStyle (kOpenCircle);
  h_rat_sub->SetMarkerColor (kBlack);
  h_rat_sub->SetLineColor (kBlack);
  h_rat_sub->Draw ("E1 same");

  l = new TLine (pTchBins[iPtZ][0], 1, pTchBins[iPtZ][nPtchBins[iPtZ]], 1);
  l->SetLineColor (kPink-8);
  l->SetLineStyle (2);
  l->Draw ("same");
  
}



#endif
