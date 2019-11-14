#ifndef __Run_C__
#define __Run_C__

#include "PhysicsAnalysis.h"
#include "FullAnalysis.h"
#include "MCAnalysis.h"
#include "MCClosureAnalysis.h"
#include "MixedMCAnalysis.h"
#include "MinbiasAnalysis.h"
#include "TruthAnalysis.h"

#include "Systematic.h"
#include "MixingVariation.h"
#include "TrackIDSystematic.h"
#include "ReweightingSystematic.h"
#include "SystematicFits.h"
#include "ElectronSystematicTable.h"

#include "PlotHybridModel.h"

const bool doSys = false;

// nominal analyses
FullAnalysis* data18 = nullptr;
//FullAnalysis* data15 = nullptr;
MinbiasAnalysis* bkg18 = nullptr;
//MinbiasAnalysis* bkg15 = nullptr;
MCAnalysis* mc = nullptr;
MixedMCAnalysis* mc_mixed = nullptr;
MCClosureAnalysis* mc_closure = nullptr;
//MCAnalysis* mc_bkg = nullptr;
MinbiasAnalysis* mc_bkg = nullptr;
TruthAnalysis* truth = nullptr;

// master systematics objects
Systematic* combSys = nullptr;
Systematic* bkgStatSys = nullptr;
Systematic* bkgMixSys = nullptr;
TrackIDSystematic* trkSys = nullptr;
ReweightingSystematic* trkEffSysWeights = nullptr;
Systematic* trkEffSys = nullptr;
ReweightingSystematic* trkPurSysWeights = nullptr;
Systematic* trkPurSys = nullptr;
Systematic* leptonRejSys = nullptr;
ElectronSystematicTable* electronPtSysTable = nullptr;
Systematic* electronPtSys = nullptr;
SystematicFits* data_muonPtUpSysFits = nullptr, *data_muonPtDownSysFits = nullptr;
Systematic* muonPtSys = nullptr;
//SystematicFits* bkg_muonPtUpSysFits = nullptr, *bkg_muonPtDownSysFits = nullptr;

// variations for systematics
MinbiasAnalysis* bkg_statUpVar = nullptr, *bkg_statDownVar = nullptr;
PhysicsAnalysis* data_bkgStatUpVar = nullptr, *data_bkgStatDownVar = nullptr;
MinbiasAnalysis* bkg_mixUpVar = nullptr, *bkg_mixDownVar = nullptr;
PhysicsAnalysis* data_mixUpVar = nullptr, *data_mixDownVar = nullptr;

PhysicsAnalysis* data_trackHItight = nullptr, *data_trkIDUpVar = nullptr, *data_trkIDDownVar = nullptr;
MinbiasAnalysis* bkg_trackHItight = nullptr, *bkg_trkIDUpVar = nullptr, *bkg_trkIDDownVar = nullptr;

PhysicsAnalysis* data_partComp = nullptr, *data_partCompUpVar = nullptr, *data_partCompDownVar = nullptr;
MinbiasAnalysis* bkg_partComp = nullptr, *bkg_partCompUpVar = nullptr, *bkg_partCompDownVar = nullptr;

PhysicsAnalysis* data_trackPurity = nullptr, *data_trkPurUpVar = nullptr, *data_trkPurDownVar = nullptr;
MinbiasAnalysis* bkg_trackPurity = nullptr, *bkg_trkPurUpVar = nullptr, *bkg_trkPurDownVar = nullptr;

PhysicsAnalysis* data_electronPtUp = nullptr, *data_electronPtDown = nullptr;
//MinbiasAnalysis* bkg_electronPtUp = nullptr, *bkg_electronPtDown = nullptr;
PhysicsAnalysis* data_muonPtUp = nullptr, *data_muonPtDown = nullptr;
//MinbiasAnalysis* bkg_muonPtUp = nullptr, *bkg_muonPtDown = nullptr;

PhysicsAnalysis* data_leptonRejVar = nullptr;


void Run () {
  data18    = new FullAnalysis ("data18");
  //data15  = new FullAnalysis ("data15");
  //data15->is2015Conds = true;
  //data15->useHijingEffs = true;
  bkg18     = new MinbiasAnalysis ();
  //bkg15->is2015Conds = true;
  //bkg15->useHijingEffs = true;
  //
  //mc      = new MCAnalysis ();
  //mc_closure = new MCClosureAnalysis ("mc_closure");
  //mc_bkg  = new MCAnalysis ("mc_bkg");
  //mc_bkg  = new MinbiasAnalysis ("bkg");
  //truth   = new TruthAnalysis ();

  if (doSys) {
    data_bkgStatUpVar       = new PhysicsAnalysis ("data_bkgStatUpVar");
    data_bkgStatDownVar     = new PhysicsAnalysis ("data_bkgStatDownVar");
    bkg_statUpVar           = new MinbiasAnalysis ("bkg_statUpVar");
    bkg_statDownVar         = new MinbiasAnalysis ("bkg_statDownVar");

    data_mixUpVar           = new PhysicsAnalysis ("data_mixUpVar");
    data_mixDownVar         = new PhysicsAnalysis ("data_mixDownVar");
    bkg_mixUpVar            = new MinbiasAnalysis ("bkg_mixUpVar");
    bkg_mixDownVar          = new MinbiasAnalysis ("bkg_mixDownVar");

    data_trackHItight       = new PhysicsAnalysis ("data_trackHITightVar");
    data_trackHItight->useHITight = true;
    //bkg_trackHItight        = new MinbiasAnalysis ("bkg_trackHITightVar");
    //bkg_trackHItight->useHITight = true;
    data_trkIDUpVar         = new PhysicsAnalysis ("data_trkIDUpVar");
    data_trkIDUpVar->useHITight = true;
    bkg_trkIDUpVar          = new MinbiasAnalysis ("bkg_trkIDUpVar");
    bkg_trkIDUpVar->useHITight = true;
    data_trkIDDownVar       = new PhysicsAnalysis ("data_trkIDDownVar");
    data_trkIDDownVar->useHITight = true;
    bkg_trkIDDownVar        = new MinbiasAnalysis ("bkg_trkIDDownVar");
    bkg_trkIDDownVar->useHITight = true;

    data_partComp        = new PhysicsAnalysis ("data_trackEffVar");
    data_partComp->doTrackEffVar = true;
    //bkg_partComp         = new MinbiasAnalysis ("bkg_trackEffVar");
    //bkg_partComp->doTrackEffVar = true;
    //data_partCompUpVar   = new PhysicsAnalysis ("data_partCompUpVar");
    //data_partCompUpVar->doTrackEffVar = true;
    //bkg_partCompUpVar    = new MinbiasAnalysis ("bkg_partCompUpVar");
    //bkg_partCompUpVar->doTrackEffVar = true;
    //data_partCompDownVar   = new PhysicsAnalysis ("data_partCompDownVar");
    //data_partCompDownVar->doTrackEffVar = true;
    //bkg_partCompDownVar    = new MinbiasAnalysis ("bkg_partCompDownVar");
    //bkg_partCompDownVar->doTrackEffVar = true;

    data_trackPurity        = new PhysicsAnalysis ("data_trackPurityVar");
    data_trackPurity->doTrackPurVar = true;
    //bkg_trackPurity         = new MinbiasAnalysis ("bkg_trackPurityVar");
    //bkg_trackPurity->doTrackPurVar = true;
    data_trkPurUpVar   = new PhysicsAnalysis ("data_trkPurUpVar");
    data_trkPurUpVar->doTrackPurVar = true;
    bkg_trkPurUpVar    = new MinbiasAnalysis ("bkg_trkPurUpVar");
    bkg_trkPurUpVar->doTrackPurVar = true;
    data_trkPurDownVar   = new PhysicsAnalysis ("data_trkPurDownVar");
    data_trkPurDownVar->doTrackPurVar = true;
    bkg_trkPurDownVar    = new MinbiasAnalysis ("bkg_trkPurDownVar");
    bkg_trkPurDownVar->doTrackPurVar = true;
    
    data_leptonRejVar       = new PhysicsAnalysis ("data_leptonRejVar");
    data_leptonRejVar->doLeptonRejVar = true;

    data_electronPtUp       = new PhysicsAnalysis ("data_electronPtUpVar");
    data_electronPtDown     = new PhysicsAnalysis ("data_electronPtDownVar");
    //bkg_electronPtUp       = new MinbiasAnalysis ("bkg");//_electronPtUpVar");
    //bkg_electronPtDown     = new MinbiasAnalysis ("bkg");//_electronPtDownVar");

    data_muonPtUp           = new PhysicsAnalysis ("data_muonPtUpVar");
    data_muonPtDown         = new PhysicsAnalysis ("data_muonPtDownVar");
    //bkg_muonPtUp           = new MinbiasAnalysis ("bkg");//_muonPtUpVar");
    //bkg_muonPtDown         = new MinbiasAnalysis ("bkg");//_muonPtDownVar");
  }

  //mc->GenerateWeights ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/PbPb18*Z*.root", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/data18hi.root", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal/eventWeightsFile.root");
  //bkg18->GenerateWeights("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/data18hi_mixed.root", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MinbiasAnalysis/Nominal/eventWeightsFile.root");
  data18->Execute   ("DataAnalysis/Nominal/data18hi.root",   "DataAnalysis/Nominal/data18hi_hists.root");
  //data15->Execute ("DataAnalysis/Nominal/data15hi.root",  "DataAnalysis/Nominal/data15hi_hists.root");

  if (doSys) {
    //data_electronPtUp->Execute      ("DataAnalysis/Variations/ElectronPtUpVariation/data18hi.root",        "DataAnalysis/Variations/ElectronPtUpVariation/data18hi_hists.root");
    //data_electronPtDown->Execute    ("DataAnalysis/Variations/ElectronPtDownVariation/data18hi.root",      "DataAnalysis/Variations/ElectronPtDownVariation/data18hi_hists.root");
    //data_muonPtUp->Execute          ("DataAnalysis/Variations/MuonPtUpVariation/data18hi.root",            "DataAnalysis/Variations/MuonPtUpVariation/data18hi_hists.root");
    //data_muonPtDown->Execute        ("DataAnalysis/Variations/MuonPtDownVariation/data18hi.root",          "DataAnalysis/Variations/MuonPtDownVariation/data18hi_hists.root");
    //data_leptonRejVar->Execute      ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/LeptonRejVariation/data18hi_hists.root");
    //data_trackHItight->Execute      ("DataAnalysis/Variations/TrackHITightWPVariation/data18hi.root",      "DataAnalysis/Variations/TrackHITightWPVariation/data18hi_hists.root");
    //data_partComp->Execute          ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackEffPionsVariation/data18hi_hists.root");
    //data_trackPurity->Execute       ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackPurityVariation/data18hi_hists.root");
    //data_trkPurUpVar->Execute       ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackPurityUpVariation/data18hi_hists.root");
    //data_trkPurDownVar->Execute     ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackPurityDownVariation/data18hi_hists.root");
  }


  data18->LoadHists   ("DataAnalysis/Nominal/v1_2_0.root");
  //data18->LoadHists   ("DataAnalysis/Nominal/data18hi_hists.root");
  //data15->LoadHists   ("DataAnalysis/Nominal/data15hi_hists.root");
  bkg18->LoadHists      ("MinbiasAnalysis/Nominal/v1_2_0.root");
  //bkg18->LoadHists      ("MinbiasAnalysis/Nominal/data18hi_hists.root");
  //bkg15->LoadHists      ("MinbiasAnalysis/Nominal/data15hi_hists.root");
  //mc->LoadHists       ("MCAnalysis/Nominal/savedHists.root");
  //mc_closure->LoadHists ("MCClosureAnalysis/Nominal/savedHists.root");
  //mc_bkg->LoadHists   ("MCAnalysis/Nominal/mc_bkg_hists.root");
  //mc_bkg->LoadHists   ("MinbiasAnalysis/Nominal/mc_bkg_hists.root");
  //mc_mixed->LoadHists   ("MixedMCAnalysis/Nominal/mc_bkg_hists.root");
  //truth->LoadHists  ("TruthAnalysis/Nominal/savedHists.root");

  /*
  data18->LoadHists   ("DataAnalysis/Nominal/data18hi_hists.root");
  //data15->LoadHists   ("DataAnalysis/Nominal/data15hi_hists.root");
  bkg18->LoadHists      ("MinbiasAnalysis/Nominal/data18hi_hists.root");
  //bkg15->LoadHists      ("MinbiasAnalysis/Nominal/data15hi_hists.root");
  mc->LoadHists       ("MCAnalysis/Nominal/savedHists.root");
  truth->LoadHists  ("TruthAnalysis/Nominal/savedHists.root");

  data18->SubtractBackground (bkg18);
  data18->CalculateIAA ();

  //data15->SubtractBackground (bkg15);
  //data15->CalculateIAA ();

  //data18->CalculateZPtDistRatio (mc);
  //data18->CalculateZEtaDistRatio (mc);
  //data18->CalculateZYDistRatio (mc);
  //data18->CalculateZMassSpectraRatio (mc);
  //data15->CalculateZMassSpectraRatio (mc);

  //mc->SubtractBackground (mc_mixed);
  //mc->SubtractBackground (mc_bkg);
  //mc->CalculateIAA ();
  //mc_closure->SubtractBackground (mc_bkg);
  //mc_closure->CalculateIAA ();
  //truth->SubtractBackground ();
  //truth->CalculateIAA ();

  if (doSys) {
    cout << "Initializing systematic objects. " << endl;
    bkgStatSys = new Systematic (data18, "bkgStatSys", "Bkg. Stat.");
    bkgStatSys->cancelIAA = false;
    bkgMixSys = new Systematic (data18, "bkgMixSys", "Mixed event");
    trkSys = new TrackIDSystematic (data18, "trkSys", "Track ID Cuts");
    trkEffSysWeights = new ReweightingSystematic ("trkEffSysWeights");
    trkEffSys = new Systematic (data18, "trkEffSys", "Particle Composition");
    trkEffSys->cancelIAA = false;
    trkPurSysWeights = new ReweightingSystematic ("trkPurSysWeights");
    trkPurSys = new Systematic (data18, "trkPurSys", "Purity Corr.");
    leptonRejSys = new Systematic (data18, "leptonRejSys", "Lepton Rejection");
    electronPtSys = new Systematic (data18, "electronPtSys", "Electron ES");
    electronPtSys->cancelIAA = false;
    data_muonPtUpSysFits = new SystematicFits ("data_muonPtUpSysFit");
    data_muonPtDownSysFits = new SystematicFits ("data_muonPtDownSysFit");
    //bkg_muonPtUpSysFits = new SystematicFits ("bkg_muonPtUpSysFit");
    //bkg_muonPtDownSysFits = new SystematicFits ("bkg_muonPtDownSysFit");
    muonPtSys = new Systematic (data18, "muonPtSys", "Muon ES");
    muonPtSys->cancelIAA = false;
    combSys = new Systematic (data18, "combSys", "Total");


    cout << "Calculating bkg. stat. systematic errors." << endl;
    data_bkgStatUpVar->CopyAnalysis (data18, false);
    data_bkgStatDownVar->CopyAnalysis (data18, false);
    bkg_statUpVar->CopyAnalysis (bkg18, false);
    bkg_statUpVar->ConvertToStatVariation (true, 1);
    bkg_statDownVar->CopyAnalysis (bkg18, false);
    bkg_statDownVar->ConvertToStatVariation (false, 1);
    data_bkgStatUpVar->SubtractBackground (bkg_statUpVar);
    data_bkgStatDownVar->SubtractBackground (bkg_statDownVar);
    data_bkgStatUpVar->CalculateIAA ();
    data_bkgStatDownVar->CalculateIAA ();
    bkgStatSys->AddVariation (data_bkgStatUpVar, true);
    bkgStatSys->AddVariation (data_bkgStatDownVar, true);
    bkgStatSys->AddVariations ();


    cout << "Calculating mixed event systematic errors." << endl;
    bkg_mixUpVar->CopyAnalysis (bkg18, false);
    bkg_mixDownVar->CopyAnalysis (bkg18, false);
    data_mixUpVar->CopyAnalysis (data18, false);
    data_mixDownVar->CopyAnalysis (data18, false);
    ApplyBkgVariation (bkg_mixUpVar, 0.005);
    ApplyBkgVariation (bkg_mixDownVar, -0.005);
    data_mixUpVar->SubtractBackground (bkg_mixUpVar);
    data_mixDownVar->SubtractBackground (bkg_mixDownVar);
    bkgMixSys->AddVariation (data_mixUpVar, true);
    bkgMixSys->AddVariation (data_mixDownVar, true);
    bkgMixSys->AddVariations ();



    cout << "Calculating track ID relative systematic errors." << endl;
    data_trackHItight->LoadHists      ("DataAnalysis/Variations/TrackHITightWPVariation/data18hi_hists.root");
    trkSys->GetRelativeVariation (data18, data_trackHItight);
    //bkg_trackHItight->LoadHists       ("MinbiasAnalysis/Variations/TrackHITightWPVariation/data18hi_hists.root");
    //data_trackHItight->SubtractBackground (bkg_trackHItight);
    cout << "Applying track ID systematic errors." << endl;
    data_trkIDUpVar->CopyAnalysis (data18, false);
    data_trkIDDownVar->CopyAnalysis (data18, false);
    bkg_trkIDUpVar->CopyAnalysis (bkg18, false);
    bkg_trkIDDownVar->CopyAnalysis (bkg18, false);
    data_trkIDUpVar->ApplyRelativeVariation (trkSys->relVar, true); // track yields go up
    data_trkIDDownVar->ApplyRelativeVariation (trkSys->relVar, false); // track yields go down
    bkg_trkIDUpVar->ApplyRelativeVariation (trkSys->relVar, true); // track yields go up
    bkg_trkIDDownVar->ApplyRelativeVariation (trkSys->relVar, false); // track yields go down
    data_trkIDUpVar->SubtractBackground (bkg_trkIDUpVar);
    data_trkIDDownVar->SubtractBackground (bkg_trkIDDownVar);
    trkSys->AddVariation (data_trkIDUpVar, true);
    trkSys->AddVariation (data_trkIDDownVar, true);
    trkSys->AddVariations ();


    cout << "Calculating particle composition systematic errors." << endl;
    data_partComp->LoadHists          ("DataAnalysis/Variations/TrackEffPionsVariation/data18hi_hists.root");
    trkEffSysWeights->GetRelativeVariations (data18, data_partComp);
    data_partCompUpVar->CopyAnalysis (data18, false);
    data_partCompDownVar->CopyAnalysis (data18, false);
    bkg_partCompUpVar->CopyAnalysis (bkg18, false);
    bkg_partCompDownVar->CopyAnalysis (bkg18, false);
    trkEffSysWeights->ApplyRelativeVariations (data_partCompUpVar, true);
    trkEffSysWeights->ApplyRelativeVariations (bkg_partCompUpVar, true);
    trkEffSysWeights->ApplyRelativeVariations (data_partCompDownVar, false);
    trkEffSysWeights->ApplyRelativeVariations (bkg_partCompDownVar, false);
    delete trkEffSysWeights;

    data_partCompUpVar->SubtractBackground (bkg_partCompUpVar);
    data_partCompDownVar->SubtractBackground (bkg_partCompDownVar);
    trkEffSys->AddVariation (data_partCompUpVar);
    trkEffSys->AddVariation (data_partCompDownVar);
    trkEffSys->AddVariations ();


    cout << "Calculating track purity systematic errors." << endl;
    data_trkPurUpVar->Execute         ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackPurityUpVariation/data18hi_hists.root");
    data_trkPurDownVar->Execute       ("DataAnalysis/Nominal/data18hi.root",                                 "DataAnalysis/Variations/TrackPurityDownVariation/data18hi_hists.root");
    data_trkPurUpVar->LoadHists         ("DataAnalysis/Variations/TrackPurityUpVariation/data18hi_hists.root");
    data_trkPurDownVar->LoadHists        ("DataAnalysis/Variations/TrackPurityDownVariation/data18hi_hists.root");
    data_trackPurity->LoadHists       ("DataAnalysis/Variations/TrackPurityVariation/data18hi_hists.root");
    //bkg_trackPurity->LoadHists        ("MinbiasAnalysis/Variations/TrackPurityVariation/data18hi_hists.root");
    //data_trackPurity->SubtractBackground (bkg_trackPurity);
    data_trkPurUpVar->LoadHists       ("DataAnalysis/Variations/TrackPurityUpVariation/data18hi_hists.root");
    trkPurSysWeights = new ReweightingSystematic ("trkPurSysWeights");
    trkPurSysWeights->GetRelativeVariations (data18, data_trkPurUpVar);
    //trkPurSysWeights->ScaleRelativeVariations (0.25);
    data_trkPurUpVar->CopyAnalysis (data18, false);
    bkg_trkPurUpVar->CopyAnalysis (bkg18, false);
    trkPurSysWeights->ApplyRelativeVariations (data_trkPurUpVar, true);
    trkPurSysWeights->ApplyRelativeVariations (bkg_trkPurUpVar, true);
    delete trkPurSysWeights;
    data_trkPurUpVar->SubtractBackground (bkg_trkPurUpVar);

    data_trkPurDownVar->LoadHists       ("DataAnalysis/Variations/TrackPurityDownVariation/data18hi_hists.root");
    trkPurSysWeights = new ReweightingSystematic ("trkPurSysWeights");
    trkPurSysWeights->GetRelativeVariations (data18, data_trkPurDownVar);
    //trkPurSysWeights->ScaleRelativeVariations (0.25);
    data_trkPurDownVar->CopyAnalysis (data18, false);
    bkg_trkPurDownVar->CopyAnalysis (bkg18, false);
    trkPurSysWeights->ApplyRelativeVariations (data_trkPurDownVar, false);
    trkPurSysWeights->ApplyRelativeVariations (bkg_trkPurDownVar, false);
    delete trkPurSysWeights;
    data_trkPurDownVar->SubtractBackground (bkg_trkPurDownVar);

    trkPurSys->AddVariation (data_trkPurUpVar);
    trkPurSys->AddVariation (data_trkPurDownVar);
    trkPurSys->AddVariations ();


    cout << "Calculating lepton rejection systematic errors." << endl;
    data_leptonRejVar->LoadHists      ("DataAnalysis/Variations/LeptonRejVariation/data18hi_hists.root");
    data_leptonRejVar->SubtractBackground (bkg18);
    leptonRejSys->AddVariation (data_leptonRejVar, -1);
    leptonRejSys->AddVariations ();


    cout << "Calculating electron ES systematic errors." << endl;
    ElectronSystematicTable* electronPtSysTable = new ElectronSystematicTable ("data18_electronPtSysTable");
    electronPtSysTable->GetRelativeVariations ("DataAnalysis/Variations/ElectronPtVariation/systematics_electron_scale.root");
    data_electronPtUp->CopyAnalysis (data18, true);
    data_electronPtDown->CopyAnalysis (data18, true);
    electronPtSysTable->ApplyRelativeVariations (data_electronPtUp, true);
    electronPtSysTable->ApplyRelativeVariations (data_electronPtDown, false);
    electronPtSys->AddVariation (data_electronPtUp);
    electronPtSys->AddVariation (data_electronPtDown);
    electronPtSys->AddVariations ();


    cout << "Calculating muon ES systematic errors." << endl;
    data_muonPtUp->LoadHists          ("DataAnalysis/Variations/MuonPtUpVariation/data18hi_hists.root");
    data_muonPtUp->SubtractBackground (bkg18);
    data_muonPtUpSysFits->GetRelativeVariations (data18, data_muonPtUp);
    data_muonPtUp->CopyAnalysis (data18, true);
    data_muonPtUpSysFits->ApplyRelativeVariations (data_muonPtUp, true);

    //bkg_muonPtUp->LoadHists          ("MinbiasAnalysis/Variations/MuonPtUpVariation/data18hi_hists.root");
    //bkg_muonPtUpSysFits->GetRelativeVariations (bkg18, bkg_muonPtUp);
    //bkg_muonPtUp->CopyAnalysis (bkg18, false);
    //bkg_muonPtDown->CopyAnalysis (bkg18, false);
    //bkg_muonPtUpSysFits->ApplyRelativeVariations (bkg_muonPtUp, true);
    //bkg_muonPtUpSysFits->ApplyRelativeVariations (bkg_muonPtDown, false);

    //data_muonPtUp->SubtractBackground (bkg_muonPtUp);
    //data_muonPtDown->SubtractBackground (bkg_muonPtDown);

    data_muonPtDown->LoadHists        ("DataAnalysis/Variations/MuonPtDownVariation/data18hi_hists.root");
    data_muonPtDown->SubtractBackground (bkg18);
    data_muonPtDownSysFits->GetRelativeVariations (data18, data_muonPtDown);
    data_muonPtDown->CopyAnalysis (data18, true);
    //data_muonPtDownSysFits->ApplyRelativeVariations (data_muonPtUp, true);
    data_muonPtDownSysFits->ApplyRelativeVariations (data_muonPtDown, false);
    delete data_muonPtDownSysFits;

    //bkg_muonPtUp->LoadHists          ("MinbiasAnalysis/Variations/MuonPtUpVariation/data18hi_hists.root");
    //bkg_muonPtUpSysFits->GetRelativeVariations (bkg18, bkg_muonPtUp);
    //bkg_muonPtUp->CopyAnalysis (bkg18, false);
    //bkg_muonPtDown->CopyAnalysis (bkg18, false);
    //bkg_muonPtUpSysFits->ApplyRelativeVariations (bkg_muonPtUp, true);
    //bkg_muonPtUpSysFits->ApplyRelativeVariations (bkg_muonPtDown, false);

    //data_muonPtUp->SubtractBackground (bkg_muonPtUp);
    //data_muonPtDown->SubtractBackground (bkg_muonPtDown);
    muonPtSys->AddVariation (data_muonPtUp);
    muonPtSys->AddVariation (data_muonPtDown);
    muonPtSys->AddVariations ();


    cout << "Adding errors in quadrature." << endl;
    combSys->AddSystematic (trkSys);
    combSys->AddSystematic (trkEffSys);
    combSys->AddSystematic (trkPurSys);
    combSys->AddSystematic (leptonRejSys);
    combSys->AddSystematic (electronPtSys);
    combSys->AddSystematic (muonPtSys);
    combSys->AddSystematic (bkgStatSys);
    combSys->AddSystematic (bkgMixSys);
    combSys->AddSystematics ();

    trkSys->SaveGraphs ("Systematics/TrackIDSys.root");
    trkEffSys->SaveGraphs ("Systematics/PartCompSys.root");
    trkPurSys->SaveGraphs ("Systematics/TrackPurSys.root");
    bkgStatSys->SaveGraphs ("Systematics/BkgStatSys.root");
    bkgMixSys->SaveGraphs ("Systematics/MixingSys.root");
    leptonRejSys->SaveGraphs ("Systematics/LeptonRejSys.root");
    electronPtSys->SaveGraphs ("Systematics/ElectronPtSys.root");
    muonPtSys->SaveGraphs ("Systematics/MuonPtSys.root");
    combSys->SaveGraphs ("Systematics/CombinedSys.root"); 

  }
  */


  SetupDirectories ("", "ZTrackAnalysis/");

}




void MakePhysicsPlots () {

  data18->PlotTrkYield (1, 0, 2, 4);
  combSys->PlotTrkYield (1, 1, 2, 4);
  bkg18->PlotTrkYield (1, 0, 2, 4);
  data18->PlotTrkYield (1, 0, 2, 4);

  max_iaa=4.4;
  combSys->PlotIAAdPhi (1, 1, 2, 4);
  data18->PlotIAAdPhi (1, 0, 2, 4);
  //max_icp=3.;
  //combSys->PlotICPdPhi (1, 1, 2, 4);
  //data18->PlotICPdPhi (1, 0, 2, 4);

  combSys->PlotTrkYieldZPt (1, 1, 2);
  data18->PlotTrkYieldZPt (1, 0, 2);

  combSys->PlotTrkYieldZPt (0, 1, 2);
  data18->PlotTrkYieldZPt (0, 0, 2);

  max_iaa = 10;
  combSys->PlotSingleIAAdPtZ (1, 1);
  data18->PlotSingleIAAdPtZ (1, 0);

  combSys->PlotSingleIAAdPtZ (0, 1);
  data18->PlotSingleIAAdPtZ (0, 0);

  max_iaa = 3.2;
  combSys->PlotIAAdPtZ (1, 1);
  data18->PlotIAAdPtZ (1, 0);

  max_iaa = 3.2;
  combSys->PlotIAAdPtZ (0, 1);
  data18->PlotIAAdPtZ (0, 0);

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
  TH1D* h = new TH1D ("h", "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
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
  h = new TH1D ("h", "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
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

  TLine* l = new TLine (ptTrkBins[iPtZ][0], 1, ptTrkBins[iPtZ][nPtTrkBins[iPtZ]], 1);
  l->SetLineColor (kPink-8);
  l->SetLineStyle (2);
  l->Draw ("same");


  ruPad->cd ();
  ruPad->SetLogx ();
  ruPad->SetLogy ();
  h = new TH1D ("h", "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
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
  h = new TH1D ("h", "", nPtTrkBins[iPtZ], ptTrkBins[iPtZ]);
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

  l = new TLine (ptTrkBins[iPtZ][0], 1, ptTrkBins[iPtZ][nPtTrkBins[iPtZ]], 1);
  l->SetLineColor (kPink-8);
  l->SetLineStyle (2);
  l->Draw ("same");
  
}




void ComparePbPbSubYields (const short iSpc = 2, const short iPtZ = nPtZBins-1) {
  TCanvas* c = new TCanvas ("c", "", 1600, 800);
  const double dPadY = 0.5;
  const double uPadY = 1. - dPadY;
  const int axisTextSize = 23;

  TH1D* h1 = nullptr, *h2 = nullptr, *h3 = nullptr, *hrat = nullptr, *eff1 = nullptr, *eff2 = nullptr, *effrat = nullptr, *purrat = nullptr;
  TGraphAsymmErrors* g = nullptr;

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    const char* uPadName = Form ("uPad_%i", iCent);
    const char* dPadName = Form ("dPad_%i", iCent);

    TPad* uPad = new TPad (uPadName, "", (1./(numCentBins))*(iCent), dPadY, (1./(numCentBins))*(iCent+1), 1);
    TPad* dPad = new TPad (dPadName, "", (1./(numCentBins))*(iCent), 0, (1./(numCentBins))*(iCent+1), dPadY);

    uPad->SetTopMargin (0.04);
    uPad->SetBottomMargin (0);
    uPad->SetLeftMargin (0.17);
    uPad->SetRightMargin (0.06);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    dPad->SetLeftMargin (0.17);
    dPad->SetRightMargin (0.06);
    uPad->Draw ();
    dPad->Draw ();


    uPad->cd ();
    uPad->SetLogx ();
    uPad->SetLogy ();

    PhysicsAnalysis* a1, *a2 = nullptr, *a3 = nullptr;

    a1 = data18;
    a2 = data_muonPtUp;
    a3 = data_muonPtDown;

    h1 = a1->h_trk_pt_ptz[2][iPtZ][iCent];
    h2 = a2->h_trk_pt_ptz[2][iPtZ][iCent];
    h3 = a3->h_trk_pt_ptz[2][iPtZ][iCent];

    //float i1 = 0, i2 = 0;
    //for (int ix = 1; ix <= h1->GetNbinsX (); ix++) {
    //  cout << h2->GetBinContent (ix) << endl;
    //  i1 += h1->GetBinContent (ix) * h1->GetBinCenter (ix) * h1->GetBinWidth (ix);
    //}
    //for (int ix = 1; ix <= h1->GetNbinsX (); ix++) {
    //  cout << h2->GetBinContent (ix) << endl;
    //  i2 += h2->GetBinContent (ix) * h2->GetBinCenter (ix) * h2->GetBinWidth (ix);
    //}

    //i1 *= pi/4;
    //i2 *= pi/4;

    //cout << "iPtZ = " << iPtZ << ", iCent = " << iCent << endl;
    //cout << "i1 = " << i1 << endl;
    //cout << "i2 = " << i2 << endl;

    //h3 = a3->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
    //h1 = (TH1D*) a1->h_z_trk_raw_pt[iSpc][iPtZ][1][iCent]->Clone ("h1");
    //h1->Add (a1->h_z_trk_raw_pt[iSpc][iPtZ][2][iCent]);
    //h2 = (TH1D*) a2->h_z_trk_raw_pt[iSpc][iPtZ][1][iCent]->Clone ("h2");
    //h2->Add (a2->h_z_trk_raw_pt[iSpc][iPtZ][2][iCent]);
    //h3 = (TH1D*) a3->h_z_trk_raw_pt[iSpc][iPtZ][1][iCent]->Clone ("h2");
    //h3->Add (a3->h_z_trk_raw_pt[iSpc][iPtZ][2][iCent]);
    //const float min = fmin (h1->GetMinimum (0), h2->GetMinimum (0));
    //const float max = fmax (h1->GetMaximum (),  h2->GetMaximum ());
    const float min = fmin (fmin (h1->GetMinimum (0), h2->GetMinimum (0)), h3->GetMinimum (0));
    const float max = fmax (fmax (h1->GetMaximum (),  h2->GetMaximum ()), h3->GetMaximum (0));

    g = a1->GetTGAE (h1);

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[2]);
    g->SetLineColor (colors[2]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    g->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    //g->GetYaxis ()->SetTitle ("N_{ch}^{total}");
    g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T}d#Delta#phi [GeV^{-1}]");
    //g->GetYaxis ()->SetTitle ("Y / Y_{bkg}");

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->Draw ("AP");


    g = a2->GetTGAE (h2);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[1]);
    g->SetLineColor (colors[1]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    g->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->Draw ("P");

    g = a3->GetTGAE (h3);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[3]);
    g->SetLineColor (colors[3]);

    g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    //g->GetXaxis ()->SetLimits (allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    g->GetYaxis ()->SetRangeUser (min, max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->Draw ("P");


    if (iCent == 0) {
      myText (0.44, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.040/uPadY);
      myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.036/uPadY);
      if (iSpc == 0)
        myText (0.44, 0.80, kBlack, "Z #rightarrow e^{+}e^{-} events", 0.036/uPadY);
      if (iSpc == 1)
        myText (0.44, 0.80, kBlack, "Z #rightarrow #mu^{+}#mu^{-} events", 0.036/uPadY);
    }
    else
      myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.036/uPadY);

    if (iCent == 1)
      myText (0.50, 0.88, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.036/uPadY);
    else if (iCent == 2) {
      if (iPtZ == nPtZBins-1)
        myText (0.50, 0.88, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.036/uPadY);
      else
        myText (0.50, 0.88, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.036/uPadY);
    }
    else if (iCent == 3) {
      //myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "Z#rightarrowee", 1.4, 0.036/uPadY);
      //myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "Z#rightarrow#mu#mu", 1.4, 0.036/uPadY);
      //myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "No correction", 1.4, 0.036/uPadY);
      //myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "Eff. corrected", 1.4, 0.036/uPadY);
      //myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "Electrons", 1.4, 0.036/uPadY);
      //myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "Muons", 1.4, 0.04/uPadY);
      myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "Central values", 1.4, 0.04/uPadY);
      myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "Muons Up", 1.4, 0.04/uPadY);
      myMarkerTextNoLine (0.50, 0.74, colors[3], kOpenSquare, "Muons Down", 1.4, 0.04/uPadY);
      //myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "#varepsilon_{Z} = 1", 1.4, 0.04/uPadY);
      //myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "#varepsilon_{Z} = #varepsilon_{trig} #times #varepsilon_{reco}", 1.4, 0.04/uPadY);
    }


    dPad->cd ();
    dPad->SetLogx ();

    hrat = (TH1D*) h2->Clone ("hrat");
    hrat->Reset ();
    //hrat->Divide (h1);
    for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
      const float y1 = h1->GetBinContent (ix);
      const float y1e = h1->GetBinError (ix);
      const float y2 = h2->GetBinContent (ix);
      const float y2e = h2->GetBinError (ix);
      hrat->SetBinContent (ix, y1/y2);
      hrat->SetBinError (ix, fabs(y1/y2)*sqrt (fabs (pow (y1e/y1,2) + pow (y2e/y2,2) - 2.*y1e*y1e/(y1*y2))));
    }
    //for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
    //  const float passes = h2->GetBinContent (ix);
    //  const float trials = h1->GetBinContent (ix);
    //  hrat->SetBinContent (ix, (passes+1) / (trials+2));
    //  hrat->SetBinError (ix, sqrt ((passes+1)*(passes+2)/((trials+2)*(trials+3)) - pow (passes+1, 2) / pow (trials+2, 2)));

    //  //hrat->SetBinError (ix, sqrt ((hrat->GetBinContent (ix)) * (1-hrat->GetBinContent (ix)) / h1->GetBinContent (ix)));
    //}

    //a1->LoadTrackingEfficiencies (true);
    //a2->LoadTrackingEfficiencies (true);
    //a1->LoadTrackingPurities (true);
    //a2->LoadTrackingPurities (true);

    //eff1 = (TH1D*) a1->h2_num_trk_effs[iCent]->ProjectionY ("1");
    //effrat = (TH1D*) a1->h2_den_trk_effs[iCent]->ProjectionY ("2");
    //eff1->Divide (effrat);
    //delete effrat;
    //eff2 = (TH1D*) a2->h2_num_trk_effs[iCent]->ProjectionY ("3");
    //effrat = (TH1D*) a2->h2_den_trk_effs[iCent]->ProjectionY ("4");
    //eff2->Divide (effrat);
    //delete effrat;

    //effrat = (TH1D*) eff2->Clone ("effrat");
    //effrat->Divide (eff1);
    //delete eff1, eff2;

    //eff1 = (TH1D*) a1->h2_num_trk_purs[iCent]->ProjectionY ("1");
    //eff1->Divide ((TH1D*) a1->h2_den_trk_purs[iCent]->ProjectionY ("2"));
    //eff2 = (TH1D*) a2->h2_num_trk_purs[iCent]->ProjectionY ("3");
    //eff2->Divide ((TH1D*) a2->h2_den_trk_purs[iCent]->ProjectionY ("4"));

    //purrat = (TH1D*) eff2->Clone ("purrat");
    //purrat->Divide (eff1);
    //delete eff1, eff2;

    //const int bin1 = effrat->FindBin (hrat->GetBinCenter (1));
    //const int bin2 = purrat->FindBin (hrat->GetBinCenter (1));
    //for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
    //  hrat->SetBinContent (ix, hrat->GetBinContent (ix) * purrat->GetBinContent (ix+bin2-1) / effrat->GetBinContent (ix+bin1-1));
    //  hrat->SetBinError (ix, hrat->GetBinError (ix) * purrat->GetBinContent (ix+bin2-1) / effrat->GetBinContent (ix+bin1-1));
    //}
    

    g = make_graph (hrat);
    delete hrat, effrat;

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[1]);
    g->SetLineColor (colors[1]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    g->GetYaxis ()->SetRangeUser (0.95, 1.05);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    //g->GetYaxis ()->SetTitle ("Z#rightarrow#mu#mu / Z#rightarrowee");
    g->GetYaxis ()->SetTitle ("Variation / Nominal");
    //g->GetYaxis ()->SetTitle ("Eff. Corrected / Uncorrected");
    //g->GetYaxis ()->SetTitle ("New ES / Old ES");
    //g->GetYaxis ()->SetTitle ("Unfolded / not unfolded");

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->GetYaxis ()->CenterTitle ();

    g->Draw ("AP");

    //TF1* fit1 = new TF1 ("fit1", "[0]", allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    //TF1* fit1 = new TF1 ("fit1", "[0]+[1]*log(x)", allXHZBins[0], allXHZBins[maxNXHZBins]);
    TF1* fit1 = new TF1 ("fit1", "[0]+[1]*log(x)", allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    fit1->SetParameter (0, 1);
    fit1->SetParameter (1, 0);
    g->Fit (fit1, "RQN0");

    fit1->SetLineColor (colors[1]);
    fit1->SetLineStyle (2);
    fit1->SetLineWidth (1);
    fit1->Draw ("same");

    //TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]", allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    //TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]-[1]*log(x)", allXHZBins[0], allXHZBins[maxNXHZBins]);
    TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]-[1]*log(x)", allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    inv_fit1->SetParameter (0, 1/fit1->GetParameter (0));
    inv_fit1->SetParameter (1, fit1->GetParameter (1)/fit1->GetParameter (0));

    inv_fit1->SetLineColor (colors[1]);
    inv_fit1->SetLineStyle (2);
    inv_fit1->SetLineWidth (1);
    inv_fit1->Draw ("same");

    cout << "chi2/ndf = " << fit1->GetChisquare () << " / " << fit1->GetNDF () << endl;


    hrat = (TH1D*) h3->Clone ("hrat");
    hrat->Reset ();
    for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
      const float y1 = h1->GetBinContent (ix);
      const float y1e = h1->GetBinError (ix);
      const float y3 = h3->GetBinContent (ix);
      const float y3e = h3->GetBinError (ix);
      hrat->SetBinContent (ix, y1/y3);
      hrat->SetBinError (ix, (y1/y3)*sqrt (fabs (pow (y1e/y1,2) + pow (y3e/y3,2) - 2.*y1e*y1e/(y1*y3))));
    }

    //for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
    //  hrat->SetBinError (ix, sqrt ((hrat->GetBinContent (ix)) * (1-hrat->GetBinContent (ix)) / h1->GetBinContent (ix)));
    //}

    //a1->LoadTrackingEfficiencies (true);
    //a2->LoadTrackingEfficiencies (true);
    //a1->LoadTrackingPurities (true);
    //a2->LoadTrackingPurities (true);

    //eff1 = (TH1D*) a1->h2_num_trk_effs[iCent]->ProjectionY ("1");
    //effrat = (TH1D*) a1->h2_den_trk_effs[iCent]->ProjectionY ("2");
    //eff1->Divide (effrat);
    //delete effrat;
    //eff2 = (TH1D*) a2->h2_num_trk_effs[iCent]->ProjectionY ("3");
    //effrat = (TH1D*) a2->h2_den_trk_effs[iCent]->ProjectionY ("4");
    //eff2->Divide (effrat);
    //delete effrat;

    //effrat = (TH1D*) eff2->Clone ("effrat");
    //effrat->Divide (eff1);
    //delete eff1, eff2;

    //eff1 = (TH1D*) a1->h2_num_trk_purs[iCent]->ProjectionY ("1");
    //eff1->Divide ((TH1D*) a1->h2_den_trk_purs[iCent]->ProjectionY ("2"));
    //eff2 = (TH1D*) a2->h2_num_trk_purs[iCent]->ProjectionY ("3");
    //eff2->Divide ((TH1D*) a2->h2_den_trk_purs[iCent]->ProjectionY ("4"));

    //purrat = (TH1D*) eff2->Clone ("purrat");
    //purrat->Divide (eff1);
    //delete eff1, eff2;

    //const int bin1 = effrat->FindBin (hrat->GetBinCenter (1));
    //const int bin2 = purrat->FindBin (hrat->GetBinCenter (1));
    //for (int iy = 1; iy <= hrat->GetNbinsX (); iy++)
    //  hrat->SetBinContent (iy, hrat->GetBinContent (iy) * purrat->GetBinContent (iy+bin2-1));// / effrat->GetBinContent (iy+bin1-1));

    g = make_graph (hrat);
    delete hrat;//, effrat;

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[3]);
    g->SetLineColor (colors[3]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    g->GetYaxis ()->SetRangeUser (0.95, 1.05);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->Draw ("P");

    //TF1* fit2 = new TF1 ("fit2", "[0]+[1]*log(x)", allXHZBins[0], allXHZBins[maxNXHZBins]);
    TF1* fit2 = new TF1 ("fit2", "[0]+[1]*log(x)", allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    fit2->SetParameter (0, 1);
    fit2->SetParameter (1, 0);
    g->Fit (fit2, "RQN0");

    fit2->SetLineColor (colors[3]);
    fit2->SetLineStyle (2);
    fit2->SetLineWidth (1);
    fit2->Draw ("same");

    cout << "chi2/ndf = " << fit2->GetChisquare () << " / " << fit2->GetNDF () << endl;

    //TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]-[1]*log(x)", allXHZBins[0], allXHZBins[maxNPtTrkBins]);
    TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]-[1]*log(x)", allPtTrkBins[0], allPtTrkBins[maxNPtTrkBins]);
    inv_fit2->SetParameter (0, 1./fit2->GetParameter (0));
    inv_fit2->SetParameter (1, fit2->GetParameter (1)/fit2->GetParameter (0));

    inv_fit2->SetLineColor (colors[3]);
    inv_fit2->SetLineStyle (2);
    inv_fit2->SetLineWidth (1);
    inv_fit2->Draw ("same");

    TLine* l = new TLine (allPtTrkBins[0], 1, allPtTrkBins[maxNPtTrkBins], 1);
    l->SetLineColor (kBlack);
    l->SetLineStyle (2);
    l->Draw ("same");

    //myText (0.22, 0.3, kBlack, Form ("Avg. = %s", FormatMeasurement (fit1->GetParameter (0), fit1->GetParError (0), 2)), 0.03/dPadY);

    delete h1, h2;
  }

}

#endif
