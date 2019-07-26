#ifndef __Run_C__
#define __Run_C__

#include "PhysicsAnalysis.h"
#include "FullAnalysis.h"
#include "DataAnalysis.h"
#include "MCAnalysis.h"
#include "MinbiasAnalysis.h"
#include "TruthAnalysis.h"
#include "MixedAnalysis.h"

#include "Systematic.h"

const bool doSys = true;

// nominal analyses
FullAnalysis* data = nullptr;
MCAnalysis* mc = nullptr;
MinbiasAnalysis* bkg = nullptr;
TruthAnalysis* truth = nullptr;

//MixedAnalysis* data_mixed = nullptr;
//PhysicsAnalysis* data_samesign = nullptr;
//FullAnalysis* data_dr10 = nullptr, *data_dr20 = nullptr, *data_dr02 = nullptr;
//MCAnalysis* mc_dr10 = nullptr, *mc_dr20 = nullptr;
//MinbiasAnalysis* bkg_nch = nullptr;

// master systematics objects
Systematic* combSys = nullptr;
Systematic* bkgSys = nullptr;
Systematic* trkSys = nullptr;
Systematic* electronPtSys = nullptr;
Systematic* muonPtSys = nullptr;
Systematic* electronLHMedSys = nullptr;
Systematic* muonLooseSys = nullptr;
//Systematic* trkEffSys = nullptr;

// variations for systematics
PhysicsAnalysis* data_trackHItight = nullptr;
MinbiasAnalysis* bkg_trackHItight = nullptr;
PhysicsAnalysis* data_electronPtUp = nullptr, *data_electronPtDown = nullptr;
PhysicsAnalysis* data_muonPtUp = nullptr, *data_muonPtDown = nullptr;
PhysicsAnalysis* data_electronLHMedium = nullptr;
PhysicsAnalysis* data_muonLoose = nullptr;

MinbiasAnalysis* bkg_statUpVar = nullptr, *bkg_statDownVar = nullptr;
PhysicsAnalysis* data_bkgStatUpVar = nullptr, *data_bkgStatDownVar = nullptr;
//PhysicsAnalysis* data_trackEffStatUpVar = nullptr, *data_trackEffStatDownVar = nullptr;

void Run () {

  //data = new FullAnalysis ("data_preleptfix", "PreLeptonFix");
  data    = new FullAnalysis ("data", "DataAnalysis/Nominal");
  mc      = new MCAnalysis ();
  bkg     = new MinbiasAnalysis ();
  truth   = new TruthAnalysis ();

  //bkg_nch = new MinbiasAnalysis ("bkg_nch", "Nominal/NchWeighted");

  //data_mixed = new MixedAnalysis ("data_mixed", "Nominal/mixed");
  //data_samesign         = new PhysicsAnalysis ("data_samesign", "Variations/SameSignLeptons");

  //data_dr10 = new DataAnalysis ("data_dr10", "Variations/TracksDR10");
  //mc_dr10   = new MCAnalysis ("mc_dr10", "Variations/TracksDR10");
  //data_dr20 = new DataAnalysis ("data_dr20", "Variations/TracksDR20");
  //mc_dr20   = new MCAnalysis ("mc_dr20", "Variations/TracksDR20");

  //data_dr02 = new DataAnalysis ("data_dr02", "Variations/TracksDR02");

  if (doSys) {
    data_trackHItight       = new PhysicsAnalysis ("data_trackHITightVar", "Variations/TrackHITightWPVariation", true);
    bkg_trackHItight        = new MinbiasAnalysis ("bkg_trackHITightVar", "Variations/TrackHITightWPVariation", true);
    data_electronPtUp       = new PhysicsAnalysis ("data_electronPtUpVar", "Variations/ElectronPtUpVariation");
    data_electronPtDown     = new PhysicsAnalysis ("data_electronPtDownVar", "Variations/ElectronPtDownVariation");
    data_muonPtUp           = new PhysicsAnalysis ("data_muonPtUpVar", "Variations/MuonPtUpVariation");
    data_muonPtDown         = new PhysicsAnalysis ("data_muonPtDownVar", "Variations/MuonPtDownVariation");
    data_electronLHMedium   = new PhysicsAnalysis ("data_electronLHMediumVar", "Variations/ElectronLHMediumWPVariation");
    data_muonLoose          = new PhysicsAnalysis ("data_muonLooseVar", "Variations/MuonLooseWPVariation");
    data_bkgStatUpVar       = new PhysicsAnalysis ("data_bkgStatUpVar", "Variations/BkgStatUpVariation");
    data_bkgStatDownVar     = new PhysicsAnalysis ("data_bkgStatDownVar", "Variations/BkgStatDownVariation");
    bkg_statUpVar           = new MinbiasAnalysis ("bkg_statUpVar");
    bkg_statDownVar         = new MinbiasAnalysis ("bkg_statDownVar");
    //data_trkEffUpStatVar    = new PhysicsAnalysis ("data_trkEffUpStatVar", "Variations/TrkEffStatUpVariation");
    //data_trkEffDownStatVar  = new PhysicsAnalysis ("data_trkEffDownStatVar", "Variations/TrkEffStatDownVariation");
  }




  //data->Execute ();
  //mc->Execute ();
  //bkg->Execute ();
  //truth->Execute ();

  //data_mixed->Execute ();
  //data_samesign->Execute ();

  //data_dr10->Execute ();
  //mc_dr10->Execute ();
  //data_dr20->Execute ();
  //mc_dr20->Execute ();
  //data_dr02->Execute ();
  //bkg_nch->Execute ();

  if (doSys) {
    //data_trackHItight->Execute ();
    //bkg_trackHItight->Execute ();
    //data_electronPtUp->Execute ();
    //data_electronPtDown->Execute ();
    //data_muonPtUp->Execute ();
    //data_muonPtDown->Execute ();
    //data_electronLHMedium->Execute ();
    //data_muonLoose->Execute ();
    //data_trkEffStatUpVar->Execute ();
    //data_trkEffStatDownVar->Execute ();
  }



  data->LoadHists ();
  mc->LoadHists ();
  bkg->LoadHists ();
  truth->LoadHists ();

  //data_mixed->LoadHists ();
  //data_samesign->LoadHists ();

  //data_dr10->LoadHists ();
  //mc_dr10->LoadHists ();
  //data_dr20->LoadHists ();
  //mc_dr20->LoadHists ();
  //data_dr02->LoadHists ();

  //bkg_nch->LoadHists ();

  if (doSys) {
    data_bkgStatUpVar->CopyAnalysis (data);
    data_bkgStatDownVar->CopyAnalysis (data);
    bkg_statUpVar->CopyAnalysis (bkg);
    bkg_statUpVar->ConvertToStatVariation (true, 1);
    bkg_statDownVar->CopyAnalysis (bkg);
    bkg_statDownVar->ConvertToStatVariation (false, 1);
    data_bkgStatUpVar->SubtractBackground (bkg_statUpVar);
    data_bkgStatDownVar->SubtractBackground (bkg_statDownVar);
    data_bkgStatUpVar->CalculateIAA ();
    data_bkgStatUpVar->CalculateICP ();
    data_bkgStatDownVar->CalculateIAA ();
    data_bkgStatDownVar->CalculateICP ();
  }

  //data->SubtractSameSigns (data_samesign);
  data->SubtractBackground (bkg);
  data->CalculateIAA ();
  data->CalculateICP ();

  data->CalculateZPtDistRatio (mc);
  data->CalculateZEtaDistRatio ();
  data->CalculateZYDistRatio ();
  data->CalculateZMassSpectraRatio (mc);

  if (doSys) {
    data_trackHItight->LoadHists ();
    bkg_trackHItight->LoadHists ();
    data_electronPtUp->LoadHists ();
    data_electronPtDown->LoadHists ();
    data_muonPtUp->LoadHists ();
    data_muonPtDown->LoadHists ();
    data_electronLHMedium->LoadHists ();
    data_muonLoose->LoadHists ();

    data_trackHItight->SubtractBackground (bkg_trackHItight);
    data_electronPtUp->SubtractBackground (bkg);
    data_electronPtDown->SubtractBackground (bkg);
    data_muonPtUp->SubtractBackground (bkg);
    data_muonPtDown->SubtractBackground (bkg);
    data_electronLHMedium->SubtractBackground (bkg);
    data_muonLoose->SubtractBackground (bkg);


    //trkEffSys = new Systematic (data, "trkEffSys", "Tracking Efficiency");
    //trkEffSys->AddVariation (data_trkEffStatUpVar, -1);
    //trkEffSys->AddVariation (data_trkEffStatDownVar, 1);

    bkgSys = new Systematic (data, "bkgSys", "Background");
    bkgSys->AddVariation (data_bkgStatUpVar, -1);
    bkgSys->AddVariation (data_bkgStatDownVar, 1);
    bkgSys->AddVariations ();

    trkSys = new Systematic (data, "trkSys", "Tracks");
    trkSys->AddVariation (data_trackHItight);
    trkSys->AddVariations ();

    electronPtSys = new Systematic (data, "electronPtSys", "Electron ES");
    electronPtSys->AddVariation (data_electronPtUp);
    electronPtSys->AddVariation (data_electronPtDown);
    electronPtSys->AddVariations ();

    muonPtSys = new Systematic (data, "muonPtSys", "Muon ES");
    muonPtSys->AddVariation (data_muonPtUp);
    muonPtSys->AddVariation (data_muonPtDown);
    muonPtSys->AddVariations ();

    electronLHMedSys = new Systematic (data, "electronLHMedium", "Electron quality");
    electronLHMedSys->AddVariation (data_electronLHMedium);
    electronLHMedSys->AddVariations ();

    muonLooseSys = new Systematic (data, "muonLoose", "Muon quality");
    muonLooseSys->AddVariation (data_muonLoose);
    muonLooseSys->AddVariations ();

    combSys = new Systematic (data, "combSys", "Total");
    //combSys->AddSystematic (trkEffSys);
    combSys->AddSystematic (trkSys);
    combSys->AddSystematic (bkgSys);
    combSys->AddSystematic (electronLHMedSys);
    combSys->AddSystematic (muonLooseSys);
    combSys->AddSystematic (electronPtSys);
    combSys->AddSystematic (muonPtSys);
    combSys->AddSystematics ();
  }

  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");

}

#endif
