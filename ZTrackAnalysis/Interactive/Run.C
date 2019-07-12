#ifndef __Run_C__
#define __Run_C__

#include "PhysicsAnalysis.h"
#include "FullAnalysis.h"
#include "DataAnalysis.h"
#include "MCAnalysis.h"
#include "MinbiasAnalysis.h"
#include "TruthAnalysis.h"

#include "Systematic.h"
//#include "WPVariation.h"

// nominal analyses
FullAnalysis* data = nullptr;//, *data_hm = nullptr;
MCAnalysis* mc = nullptr;
MinbiasAnalysis* bkg = nullptr;
TruthAnalysis* truth = nullptr;

// master systematics objects
Systematic* combSys = nullptr;
Systematic* bkgSys = nullptr;
Systematic* trkSys = nullptr;
Systematic* electronPtSys = nullptr;
Systematic* muonPtSys = nullptr;
Systematic* electronLHMedSys = nullptr;
Systematic* muonLooseSys = nullptr;

FullAnalysis* data_dr10 = nullptr, *data_dr20 = nullptr, *data_dr02 = nullptr;
MCAnalysis* mc_dr10 = nullptr, *mc_dr20 = nullptr;

// variations for systematics
PhysicsAnalysis* data_trackHItight = nullptr;
MinbiasAnalysis* bkg_trackHItight = nullptr;
PhysicsAnalysis* data_electronPtUp = nullptr, *data_electronPtDown = nullptr;
PhysicsAnalysis* data_muonPtUp = nullptr, *data_muonPtDown = nullptr;
PhysicsAnalysis* data_electronLHMedium = nullptr;
PhysicsAnalysis* data_muonLoose = nullptr;

MinbiasAnalysis* bkg_statUpVar = nullptr, *bkg_statDownVar = nullptr;
PhysicsAnalysis* data_bkgStatUpVar = nullptr, *data_bkgStatDownVar = nullptr;

void Run () {

  //data = new FullAnalysis ("data_preleptfix", "PreLeptonFix");
  data    = new FullAnalysis ("data", "DataAnalysis/Nominal");
  //data_hm = new DataAnalysis ("data_hm");
  mc      = new MCAnalysis ();
  bkg     = new MinbiasAnalysis ();
  truth   = new TruthAnalysis ();

  data_dr10 = new DataAnalysis ("data_dr10", "Variations/TracksDR10");
  mc_dr10   = new MCAnalysis ("mc_dr10", "Variations/TracksDR10");
  data_dr20 = new DataAnalysis ("data_dr20", "Variations/TracksDR20");
  mc_dr20   = new MCAnalysis ("mc_dr20", "Variations/TracksDR20");

  data_dr02 = new DataAnalysis ("data_dr02", "Variations/TracksDR02");

  data_trackHItight     = new PhysicsAnalysis ("data_trackHITightVar", "Variations/TrackHITightWPVariation", true);
  bkg_trackHItight      = new MinbiasAnalysis ("bkg_trackHITightVar", "Variations/TrackHITightWPVariation", true);
  data_electronPtUp     = new PhysicsAnalysis ("data_electronPtUpVar", "Variations/ElectronPtUpVariation");
  data_electronPtDown   = new PhysicsAnalysis ("data_electronPtDownVar", "Variations/ElectronPtDownVariation");
  data_muonPtUp         = new PhysicsAnalysis ("data_muonPtUpVar", "Variations/MuonPtUpVariation");
  data_muonPtDown       = new PhysicsAnalysis ("data_muonPtDownVar", "Variations/MuonPtDownVariation");
  //data_electronLHMedium = new PhysicsAnalysis ("data_electronLHMediumVar", "Variations/ElectronLHMediumVariation");
  //data_muonLoose        = new PhysicsAnalysis ("data_muonLooseVar", "Variations/MuonLooseVariation");
  data_bkgStatUpVar     = new PhysicsAnalysis ("data_bkgStatUpVar", "Variations/BkgStatUpVariation");
  data_bkgStatDownVar   = new PhysicsAnalysis ("data_bkgStatDownVar", "Variations/BkgStatDownVariation");
  bkg_statUpVar         = new MinbiasAnalysis ("bkg_statUpVar");
  bkg_statDownVar       = new MinbiasAnalysis ("bkg_statDownVar");



  //data->Execute ();
  //mc->Execute ();
  //bkg->Execute ();
  //truth->Execute ();

  //data_dr10->Execute ();
  //mc_dr10->Execute ();
  //data_dr20->Execute ();
  //mc_dr20->Execute ();

  //data_dr02->Execute ();

  data_trackHItight->Execute ();
  //bkg_trackHItight->Execute ();
  //data_electronPtUp->Execute ();
  //data_electronPtDown->Execute ();
  //data_muonPtUp->Execute ();
  //data_muonPtDown->Execute ();
  //data_electronLHMedium->Execute ();
  //data_muonLoose->Execute ();



  data->LoadHists ();
  //data_hm->CopyAnalysis (data);
  mc->LoadHists ();
  bkg->LoadHists ();
  truth->LoadHists ();

  //data_hm->SubtractBackground ();
  //data_hm->CalculateIAA ();
  //data_hm->CalculateICP ();
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


  data->SubtractBackground (bkg);
  data->CalculateIAA ();
  data->CalculateICP ();

  data->CalculateZPtDistRatio (mc);
  data->CalculateZEtaDistRatio ();
  data->CalculateZYDistRatio ();
  data->CalculateZMassSpectraRatio (mc);

  data_dr10->LoadHists ();
  mc_dr10->LoadHists ();
  data_dr20->LoadHists ();
  mc_dr20->LoadHists ();

  data_dr02->LoadHists ();

  data_trackHItight->LoadHists ();
  bkg_trackHItight->LoadHists ();
  data_electronPtUp->LoadHists ();
  data_electronPtDown->LoadHists ();
  data_muonPtUp->LoadHists ();
  data_muonPtDown->LoadHists ();
  //data_electronLHMedium->LoadHists ();
  //data_muonLoose->LoadHists ();

  data_trackHItight->SubtractBackground (bkg_trackHItight);
  data_electronPtUp->SubtractBackground (bkg);
  data_electronPtDown->SubtractBackground (bkg);
  data_muonPtUp->SubtractBackground (bkg);
  data_muonPtDown->SubtractBackground (bkg);
  //data_electronLHMedium->SubtractBackground (bkg);
  //data_muonLoose->SubtractBackground (bkg);


  bkgSys = new Systematic (data, "bkgSys", "Background");
  //bkgSys->AddVariation (data_hm);
  bkgSys->AddVariation (data_bkgStatUpVar);
  bkgSys->AddVariation (data_bkgStatDownVar);
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

  combSys = new Systematic (data, "combSys", "Total");
  combSys->AddSystematic (bkgSys);
  combSys->AddSystematic (trkSys);
  combSys->AddSystematic (electronPtSys);
  combSys->AddSystematic (muonPtSys);
  combSys->AddSystematics ();

  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");

}



//void PlotDRVariation () {
//  TCanvas* c = new TCanvas ("c_DRVar", "", 800, 600);
//
//  TH1D* h1 = data->h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent];;
//  TH1D* h2 = data_
//  TH1D* h2 = data_
//}

#endif
