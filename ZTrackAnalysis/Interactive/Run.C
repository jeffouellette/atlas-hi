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
FullAnalysis* data = nullptr, *data_hm = nullptr;
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


// variations for systematics
PhysicsAnalysis* data_trackHItight = nullptr;
MinbiasAnalysis* bkg_trackHItight = nullptr;
PhysicsAnalysis* data_electronPtUp = nullptr;
PhysicsAnalysis* data_electronPtDown = nullptr;
PhysicsAnalysis* data_muonPtUp = nullptr;
PhysicsAnalysis* data_muonPtDown = nullptr;
PhysicsAnalysis* data_electronLHMedium = nullptr;
PhysicsAnalysis* data_muonLoose = nullptr;

void Run () {

  //data = new FullAnalysis ("data_preleptfix", "PreLeptonFix");
  data    = new FullAnalysis ("data", "DataAnalysis/Nominal");
  data_hm = new DataAnalysis ("data_hm");
  mc      = new MCAnalysis ();
  bkg     = new MinbiasAnalysis ();
  truth   = new TruthAnalysis ();

  data_trackHItight     = new PhysicsAnalysis ("data_trackHITightVar", "Variations/TrackHITightWPVariation");
  bkg_trackHItight      = new MinbiasAnalysis ("bkg_trackHITightVar", "Variations/TrackHITightWPVariation");
  data_electronPtUp     = new PhysicsAnalysis ("data_electronPtUpVar", "Variations/ElectronPtUpVariation");
  data_electronPtDown   = new PhysicsAnalysis ("data_electronPtDownVar", "Variations/ElectronPtDownVariation");
  data_muonPtUp         = new PhysicsAnalysis ("data_muonPtUpVar", "Variations/MuonPtUpVariation");
  data_muonPtDown       = new PhysicsAnalysis ("data_muonPtDownVar", "Variations/MuonPtDownVariation");
  data_electronLHMedium = new PhysicsAnalysis ("data_electronLHMediumVar", "Variations/ElectronLHMediumVariation");
  data_muonLoose        = new PhysicsAnalysis ("data_muonLooseVar", "Variations/MuonLooseVariation");



  //data->Execute ();
  //mc->Execute ();
  //bkg->Execute ();
  //truth->Execute ();

  //data_trackHItight->Execute ();
  //bkg_trackHItight->Execute ();
  //data_electronPtUp->Execute ();
  //data_electronPtDown->Execute ();
  //data_muonPtUp->Execute ();
  //data_muonPtDown->Execute ();
  //data_electronLHMedium->Execute ();
  //data_muonLoose->Execute ();



  data->LoadHists ();
  data_hm->CopyAnalysis (data);
  mc->LoadHists ();
  bkg->LoadHists ();
  //truth->LoadHists ();

  data->SubtractBackground (bkg);

  data->CalculateZPtDistRatio (mc);
  data->CalculateZEtaDistRatio ();
  data->CalculateZYDistRatio ();
  data->CalculateZMassSpectraRatio (mc);
  data_hm->SubtractBackground (data_hm);



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


  bkgSys = new Systematic (data, "bkgSys");
  bkgSys->AddVariation (data_hm);
  bkgSys->AddVariations ();

  trkSys = new Systematic (data, "trkSys");
  trkSys->AddVariation (data_trackHItight);
  trkSys->AddVariations ();

  electronPtSys = new Systematic (data, "electronPtSys");
  electronPtSys->AddVariation (data_electronPtUp);
  electronPtSys->AddVariation (data_electronPtDown);
  electronPtSys->AddVariations ();

  muonPtSys = new Systematic (data, "muonPtSys");
  muonPtSys->AddVariation (data_muonPtUp);
  muonPtSys->AddVariation (data_muonPtDown);
  muonPtSys->AddVariations ();

  combSys = new Systematic (data, "combSys");
  combSys->AddSystematic (bkgSys);
  combSys->AddSystematic (trkSys);
  combSys->AddSystematic (electronPtSys);
  combSys->AddSystematic (muonPtSys);
  combSys->AddSystematics ();

  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");

}

#endif
