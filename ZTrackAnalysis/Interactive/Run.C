#ifndef __Run_C__
#define __Run_C__

#include "Analysis.h"
#include "DataAnalysis.h"
#include "MCAnalysis.h"
#include "MinbiasAnalysis.h"
#include "TruthAnalysis.h"

#include "Systematics.h"
//#include "WPVariation.h"

// nominal analyses
DataAnalysis* data = nullptr, *data_hm = nullptr;
MCAnalysis* mc = nullptr;
MinbiasAnalysis* bkg = nullptr;
TruthAnalysis* truth = nullptr;

// master systematics object
Systematics* combSys = nullptr;
Systematics* bkgSys = nullptr;
Systematics* trkSys = nullptr;

// variations for systematics
DataAnalysis* data_trackHIloose = nullptr;
MinbiasAnalysis* bkg_trackHIloose = nullptr;

void Run () {

  data = new DataAnalysis ();
  data_hm = new DataAnalysis ("data_hm");
  mc = new MCAnalysis ();
  bkg = new MinbiasAnalysis ();
  truth = new TruthAnalysis ();

  data_trackHIloose = new DataAnalysis ("data_trackHILooseVar", "Variations/TrackHILooseWPVariation");
  bkg_trackHIloose = new MinbiasAnalysis ("bkg_trackHILooseVar", "Variations/TrackHILooseWPVariation");



  //data->Execute ();
  //mc->Execute ();
  //bkg->Execute ();
  //truth->Execute ();

  //data_trackHIloose->Execute ();
  //bkg_trackHIloose->Execute ();



  data->LoadHists ();
  data_hm->CopyAnalysis (data);
  mc->LoadHists ();
  bkg->LoadHists ();
  truth->LoadHists ();

  data->SubtractBackground (bkg);
  data->CalculateZMassSpectraRatio (mc);
  data_hm->SubtractBackground (data_hm);



  data_trackHIloose->LoadHists ();
  bkg_trackHIloose->LoadHists ();

  data_trackHIloose->SubtractBackground (bkg_trackHIloose);


  bkgSys = new Systematics (data);
  bkgSys->AddVariation (data_hm);
  bkgSys->CombineErrors ();

  trkSys = new Systematics (data);
  trkSys->AddVariation (data_trackHIloose);
  trkSys->CombineErrors ();

  combSys = new Systematics (data);
  combSys->AddSystematic (bkgSys);
  combSys->AddSystematic (trkSys);

  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");

}

#endif
