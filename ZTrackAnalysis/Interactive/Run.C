#ifndef __Run_C__
#define __Run_C__

#include "Analysis.h"
#include "DataAnalysis.h"
#include "MCAnalysis.h"
#include "MinbiasAnalysis.h"
#include "TruthAnalysis.h"

#include "Systematics.h"
#include "WPVariation.h"

// nominal analyses
DataAnalysis* data = nullptr;
MCAnalysis* mc = nullptr;
MinbiasAnalysis* bkg = nullptr;
TruthAnalysis* truth = nullptr;

// master systematics object
Systematics* combSystematics = nullptr;

// variations for systematics
DataAnalysis* data_trackminbias = nullptr;
MinbiasAnalysis* bkg_trackminbias = nullptr;

void Run () {

  data = new DataAnalysis ();
  mc = new MCAnalysis ();
  bkg = new MinbiasAnalysis ();
  truth = new TruthAnalysis ();

  data_trackminbias = new DataAnalysis ("data_trackMinbiasVar", "Variations/TrackMinbiasWPVariation");
  bkg_trackminbias = new MinbiasAnalysis ("bkg_trackMinbiasVar", "Variations/TrackMinbiasWPVariation");



  //data->Execute ();
  //mc->Execute ();
  //bkg->Execute ();
  //truth->Execute ();

  //data_trackminbias->Execute ();
  //bkg_trackminbias->Execute ();



  data->LoadHists ();
  mc->LoadHists ();
  bkg->LoadHists ();
  truth->LoadHists ();

  data->SubtractBackground (bkg);



  data_trackminbias->LoadHists ();
  bkg_trackminbias->LoadHists ();

  data_trackminbias->SubtractBackground (bkg_trackminbias);



  combSystematics = new Systematics (data);
  combSystematics->AddVariation (data_trackminbias);
  combSystematics->CombineErrors ();

  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");

}

#endif
