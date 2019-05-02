#ifndef __Run_C__
#define __Run_C__

#include "Analysis.h"
#include "DataAnalysis.h"
#include "MCAnalysis.h"
#include "MinbiasAnalysis.h"
#include "TruthAnalysis.h"
#include "Signal.h"

Signal* sig = nullptr;
DataAnalysis* data = nullptr;
MCAnalysis* mc = nullptr;
MinbiasAnalysis* bkg = nullptr;
TruthAnalysis* truth = nullptr;

void Run () {

  data = new DataAnalysis ();
  mc = new MCAnalysis ();
  bkg = new MinbiasAnalysis ();
  truth = new TruthAnalysis ();

  data->Execute ();
  mc->Execute ();
  bkg->Execute ();
  truth->Execute ();

  data->LoadHists ();
  mc->LoadHists ();
  bkg->LoadHists ();
  truth->LoadHists ();

  sig = new Signal (data, bkg);

  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");

}

#endif
