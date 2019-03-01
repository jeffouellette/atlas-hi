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
  //data->Execute ();
  data->LoadHists ();

  mc = new MCAnalysis ();
  //mc->Execute ();
  mc->LoadHists ();

  //bkg = new MinbiasAnalysis ();
  //bkg->Execute ();
  //bkg->LoadHists ();

  //truth = new TruthAnalysis ();
  //truth->Execute ();
  //truth->LoadHists ();

  //sig = new Signal (data, bkg);
  //sig->GenerateHistograms ();

}

#endif
