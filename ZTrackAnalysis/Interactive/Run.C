#ifndef __Run_C__
#define __Run_C__

#include "Analysis.h"
#include "ZTrackAnalysis.h"
#include "MinbiasAnalysis.h"
#include "Signal.h"

Signal* sig = nullptr;
ZTrackAnalysis* obs = nullptr;
MinbiasAnalysis* bkg = nullptr;

void Run () {

  obs = new ZTrackAnalysis ();
  obs->Execute ();
  obs->LoadHists ("ztrack");
  ResetDirectories ();

  bkg = new MinbiasAnalysis ();
  bkg->Execute ();
  bkg->LoadHists ("minbias");
  ResetDirectories ();

  sig = new Signal (obs, bkg);
  sig->GenerateHistograms ();

}

#endif
