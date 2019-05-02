#ifndef __Systematics_h__
#define __Systematics_h__

#include "Params.h"
#include "Analysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

virtual class Systematics {

  private:
  Analysis* upVar = nullptr;
  Analysis* downVar = nullptr;

  protected:
  const char* name;
  const char* directory;

  public:
  Systematics () {
    name = "data";
    directory = "Systematics/";
    SetupDirectories (directory, "ZTrackAnalysis/");
  }

  void Execute ();
};


////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics::Execute () {

  upVar->Execute ();
  downVar->Execute ();

}

#endif
