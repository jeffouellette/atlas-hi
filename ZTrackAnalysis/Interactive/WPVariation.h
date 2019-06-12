#ifndef __WPVariation_h__
#define __WPVariation_h__

#include "Params.h"
#include "Analysis.h"

#include <GlobalParams.h>
#include <Utilities.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class WPVariation : public Analysis {

  private:
  string analysisType; // Data, Minbias, etc.
  string wpName; // TrackMinbias, LHTightElectron, etc.

  public:
  WPVariation (string _wpName, string _analysisType) : Analysis () {
    wpName = _wpName;
    analysisType = _analysisType;
    name = Form ("%s_%sWPVariation", analysisType.c_str (), wpName.c_str ());
    directory = Form ("%s/Variations/%sWPVariation/", analysisType.c_str (), wpName.c_str ());
    plotFill = false;
    LoadTrackingEfficiencies ();
    SetupDirectories (directory, "ZTrackAnalysis/");
  }
};


#endif
