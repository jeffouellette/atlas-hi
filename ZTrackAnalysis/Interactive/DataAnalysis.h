#ifndef __DataAnalysis_h__
#define __DataAnalysis_h__

#include "Params.h"
#include "FullAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class DataAnalysis : public FullAnalysis {

  public:
  DataAnalysis (const char* _name = "data", const char* subDir = "Nominal", const bool _useHITight = false) : FullAnalysis () {
    name = _name;
    directory = Form ("DataAnalysis/%s/", subDir);
    plotFill = false;
    useHITight = _useHITight;
    LoadTrackingEfficiencies ();
    SetupDirectories (directory, "ZTrackAnalysis/");
  }

};

#endif
