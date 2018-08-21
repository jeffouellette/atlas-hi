#ifndef __ElectronContaminationStudyHist_h__
#define __ElectronContaminationStudyHist_h__

#include "Params.h"
#include <Initialization.h>

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {

static const double pebins[18] = {20, 25, 35, 45, 55, 65, 75, 85, 105, 125, 150, 175, 200, 250, 300, 350, 400, 550};
static const short numpebins = sizeof(pebins)/sizeof(pebins[0]) - 1;

static const double eetabins[6] = { -2.37, -1.56, -1.37, 1.37, 1.56, 2.37};
static const short numeetabins = sizeof(eetabins)/sizeof(eetabins[0]) - 1;

void ElectronContaminationStudyHist ();

} // end namespace

#endif
