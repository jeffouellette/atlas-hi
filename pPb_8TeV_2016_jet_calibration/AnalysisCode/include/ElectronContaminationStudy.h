#ifndef __ElectronContaminationStudy_h__
#define __ElectronContaminationStudy_h__

#include <TString.h>

namespace pPb8TeV2016JetCalibration {

/**
 * Primary macro.
 * isPeriodA: flag that is raised for MC (meaningless if isMC is false)
 * inFileName: Input root file name where tree is stored; if == "" code will try to guess file name based on other info
 */
void ElectronContaminationStudy (const bool isPeriodA = false, 
                                 const TString inFileName = "");

} // end namespace

#endif
