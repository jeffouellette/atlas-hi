#ifndef __EMTopoComparison_h__
#define __EMTopoComparison_h__

#include "Params.h"
#include <Initialization.h>
#include "TreeVariables.h"

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {

//static const short electronTrigLength = 1;
//static const char* electronTriggerNames[electronTrigLength] = {
// "HLT_e15_lhloose_nod0"
//};
//static const float electronTriggerMinPtCuts[electronTrigLength] = {15};
//static const float electronTriggerMaxPtCuts[electronTrigLength] = {100000};
//
//static const short muonTrigLength = 1;
//static const char* muonTriggerNames[muonTrigLength] = {
// "HLT_mu15",
//};
//static const float muonTriggerMinPtCuts[muonTrigLength] = {15};
//static const float muonTriggerMaxPtCuts[muonTrigLength] = {100000};
//
//static const short photonTrigLength = 6;
//static const char* photonTriggerNames[photonTrigLength] = {
// "HLT_g10_loose",
// "HLT_g15_loose",
// "HLT_g20_loose",
// "HLT_g25_loose",
// "HLT_g30_loose",
// "HLT_g35_loose"
//};
//static const float photonTriggerMinPtCuts[photonTrigLength] = {15, 20, 25, 30, 35, 40};
//static const float photonTriggerMaxPtCuts[photonTrigLength] = {20, 25, 30, 35, 40, 100000};


/**
 * Calculates the original systematic error on this jet from the cross calib.
 * jpt: pt of the jet
 * jeta: eta of the jet
 */
double GetXCalibSystematicError (const double jpt, const double jeta);


/**
 * Primary macro.
 * dataSet: Data set identifier. This should be a run number for data or some other identifier for MC (e.g., slice number).
 * luminosity: Integrated luminosity for this run. Presumed constant over the run period. Meaningless for MC.
 * isMC: is data/MC flag.
 * isPeriodA: flag that is raised for MC (meaningless if isMC is false)
 * inFileName: Input root file name where tree is stored; if == "" code will try to guess file name based on other info
 */
void EMTopoComparison (const int dataSet,
                       const double luminosity = 0, 
                       const bool isMC = false,
                       const bool isPeriodA = false, 
                       const TString inFileName = "",
                       const double crossSection_microbarns = 0,
                       const double filterEfficiency = 0,
                       const int numberEvents = 0);

} // end namespace

#endif
