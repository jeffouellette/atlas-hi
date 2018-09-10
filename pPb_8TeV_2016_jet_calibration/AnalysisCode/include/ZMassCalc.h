#ifndef __ZMassCalc_h__
#define __ZMassCalc_h__

#include "Params.h"
#include <Initialization.h>
#include "TreeVariables.h"

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {

static const short electronTrigLength = 1;
static const char* electronTriggerNames[electronTrigLength] = {
 "HLT_e15_lhloose_nod0"
};
static const float electronTriggerMinPtCuts[electronTrigLength] = {15};
static const float electronTriggerMaxPtCuts[electronTrigLength] = {100000};

static const short muonTrigLength = 1;
static const char* muonTriggerNames[muonTrigLength] = {
 "HLT_mu15",
};
static const float muonTriggerMinPtCuts[muonTrigLength] = {8};
static const float muonTriggerMaxPtCuts[muonTrigLength] = {100000};

//static const short electronTrigLength = 3;
//static const char* electronTriggerNames[electronTrigLength] = {
// "HLT_e20_lhloose",
// "HLT_e22_lhloose",
// "HLT_e24_lhloose"
//};
//static const float electronTriggerMinPtCuts[electronTrigLength] = {20, 22, 24};
//static const float electronTriggerMaxPtCuts[electronTrigLength] = {100000, 100000, 100000};
//
//static const short muonTrigLength = 3;
//static const char* muonTriggerNames[muonTrigLength] = {
// "HLT_mu15",
// "HLT_mu18",
// //"HLT_mu20",
// "HLT_mu20_L1MU15"
//};
//static const float muonTriggerMinPtCuts[muonTrigLength] = {15, 18, 20};
//static const float muonTriggerMaxPtCuts[muonTrigLength] = {100000, 100000, 100000};


/**
 * Primary macro.
 * dataSet: Data set identifier. This should be a run number for data or some other identifier for MC (e.g., slice number).
 * luminosity: Integrated luminosity for this run. Presumed constant over the run period. Meaningless for MC.
 * isMC: is data/MC flag.
 * isPeriodA: flag that is raised for MC (meaningless if isMC is false)
 * inFileName: Input root file name where tree is stored; if == "" code will try to guess file name based on other info
 */
void ZMassCalc (const int dataSet,
                const double luminosity = 0, 
                const bool isMC = false,
                const bool isPeriodA = false, 
                const TString inFileName = "");

} // end namespace

#endif
