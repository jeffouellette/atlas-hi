#ifndef __FCalDistribution_h__
#define __FCalDistribution_h__

#include <TString.h>

namespace pPb8TeV2016JetCalibration {

static const short electronTrigLength = 1;
static const char* electronTriggerNames[electronTrigLength] = {
 "HLT_e15_loose_ion_L1EM12"
};
static const float electronTriggerMinPtCuts[electronTrigLength] = {15};
static const float electronTriggerMaxPtCuts[electronTrigLength] = {100000};

static const short muonTrigLength = 1;
static const char* muonTriggerNames[muonTrigLength] = {
 "HLT_mu8",
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

static const short photonTrigLength = 6;//7;
static const char* photonTriggerNames[photonTrigLength] = {
 "HLT_g10_loose",
 "HLT_g15_loose",
 "HLT_g20_loose",
 "HLT_g25_loose",
 "HLT_g30_loose",
 "HLT_g35_loose",
 //"HLT_g60_loose"
};
static const float photonTriggerMinPtCuts[photonTrigLength] = {15, 20, 25, 30, 35, 40};//40, 65};
static const float photonTriggerMaxPtCuts[photonTrigLength] = {20, 25, 30, 35, 40, 100000};//65, 100000};

TString GetIdentifier (const int dataSet, const bool isMC, const bool isValidationSample, const bool periodA);

void FCalDistribution (const int dataSet,
                       const double luminosity = 0, 
                       const bool isMC = false,
                       const bool isMCperiodAflag = false, 
                       const TString inFileName = "");

}

#endif
