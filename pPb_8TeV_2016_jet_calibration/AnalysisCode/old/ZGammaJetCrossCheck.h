#include "../Params.h"
#include <Initialization.h>
#include "../TreeVariables.h"

#ifndef __ZGammaJetCrossCheck_h__
#define __ZGammaJetCrossCheck_h__

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
 "HLT_g35_loose"
// "HLT_g60_loose"
};
static const float photonTriggerMinPtCuts[photonTrigLength] = {15, 20, 25, 30, 35, 40};//40, 65};
static const float photonTriggerMaxPtCuts[photonTrigLength] = {20, 25, 30, 35, 40, 100000};//65, 10000};


/**
 * Calculates the original systematic error on this jet from the cross calib.
 * jpt: pt of the jet
 * jeta: eta of the jet
 */
double GetXCalibSystematicError(const double jpt, const double jeta);


/**
 * Calculates the additional systematic on jet pt given the jet eta and the
 * reference vector boson pt^ref.
 * jeta: eta of the jet
 * refpt: reference pt of the vector boson
 */
double GetNewXCalibSystematicError(TFile* file, const double jeta, const double refpt, const bool periodA);


/**
 * Primary macro.
 * dataSet: Data set identifier. This should be a run number for data or some other identifier for MC (e.g., slice number).
 * luminosity: Integrated luminosity for this run. Presumed constant over the run period. Meaningless for MC.
 * isMC: is data/MC flag.
 * isMCperiodAflag: flag that is raised for MC (meaningless if isMC is false)
 * inFileName: Input root file name where tree is stored; if == "" code will try to guess file name based on other info
 */
void ZGammaJetCrossCheck (const int dataSet,
                          const double luminosity = 0, 
                          const bool isMC = false,
                          const bool isMCperiodAflag = false, 
                          const TString inFileName = "");

#endif
