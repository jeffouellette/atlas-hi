#ifndef __Params_h__
#define __Params_h__

#include <vector>
#include <GlobalParams.h>
#include <Utilities.h>
#include <TString.h>

using namespace std;
using namespace atlashi;

namespace JetTrackAnalysis {

const bool useGaussian = false; // whether to use Gaussian fits when creating custom TProfiles
const bool ZYAM = false;

const short rebinFactor = 20;

const double etabins[9] = {-4.5, -3.2, -2, -1, 0, 1, 2, 3.2, 4.5}; 
const short numetabins = sizeof (etabins)/sizeof (etabins[0]) - 1;

const short numphibins = 256;
const double* phibins = linspace (-TMath::Pi ()/2., 3*TMath::Pi ()/2, numphibins);

const float jet_pt_cut = 20;
const float trk_pt_cut = 0.5;

const float trijetMaxPtRatio = 0.4;

//// 0%   10%   20%   40%   60%   90%
//const float centCuts[6] = {1000, 63.6967, 47.5088, 28.6388, 15.9457, 3.04914}; // cuts in FCal
//const float centBins[6] = {0, 10, 20, 40, 60, 90};

//// 0%   5%   10%   15%   20%   40%   60%   90%
//const float centCuts[8] = {1000, 77.5874, 63.6967, 54.2686, 47.5088, 28.6388, 15.9457, 3.04914};
//const float centBins[8] = {0, 5, 10, 15, 20, 40, 60, 90};

// 0%   15%   40%   90%
const float centCuts[4] = {1000, 54.2686, 28.6388, 3.04914};
const float centBins[4] = {0, 15, 40, 90};
const int numCentBins = sizeof (centCuts) / sizeof (centCuts[0]) - 1;

// leading jet pT selections
const double ptBins[6] = {20, 30, 40, 60, 80, 100};
const short numPtBins = sizeof (ptBins) / sizeof (ptBins[0]) - 1;

// z vertex bins
const short numVertZBins = 30;
const double* vertZBins = linspace (-15, 15, numVertZBins);

const int nJetTrigIon = 1;
const TString jetTrigNamesIon[nJetTrigIon] = {
//  "HLT_j30_ion_0eta490_L1TE10",
//  "HLT_j15_ion_n320eta490_L1MBTS_1_1",
//  "HLT_j25_ion_n320eta490_L1TE5",
//  "HLT_j35_ion_n320eta490_L1TE10"
  "HLT_j40_ion_L1J5"
};
const float jetTrigMinPtCutsIon[nJetTrigIon] = {
//  30,
//  15,
//  25,
//  35,
  40
};
const float jetTrigMaxPtCutsIon[nJetTrigIon] = {
//  10000,
//  25,
//  15,
//  10000,
  10000
};
const float jetTrigMinEtaCutsIon[nJetTrigIon] = {
//  -3.2,
//  -4.9,
//  -4.9,
//  -4.9
  -3.2
};
const float jetTrigMaxEtaCutsIon[nJetTrigIon] = {
//   4.9,
//  -3.2,
//  -3.2,
//  -3.2
  3.2
};

const int nJetTrigPP = 1;
const TString jetTrigNamesPP[nJetTrigPP] = {
//  "HLT_j30_0eta490_L1TE10",
//  "HLT_j15_p320eta490_L1MBTS_1_1",
//  "HLT_j25_p320eta490_L1TE5",
//  "HLT_j35_p320eta490_L1TE10",
  "HLT_j40_L1J5"
};
const float jetTrigMinPtCutsPP[nJetTrigPP] = {
//  30,
//  15,
//  25,
//  35,
  40
};
const float jetTrigMaxPtCutsPP[nJetTrigPP] = {
//  10000,
//  25,
//  35,
//  10000,
  10000
};
const float jetTrigMinEtaCutsPP[nJetTrigPP] = {
//  -4.9,
//   3.2,
//   3.2,
//   3.2,
 -3.2
};
const float jetTrigMaxEtaCutsPP[nJetTrigPP] = {
//   3.2,
//   4.9,
//   4.9,
//   4.9,
  3.2
};

const int nJetTrig = 1;
const TString jetTrigNames [nJetTrig] = {
  "HLT_mb_mbts_L1MBTS_1"
};
const float jetTrigMinPtCuts [nJetTrig] = {
  20
};
const float jetTrigMaxPtCuts [nJetTrig] = {
  10000
};
const float jetTrigMinEtaCuts [nJetTrig] = {
  -4.9
};
const float jetTrigMaxEtaCuts [nJetTrig] = {
   4.9
};

} // end namespace

#endif
