#ifndef __Params_h__
#define __Params_h__

#include <string>

using namespace std;

namespace offlineAnalyses {

const int electronTrigN = 22;
const string electronTrigNames[electronTrigN] = {
  "HLT_e13_etcut_ion",
  "HLT_e15_etcut_ion",
  "HLT_e15_lhloose_ion_L1EM12",
  "HLT_e15_lhmedium_ion_L1EM12",
  "HLT_e15_loose_ion",
  "HLT_e15_loose_ion_idperf",
  "HLT_e15_medium_ion",
  "HLT_e18_etcut_ion",
  "HLT_e18_lhloose_ion",
  "HLT_e18_lhmedium_ion",
  "HLT_e18_loose_ion",
  "HLT_e18_loose_ion_idperf",
  "HLT_e18_medium_ion",
  "HLT_e20_etcut_ion",
  "HLT_e20_lhloose_ion",
  "HLT_e20_lhmedium_ion",
  "HLT_e20_loose_ion",
  "HLT_e20_loose_ion_idperf",
  "HLT_e20_medium_ion",
  "HLT_e30_etcut_ion",
  "HLT_e50_etcut_ion",
  "HLT_2e20_loose_ion"
};

const int photonTrigN = 10;
const string photonTrigNames[photonTrigN] = {
  "HLT_g13_etcut_ion",
  "HLT_g15_loose_ion",
  "HLT_g18_etcut",
  "HLT_g18_etcut_ion",
  "HLT_g20_loose",
  "HLT_g20_loose_ion",
  "HLT_g28_etcut_ion",
  "HLT_g30_loose_ion",
  "HLT_g50_loose",
  "HLT_g50_loose_ion"
};
const float photonTrigMinPtCuts[photonTrigN] = {
  13,
  15,
  18,
  18,
  20,
  20,
  28,
  30,
  50,
  50
};

const int jetTrigN = 15;
const string jetTrigNames[jetTrigN] = {
  "HLT_j50_ion_L1J12",
  "HLT_j60_ion_L1J20",
  "HLT_j75_ion_L1J30",
  "HLT_j85_ion_L1J30",
  "HLT_j100_ion_L1J30",
  "HLT_j110_ion_L1J50",
  "HLT_j120_ion_L1J50",
  "HLT_j150_ion_L1J50",
  "HLT_j180_ion_L1J50",
  "HLT_j200_ion_L1J50",
  "HLT_j50_320eta490_ion",
  "HLT_j60_320eta490_ion",
  "HLT_j70_320eta490_ion",
  "HLT_j100_320eta490_ion",
  "HLT_j170_320eta490_ion"
};

const double pbins[11] = {0, 20, 30, 40, 50, 60, 80, 100, 125, 150, 200};
const int numpbins = sizeof (pbins) / sizeof (pbins[0]) - 1;

const double etabins[4] = {0, 1.37, 1.52, 2.37};
const int numetabins = sizeof (etabins) / sizeof (etabins[0]) - 1;

const float showerShapeBinsLow[12] = {-0.02, -0.04, 0, 0.7, 0.5, 0.2, 0.005, 0.5, -0.2, 0, -10, 0.6};
const float showerShapeBinsHigh[12] = {0.02, 0.04, 100000, 1, 1.1, 1, 0.02, 4, 0.6, 0.9, 350, 1};

const int showerShapeMCRebinFactors[12] = {4, 4, 8, 8, 4, 8, 8, 4, 16, 8, 4, 4};

const int runNumbers[15] = {365498, 365502, 365512, 365573, 365602, 365627, 365678, 365681, 365709, 365752, 365763, 365768, 365834, 365932, 366011};
const int numRunNumbers = sizeof (runNumbers) / sizeof (runNumbers[0]);

} // end namespace

#endif
