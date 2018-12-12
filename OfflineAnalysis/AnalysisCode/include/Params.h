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

const float showerShapeBinsLow[13] = {-0.02, -0.04, 0, 0.7, 0.5, 0.2, 0.005, 0.5, -0.2, -0.02, 0, -10, 0.6};
const float showerShapeBinsHigh[13] = {0.02, 0.04, 100000, 1, 1.1, 1, 0.02, 4, 0.6, 0.08, 0.9, 350, 1};

const int showerShapeMCRebinFactors[13] = {4, 4, 8, 8, 4, 8, 8, 4, 16, 16, 8, 4, 4};

const int runNumbers[48] = {365498, 365502, 365512, 365573, 365602, 365627, 365678, 365681, 365709, 365752, 365763, 365768, 365834, 365914, 365932, 366011, 366029, 366092, 366142, 366268, 366337, 366383, 366413, 366476, 366526, 366528, 366627, 366691, 366754, 366805, 366860, 366878, 366919, 366931, 366994, 367023, 367099, 367134, 367165, 367170, 367233, 367273, 367318, 367321, 367363, 367364, 367365, 367384};
const int numRunNumbers = sizeof (runNumbers) / sizeof (runNumbers[0]);

const float lumis[48] = {0.2674, 1.018, 0.9223, 7.852, 9.897, 25.151, 16.219, 24.847, 35.472, 36.751, 5.484, 7.659, 5.651, 6.595, 41.822, 43.806, 42.881, 45.336, 41.946, 54.532, 44.538, 0.1853, 40.729, 42.383, 1.034, 33.016, 36.795, 55.241, 61.82, 72.872, 65.283, 64.03, 73.819, 57.678, 9.506, 64.851, 65.985, 72.498, 8.648, 47.998, 55.045, 53.962, 64.46, 65.577, 13.058, 66.662, 51.885, 7.957};

} // end namespace

#endif
