#ifndef __Params_h__
#define __Params_h__

#include <string>

using namespace std;

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

const double pbins[22] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120, 150, 180, 220, 240, 300, 350};
const int numpbins = sizeof (pbins) / sizeof (pbins[0]) - 1;

const double etabins[4] = {0, 1.37, 1.52, 2.37};
const int numetabins = sizeof (etabins) / sizeof (etabins[0]) - 1;

#endif
