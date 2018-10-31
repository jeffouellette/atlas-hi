#ifndef __Params_h__
#define __Params_h__

#include <GlobalParams.h>
#include <Utils.h>
#include <TFile.h>

using namespace atlashi;

namespace JetCalibration {

/** General parameters **/

//TFile* xCalibSystematicsFile = NULL;
//TFile* dataOverMCFile = NULL;

const bool runValidation = true; // Use validation (signal only) gamma+jet sample instead of data overlay

// pT and analysis cuts
const double trk_pt_cut = 0.5; // Cut on track pt in GeV
const double jet_pt_cut = 20; // Cut on jet pt
const double photon_pt_cut = 20; // Cut on photon pt
const double electron_pt_cut = 20; // Cut on electron pt
const double muon_pt_cut = 20; // Cut on muon pt
const double Z_pt_cut = 0; // Cut on Z pt

const double Z_mass_lower_cut = 25; // Cuts on invariant mass Z_mass +/- Z_mass_cut
const double Z_mass_upper_cut = 15;

const double isolationEnergyIntercept = 4.8; // Cut on photon isolation energy [GeV]
const double isolationEnergySlope = 0.0042; // Slope of photon isolation energy cut

// plotting variables
const Color_t dataColor = kBlack; // plot color for data
const Color_t mcOverlayColor = kRed; // plot color for MC overlay
const Color_t mcSignalColor = kBlue; // plot color for MC signal

const double Z_mass_fitNsigma = 1.5; // Number of sigma around the Z mass to fit invariant mass peak

const bool skipOldInsitu = true; // whether to skip plotting results with old insitu factors (2015)
const bool skipSignalMC = true; // whether to skip plotting results with signal only photon+jet MC

const bool useGaussian = false; // whether to use Gaussian fits when creating custom TProfiles

const bool plot_xjref = false; // whether to plot xjref distributions

const bool calcPtClosure = true; // whether the energy scale calculations will use pT (true) or E (false)

const bool subtractBackground = true; // whether to subtract the calculated background

const bool exclusive = false; // whether to calculate xjref outside of eta bounds (vs. inside)


/** binning parameters **/

const short rebinFactor = 20;
const short numxjrefbins = 400;
const double* xjrefbins = linspace (0, 2.0, numxjrefbins);

const short numrtrkbins = 400;
const double* rtrkbins = linspace (0, 2.0, numrtrkbins);

const double dpbins[7] = {17, 35, 50, 70, 140, 280, 500};
//const double dpbins[2] = {60, 140};
const short numdpbins = sizeof (dpbins)/sizeof (dpbins[0]) - 1;

//const double pbins[31] = {15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 90., 100., 110., 120., 130., 140., 150., 165., 180., 200., 225., 250., 275., 300., 350., 420., 500.};
const double pbins[16] = {20., 25., 30., 35., 40., 50., 60., 70., 90., 110., 140., 180., 220., 280., 350., 500.};
const short numpbins = sizeof (pbins)/sizeof (pbins[0]) - 1;

//const double pzbins[12] = {20., 30., 40., 50., 60., 75., 90., 110., 140., 180., 220., 300.};
//const short numpzbins = sizeof (pzbins)/sizeof (pzbins[0]) - 1;

//const double etabins[7] = {-2.5, -1.3, -0.5, 0, 0.5, 1.3, 2.5};
const double etabins[15] = {-4.4, -3.6, -2.8, -2.1, -1.2, -0.8, -0.3, 0, 0.3, 0.8, 1.2, 2.1, 2.8, 3.6, 4.4};
//const double etabins[37] = {-4.4, -3.6, -3.2, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.6, 4.4};
const short numetabins = sizeof (etabins)/sizeof (etabins[0]) - 1;

const short numphibins = 48;
const double* phibins = linspace (-pi, pi, numphibins);

//const short numetabins = 98;
//const double* etabins = linspace (-4.9, 4.9, numetabins);

const double zetabins[7] = {-2.4, -1.52, -1.37, 0, 1.37, 1.52, 2.4};
const short numzetabins = sizeof (zetabins) / sizeof (zetabins[0]) - 1;

const double xcalibEtabins[8] = {0, 0.3, 0.8, 1.2, 2.1, 2.8, 3.6, 4.4};
const short numXCalibEtabins = sizeof (xcalibEtabins)/sizeof (xcalibEtabins[0]) - 1;

const short numclosurebins = 100;
const double* closurebins = linspace (0, 2.0, numclosurebins);

const double maxSigma = 0.40;
const short numSigmaBins = 80;

// bin min and max for combined results
const int eta_lo_comb = 5;
const int eta_hi_comb = 10;
const int p_lo_comb = 8;
const int p_hi_comb = numpbins;


/** Trigger declarations **/

const short electronTrigLength = 1;
const char* electronTriggerNames[electronTrigLength] = {
 "HLT_e15_lhloose_nod0"
 //"HLT_e15_loose_ion_L1EM12" // old - not used in HION5
};
const float electronTriggerMinPtCuts[electronTrigLength] = {15};
const float electronTriggerMaxPtCuts[electronTrigLength] = {100000};

const short muonTrigLength = 1;
const char* muonTriggerNames[muonTrigLength] = {
 "HLT_mu15"
 //"HLT_mu8" // old - not used in HION5
};
const float muonTriggerMinPtCuts[muonTrigLength] = {15};
const float muonTriggerMaxPtCuts[muonTrigLength] = {100000};

const short photonTrigLength = 6;
const char* photonTriggerNames[photonTrigLength] = {
 "HLT_g10_loose",
 "HLT_g15_loose",
 "HLT_g20_loose",
 "HLT_g25_loose",
 "HLT_g30_loose",
 "HLT_g35_loose"
};
const float photonTriggerMinPtCuts[photonTrigLength] = {15, 20, 25, 30, 35, 40};
const float photonTriggerMaxPtCuts[photonTrigLength] = {20, 25, 30, 35, 40, 100000};

const short jetTrigLength = 1;
const char* jetTriggerNames[jetTrigLength] = {
// "HLT_j30_ion_0eta490_L1TE10",
// "HLT_j35_ion_n320eta490_L1TE10",
// "HLT_j40_ion_L1J5",
 "HLT_j50_ion_L1J10",
// "HLT_j60_ion_L1J20",
// "HLT_j90_ion_L1J20",
// "HLT_j100_ion_L1J20"
};

const float jetTriggerMinPtCuts[jetTrigLength] = {60};//40, 45, 55, 60, 70, 100, 110};
const float jetTriggerMaxPtCuts[jetTrigLength] = {100000};//45, 50, 60, 70, 100, 110, 100000};

} // end namespace

#endif
