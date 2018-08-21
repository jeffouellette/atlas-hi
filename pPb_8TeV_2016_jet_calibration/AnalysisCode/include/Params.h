#ifndef __Params_h__
#define __Params_h__

#include <GlobalParams.h>

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {

/** User defined parameters **/

const bool runPeriodA = true; // Analyze period A data
const bool runPeriodB = true; // Analyze period B data
const bool runValidation = false; // Use validation (signal only) gamma+jet sample instead of data overlay

const double trk_pt_cut = 0.9;
const double jet_pt_cut = 20; // Cut on jet pt
const double photon_pt_cut = 20; // Cut on photon pt
const double electron_pt_cut = 20; // Cut on electron pt
const double muon_pt_cut = 20; // Cut on muon pt
const double Z_pt_cut = 0; // Cut on Z pt
const double Z_mass_lower_cut = 25; // Cuts on invariant mass Z_mass +/- Z_mass_cut
const double Z_mass_upper_cut = 15;
const double isolationEnergyIntercept = 4.8; // Cut on photon isolation energy [GeV]
const double isolationEnergySlope = 0.0042; // Slope of photon isolation energy cut

const Color_t data_color = kBlack; // plot color for data
const Color_t mc_color = kRed; // plot color for MC
const double Z_mass_fitNsigma = 1.5; // Number of sigma around the Z mass to fit invariant mass peak

const bool useGaussian = false; // 

const bool plot_xjref = false; // whether to plot xjref distributions

/** End user defined parameters **/


/** General (non-user defined) paramters **/

const short rebinFactor = 20;
const short numxjrefbins = 400;
const double* xjrefbins = linspace(0, 2.0, numxjrefbins);

const short numrtrkbins = 400;
const double* rtrkbins = linspace(0, 2.0, numrtrkbins);

//const double pbins[31] = {15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 90., 100., 110., 120., 130., 140., 150., 165., 180., 200., 225., 250., 275., 300., 350., 420., 500.};
const double pbins[14] = {20., 30., 40., 50., 60., 75., 90., 110., 140., 180., 220., 300., 400., 500.};
const short numpbins = sizeof(pbins)/sizeof(pbins[0]) - 1;

const double pzbins[12] = {20., 30., 40., 50., 60., 75., 90., 110., 140., 180., 220., 300.};
const short numpzbins = sizeof(pzbins)/sizeof(pzbins[0]) - 1;

//const double etabins[7] = {-2.5, -1.3, -0.5, 0, 0.5, 1.3, 2.5};
//const short numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;
const double etabins[15] = {-4.4, -3.6, -2.8, -2.1, -1.2, -0.8, -0.3, 0, 0.3, 0.8, 1.2, 2.1, 2.8, 3.6, 4.4};
const short numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;

const double xcalibEtabins[8] = {0, 0.3, 0.8, 1.2, 2.1, 2.8, 3.6, 4.4};
const short numXCalibEtabins = sizeof(xcalibEtabins)/sizeof(xcalibEtabins[0]) - 1;

const double maxSigma = 0.40;
const short numSigmaBins = 80;

/** End general parameters **/

} // end namespace

#endif
