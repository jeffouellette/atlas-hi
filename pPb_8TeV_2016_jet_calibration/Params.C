#include "../GlobalParams.C"

/** User defined parameters **/

static const bool runPeriodA = true; // Analyze period A data
static const bool runPeriodB = true; // Analyze period B data

static const double electron_pt_cut = 20; // Cut on electron pt
static const double muon_pt_cut = 20; // Cut on muon pt
static const double Z_mass_lower_cut = 25; // Cuts on invariant mass Z_mass +/- Z_mass_cut
static const double Z_mass_upper_cut = 15;
static const double isolationEnergyCut = 5; // Cut on photon isolation energy

static const Color_t data_color = kBlack;
static const Color_t mc_color = kRed;
static const double Z_mass_fitNsigma = 1.3;

/** End user defined parameters **/


/** General (non-user defined) paramters **/

//static const double xjrefbins[181] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.30, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57
static const short numxjrefbins = 180;
static const double* xjrefbins = linspace(0, 1.8, numxjrefbins);
static const double pbins[15] = {20., 30., 40., 50., 60., 75., 90., 105., 120., 140., 160., 180., 200., 250., 300.};
static const short numpbins = sizeof(pbins)/sizeof(pbins[0]) - 1;
//static const double etabins[7] = {-1.3, -0.8, -0.3, 0, 0.3, 0.8, 1.3};
static const double etabins[7] = {-2.5, -1.3, -0.5, 0, 0.5, 1.3, 2.5};
static const short numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;
static const double xcalibEtabins[8] = {0, 0.3, 0.8, 1.2, 2.1, 2.8, 3.6, 4.4};
static const short numXCalibEtabins = sizeof(xcalibEtabins)/sizeof(xcalibEtabins[0]) - 1;

/** End general parameters **/
