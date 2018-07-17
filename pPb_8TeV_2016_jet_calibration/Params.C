#include "../GlobalParams.C"

/** User defined parameters **/

static const bool runPeriodA = true; // Analyze period A data
static const bool runPeriodB = true; // Analyze period B data
static const bool runValidation = false; // Use validation (signal only) gamma+jet sample instead of data overlay

static const double electron_pt_cut = 20; // Cut on electron pt
static const double muon_pt_cut = 20; // Cut on muon pt
static const double Z_pt_cut = 20; // Cut on Z pt
static const double Z_mass_lower_cut = 25; // Cuts on invariant mass Z_mass +/- Z_mass_cut
static const double Z_mass_upper_cut = 15;
static const double isolationEnergyIntercept = 4.8; // Cut on photon isolation energy [GeV]
static const double isolationEnergySlope = 0.0042; // Slope of photon isolation energy cut

static const Color_t data_color = kBlack; // plot color for data
static const Color_t mc_color = kRed; // plot color for MC
static const double Z_mass_fitNsigma = 1.5; // Number of sigma around the Z mass to fit invariant mass peak

static const bool useGaussian = false; // 

/** End user defined parameters **/


/** General (non-user defined) paramters **/

static const short rebinFactor = 20;
static const short numxjrefbins = 400;
static const double* xjrefbins = linspace(0, 2.0, numxjrefbins);
static const double pgammabins[29] = {15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 90., 100., 110., 120., 130., 140., 150., 165., 180., 200., 225., 250., 275., 300., 350.};
//static const double pgammabins[12] = {20., 30., 40., 50., 60., 75., 90., 110., 140., 180., 220., 300.};
static const short numpgammabins = sizeof(pgammabins)/sizeof(pgammabins[0]) - 1;

static const double pzbins[12] = {20., 30., 40., 50., 60., 75., 90., 110., 140., 180., 220., 300.};
static const short numpzbins = sizeof(pzbins)/sizeof(pzbins[0]) - 1;
//static const double etabins[7] = {-1.3, -0.8, -0.3, 0, 0.3, 0.8, 1.3};
static const double etabins[7] = {-2.5, -1.3, -0.5, 0, 0.5, 1.3, 2.5};
static const short numetabins = sizeof(etabins)/sizeof(etabins[0]) - 1;
static const double xcalibEtabins[8] = {0, 0.3, 0.8, 1.2, 2.1, 2.8, 3.6, 4.4};
static const short numXCalibEtabins = sizeof(xcalibEtabins)/sizeof(xcalibEtabins[0]) - 1;
static const double maxSigma = 0.40;
static const short numSigmaBins = 80;

/** End general parameters **/
