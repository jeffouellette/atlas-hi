#ifndef __Params_h__
#define __Params_h__

#include <GlobalParams.h>

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {

/** User defined parameters **/

extern const bool runPeriodA; // Analyze period A data
extern const bool runPeriodB; // Analyze period B data
extern const bool runValidation; // Use validation (signal only) gamma+jet sample instead of data overlay

extern const double trk_pt_cut;
extern const double jet_pt_cut; // Cut on jet pt
extern const double photon_pt_cut; // Cut on photon pt
extern const double electron_pt_cut; // Cut on electron pt
extern const double muon_pt_cut; // Cut on muon pt
extern const double Z_pt_cut; // Cut on Z pt
extern const double Z_mass_lower_cut; // Cuts on invariant mass Z_mass +/- Z_mass_cut
extern const double Z_mass_upper_cut;
extern const double isolationEnergyIntercept; // Cut on photon isolation energy [GeV]
extern const double isolationEnergySlope; // Slope of photon isolation energy cut

extern const Color_t dataColor; // plot color for data
extern const Color_t mcOverlayColor; // plot color for MC overlay
extern const Color_t mcSignalColor; // plot color for MC signal

extern const double Z_mass_fitNsigma; // Number of sigma around the Z mass to fit invariant mass peak

extern const bool useGaussian; // 

extern const bool plot_xjref; // whether to plot xjref distributions

/** End user defined parameters **/


/** General (non-user defined) paramters **/

extern const short rebinFactor;
extern const short numxjrefbins;
extern const double* xjrefbins;

extern const short numrtrkbins;
extern const double* rtrkbins;

//extern const double pbins[31];
extern const double pbins[14];
extern const short numpbins;

extern const double pzbins[12];
extern const short numpzbins;

//extern const double etabins[7];
extern const double etabins[15];
//extern const double etabins[37];
extern const short numetabins;

extern const double xcalibEtabins[8];
extern const short numXCalibEtabins;

extern const double maxSigma;
extern const short numSigmaBins;

/** End general parameters **/

} // end namespace

#endif
