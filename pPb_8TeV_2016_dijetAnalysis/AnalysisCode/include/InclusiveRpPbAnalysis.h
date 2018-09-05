#ifndef __InclusiveRpPbAnalysis_h__
#define __InclusiveRpPbAnalysis_h__

#include "../Util.C"

namespace pPb8TeV2016DijetAnalysis {

extern const double etaCoM = -0.465; // boost into CoM frame in period B kinematics

extern const int numppPtbins = 36;
extern const int numppYbins = 12;
extern double ppPtbins[(numppPtbins+1)*numppYbins];// = {}; //= new double[(numppPtbins+1)*numppYbins];
extern const double ppYbins[] = {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3};

//extern const int numpPbYbins = 18;
extern const double pPbYbins[] = {-3.465, -3.2, -2.965, -2.465, -2, -1.965, -1.465, -1, -0.965, -0.465, 0, 0.035, 0.535, 1, 1.035, 1.535, 2, 2.035, 2.535};
extern const int numpPbYbins = sizeof(pPbYbins)/sizeof(pPbYbins[0]) - 1;

int getppPbin (double pt);


int getppYbin (double eta);


int getpPbYbin (double eta);

/**
 * Loads binning information from the 8TeV pp reference data.
 */
TH1D** setupPPConfiguration();


void InclusiveRpPbAnalysis(const int thisRunNumber, // Run number identifier.
                           double luminosity); // Integrated luminosity for this run. Presumed constant over the run period.

} // end namespace

#endif
