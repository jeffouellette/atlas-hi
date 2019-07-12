#ifndef __Utilities_h__
#define __Utilities_h__

#include <TFile.h>
#include <TString.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraphAsymmErrors.h>

namespace atlashi {

/**
 * Safely deletes an object so that its pointer is also deleted.
 */
template <typename T> inline void SaferDelete (T* &t) {
  if (t) { delete t; t = nullptr; }
}

/**
 * Returns a TString summarizing a measurement.
 * By default, there is only 1 significant digit in the error.
 * E.g., FormatMeasurement (40.58, 1.29) returns "40#pm1".
 * Or, FormatMeasurement (40.58, 1.29, 2) returns "40.6#pm1.3".
 */
const char* FormatMeasurement (double val, double err, const int n=1);


/**
 * Modifies the directory strings to point to the correct locations.
 */
void SetupDirectories (const TString dataSubDir, const TString thisWorkPath);


/**
 * Clears sub-directory information from the directory strings
 */
void ResetDirectories ();


/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
double* linspace (double lo, double hi, int num);


/**
 * Returns a logarithmically spaced array, where the 0th element is lo and the num-th element is hi.
 */
double* logspace (double lo, double hi, int num);


/**
 * Returns the equivalent angle in the range 0 to 2pi.
 */
double InTwoPi (double phi);


/**
 * Returns the difference between two angles in 0 to pi.
 * If sign is true, will return a signed version such that phi2 = phi1 + dphi
 */
double DeltaPhi (double phi1, double phi2, const bool sign = false);


/**
 * Returns dR between two eta, phi coordinates.
 */
double DeltaR (const double eta1, const double eta2, const double phi1, const double phi2 );


/**
 * Returns true iff this eta, phi coordinate lies in the disabled HEC region.
 */
bool InDisabledHEC (const double eta, double phi, const double dr = 0.4);


/**
 * Returns true iff this eta lies within the EMCal.
 */
bool InEMCal (const float eta);


/**
 * Returns true iff this object is within a given radius in the HCal.
 */
bool InHadCal (const float eta, const float R = 0.4);


/**
 * Sets all the errors in this histogram to 0.
 */
void ResetHistErrors (TH1D* h);


/**
 * Adds nSigma statistical error variations to this histogram
 */
void AddStatVar (TH1D* h, const bool upvar = true, const float nSigma = 1);


/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 */
void AddErrorsInQuadrature (TH1D* master, TH1D* sys);


/**
 * Calculates simple systematics as maximum variations on the nominal.
 * Intended for combining up/down variations in an expandable way.
 */
void CalcSystematics (TH1D* sys, TH1D* var);


/**
 * Calculates the systematic errors on optimal, storing the results in graph.
 */
void CalcSystematics (TGraphAsymmErrors* graph, const TH1* optimal, const TH1* sys_hi, const TH1* sys_lo);


/**
 * Calculates the systematic errors on optimal, storing the results in graph.
 */
void CalcSystematics (TGraphAsymmErrors* graph, const TGraphAsymmErrors* optimal, const TGraph* sys_hi, const TGraph* sys_lo);


/**
 * Sets the bin contents in target as the error / central values in centralValues
 */
void SaveRelativeErrors (TH1D* errors, TH1D* centralValues);


/**
 * Creates a projection of a TH3 with the specified axes by integrating between min and max on the 3rd axis of the TH3.
 */
TH2D* Project2D (TString name, TH3D* h3, const TString xaxis, const TString yaxis, const int min, const int max, const bool exclusive = false);


/**
 * Separates each point on a TGraphAsymmErrors by delta along the x axis, so that the errors don't overlap.
 */
void deltaize (TGraphAsymmErrors* tg, const double delta = 0, const bool logx = false);


/**
 * Makes a TGraphAsymmErrors from the input histogram.
 */
TGraphAsymmErrors* make_graph (TH1* h, const float cutoff = -1);


/**
 * Recenters a TGraphAsymmErrors point for a log scale.
 */
void RecenterGraph (TGraphAsymmErrors* g);


/**
 * Returns the appropriate file in the given directory.
 * For MC, inFileName MUST be specified.
 */
TFile* GetFile (const char* directory, const int dataSet, const bool isMC, const char* inFileName = "");


/**
 * Returns an abbreviated, unique identifier for a given dataset.
 */
TString GetIdentifier (const int dataSet, const char* inFileName, const bool isMC, const bool isSignalOnlySample, const bool periodA);


/**
 * Returns the TProfile of an input histogram along the x axis. Can use either statistical mean or gaussian mean.
 */
TH1D* GetProfileX (const TString name, TH2D* hist, const int nbinsx, const double* xbins, const bool useFit = false, const double* xlo = NULL, const double* xhi = NULL);


/**
 * Returns the TProfile of an input histogram along the y axis. Can use either statistical mean or gaussian mean.
 */
TH1D* GetProfileY (const TString name, TH2D* hist, const int nbinsy, const double* ybins, const bool useFit = false, const double* ylo = NULL, const double* yhi = NULL);


/**
 * Returns a histogram with the TProfile of data over the TProfile of MC along either the x or y axes. Can use either the statistical or gaussian mean.
 */
TH1D* GetDataOverMC (const TString name, TH2D* data, TH2D* mc, const int numbins, const double* bins, const bool useFit = false, const TString axis = "x", const double* lo = NULL, const double* hi = NULL);


/**
 * Converts a TProfile to a TH1D.
 */
TH1D* TProfile2TH1D (const char* name, TProfile* p, const int nx, const double* x);


/**
 * Reflects the contents of h around the n-th bin in x.
 */
void GetReflectionX (TH1D* h, const int n);

} // end namespace

#endif
