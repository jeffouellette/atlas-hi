//
//   @file    AtlasUtils.h         
//   
//
//   @author M.Sutton
// 
//   Copyright (C) 2010 Atlas Collaboration
//
//   $Id: AtlasUtils.h, v0.0   Thu 25 Mar 2010 10:34:20 CET $


#ifndef __AtlasUtils_h__
#define __AtlasUtils_h__

#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TProfile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TEfficiency.h>

typedef TGraphAsymmErrors TGAE;

void SetStyle ();

void ATLAS_LABEL (double x, double y, Color_t color=1); 

TGraphErrors* myTGraphErrorsDivide (TGraphErrors* g1, TGraphErrors* g2);

TGraphAsymmErrors* myTGraphErrorsDivide (TGraphAsymmErrors* g1, TGraphAsymmErrors* g2);

TGraphAsymmErrors* myMakeBand (TGraphErrors* g0, TGraphErrors* g1, TGraphErrors* g2);

void myAddtoBand (TGraphErrors* g1, TGraphAsymmErrors* g2);

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


///**
// * Modifies the directory strings to point to the correct locations.
// */
//void SetupDirectories (const TString dataSubDir, const TString thisWorkPath);
//
//
///**
// * Clears sub-directory information from the directory strings
// */
//void ResetDirectories ();



///**
// * Returns the appropriate file in the given directory.
// * For MC, inFileName MUST be specified.
// */
//TFile* GetFile (const char* directory, const int dataSet, const bool isMC, const char* inFileName = "");


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


#endif // __AtlasUtils_h__
