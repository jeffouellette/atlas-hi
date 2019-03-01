#ifndef __Params_h__
#define __Params_h__

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

const Color_t colors[10] =    {kBlack, kRed, kBlue, kGreen+2, kMagenta, 8, kCyan+1, kOrange, kViolet, kGray};
const Color_t fillColors[10] = {kGray, kRed-10, kBlue-10, kGreen-10, kMagenta-10, kGreen-2, kCyan-6, kOrange, kViolet-2, kGray};
const float fillAlpha = 1;

const int nXZTrkBins = 75;
const double* xZTrkBins = logspace (2e-3, 2e1, nXZTrkBins);

const int nPtTrkBins = 6;
const double* ptTrkBins = logspace (1, 80, nPtTrkBins);

double* doubleLogSpace (const double lo, const double hi, const int num) {
  double* arr = Get1DArray <double> (2*num+1);
  double* logarr = logspace (lo, hi, num-1);
  for (int i = 0; i < num; i++) {
    arr[num+i+1] = logarr[i];
    arr[i] = -logarr[num-1-i];
  }
  arr[num] = 0;
  delete logarr;
  return arr;
}

const int nFCalBins = 300;
const double* fcalBins = linspace (0, 6000, nFCalBins);

const int numZMissingPtBins = 16;
const double* zMissingPtBins = doubleLogSpace (1, 150, 8);
//const double* zMissingPtBins = linspace (-150, 150, numZMissingPtBins);
const double phiTrkBins[3] = {0, pi/8, pi/2};
const int numPhiTrkBins = sizeof (phiTrkBins) / sizeof (phiTrkBins[0]) - 1;

// Centrality cuts in GeV:  80%  -  30%  -   0%
//const double centBins[3] = {63.719, 1368.75, 5000};
//const int centCuts[3] =    {80,     30,      0};

// Centrality cuts in GeV:  80%  -  40%  -  15%  -   0%
//const double centBins[4] = {63.719, 875.41, 2476.58, 5000};
//const int centCuts[4] =    {80,     40,     15,      0};

// Centrality cuts in GeV:  80%  -  30%  -  10%  -   0%
const double centBins[4] = {63.719, 1368.75, 2989.31, 5000};
const int centCuts[4] =    {80,     30,     10,      0};

// Centrality cuts in GeV:    80%   -   70%   -   60%  -   50%  -   40%  -   30%  -   20%  -   10%  -   0%
//const double centBins[9]   = {63.719,   144.14,   289.595, 525.092, 875.41,  1368.75, 2046.51, 2989.31, 5000};
//const double centCuts[9]   = {80,       70,       60,      50,      40,      30,      20,      10,      0};
//const double centThicks[8] = {0.223102, 0.564885, 1.2811,  2.63435, 4.94611, 8.63844, 14.3345, 23.3523};
const int numCentBins = sizeof (centBins) / sizeof (centBins[0]); // no minus 1 to include pp bin

const double phiLowBins[3] = {0, pi/2, 7*pi/8};
const double phiHighBins[3] = {pi/2, 7*pi/8, pi};
const int numPhiBins = sizeof (phiLowBins) / sizeof (phiLowBins[0]);

const double zPtBins[3] = {0, 5, 10000};
const int nPtZBins = sizeof (zPtBins) / sizeof (zPtBins[0]) - 1;

#endif
