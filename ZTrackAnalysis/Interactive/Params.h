#ifndef __Params_h__
#define __Params_h__

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

typedef TGraphAsymmErrors TGAE;

//const Color_t colors[10] =    {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta, kViolet-3, kCyan+1, kOrange+1, kGreen-7, kAzure+7};
const Color_t colors[10] =    {kBlack, kRed+1, kAzure+2, kGreen+2, kOrange+1, kViolet-3, kCyan+1, kMagenta, kGreen-7, kBlue+1};
const Color_t fillColors[10] = {kGray, kRed-10, kBlue-10, kGreen-10, kOrange, kGreen-2, kCyan-6, kMagenta-10, kViolet-2, kGray};
//const Color_t fillColors[10] = {kGray, kRed-10, kAzure+10, kGreen-10, kMagenta-10, kGreen-2, kCyan-6, kOrange, kViolet-2, kGray};
const float fillAlpha = 1;

//enum XAxisVar { trkPt = 0, xZh = 1}
//enum BinVar { trkPt = 0, dphi = 1, xZh = 2, zPt = 3}

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

const int nFCalBins = 250;
const double* fcalBins = linspace (0, 5000, nFCalBins);

const int numZMissingPtBins = 16;
const double* zMissingPtBins = doubleLogSpace (1, 150, 8);
//const double* zMissingPtBins = linspace (-150, 150, numZMissingPtBins);
const double phiTrkBins[3] = {0, pi/8, pi/2};
const int numPhiTrkBins = sizeof (phiTrkBins) / sizeof (phiTrkBins[0]) - 1;
const double etaTrkBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
const int numEtaTrkBins = sizeof (etaTrkBins) / sizeof (etaTrkBins[0]) - 1;
const int numFinerEtaTrkBins = 40;
const double* finerEtaTrkBins = linspace (-2.5, 2.5, numFinerEtaTrkBins);

// Centrality cuts in GeV:  80%  -  40%  -  15%  -   0%
//const double centBins[4] = {63.719, 875.41, 2476.58, 5000};
//const int centCuts[4] =    {80,     40,     15,      0};

// Centrality cuts in GeV:    80%   -   70%   -   60%  -   50%  -   40%  -   30%  -   20%  -   10%  -   0%
//const double centBins[9]   = {63.719,   144.14,   289.595, 525.092, 875.41,  1368.75, 2046.51, 2989.31, 5000};
//const double centCuts[9]   = {80,       70,       60,      50,      40,      30,      20,      10,      0};
//const double centThicks[8] = {0.223102, 0.564885, 1.2811,  2.63435, 4.94611, 8.63844, 14.3345, 23.3523};

// Centrality cuts in GeV:  80%  -  30%  -  10%  -   0%
//const double centBins[4] = {63.719, 1368.75, 2989.31, 5000}; // 2015 recommendations
const double centBins[4] = {66.402, 1378.92, 2995.94, 5000}; // updated 2015 recommendations
const int centCuts[4]    = {80,     30,     10,      0};

const double finerCentBins[10] = {66.402, 148.625, 296.17, 533.608, 885.172, 1378.92, 2055.77, 2995.94, 3622.6, 5000};
const int finerCentCuts[10]    = {80,     70,      60,     50,      40,      30,      20,      10,      5,      0};

const int numCentBins = sizeof (centBins) / sizeof (centBins[0]); // no minus 1 to include pp bin
const int numFinerCentBins = sizeof (finerCentBins) / sizeof (finerCentBins[0]);

const double phiLowBins[3] = {0,      3*pi/4,     15*pi/16};
const double phiHighBins[3] = {pi/2,  15*pi/16,   pi};
const int numPhiBins = sizeof (phiLowBins) / sizeof (phiLowBins[0]);

//const double zPtBins[5] = {0, 5, 20, 40, 10000};
const double zPtBins[6] = {0, 5, 15, 30, 60, 10000};
const int nPtZBins = sizeof (zPtBins) / sizeof (zPtBins[0]) - 1;

const int nZHBins = 6;
const double* zHBins = logspace (1./15., 1, nZHBins);

const double trk_min_pt = 1;
const double trk_max_pt = 60;
const int nPtTrkBins = 7;

double** init_ptTrkBins () {
  double** _ptTrkBins = Get1DArray <double*> (nPtZBins);

  for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    if (trk_min_pt > zPtBins[iPtZ]);
      _ptTrkBins[iPtZ] = logspace (trk_min_pt, trk_max_pt, nPtTrkBins);
    else
      _ptTrkBins[iPtZ] = logspace (trk_min_pt, fmin (trk_max_pt, zPtBins[iPtZ]), nPtTrkBins);
  }
  return _ptTrkBins;
}
double** ptTrkBins = init_ptTrkBins (); // iPtZ, iPtTrk
//logspace (trk_min_pt, 65, nPtTrkBins);


void PrintPtBins () {
  for (int i = 0; i <= nPtTrkBins; i++) {
    cout << ptTrkBins[i] << endl;
  }
}


long gcd(long a, long b)
{
    if (a == 0)
        return b;
    else if (b == 0)
        return a;

    if (a < b)
        return gcd(a, b % a);
    else
        return gcd(b, a % b);
}


const char* GetPiString (double phi) {
  phi = phi / pi;
  while (phi < 0)
    phi += 2;
  while (phi >= 2)
    phi -= 2;

  if (phi == 0)
    return "0";

  int denom, numer;

  double integerFloor = floor (phi);
  double fracPart = phi - integerFloor;

  if (integerFloor == phi) {
    denom = 1;
    numer = phi;
  }
  else {
    const long precision = 1000000000;
    long gcd_ = gcd (round (fracPart * precision), precision);

    denom = (int)(precision / gcd_);
    numer = (int)(round (fracPart * precision) / gcd_);
  }

  if (numer == 1 && denom == 1)
    return "#pi";
  else if (numer == 1)
    return Form ("#pi/%i", denom);
  else if (denom == 1)
    return Form ("%i#pi", numer);
  else
    return Form ("%i#pi/%i", numer, denom);
}


short GetidPhi (const float dphi) {
  short idPhi = 0;
  while (idPhi < numPhiBins) {
    if (phiLowBins[idPhi] < dphi && dphi < phiHighBins[idPhi])
      break;
    else
      idPhi++;
  }
  return idPhi;
}


short GetiZH (const float zH) {
  short iZH = 0;
  while (iZH < nZHBins) {
    if (zHBins[iZH+1] < zH)
      iZH++;
    else
      break;
  }
  return iZH;
}

#endif
