#ifndef __Params_h__
#define __Params_h__

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

typedef TGraphAsymmErrors TGAE;

const Style_t markerStyles[6] = {kFullCircle, kFullSquare, kFullCrossX, kOpenCircle, kOpenSquare, kOpenCrossX};
//const Color_t colors[10] =    {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta, kViolet-3, kCyan+1, kOrange+1, kGreen-7, kAzure+7};
const Color_t colors[10] =    {kBlack, kRed+1, kAzure+2, kGreen+2, kOrange-3, kViolet-3, kCyan+1, kMagenta, kGreen-7, kBlue+1};
const Color_t fillColors[10] = {kGray, kRed-10, kBlue-10, kGreen-10, kOrange, kViolet-2, kCyan-6, kMagenta-10, kGreen-2, kGray};
//const Color_t fillColors[10] = {kGray, kRed-10, kAzure+10, kGreen-10, kMagenta-10, kGreen-2, kCyan-6, kOrange, kViolet-2, kGray};
const float fillAlpha = 1;

//const int GroupA[8] = {365502, 365512, 365573, 365602, 365627, 365678, 365681, 365709};
//const int GroupB[7] = {365752, 365834, 365914, 365932};
//const int GroupC[3] = {366011, 366029, 366092};
//const int GroupD[6] = {366142, 366268, 366337, 366383, 366413, 366476};
//const int GroupE[9] = {366526, 366528, 366627, 366691, 366754, 366805};
//const int GroupF[6] = {366860, 366878, 366919, 366931, 366994, 367023};
//const int GroupG[6] = {367099, 367134, 367165, 367170, 367233, 367273};
//const int GroupH[6] = {367318, 367321, 367363, 367364, 367365, 367384};

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

float max_iaa = 1.4;
float max_icp = 1.4;
float max_rel_sys = 0.4;

int mixingFraction = 40;

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


// Centrality cuts in GeV:  80%  -  30%  -  10%  -   0%
//const double centBins[4] = {63.719, 1368.75, 2989.31, 5000}; // 2015 recommendations
const double centBins[4] = {66.402, 1378.92, 2995.94, 5000}; // updated 2015 recommendations
const int centCuts[4]    = {80,     30,     10,      0};
const int numCentBins = sizeof (centBins) / sizeof (centBins[0]); // no minus 1 to include pp bin

short GetCentBin (const float fcal_et) {
  short i = 0;
  while (i < numCentBins) {
    if (fcal_et < centBins[i])
      break;
    i++;
  }
  return i;
}

const double finerCentBins[10] = {66.402, 148.625, 296.17, 533.608, 885.172, 1378.92, 2055.77, 2995.94, 3622.6, 5000};
const int finerCentCuts[10]    = {80,     70,      60,     50,      40,      30,      20,      10,      5,      0};
const int numFinerCentBins = sizeof (finerCentBins) / sizeof (finerCentBins[0]);
short GetFinerCentBin (const float fcal_et) {
  short i = 0;
  while (i < numFinerCentBins) {
    if (fcal_et < finerCentBins[i])
      break;
    i++;
  }
  return i;
}


const double phiLowBins[3] = {0,      3*pi/4,     15*pi/16};
const double phiHighBins[3] = {pi/2,  15*pi/16,   pi};
const int numPhiBins = sizeof (phiLowBins) / sizeof (phiLowBins[0]);


//const double zPtBins[5] = {0, 5, 20, 40, 10000};
const double zPtBins[6] = {0, 5, 15, 30, 60, 10000};
const int nPtZBins = sizeof (zPtBins) / sizeof (zPtBins[0]) - 1;
short GetPtZBin (const float zPt) {
  short iPtZ = 0;
  while (iPtZ < nPtZBins) {
    if (zPt < zPtBins[iPtZ+1])
      break;
    else
      iPtZ++;
  }
  return iPtZ;
}


const double trk_min_pt = 1;
const double trk_max_pt = 60;


const int maxNPtTrkBins = 6;
double allPtTrkBins[7] = {1, 2, 4, 8, 15, 30, 60};
int* nPtTrkBins = Get1DArray <int> (nPtZBins);

double** init_ptTrkBins () {
  double** _ptTrkBins = Get1DArray <double*> (nPtZBins);
  for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    nPtTrkBins[iPtZ] = 0;
    while (nPtTrkBins[iPtZ] < maxNPtTrkBins && allPtTrkBins[nPtTrkBins[iPtZ]] < zPtBins[iPtZ])
      nPtTrkBins[iPtZ]++;

    if (nPtTrkBins[iPtZ] == 0)
      nPtTrkBins[iPtZ] = 1;
  }

  for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    _ptTrkBins[iPtZ] = Get1DArray<double> (nPtTrkBins[iPtZ]+1);
    for (int iPtTrk = 0; iPtTrk <= nPtTrkBins[iPtZ]; iPtTrk++) {
      _ptTrkBins[iPtZ][iPtTrk] = allPtTrkBins[iPtTrk];
    }
  }

  return _ptTrkBins;
}
double** ptTrkBins = init_ptTrkBins (); // iPtZ, iPtTrk


const int maxNXHZBins = 6;
double allXHZBins[7] = {1./60., 1./30., 1./15., 1./8., 1./4., 1./2., 1.};
int* nXHZBins = Get1DArray <int> (nPtZBins);

double** init_xHZBins () {
  double** _xHZBins = Get1DArray <double*> (nPtZBins);

  for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    nXHZBins[iPtZ] = maxNXHZBins;
    while (nXHZBins[iPtZ] > 0 && zPtBins[iPtZ] > 0 && allXHZBins[maxNXHZBins-nXHZBins[iPtZ]] < 1./zPtBins[iPtZ])
      nXHZBins[iPtZ]--;

    if (nXHZBins[iPtZ] == 0)
      nXHZBins[iPtZ] = 1;
  }

  for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    _xHZBins[iPtZ] = Get1DArray<double> (nXHZBins[iPtZ]+1);
    for (int iXHZ = 0; iXHZ <= nXHZBins[iPtZ]; iXHZ++) {
      _xHZBins[iPtZ][iXHZ] = allXHZBins[iXHZ+(maxNXHZBins-nXHZBins[iPtZ])];
    }
  }

  return _xHZBins;
}
double** xHZBins = init_xHZBins ();


void PrintPtBins (const short iPtZ) {
  for (int i = 0; i <= nPtTrkBins[iPtZ]; i++) {
    cout << ptTrkBins[iPtZ][i] << endl;
  }
}

void PrintXHZBins (const short iPtZ) {
  for (int i = 0; i <= nXHZBins[iPtZ]; i++) {
    cout << xHZBins[iPtZ][i] << endl;
  }
}

void LatexXHZBins () {
  cout << "\\hline" << endl;
  for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    cout << "\\multicolumn{2}{";
    if (iPtZ == 2)
      cout << "|";
    cout << " c |";
    if (iPtZ != nPtZBins-1)
      cout << "|";
    cout << "}{";
    cout << "\\pt^\\mathrm{Z} = \\SI{" << zPtBins[iPtZ] << "}{\\GeV}} & ";
  }
  cout << "\\\\ \\hline" << endl;
  for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    cout << "$x_\\mathrm{hZ}$ & $\\pt^\\mathrm{ch}$ [\\GeV] ";
    if (iPtZ != nPtZBins-1)
      cout << "& ";
  }
  cout << "\\\\ \\hline" << endl;
  for (int i = 0; i <= maxNXHZBins; i++) {
    for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      if (i < nXHZBins[iPtZ])
        cout << xHZBins[iPtZ][i] << " & " << xHZBins[iPtZ][i] * zPtBins[iPtZ] << " ";
      if (iPtZ != nPtZBins-1)
        cout << "& ";
    }
    //if (i != nXHZBins[iPtZ])
    cout << "\\\\ \\hline";
    cout << endl;
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


short GetiZH (const float zH, const int iPtZ) {
  short iZH = 0;
  while (iZH < nXHZBins[iPtZ]) {
    if (xHZBins[iPtZ][iZH+1] < zH)
      iZH++;
    else
      break;
  }
  return iZH;
}


const float superFineCentBins[66] = {
    66.402, // 80%
  //  72.411, // 79%
    78.834, // 78%
  //  85.729, // 77%
    93.071, // 76%
  // 100.928, // 75%
   109.333, // 74%
  // 118.249, // 73%
   127.745, // 72%
  // 137.874, // 71%
   148.625, // 70%
  // 160.033, // 69%
   172.097, // 68%
  // 184.922, // 67%
   198.428, // 66%
  // 212.684, // 65%
   227.751, // 64%
  // 243.588, // 63%
   260.219, // 62%
  // 277.742, // 61%
   296.17,  // 60%
  // 315.523, // 59%
   335.738, // 58%
  // 356.885, // 57%
   378.968, // 56%
  // 402.144, // 55%
   426.354, // 54%
  // 451.509, // 53%
   477.734, // 52%
  // 505.085, // 51%
   533.608, // 50%
  // 563.263, // 49%
   594.07,  // 48%
  // 626.047, // 47%
   659.269, // 46%
  // 693.606, // 45%
   729.251, // 44%
  // 766.305, // 43%
   804.607, // 42%
  // 844.192, // 41%
   885.172, // 40%
  // 927.582, // 39%
   971.487, // 38%
  //1016.86,  // 37%
  1063.71,  // 36%
  //1112.2,   // 35%
  1162.33,  // 34%
  //1213.91,  // 33%
  1267.07,  // 32%
  //1322.12,  // 31%
  1378.92,  // 30%
  1437.29,  // 29%
  1497.54,  // 28%
  1560.05,  // 27%
  1624.34,  // 26%
  1690.47,  // 25%
  1759.06,  // 24%
  1829.74,  // 23%
  1902.73,  // 22%
  1978.02,  // 21%
  2055.77,  // 20%
  2136.17,  // 19%
  2218.88,  // 18%
  2304.47,  // 17%
  2393.11,  // 16%
  2484.75,  // 15%
  2579.56,  // 14%
  2677.6,   // 13%
  2779.68,  // 12%
  2885.71,  // 11%
  2995.94,  // 10%
  0.5*(2995.94+3110.27), // middle bin
  3110.27,  //  9%
  0.5*(3110.27+3229.67), // middle bin
  3229.67,  //  8%
  0.5*(3229.67+3354.66), // middle bin
  3354.66,  //  7%
  0.5*(3354.66+3485.57), // middle bin
  3485.57,  //  6%
  0.5*(3485.57+3622.6),   // middle bin
  3622.6,   //  5%
  0.5*(3622.6+3767.),     // middle bin
  3767.,    //  4%
  0.5*(3767.+3920.41),    // middle bin
  3920.41,  //  3%
  0.5*(3920.41+4083.38),  // middle bin
  4083.38,  //  2%
  0.5*(4083.38+4263.72),  // middle bin
  4263.72,  //  1%
  0.5*(4263.72+5000),     // middle bin
  5000      //  0%,   entry in array is numSuperFineCentBins-1
};
const short numSuperFineCentBins = sizeof (superFineCentBins) / sizeof (superFineCentBins[0]);

short GetSuperFineCentBin (const float fcal_et) {
  short i = 0;
  while (i < numSuperFineCentBins) {
    if (fcal_et < superFineCentBins[i])
      break;
    i++;
  }
  return i;
}


#endif
