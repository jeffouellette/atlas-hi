#include "GlobalParams.h"
#include <iostream>
#include <iomanip>

/**
 * Contains useful variables and directory information for 2016 pPb data analyses.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

namespace atlashi {

TString workPath = ""; // Home analysis directory, should be modified in code outside this path structure
TString externalWorkPath = ""; // External drive storage directory, should be modified in code below
TString rootPath = ""; // Where analyzed *.root files are stored. Different analysis modules have different subdirectories here.
TString dataPath = ""; // Where the *.root raw data files (from the CERN grid) are stored.
TString plotPath = ""; // Where plots are stored.
TString ptPath = ""; // Where the pt analysis module output is stored.
TString trigPath = "";  // Where the trigger fire count module output is stored.
TString effPath = ""; // Where the trigger efficiency module output is stored.
TString xPath = ""; // Where the xa/xp module output is stored.
TString RpPbPath = ""; // Where the R_pPb module output is stored.

int numtrigs = 0; // number of triggers that have been initialized

const char* FormatMeasurement (double val, double err, const int n) {
  assert (n > 0); // sanity check

  TString out = "";

  std::string valStr = Form ("%g", val);
  std::string errStr = Form ("%g", err);

  if (err == 0) {
    short valDec = 0;
    while (valDec < (short)valStr.length() && valStr[valDec] != '.')
     valDec++;
    if (valDec == (short)valStr.length()) {
     valStr = valStr + ".";
    }

    if (valDec > n) { // if decimal is after least significant digit
      const double factorOfTen = pow (10, n-valDec); // e.g., for "1520" and n=2, get 0.01
      val = floor (factorOfTen * val + 0.5) / factorOfTen;
    }
    else { // if decimal is before least significant digit, e.g. for "15.24" and n=3.
      const double factorOfTen = pow (10, n-valDec+1); // e.g. for 15.24 and n=3 get 10
      val = floor (factorOfTen * val + 0.5) / factorOfTen;
    }
    return Form ("%g", val);
  }

  if (err < 1) {
   // find the first significant digit
   int errStart = 0;
   while (errStr[errStart] == '0' || errStr[errStart] == '.')
    errStart++;
   errStart = errStart - 2; // first 2 characters are necessarly "0."

   // round the value and error to the appropriate decimal place
   const double factorOfTen = pow(10, errStart+n);
   val = floor (factorOfTen * val + 0.5) / factorOfTen;
   err = floor (factorOfTen * err + 0.5) / factorOfTen;

   // recast to string
   valStr = Form ("%g", val);
   errStr = Form ("%g", err);

   // find where the decimal place is, append it if it's not present
   short valDec = 0;
   while (valDec < (short)valStr.length() && valStr[valDec] != '.')
    valDec++;
   if (valDec == (short)valStr.length()) {
    valStr = valStr + ".";
   }

   // pad with zeroes
   while ((short)valStr.length() < valDec + errStart + 1 + n)
    valStr = valStr + "0";
   while ((short)errStr.length() < errStart + 2 + n)
    errStr = errStr + "0";

   //// find where we are truncating the value
   //const short valCut = valDec + errStart + 2 + n;

   //// now truncate
   //valStr = valStr.substr (0, valCut);
   //errStr = errStr.substr (0, errStart + 2 + n);
  }
  else { // now if err>1
   // find the decimal place
   short errDec = 0;
   while (errDec < (short)errStr.length() && errStr[errDec] != '.')
    errDec++;

   // round the value and error to the appropriate decimal place
   if (errDec < (short)errStr.length() - 1) {
    const double factorOfTen = pow (10., n - errDec);
    if (n - errDec >= 0) {
     val = floor (factorOfTen * val + 0.5) / factorOfTen;
     err = floor (factorOfTen * err + 0.5) / factorOfTen;
    } else {
     val = floor (factorOfTen * val) / factorOfTen;
     err = floor (factorOfTen * err) / factorOfTen;
    }
   }

   // recast to string
   valStr = Form ("%g", val);
   errStr = Form ("%g", err);
  }

  // now save the value string to the output
  out = valStr + " #pm " + errStr + "";

  return out.Data();
}


void SetupDirectories (const TString dataSubDir, const TString thisWorkPath) {

  if (thisWorkPath[thisWorkPath.Length()-1] != '/') {
    workPath = homePath + thisWorkPath + "/";
    externalWorkPath = drivePath + thisWorkPath + "/";
  } else {
    workPath = homePath + thisWorkPath;
    externalWorkPath = drivePath + thisWorkPath;
  }

  rootPath = workPath + "rootFiles/" + dataSubDir;
  dataPath = externalWorkPath + "data/" + dataSubDir;
  plotPath = workPath + "Plots/" + dataSubDir;
  ptPath = rootPath + "ptData/";
  trigPath = rootPath + "trigData/";
  effPath = rootPath + "effData/";
  xPath = rootPath + "xData/";
  RpPbPath = rootPath + "RpPbData/";

  return;
}


double* linspace (double lo, double hi, int num) {
  double* arr = new double[num+1];
  double delta = ((double)(hi)-(double)(lo))/(double)(num);
  for (int i = 0; i <= num; i++) {
    arr[i] = lo + i * delta;
  }
  return arr;
}


double* logspace (double lo, double hi, int num) {
  double loghi = TMath::Log2(hi);
  if (lo == 0) {
    double* arr = linspace(TMath::Log2(hi/(100*num)), loghi, num);
    for (int i = 0; i <= num; i++) {
      arr[i] = TMath::Power(2, arr[i]);
    }
    return arr;
  } else {
    double loglo = TMath::Log2(lo);
    double* arr = linspace(loglo, loghi, num);
    for (int i = 0; i <= num; i++) {
      arr[i] = TMath::Power(2, arr[i]);
    }
    return arr;
  }
}


double InTwoPi (double phi) {
  while (phi < 0 || 2*pi <= phi) {
   if (phi < 0) phi += 2*pi;
   else phi -= 2*pi;
  }
  return phi;
}


double DeltaPhi (double phi1, double phi2) {
  phi1 = InTwoPi(phi1);
  phi2 = InTwoPi(phi2);
  double dphi = abs(phi1 - phi2);
  while (dphi > pi) dphi = abs(dphi - 2*pi);
  return dphi;
}


double DeltaR (const double eta1, const double eta2, const double phi1, const double phi2 ) {
 const double deta = eta1 - eta2;
 const double dphi = DeltaPhi(phi1, phi2);
 return sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
}


bool InDisabledHEC (const double eta, double phi) {
  phi = InTwoPi (phi);
  return lowerEtaCut < eta &&
         eta < upperEtaCut &&
         lowerPhiCut < phi &&
         phi < upperPhiCut;
}


bool InEMCal (const float eta) {
  return TMath::Abs(eta) < 1.37 || (1.56 < TMath::Abs(eta) && TMath::Abs(eta) < 2.47);
}


bool InHadCal (const float eta, const float R) {
  return TMath::Abs(eta) <= 4.9 - R;
}


void CalcSystematics (TGraphAsymmErrors* graph, TH1* optimal, TH1* sys_hi, TH1* sys_lo) {
  for (int xbin = 1; xbin <= optimal->GetNbinsX(); xbin++) {
   const double content = optimal->GetBinContent(xbin);
   const double lo = sys_lo->GetBinContent(xbin);
   const double hi = sys_hi->GetBinContent(xbin);

   const double err_lo = content - lo;
   const double err_hi = hi - content;

   graph->SetPoint(xbin-1, optimal->GetBinCenter(xbin), content);
   graph->SetPointEXlow(xbin-1, optimal->GetBinCenter(xbin) - optimal->GetBinLowEdge(xbin));
   graph->SetPointEXhigh(xbin-1, optimal->GetBinLowEdge(xbin+1) - optimal->GetBinCenter(xbin));

   if (err_lo < 0 && err_hi < 0) {
    graph->SetPointEYlow(xbin-1, -err_hi);
    graph->SetPointEYhigh(xbin-1, -err_lo);
   }
   else {
    graph->SetPointEYlow(xbin-1, err_lo);
    graph->SetPointEYhigh(xbin-1, err_hi);
   }
  }
  return;
}

} // end namespace
