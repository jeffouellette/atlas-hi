#include "Utilities.h"
#include "GlobalParams.h"

#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TF1.h>

using namespace std;

namespace atlashi {

/**
 * Returns a TString summarizing a measurement.
 * By default, there is only 1 significant digit in the error.
 * E.g., FormatMeasurement (40.58, 1.29) returns "40#pm1".
 * Or, FormatMeasurement (40.58, 1.29, 2) returns "40.6#pm1.3".
 */
const char* FormatMeasurement (double val, double err, const int n) {
  assert (n > 0); // sanity check

  TString out = "";

  string valStr = Form ("%g", val);
  string errStr = Form ("%g", err);

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
   const double factorOfTen = pow (10, errStart+n);
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


/**
 * Modifies the directory strings to point to the correct locations.
 */
void SetupDirectories (const TString dataSubDir, const TString thisWorkPath) {
  ResetDirectories ();

  workPath = homePath + "/" + thisWorkPath;
  externalWorkPath = drivePath + "/" + thisWorkPath;
  intWorkPath = intPath + "/" + thisWorkPath;

  rootPath = intWorkPath + "/rootFiles/" + dataSubDir;
  dataPath = externalWorkPath + "/data/" + dataSubDir;
  plotPath = workPath + "/Plots/" + dataSubDir;
  ptPath = rootPath + "/ptData/";
  trigPath = rootPath + "/trigData/";
  effPath = rootPath + "/effData/";
  xPath = rootPath + "/xData/";
  RpPbPath = rootPath + "/RpPbData/";

  return;
}


/**
 * Clears sub-directory information from the directory strings
 */
void ResetDirectories () {
  workPath = ""; // Home analysis directory, should be modified in code outside this path structure
  externalWorkPath = ""; // External drive storage directory, should be modified in code below
  intWorkPath = ""; // Base directory for intermediate root files, should be modified in code outside this path structure
  rootPath = ""; // Where analyzed *.root files are stored. Different analysis modules have different subdirectories here.
  dataPath = ""; // Where the *.root raw data files (from the CERN grid) are stored.
  plotPath = ""; // Where plots are stored.
  ptPath = ""; // Where the pt analysis module output is stored.
  trigPath = "";  // Where the trigger fire count module output is stored.
  effPath = ""; // Where the trigger efficiency module output is stored.
  xPath = ""; // Where the xa/xp module output is stored.
  RpPbPath = ""; // Where the R_pPb module output is stored.
}


/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
double* linspace (double lo, double hi, int num) {
  double* arr = new double[num+1];
  double delta = ((double)(hi)-(double)(lo))/(double)(num);
  for (int i = 0; i <= num; i++) {
    arr[i] = lo + i * delta;
  }
  return arr;
}


/**
 * Returns a logarithmically spaced array, where the 0th element is lo and the num-th element is hi.
 */
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


/**
 * Returns the equivalent angle in the range 0 to 2pi.
 */
double InTwoPi (double phi) {
  while (phi < 0 || 2*pi <= phi) {
   if (phi < 0) phi += 2*pi;
   else phi -= 2*pi;
  }
  return phi;
}


/**
 * Returns the difference between two angles in 0 to pi.
 * If sign is true, will return a signed version such that phi2 = phi1 + dphi
 */
double DeltaPhi (double phi1, double phi2, const bool sign) {
  phi1 = InTwoPi(phi1);
  phi2 = InTwoPi(phi2);
  double dphi = abs(phi1 - phi2);
  while (dphi > pi) dphi = abs (dphi - 2*pi);

  if (sign && InTwoPi (phi2 + dphi) == phi1)
     dphi *= -1;

  return dphi;
}


/**
 * Returns dR between two eta, phi coordinates.
 */
double DeltaR (const double eta1, const double eta2, const double phi1, const double phi2 ) {
 const double deta = eta1 - eta2;
 const double dphi = DeltaPhi (phi1, phi2, false);
 return sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
}


/**
 * Returns true iff this eta, phi coordinate lies in the disabled HEC region.
 */
bool InDisabledHEC (const double eta, double phi, const double dr) {
  phi = InTwoPi (phi);
  return 1.5-dr < eta && eta < 3.2+dr &&
         pi-dr < phi && phi < 3*pi/2+dr;
}


/**
 * Returns true iff this eta lies within the EMCal.
 */
bool InEMCal (const float eta) {
  return TMath::Abs(eta) < 1.37 || (1.56 < TMath::Abs(eta) && TMath::Abs(eta) < 2.47);
}


/**
 * Returns true iff this object is within a given radius in the HCal.
 */
bool InHadCal (const float eta, const float R) {
  return TMath::Abs(eta) <= 4.9 - R;
}


/**
 * Sets all the errors in this histogram to 0.
 */
void ResetHistErrors (TH1D* h) {
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    h->SetBinError (ix, 0);
  }
}


/**
 * Sets all the errors in this TGAE to 0.
 */
void ResetTGAEErrors (TGAE* g) {
  for (int ix = 0; ix < g->GetN (); ix++) {
    g->SetPointEYhigh (ix, 0);
    g->SetPointEYlow (ix, 0);
  }
}


/**
 * Sets all the x errors in this TGAE to 0.
 */
void ResetXErrors (TGAE* tg) {
  for (int ix = 0; ix < tg->GetN (); ix++) {
    tg->SetPointEXlow (ix, 0);
    tg->SetPointEXhigh (ix, 0);
  }
  return;
}


/**
 * Adds nSigma statistical error variations to this histogram
 */
void AddStatVar (TH1D* h, const bool upvar, const float nSigma) {
  if (!h)
    return;
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    h->SetBinContent (ix, h->GetBinContent (ix) + (upvar ? 1.:-1.) * nSigma * h->GetBinError (ix));
  }
  return;
}



/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 */
void AddErrorsInQuadrature (TH1D* master, TH1D* sys) {
  for (int ix = 1; ix <= master->GetNbinsX (); ix++) {
    if (master->GetBinContent (ix) != sys->GetBinContent (ix))
      cout << "Warning: In Utilities.cxx::AddErrorsInQuadrature: unequal nominal bin contents between systematics! Results will be skewed." << endl;

    const float newErr = sqrt (pow (master->GetBinError (ix), 2) + pow (sys->GetBinError (ix), 2));
    master->SetBinError (ix, newErr);
  }
}


/**
 * Adds independent systematic errors in quadrature, storing the sum in master
 */
void AddErrorsInQuadrature (TGAE* master, TGAE* sys, const bool doXErrs) {
  for (int ix = 0; ix < master->GetN (); ix++) {
    //double xm=0, ym=0, xs=0, ys=0;
    //master->GetPoint (ix, xm, ym);
    //sys->GetPoint (ix, xs, ys);

    //if (xm != xs || ym != ys)
    //  cout << "Warning: In Utilities.cxx::AddErrorsInQuadrature: Bins don't match!" << endl;

    float newErr = sqrt (pow (master->GetErrorYhigh (ix), 2) + pow (sys->GetErrorYhigh (ix), 2));
    master->SetPointEYhigh (ix, newErr);
    newErr = sqrt (pow (master->GetErrorYlow (ix), 2) + pow (sys->GetErrorYlow (ix), 2));
    master->SetPointEYlow (ix, newErr);

    if (doXErrs) {
      newErr = sqrt (pow (master->GetErrorXhigh (ix), 2) + pow (sys->GetErrorXhigh (ix), 2));
      master->SetPointEXhigh (ix, newErr);
      newErr = sqrt (pow (master->GetErrorXlow (ix), 2) + pow (sys->GetErrorXlow (ix), 2));
      master->SetPointEXlow (ix, newErr);
    }
  }
}


/**
 * Calculates simple systematics as maximum variations on the nominal.
 * Intended for combining up/down variations in an expandable way.
 */
void CalcSystematics (TH1D* sys, TH1D* var) {
  for (int ix = 1; ix <= sys->GetNbinsX (); ix++) {
    const float newErr = fabs (var->GetBinContent (ix) - sys->GetBinContent (ix));
    sys->SetBinError (ix, fmax (newErr, sys->GetBinError (ix)));
  }
}


/**
 * Calculates simple systematics as maximum variations on the nominal.
 * Intended for combining up/down variations in an expandable way.
 * dir represents the variation direction, -1 = down, 0 = either, 1 = up
 */
void CalcSystematics (TGAE* sys, TH1D* var, const bool applyBothWays) {
  for (int ix = 0; ix < sys->GetN (); ix++) {
    double x, y;
    sys->GetPoint (ix, x, y);
    const float newErr = var->GetBinContent (ix+1) - y;

    if (applyBothWays || newErr > 0) {
      sys->SetPointEYhigh (ix, fmax (fabs (newErr), sys->GetErrorYhigh (ix)));
    }
    if (applyBothWays || newErr < 0) {
      sys->SetPointEYlow (ix, fmax (fabs (newErr), sys->GetErrorYlow (ix)));
    }
  }
}


void CalcSystematics (TGAE* graph, const TH1* optimal, const TH1* sys_hi, const TH1* sys_lo) {
  for (int ix = 1; ix <= optimal->GetNbinsX(); ix++) {
    const double content = optimal->GetBinContent (ix);
    const double lo = sys_lo->GetBinContent (ix);
    const double hi = sys_hi->GetBinContent (ix);

    const double err_lo = content - lo;
    const double err_hi = hi - content;

  

    graph->SetPoint (ix-1, optimal->GetBinCenter (ix), content);
    graph->SetPointEXlow (ix-1, optimal->GetBinCenter (ix) - optimal->GetBinLowEdge (ix));
    graph->SetPointEXhigh (ix-1, optimal->GetBinLowEdge (ix+1) - optimal->GetBinCenter (ix));

    if (err_lo < 0 && err_hi < 0) {
      graph->SetPointEYlow (ix-1, -err_hi);
      graph->SetPointEYhigh (ix-1, -err_lo);
    }
    else if (err_lo >= 0 && err_hi >= 0) {
      graph->SetPointEYlow (ix-1, err_lo);
      graph->SetPointEYhigh (ix-1, err_hi);
    }
    else {
      graph->SetPointEYlow (ix-1, sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
      graph->SetPointEYhigh (ix-1, sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
    }

  }
  return;
}


/**
 * Calculates the systematic errors on optimal, storing the results in graph.
 */
void CalcSystematics (TGAE* graph, const TGAE* optimal, const TGraph* sys_hi, const TGraph* sys_lo) {
  for (int ix = 0; ix < optimal->GetN(); ix++) {
    double x, y, xl, yl, xh, yh;
    optimal->GetPoint (ix, x, y);
    sys_lo->GetPoint (ix, xl, yl);
    sys_hi->GetPoint (ix, xh, yh);

    const double err_lo = y - yl;
    const double err_hi = yh - y;

    graph->SetPoint (ix, x, y);
    graph->SetPointEXlow (ix, optimal->GetErrorXlow (ix));
    graph->SetPointEXhigh (ix, optimal->GetErrorXhigh (ix));

    if (err_lo < 0 && err_hi < 0) {
      graph->SetPointEYlow (ix, -err_hi);
      graph->SetPointEYhigh (ix, -err_lo);
    }
    else if (err_lo >= 0 && err_hi >= 0) {
      graph->SetPointEYlow (ix, err_lo);
      graph->SetPointEYhigh (ix, err_hi);
    }
    else {
      graph->SetPointEYlow (ix, sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
      graph->SetPointEYhigh (ix, sqrt (0.5*(err_lo*err_lo + err_hi*err_hi)));
    }
  }
  return;
}


/**
 * Sets the bin contents in target as the error in errors / central values in centralValues
 */
void SaveRelativeErrors (TH1D* errors, TH1D* centralValues) {
  for (int ix = 1; ix <= errors->GetNbinsX (); ix++) {
    if (centralValues->GetBinContent (ix) != 0)
      errors->SetBinContent (ix, errors->GetBinError (ix) / centralValues->GetBinContent (ix));
    errors->SetBinError (ix, 0);
  }
}


/**
 * Sets the bin contents in highs and lows as the respective errors / central values in centralValues
 */
void SaveRelativeErrors (TGAE* errors, TGAE* centralValues, TH1D* highs, TH1D* lows) {
  for (int ix = 0; ix < centralValues->GetN (); ix++) {
    double x = 0, y = 0, eyhi = 0, eylo = 0;
    centralValues->GetPoint (ix, x, y);
    eyhi = errors->GetErrorYhigh (ix);
    eylo = errors->GetErrorYlow (ix);

    if (y != 0) {
      highs->SetBinContent (ix+1, eyhi / y);
      lows->SetBinContent (ix+1, -eylo / y);
    }
    highs->SetBinError (ix+1, 0);
    lows->SetBinError (ix+1, 0);
  }
}


TH2D* Project2D (TString name, TH3D* h3, const TString xaxis, const TString yaxis, const int min, const int max, const bool exclusive) {
  int nx, ny, nz;
  double* xbins = NULL;
  double* ybins = NULL;

  if (xaxis == yaxis || (!(xaxis == "x" || xaxis == "y" || xaxis == "z") || !(yaxis == "x" || yaxis == "y" || yaxis == "z")))
   return NULL;

  // determine which axis in the TH3 is the xaxis in the TH2
  if      (xaxis == "x") { nx = h3->GetNbinsX(); xbins = (double*)h3->GetXaxis()->GetXbins()->GetArray(); }
  else if (xaxis == "y") { nx = h3->GetNbinsY(); xbins = (double*)h3->GetYaxis()->GetXbins()->GetArray(); }
  else if (xaxis == "z") { nx = h3->GetNbinsZ(); xbins = (double*)h3->GetZaxis()->GetXbins()->GetArray(); }

  // determine which axis in the TH3 is the yaxis in the TH2
  if      (yaxis == "x") { ny = h3->GetNbinsX(); ybins = (double*)h3->GetXaxis()->GetXbins()->GetArray(); }
  else if (yaxis == "y") { ny = h3->GetNbinsY(); ybins = (double*)h3->GetYaxis()->GetXbins()->GetArray(); }
  else if (yaxis == "z") { ny = h3->GetNbinsZ(); ybins = (double*)h3->GetZaxis()->GetXbins()->GetArray(); }

  if (!xbins || !ybins) return NULL;

  // determine which axis in TH3 is the extra "z" axis
  if      (xaxis != "x" && yaxis != "x") nz = h3->GetNbinsX();
  else if (xaxis != "y" && yaxis != "y") nz = h3->GetNbinsY();
  else if (xaxis != "z" && yaxis != "z") nz = h3->GetNbinsZ();

  if (name == "") name = TString(h3->GetName()) + "_" + xaxis + yaxis; 
  TH2D* h2 = new TH2D (name, "", nx, xbins, ny, ybins);

  for (int ix = 1; ix <= nx; ix++) {
   for (int iy = 1; iy <= ny; iy++) {
    double content = 0;
    double var = 0;

    int h3ix, h3iy, h3iz;
    for (int iz = 1; iz <= nz; iz++) { // loop over extra "z" axis (not in TH2)
     if (!exclusive && (iz < min || max < iz))
      continue;
     else if (exclusive && min <= iz && iz <= max)
      continue;

     if      (xaxis == "x" && yaxis == "y") { h3ix = ix; h3iy = iy; h3iz = iz; }
     else if (xaxis == "x" && yaxis == "z") { h3ix = ix; h3iz = iy; h3iy = iz; }
     else if (xaxis == "y" && yaxis == "x") { h3iy = ix; h3ix = iy; h3iz = iz; }
     else if (xaxis == "y" && yaxis == "z") { h3iy = ix; h3iz = iy; h3ix = iz; }
     else if (xaxis == "z" && yaxis == "x") { h3iz = ix; h3ix = iy; h3iy = iz; }
     else if (xaxis == "z" && yaxis == "y") { h3iz = ix; h3iy = iy; h3ix = iz; }

     content += h3->GetBinContent (h3ix, h3iy, h3iz);
     var += pow (h3->GetBinError (h3ix, h3iy, h3iz), 2);
     
    } // end loop over extra axis

    h2->SetBinContent (ix, iy, content);
    h2->SetBinError (ix, iy, sqrt (var));
   }
  }

  if      (xaxis == "x") h2->GetXaxis()->SetTitle (h3->GetXaxis()->GetTitle());
  else if (xaxis == "y") h2->GetXaxis()->SetTitle (h3->GetYaxis()->GetTitle());
  else if (xaxis == "z") h2->GetXaxis()->SetTitle (h3->GetZaxis()->GetTitle());

  if      (yaxis == "x") h2->GetYaxis()->SetTitle (h3->GetXaxis()->GetTitle());
  else if (yaxis == "y") h2->GetYaxis()->SetTitle (h3->GetYaxis()->GetTitle());
  else if (yaxis == "z") h2->GetYaxis()->SetTitle (h3->GetZaxis()->GetTitle());

  return h2;
}


/**
 * Returns the appropriate file in the given directory.
 * For MC, inFileName MUST be specified.
 */
TFile* GetFile (const char* directory, const int dataSet, const bool isMC, const char* inFileName) {
  TFile* file = NULL;

  // First figure out the file we are looking for
  TString fileIdentifier;
  if (TString (inFileName) == "") {
   if (!isMC) fileIdentifier = to_string (dataSet);
   else {
    cout << "Error: In Utilities.cxx: Cannot identify this MC file! Will return null!" << endl;
    return NULL;
   }
  }
  else fileIdentifier = inFileName;

  // Now get the list of files
  const TString dataPathTemp = dataPath + "/" + directory + "/";
  TSystemDirectory dir (dataPathTemp.Data (), dataPathTemp.Data ());
  TList* sysfiles = dir.GetListOfFiles ();
  if (!sysfiles) {
   cout << "Error: In Utilities.cxx: Cannot get list of files! Will return null!" << endl;
   return NULL;
  }
  TSystemFile* sysfile;
  TString fname;
  TIter next (sysfiles);

  while ( (sysfile = (TSystemFile*)next ())) {
   fname = sysfile->GetName ();
   if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
    if (debugStatements) cout << "Status: In Utilities.cxx: Found " << fname.Data () << endl;
    
    if (fname.Contains (fileIdentifier)) {
     file = new TFile (dataPathTemp+fname, "READ");
     break;
    }
   }
  }

  if (!file) {
   cout << "Error: In Utilities.cxx: TFile not obtained for given data set. Will return null!" << endl;
   return NULL;
  }
  else return file;
}
  

/**
 * Returns an abbreviated, unique identifier for a given dataset.
 */
TString GetIdentifier (const int dataSet, const char* inFileName, const bool isMC, const bool isSignalOnlySample, const bool periodA) {
  if (!isMC) return to_string (dataSet);

  TString id = (periodA ? "pPb_" : "Pbp_");

  id = id + (isSignalOnlySample ? "Signal_" : "Overlay_");

  if (TString (inFileName).Contains ("42310") && TString (inFileName).Contains ("Slice")) { // gamma+jet
   if (dataSet <= 0) return "";
   id = id + "GammaJet_Slice" + to_string (dataSet);
  }
  else if (TString (inFileName).Contains ("ZeeJet")) { // Zee+jet
   if (dataSet < 0) return "";
   id = id + "ZeeJet" + (dataSet == 0 ? "" : "_Slice" + to_string (dataSet));
  }
  else if (TString (inFileName).Contains ("ZmumuJet")) { // Zmumu+jet
   if (dataSet != 0) return "";
   id = id + "ZmumuJet";
  }
  else if (TString (inFileName).Contains ("Dstar")) { // D* samples
   id = id + "Dstar";
   if (TString (inFileName).Contains ("Minus"))
    id = id + "Minus";
   else if (TString (inFileName).Contains ("Plus"))
    id = id + "Plus";

   if (TString (inFileName).Contains ("JZ1WA"))
    id = id + "_JZ1WA";
   else if (TString (inFileName).Contains ("JZRW1B"))
    id = id + "_JZRW1B";
  }
  else if (TString (inFileName).Contains ("jetjet")) { // dijet
   if (dataSet <= 0) return "";
   id = id + "Dijet_Slice" + to_string (dataSet);
  }

  return id;
}


/**
 * Separates each point on a TGAE by delta along the x axis, so that the errors don't overlap.
 */
void deltaize (TGAE* tg, const double delta, const bool logx) {
  double x, y, exh, exl;
  for (int n = 0; n < tg->GetN (); n++) {
    tg->GetPoint (n, x, y);
    exh = x + tg->GetErrorXhigh (n);
    exl = x - tg->GetErrorXhigh (n);
    
    if (logx) tg->SetPoint (n, x*delta, y);
    else tg->SetPoint (n, x + delta, y);

    tg->GetPoint (n, x, y);
    exh = exh - x;
    exl = x - exl;

    tg->SetPointEXhigh (n, exh);
    tg->SetPointEXlow (n, exl);
  }
}


/**
 * Makes a TGAE from the input histogram.
 */
TGAE* make_graph (TH1* h, const float cutoff) {
  TGAE* tg = new TGAE ();

  const float xlo = h->GetBinLowEdge (1);
  const float xhi = h->GetBinLowEdge (h->GetNbinsX ()) + h->GetBinWidth (h->GetNbinsX ());

  for (int n = 0; n < h->GetNbinsX (); n++) {
    if (cutoff >= 0 && h->GetBinContent (n+1) <= cutoff) {
      continue;
    }
    else {
      tg->SetPoint (tg->GetN (), h->GetBinCenter (n+1), h->GetBinContent (n+1));
      tg->SetPointError (tg->GetN ()-1, h->GetBinWidth (n+1) / 2, h->GetBinWidth (n+1) / 2, h->GetBinError (n+1), h->GetBinError (n+1));
    }
  }

  tg->GetXaxis ()->SetRangeUser (xlo, xhi);
  return tg;
}


/**
 * Converts a TEfficiency to a TGAE
 */
TGAE* TEff2TGAE (TEfficiency* e) {
  const int nx = e->GetTotalHistogram ()->GetNbinsX ();
  TGAE* g = new TGAE (nx);
  for (int ix = 0; ix < nx; ix++) {
    g->SetPoint (ix, e->GetTotalHistogram ()->GetBinCenter (ix+1), e->GetEfficiency (ix+1));
    g->SetPointEXhigh (ix, 0.5*e->GetTotalHistogram ()->GetBinWidth (ix+1));
    g->SetPointEXlow (ix, 0.5*e->GetTotalHistogram ()->GetBinWidth (ix+1));
    g->SetPointEYhigh (ix, e->GetEfficiencyErrorUp (ix+1));
    g->SetPointEYlow (ix, e->GetEfficiencyErrorLow (ix+1));
  }
  return g;
}


/**
 * Recenters a TGAE point for a log scale.
 */
void RecenterGraph (TGAE* g) {
  double x, y, newx, xlo, xhi, logDelta;
  for (int n = 0; n < g->GetN (); n++) {
    g->GetPoint (n, x, y);
    xlo = x - g->GetErrorXlow (n);
    xhi = x + g->GetErrorXhigh (n);
    logDelta = log10 (xhi) - log10 (xlo);
    newx = pow (10, log10 (xlo) + 0.5*logDelta);
    g->SetPoint (n, newx, y);
    g->SetPointEXlow (n, newx-xlo);
    g->SetPointEXhigh (n, xhi-newx);
  }
}


/**
 * Returns the TProfile of an input histogram along the x axis. Can use either statistical mean or gaussian mean.
 */
TH1D* GetProfileX (const TString name, TH2D* hist, const int nbinsx, const double* xbins, const bool useFit, const double* xlos, const double* xhis) {

  TH1D* prof = new TH1D (name, "", nbinsx, xbins);

  const bool useCustomBounds = (xlos && xhis);

  for (int xbin = 1; xbin <= nbinsx; xbin++) {
    TH1D* proj = hist->ProjectionY ("proj", xbin, xbin);
    if (proj->Integral () != 0)
      proj->Scale (1.0 / proj->Integral ());

    double mean, mean_err;

    // Calculate gaussian mean
    if (useFit) {
      TF1* gaus;
      if (useCustomBounds)
        gaus = new TF1 ("gaus", "gaus(0)", xlos[xbin-1], xhis[xbin-1]);
      else
        gaus = new TF1 ("gaus", "gaus(0)", proj->GetMean ()-2.8*proj->GetStdDev (), proj->GetMean ()+2.8*proj->GetStdDev ());
      gaus->SetParameter (1, proj->GetMean ());
      gaus->SetParameter (2, proj->GetStdDev ());
      proj->Fit (gaus, "Q0R");
      mean = gaus->GetParameter (1);
      mean_err = gaus->GetParError (1);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // Calculate statistical mean
    else {
      mean = proj->GetMean ();
      mean_err = proj->GetMeanError ();
    }

    prof->SetBinContent (xbin, mean);
    prof->SetBinError (xbin, mean_err);
    if (proj) { delete proj; proj = NULL; }
  }

  return prof;
}


/**
 * Returns the TProfile of an input histogram along the y axis. Can use either statistical mean or gaussian mean.
 */
TH1D* GetProfileY (const TString name, TH2D* hist, const int nbinsy, const double* ybins, const bool useFit, const double* ylos, const double* yhis) {

  TH1D* prof = new TH1D (name, "", nbinsy, ybins);

  const bool useCustomBounds = (ylos && yhis);

  for (int ybin = 1; ybin <= nbinsy; ybin++) {
    TH1D* proj = hist->ProjectionX ("proj", ybin, ybin);
    if (proj->Integral () != 0)
      proj->Scale (1.0 / proj->Integral ());

    double mean, mean_err;

    // Calculate gaussian mean
    if (useFit) {
      TF1* gaus;
      if (useCustomBounds)
        gaus = new TF1 ("gaus", "gaus(0)", ylos[ybin-1], yhis[ybin-1]);
      else
        gaus = new TF1 ("gaus", "gaus(0)", proj->GetMean ()-2.8*proj->GetStdDev (), proj->GetMean ()+2.8*proj->GetStdDev ());
      gaus->SetParameter (1, proj->GetMean ());
      gaus->SetParameter (2, proj->GetStdDev ());
      proj->Fit (gaus, "Q0R");
      mean = gaus->GetParameter (1);
      mean_err = gaus->GetParError (1);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // Calculate statistical mean
    else {
      mean = proj->GetMean ();
      mean_err = proj->GetMeanError ();
    }

    prof->SetBinContent (ybin, mean);
    prof->SetBinError (ybin, mean_err);
    if (proj) { delete proj; proj = NULL; }
  }

  return prof;
}


/**
 * Returns a histogram with the TProfile of data over the TProfile of MC along either the x or y axes. Can use either the statistical or gaussian mean.
 */
TH1D* GetDataOverMC (const TString name, TH2D* data, TH2D* mc, const int numbins, const double* bins, const bool useFit, const TString axis, const double* los, const double* his) {

  // figure out which axis to use
  if (axis != "x" && axis != "y" && axis != "X" && axis != "Y") {
    cout << "JetCalibration::GetDataOverMC: Invalid axis specified!" << endl;
    return NULL;
  }
  const bool useXaxis = (axis == "x" || axis == "X");
  const bool useCustomBounds = (los && his);

  TH1D* dataOverMC = new TH1D (name, "", numbins, bins);

  TH1D* proj = NULL;

  double dataAvg, dataErr, mcAvg, mcErr;

  for (int bin = 1; bin <= numbins; bin++) {

    // first calculate the data value (numerator)
    if (useXaxis) proj = data->ProjectionY (name + Form ("data_bin%i", bin), bin, bin);
    else proj = data->ProjectionX (name + Form ("data_bin%i", bin), bin, bin);
    if (proj->Integral () != 0)
      proj->Scale (1.0 / proj->Integral ());

    // Calculate gaussian mean
    if (useFit) {
      TF1* gaus;
      if (useCustomBounds)
        gaus = new TF1 ("gaus", "gaus(0)", los[bin-1], his[bin-1]);
      else
        gaus = new TF1 ("gaus", "gaus(0)", proj->GetMean ()-2.8*proj->GetStdDev (), proj->GetMean ()+2.8*proj->GetStdDev ());
      gaus->SetParameter (1, proj->GetMean ());
      gaus->SetParameter (2, proj->GetStdDev ());
      proj->Fit (gaus, "Q0R");
      dataAvg = gaus->GetParameter (1);
      dataErr = gaus->GetParError (1);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // Calculate statistical mean
    if (!useFit) {
      dataAvg = proj->GetMean ();
      dataErr = proj->GetMeanError ();
    }

    if (proj) { delete proj; proj = NULL; }

    // next calculate the MC value (denominator)
    if (useXaxis) proj = mc->ProjectionY (name + Form ("mc_bin%i", bin), bin, bin);
    else proj = mc->ProjectionX (name + Form ("mc_bin%i", bin), bin, bin);
    if (proj->Integral () != 0)
      proj->Scale (1.0 / proj->Integral ());

    // Calculate gaussian mean
    if (useFit) {
      TF1* gaus;
      if (useCustomBounds)
        gaus = new TF1 ("gaus", "gaus(0)", los[bin-1], his[bin-1]);
      else
        gaus = new TF1 ("gaus", "gaus(0)", proj->GetMean ()-2.8*proj->GetStdDev (), proj->GetMean ()+2.8*proj->GetStdDev ());
      gaus->SetParameter (1, proj->GetMean ());
      gaus->SetParameter (2, proj->GetStdDev ());
      proj->Fit (gaus, "Q0R");
      mcAvg = gaus->GetParameter (1);
      mcErr = gaus->GetParError (1);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // Calculate statistical mean
    else {
      mcAvg = proj->GetMean ();
      mcErr = proj->GetMeanError ();
    }

    if (proj) { delete proj; proj = NULL; }

    // now set the nominal value to the numerator over denominator
    if (!(dataAvg == 0 || isnan (dataAvg) || mcAvg == 0 || isnan (mcAvg))) {
      const double dataOverMCavg = dataAvg/mcAvg;
      const double dataOverMCerr = sqrt (pow (dataErr/mcAvg, 2) + pow (dataOverMCavg * mcErr/mcAvg, 2));
      dataOverMC->SetBinContent (bin, dataOverMCavg);
      dataOverMC->SetBinError (bin, dataOverMCerr);
    }
  }

  if (useXaxis) dataOverMC->GetXaxis ()->SetTitle (data->GetXaxis ()->GetTitle ()); // copy x axis title
  else dataOverMC->GetXaxis ()->SetTitle (data->GetYaxis ()->GetTitle ()); // copy y axis title
  return dataOverMC;
}


/**
 * Converts a TProfile to a TH1D.
 */
TH1D* TProfile2TH1D (const char* name, TProfile* p, const int nx, const double* x) {
  TH1D* h = new TH1D (name, "", nx, x);
  h->Sumw2 ();

  for (int ix = 1; ix <= nx; ix++) {
    h->SetBinContent (ix, p->GetBinContent (ix));
    h->SetBinError (ix, p->GetBinError (ix));
  }
  return h;
}


/**
 * Reflects the contents of h around the n-th bin in x.
 */
void GetReflectionX (TH1D* h, const int n) {
  if (n < 1 || h->GetNbinsX () < n) {
    cout << "Invalid choice of n! Exiting gracefully." << endl;
    return;
  }

  const int d = std::max (n, h->GetNbinsX () - n);
  int start = n - d + (h->GetNbinsX () % 2) + 1;
  int end = n + d;
  while (start < end) {
    if (1 <= start && end <= h->GetNbinsX ()) {
      const double temp = h->GetBinContent (start); 
      const double temperr = h->GetBinError (start);

      h->SetBinContent (start, h->GetBinContent (end)); 
      h->SetBinError (start, h->GetBinError (end));

      h->SetBinContent (end, temp);
      h->SetBinError (end, temperr);
    }
    else {
      if (start < 1) {
        h->SetBinContent (end, 0);
        h->SetBinError (end, 0);
      }
      if (end > h->GetNbinsX ()) {
        h->SetBinContent (start, 0);
        h->SetBinError (start, 0);
      }
    }
    start++; 
    end--; 
  } 

  return;
}


/**
 * Applies new binning to a histogram
 * BE CAREFUL: if bins edges don't overlap, this can lead to unexpected behavior!
 */ 
void RebinSomeBins (TH1D* &h, int nbins, double* bins) {
  if (nbins >= h->GetNbinsX ()) {
    cout << "More new bins than old bins, returning." << endl;
    return;
  }

  const TString name = h->GetName ();
  h->SetName ("temp");

  int noldbins = h->GetNbinsX ();
  double* oldbins = new double[noldbins+1];
  for (int ix = 1; ix <= noldbins; ix++) {
    oldbins[ix-1] = h->GetBinLowEdge (ix);
  }
  oldbins[noldbins] = h->GetBinLowEdge (noldbins) + h->GetBinWidth (noldbins);

  TH1D* hnew = new TH1D (name, "", nbins, bins);
  hnew->Sumw2 ();
  for (int ix = 0; ix < nbins; ix++) {
    for (int ixprime = 0; ixprime < noldbins; ixprime++) {
      if (bins[ix] <= oldbins[ixprime] && oldbins[ixprime+1] <= bins[ix+1]) {
        hnew->SetBinContent (ix+1, hnew->GetBinContent (ix+1) + h->GetBinContent (ixprime+1));
        hnew->SetBinError (ix+1, hnew->GetBinError (ix+1) + pow (h->GetBinError (ixprime+1), 2));
      }
      else if (oldbins[ixprime] <= bins[ix] && bins[ix+1] <= oldbins[ixprime+1]) {
        hnew->SetBinContent (ix+1, hnew->GetBinContent (ix+1) + h->GetBinContent (ixprime+1));
        hnew->SetBinError (ix+1, hnew->GetBinError (ix+1) + pow (h->GetBinError (ixprime+1), 2));
      }
    }
    hnew->SetBinError (ix+1, sqrt (hnew->GetBinError (ix+1)));
  }
  delete[] oldbins;

  delete h;
  h = hnew;
}


/**
 * Adds a to h without propagating errors (e.g. for subtracting a background)
 */
void AddNoErrors (TH1D* h, TH1D* a, const float sf) {
  assert (h->GetNbinsX () == a->GetNbinsX ());
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    h->SetBinContent (ix, h->GetBinContent (ix) + sf * a->GetBinContent (ix));
  }
}


/**
 * Multiplies h by a raised to power without propagating errors (e.g. for subtracting a background)
 */
void MultiplyNoErrors (TH1D* h, TH1D* a, const float power) {
  for (int ix = 1; ix <= h->GetNbinsX (); ix++) {
    h->SetBinContent (ix, h->GetBinContent (ix) * pow (a->GetBinContent (ix), power));
  }
}


} // end namespace
