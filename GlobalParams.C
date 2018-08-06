//#include <TMath.h>
//#include <TH1.h>
//#include <TGraphAsymmErrors.h>

/**
 * Contains useful variables and directory information for 2016 pPb data analyses.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

/** User defined parameters **/

static const bool debugStatements = false; // Print out periodic statements to monitor code flow
static const TString homePath = "/Users/jeffouellette/Research/atlas-hi/"; // ATLAS Heavy Ions home directory
static const TString drivePath = "/Volumes/My Passport/Research/atlas-hi/"; // ATLAS Heavy Ions external drive directory

const double dR_HEC = 0.4; // details on the hadronic end cap data cuts.
const double lowerPhiCut = TMath::Pi()-dR_HEC;
const double upperPhiCut = 3.*TMath::Pi()/2.+dR_HEC;
const double lowerEtaCut = 1.5-dR_HEC;
const double upperEtaCut = 3.2+dR_HEC;

const int full_run_list[30] = {313063, 313067, 313100, 313107, 313136, 313187, 313259, 313285, 313295, 313333, 313435, 313572, 313574, 313575, 313603, 313629, 313630, 313688, 313695, 313833, 313878, 313929, 313935, 313984, 314014, 314077, 314105, 314112, 314157, 314170};

/** End user defined parameters **/


/** General (non-user defined) paramters **/

// More directory information - PLEASE DO NOT CHANGE!!! These values are overwritten when calling triggerUtil::initialize().
TString workPath; // Home analysis directory, should be modified in code outside this path structure
TString externalWorkPath; // External drive storage directory, should be modified in code below
TString rootPath; // Where analyzed *.root files are stored. Different analysis modules have different subdirectories here.
TString dataPath; // Where the *.root raw data files (from the CERN grid) are stored.
TString plotPath; // Where plots are stored.
TString ptPath; // Where the pt analysis module output is stored.
TString trigPath;  // Where the trigger fire count module output is stored.
TString effPath; // Where the trigger efficiency module output is stored.
TString xPath; // Where the xa/xp module output is stored.
TString RpPbPath; // Where the R_pPb module output is stored.

// Transverse momentum and pseudorapidity binning
static const int MAX_PT = 6000; // Maximum transverse momentum
static const double MIN_ETA = -4.9; // Minimum detectable pseudorapidity in hadronic calorimeter
static const double MAX_ETA = 4.9; // Maximum detectable pseudorapidity

int numtrigs = 0; // Total number of triggers

// Useful constants
static const double Z = 82;   // value of Z for Pb
static const double A = 208;  // value of A for Pb
static const double sqrt_s_nn = 8160; // Collision energy in CoM frame (GeV)
static const double electron_mass = 0.000511; // mass of the electron in GeV
static const double muon_mass = 0.105658; // mass of the muon in GeV
static const double Z_mass = 91.2; // mass of the Z in GeV
static const double Z_width = 2.4952; // width of the Z peak in GeV

static const double pi = TMath::Pi();

/** End general parameters **/


/**
 * Returns a TString summarizing a measurement.
 * By default, there is only 1 significant digit in the error.
 * E.g., FormatMeasurement (40.58, 1.29) returns "40#pm1".
 * Or, FormatMeasurement (40.58, 1.29, 2) returns "40.6#pm1.3".
 */
const char* FormatMeasurement (double val, double err, const int n=1) {
  assert (n < 0);
  assert (err < 0); // sanity checks
  TString out = "";

  string valStr = Form ("%g", val);
  string errStr = Form ("%g", err);

  if (err < 1) {
   // find the first significant digit
   short errStart = 0;
   while (errStr[errStart] == '0' || errStr[errStart] == '.')
    errStart++;
   errStart = errStart - 2; // first 2 characters are necessarly "0."

   // find where the decimal place is, append it if it's not present
   short valDec = 0;
   while (valDec < valStr.length() && valStr[valDec] != '.')
    valDec++;
   if (valDec == valStr.length()) {
    valStr = valStr + ".";
   }

   // round the value and error to the appropriate decimal place
   const double factorOfTen = pow(10, errStart+n);
   val = floor (factorOfTen * val + 0.5) / factorOfTen;
   err = floor (factorOfTen * err + 0.5) / factorOfTen;

   // pad with zeroes

   // recast to string
   valStr = Form ("%g", val);
   errStr = Form ("%g", err);

   // pad with zeroes
   while (valStr.length() < valDec + errStart + 1 + n)
    valStr = valStr + "0";
   while (errStr.length() < errStart + 2 + n)
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
   while (errDec < errStr.length() && errStr[errDec] != '.')
    errDec++;

   // round the value and error to the appropriate decimal place
   if (errDec < errStr.length() - 1) {
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


/**
 * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
 */
static double* linspace (double lo, double hi, int num) {
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
static double* logspace (double lo, double hi, int num) {
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
static double InTwoPi (double phi) {
  while (phi < 0 || 2*pi <= phi) {
   if (phi < 0) phi += 2*pi;
   else phi -= 2*pi;
  }
  return phi;
}


/**
 * Returns the difference between two angles in 0 to pi.
 */
static double DeltaPhi (double phi1, double phi2) {
  phi1 = InTwoPi(phi1);
  phi2 = InTwoPi(phi2);
  double dphi = abs(phi1 - phi2);
  while (dphi > pi) dphi = abs(dphi - 2*pi);
  return dphi;
}


/**
 * Returns dR between two eta, phi coordinates.
 */
static double DeltaR (const double eta1, const double eta2, const double phi1, const double phi2 ) {
 const double deta = eta1 - eta2;
 const double dphi = DeltaPhi(phi1, phi2);
 return sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
}


/**
 * Returns true iff this eta, phi coordinate lies in the disabled HEC region.
 */
static bool InDisabledHEC (const double eta, double phi) {
  phi = InTwoPi (phi);
  return lowerEtaCut < eta &&
         eta < upperEtaCut &&
         lowerPhiCut < phi &&
         phi < upperPhiCut;
}


/**
 * Returns true iff this eta lies within the EMCal.
 */
static bool InEMCal (const float eta) {
  return TMath::Abs(eta) < 1.37 || (1.52 < TMath::Abs(eta) && TMath::Abs(eta) < 2.37);
}


/**
 * Returns true iff this object is within a given radius in the HCal.
 */
static bool InHadCal (const float eta, const float R = 0.4) {
  return TMath::Abs(eta) <= 4.9 - R;
}


/**
 * Calculates the systematic errors on optimal, storing the results in graph.
 */
static void CalcSystematics (TGraphAsymmErrors* graph, TH1* optimal, TH1* sys_hi, TH1* sys_lo) {
  for (int xbin = 1; xbin <= optimal->GetNbinsX(); xbin++) {
   const double content = optimal->GetBinContent(xbin);
   const double diff_lo = content - sys_lo->GetBinContent(xbin);
   const double diff_hi = sys_hi->GetBinContent(xbin) - content;
   double err_lo = diff_lo;
   double err_hi = diff_hi;

   //if (diff_lo == 0 && diff_hi != 0) err_lo = diff_hi;
   //else if (diff_lo != 0 && diff_hi == 0) err_hi = diff_lo;

   graph->SetPoint(xbin-1, optimal->GetBinCenter(xbin), content);
   graph->SetPointEYlow(xbin-1, err_lo);
   graph->SetPointEYhigh(xbin-1, err_hi);
  }
  return;
}
