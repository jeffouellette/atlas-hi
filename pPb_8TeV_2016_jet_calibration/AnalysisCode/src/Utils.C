#include "Utils.h"
#include "Params.h"

#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

namespace pPb8TeV2016JetCalibration {

TFile* xCalibSystematicsFile = NULL;
TFile* dataOverMCFile = NULL;

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
    cout << "Error: In Utils.C: Cannot identify this MC file! Will return null!" << endl;
    return NULL;
   }
  }
  else fileIdentifier = inFileName;

  // Now get the list of files
  const TString dataPathTemp = dataPath + "/" + directory + "/";
  TSystemDirectory dir (dataPathTemp.Data (), dataPathTemp.Data ());
  TList* sysfiles = dir.GetListOfFiles ();
  if (!sysfiles) {
   cout << "Error: In Utils.C: Cannot get list of files! Will return null!" << endl;
   return NULL;
  }
  TSystemFile* sysfile;
  TString fname;
  TIter next (sysfiles);

  while ( (sysfile = (TSystemFile*)next ())) {
   fname = sysfile->GetName ();
   if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
    if (debugStatements) cout << "Status: In Utils.C: Found " << fname.Data () << endl;
    
    if (fname.Contains (fileIdentifier)) {
     file = new TFile (dataPathTemp+fname, "READ");
     break;
    }
   }
  }

  if (!file) {
   cout << "Error: In Utils.C: TFile not obtained for given data set. Will return null!" << endl;
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

  if (TString (inFileName).Contains ("jetjet")) { // dijet
   if (dataSet <= 0) return "";
   id = id + "Dijet_Slice" + to_string (dataSet);
  }
  else if (TString (inFileName).Contains ("42310") && TString (inFileName).Contains ("Slice")) { // gamma+jet
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

  return id;
}


/**
 * Returns the initial systematic error on the 2015 cross-calibration as a function of jet pT and eta.
 * Requires xCalibSystematicsFile to be defined and open, else will return 0.
 */
double GetXCalibSystematicError (const double jpt, const double jeta) {
  TFile* file = xCalibSystematicsFile;

  if (!file || !file->IsOpen ()) {
   cout << "Warning: In Utils.C: Cross calibration systematics file not open! Will return 0." << endl;
   return 0;
  }

  if (TMath::Abs (jeta) < xcalibEtabins[0] ||
      xcalibEtabins[sizeof (xcalibEtabins)/sizeof (xcalibEtabins[0]) -1] < TMath::Abs (jeta)) {
   return 0;
  }

  short iEta = 0;
  while (xcalibEtabins[iEta] < TMath::Abs (jeta)) iEta++;
  iEta--;

  const TString hname = TString ("fsys_rel_") + Form ("%i", iEta);
  TH1D* fsys_rel = (TH1D*)file->Get (hname.Data ());

  return TMath::Abs (fsys_rel->GetBinContent (fsys_rel->FindBin (jpt)) - 1) * jpt;
}


/**
 * Returns the additional systematics associated with applying the 2015 cross-calibration to the 2016 pPb, as a function of jet pT and eta.
 * Requires dataOverMCFile to be defined and open, else will return 0.
 */
double GetNewXCalibSystematicError (const double jeta, const double refpt, const bool periodA) {
  return 0;
  TFile* file = dataOverMCFile;
  if (!file || !file->IsOpen ())
   return 0;

  if (jeta < etabins[0] ||
      etabins[numetabins] < jeta) {
   return 0;
  }

  short bin = 0;
  while (bin <= numetabins && etabins[bin] < jeta) bin++;
  bin--;

  const char* period = (periodA ?  "periodA" : "periodB");
  const TString hname = TString (Form ("gJetPtRatio_diff%i_stat_%s", bin, period));
  TH1D* hist = (TH1D*)file->Get (hname.Data ());

  return hist->GetBinContent (hist->FindBin (refpt)) * refpt;
}


/**
 * Separates each point on a TGraphAsymmErrors by delta along the x axis, so that the errors don't overlap.
 */
void deltaize (TGraphAsymmErrors* tg, const double delta, const bool logx) {
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
 * Makes a TGraphAsymmErrors from the input histogram.
 */
TGraphAsymmErrors* make_graph (TH1* h, const float cutoff) {
  TGraphAsymmErrors* tg = new TGraphAsymmErrors (TString (h->GetName ()) + "_graph");

  for (int n = 0; n < h->GetNbinsX (); n++) {
    if (cutoff >= 0 && h->GetBinContent (n+1) <= cutoff) {
      continue;
    }
    else {
      tg->SetPoint (tg->GetN (), h->GetBinCenter (n+1), h->GetBinContent (n+1));
      tg->SetPointError (tg->GetN ()-1, h->GetBinWidth (n+1) / 2, h->GetBinWidth (n+1) / 2, h->GetBinError (n+1), h->GetBinError (n+1));
    }
  }
  return tg;
}


/**
 * Returns the TProfile of an input histogram along the x axis. Can use either statistical mean or gaussian mean.
 */
TH1D* GetProfileX (const TString name, TH2D* hist, const int nbinsx, const double* xbins, const bool useFit) {

  TH1D* prof = new TH1D (name, "", nbinsx, xbins);

  for (int xbin = 1; xbin <= nbinsx; xbin++) {
    TH1D* projy = hist->ProjectionY ("projy", xbin, xbin);
    //projy->Rebin (rebinFactor);

    //projy->GetXaxis ()->SetLimits (0, 2.0);
    double mean, mean_err;
    double chi_square = 0;
    int numNonzeroBins = 0;
    for (int xbinprime = 1; xbinprime <= projy->GetNbinsX (); xbinprime++)
      if (projy->GetBinContent (xbinprime) > 0) numNonzeroBins++;

    // Calculate gaussian mean
    if (useFit && useGaussian && numNonzeroBins > 4) {
      TF1* gaus = new TF1 ("gaus", "gaus (0)", projy->GetXaxis ()->GetBinLowEdge (1), projy->GetXaxis ()->GetBinUpEdge (projy->GetNbinsX ()));
      projy->Fit (gaus, "Q0R");
      mean = gaus->GetParameter (1);
      mean_err = gaus->GetParError (1);
      chi_square = gaus->GetChisquare () / (projy->GetNbinsX () - 3);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // If statistics are too poor for gaussian mean or it is not desired, use the statistical mean.
    if (!useGaussian || !useFit || chi_square > 1.0 || numNonzeroBins <= 4) {
      mean = projy->GetMean ();
      mean_err = projy->GetMeanError ();
    }

    prof->SetBinContent (xbin, mean);
    prof->SetBinError (xbin, mean_err);
    if (projy) { delete projy; projy = NULL; }
  }

  return prof;
}


/**
 * Returns the TProfile of an input histogram along the y axis. Can use either statistical mean or gaussian mean.
 */
TH1D* GetProfileY (const TString name, TH2D* hist, const int nbinsy, const double* ybins, const bool useFit) {

  TH1D* prof = new TH1D (name, "", nbinsy, ybins);

  for (int ybin = 1; ybin <= nbinsy; ybin++) {
    TH1D* projx = hist->ProjectionX ("projx", ybin, ybin);
    //projx->Rebin (rebinFactor);
    //projx->GetXaxis ()->SetLimits (0, 2.0);

    double mean, mean_err;
    double chi_square = 0;
    int numNonzeroBins = 0;
    for (int ybin = 1; ybin <= projx->GetNbinsX (); ybin++)
      if (projx->GetBinContent (ybin) > 0) numNonzeroBins++;

    // Calculate gaussian mean
    if (useFit && useGaussian && numNonzeroBins > 4) {
      TF1* gaus = new TF1 ("gaus", "gaus (0)", projx->GetXaxis ()->GetBinLowEdge (1), projx->GetXaxis ()->GetBinUpEdge (projx->GetNbinsX ()));
      projx->Fit (gaus, "Q0R");
      mean = gaus->GetParameter (1);
      mean_err = gaus->GetParError (1);
      chi_square = gaus->GetChisquare () / (projx->GetNbinsX () - 3);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // If statistics are too poor for gaussian mean or it is not desired, use the statistical mean.
    if (!useGaussian || !useFit || chi_square > 1.0 || numNonzeroBins <= 4) {
      mean = projx->GetMean ();
      mean_err = projx->GetMeanError ();
    }

    prof->SetBinContent (ybin, mean);
    prof->SetBinError (ybin, mean_err);
    if (projx) { delete projx; projx = NULL; }
  }

  return prof;
}


/**
 * Returns a histogram with the TProfile of data over the TProfile of MC along either the x or y axes. Can use either the statistical or gaussian mean.
 */
TH1D* GetDataOverMC (const TString name, TH2D* data, TH2D* mc, const int numbins, const double* bins, const bool useFit, const TString axis) {

  // figure out which axis to use
  if (axis != "x" && axis != "y" && axis != "X" && axis != "Y") {
    cout << "pPb8TeV2016JetCalibration::GetDataOverMC: Invalid axis specified!" << endl;
    return NULL;
  }
  const bool useXaxis = (axis == "x" || axis == "X");

  TH1D* dataOverMC = new TH1D (name, "", numbins, bins);

  TH1D* proj = NULL;

  double dataAvg, dataErr, mcAvg, mcErr, chi_square;
  int numNonzeroBins;

  for (int bin = 1; bin <= numbins; bin++) {
    // first calculate the data value (numerator)
    if (useXaxis) proj = data->ProjectionY (name + Form ("data_bin%i", bin), bin, bin);
    else proj = data->ProjectionX (name + Form ("data_bin%i", bin), bin, bin);
    proj->Rebin (rebinFactor);

    chi_square = 0;
    numNonzeroBins = 0;
    for (int binprime = 1; binprime <= proj->GetNbinsX (); binprime++)
      if (proj->GetBinContent (bin) > 0) numNonzeroBins++;

    // Calculate gaussian mean
    if (useFit && useGaussian && numNonzeroBins > 4) {
      TF1* gaus = new TF1 ("gaus", "gaus(0)", proj->GetXaxis ()->GetBinLowEdge (1), proj->GetXaxis ()->GetBinUpEdge (proj->GetNbinsX ()));
      proj->Fit (gaus, "Q0R");
      dataAvg = gaus->GetParameter (1);
      dataErr = gaus->GetParError (1);
      chi_square = gaus->GetChisquare () / (proj->GetNbinsX () - 3);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // If statistics are too poor for gaussian mean or it is not desired, use the statistical mean.
    if (!useGaussian || !useFit || chi_square > 1.0 || numNonzeroBins <= 4) {
      dataAvg = proj->GetMean ();
      dataErr = proj->GetMeanError ();
    }

    if (proj) { delete proj; proj = NULL; }

    // next calculate the MC value (denominator)
    if (useXaxis) proj = mc->ProjectionY (name + Form ("mc_bin%i", bin), bin, bin);
    else proj = mc->ProjectionY (name + Form ("mc_bin%i", bin), bin, bin);
    proj->Rebin (rebinFactor);

    chi_square = 0;
    numNonzeroBins = 0;

    for (int binprime = 1; binprime <= proj->GetNbinsX (); binprime++)
      if (proj->GetBinContent (binprime) != 0) numNonzeroBins++;

    // Calculate gaussian mean
    if (useFit && useGaussian && numNonzeroBins > 4) {
      TF1* gaus = new TF1 ("gaus", "gaus(0)", proj->GetXaxis ()->GetBinLowEdge (1), proj->GetXaxis ()->GetBinUpEdge (proj->GetNbinsX ()));
      proj->Fit (gaus, "Q0R");
      mcAvg = gaus->GetParameter (1);
      mcErr = gaus->GetParError (1);
      chi_square = gaus->GetChisquare () / (proj->GetNbinsX () - 3);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    // If statistics are too poor for gaussian mean or it is not desired, use the statistical mean.
    if (!useGaussian || !useFit || chi_square > 1.0 || numNonzeroBins <= 4) {
      mcAvg = proj->GetMean ();
      mcErr = proj->GetMeanError ();
    }

    if (proj) { delete proj; proj = NULL; }

    // now set the nominal value to the numerator over denominator
    if (!(dataAvg == 0 || isnan (dataAvg) || mcAvg == 0 || isnan (mcAvg))) {
      const double dataOverMCavg = dataAvg/mcAvg;
      const double dataOverMCerr = dataOverMCavg * TMath::Sqrt (TMath::Power (dataErr/dataAvg, 2) + TMath::Power (mcErr/mcAvg, 2));
      dataOverMC->SetBinContent (bin, dataOverMCavg);
      dataOverMC->SetBinError (bin, dataOverMCerr);
    }
  }

  if (useXaxis) dataOverMC->GetXaxis ()->SetTitle (data->GetXaxis ()->GetTitle ()); // copy x axis title
  else dataOverMC->GetXaxis ()->SetTitle (data->GetYaxis ()->GetTitle ()); // copy y axis title
  return dataOverMC;
}

} // end namespace
