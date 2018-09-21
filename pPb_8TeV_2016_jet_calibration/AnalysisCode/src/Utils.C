#include "Utils.h"
#include "Params.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TF1.h>

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {

TH1D* GetProfileX (const TString name, TH2D* hist, const int nbinsx, const double* xbins, const bool useFit) {

  TH1D* prof = new TH1D (name, "", nbinsx, xbins);

  for (int xbin = 1; xbin <= nbinsx; xbin++) {
    TH1D* projy = hist->ProjectionY ("projy", xbin, xbin);
    projy->Rebin (rebinFactor);

    //projy->GetXaxis ()->SetLimits (0, 2.0);
    double mean, mean_err;
    double chi_square = 0;
    int numNonzeroBins = 0;
    for (int xbinprime = 1; xbinprime <= projy->GetNbinsX (); xbinprime++)
      if (projy->GetBinContent (xbinprime) > 0) numNonzeroBins++;

    if (useFit && useGaussian && numNonzeroBins > 4) {
      TF1* gaus = new TF1 ("gaus", "gaus (0)", projy->GetXaxis ()->GetBinLowEdge (1), projy->GetXaxis ()->GetBinUpEdge (projy->GetNbinsX ()));
      projy->Fit (gaus, "Q0R");
      mean = gaus->GetParameter (1);
      mean_err = gaus->GetParError (1);
      chi_square = gaus->GetChisquare () / (projy->GetNbinsX () - 3);
      if (gaus) { delete gaus; gaus = NULL; }
    }

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


TH1D* GetProfileY (const TString name, TH2D* hist, const int nbinsy, const double* ybins, const bool useFit) {

  TH1D* prof = new TH1D (name, "", nbinsy, ybins);

  for (int ybin = 1; ybin <= nbinsy; ybin++) {
    TH1D* projx = hist->ProjectionX ("projx", ybin, ybin);
    projx->Rebin (rebinFactor);
    //projx->GetXaxis ()->SetLimits (0, 2.0);

    double mean, mean_err;
    double chi_square = 0;
    int numNonzeroBins = 0;
    for (int ybin = 1; ybin <= projx->GetNbinsX (); ybin++)
      if (projx->GetBinContent (ybin) > 0) numNonzeroBins++;

    if (useFit && useGaussian && numNonzeroBins > 4) {
      TF1* gaus = new TF1 ("gaus", "gaus (0)", projx->GetXaxis ()->GetBinLowEdge (1), projx->GetXaxis ()->GetBinUpEdge (projx->GetNbinsX ()));
      projx->Fit (gaus, "Q0R");
      mean = gaus->GetParameter (1);
      mean_err = gaus->GetParError (1);
      chi_square = gaus->GetChisquare () / (projx->GetNbinsX () - 3);
      if (gaus) { delete gaus; gaus = NULL; }
    }

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
    if (!useXaxis) proj = data->ProjectionY (name + Form ("data_bin%i", bin), bin, bin);
    else proj = data->ProjectionX (name + Form ("data_bin%i", bin), bin, bin);
    proj->Rebin (rebinFactor);

    chi_square = 0;
    numNonzeroBins = 0;
    for (int binprime = 1; binprime <= proj->GetNbinsX (); binprime++)
      if (proj->GetBinContent (bin) > 0) numNonzeroBins++;

    if (useFit && useGaussian && numNonzeroBins > 4) {
      TF1* gaus = new TF1 ("gaus", "gaus (0)", proj->GetXaxis ()->GetBinLowEdge (1), 2.0);//proj->GetXaxis ()->GetBinUpEdge (proj->GetNbinsX ()));
      proj->Fit (gaus, "Q0R");
      dataAvg = gaus->GetParameter (1);
      dataErr = gaus->GetParError (1);
      chi_square = gaus->GetChisquare () / (proj->GetNbinsX () - 3);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    if (!useGaussian || !useFit || chi_square > 1.0 || numNonzeroBins <= 4) {
      dataAvg = proj->GetMean ();
      dataErr = proj->GetMeanError ();
    }

    if (proj) { delete proj; proj = NULL; }

    // next calculate the MC value (denominator)
    if (!useXaxis) proj = mc->ProjectionY (name + Form ("mc_bin%i", bin), bin, bin);
    else proj = mc->ProjectionY (name + Form ("mc_bin%i", bin), bin, bin);
    proj->Rebin (rebinFactor);

    chi_square = 0;
    numNonzeroBins = 0;

    for (int binprime = 1; binprime <= proj->GetNbinsX (); binprime++)
      if (proj->GetBinContent (binprime) != 0) numNonzeroBins++;

    if (useFit && useGaussian && numNonzeroBins > 4) {
      TF1* gaus = new TF1 ("gaus", "gaus (0)", proj->GetXaxis ()->GetBinLowEdge (1), 2.0);//proj->GetXaxis ()->GetBinUpEdge (proj->GetNbinsX ()));
      proj->Fit (gaus, "Q0R");
      mcAvg = gaus->GetParameter (1);
      mcErr = gaus->GetParError (1);
      chi_square = gaus->GetChisquare () / (proj->GetNbinsX () - 3);
      if (gaus) { delete gaus; gaus = NULL; }
    }

    if (!useGaussian || !useFit || chi_square > 1.0 || numNonzeroBins <= 4) {
      mcAvg = proj->GetMean ();
      mcErr = proj->GetMeanError ();
    }

    if (proj) { delete proj; proj = NULL; }

    // now set the nominal value to the numerator over denominator
    const double dataOverMCavg = dataAvg/mcAvg;
    const double dataOverMCerr = dataOverMCavg * TMath::Sqrt (TMath::Power (dataErr/dataAvg, 2) + TMath::Power (mcErr/mcAvg, 2));
    if (!isnan (dataOverMCavg) && !isnan (dataOverMCerr)) {
      dataOverMC->SetBinContent (bin, dataOverMCavg);
      dataOverMC->SetBinError (bin, dataOverMCerr);
    }
  }

  dataOverMC->GetXaxis ()->SetTitle (data->GetXaxis ()->GetTitle ());
  return dataOverMC;
}

} // end namespace
