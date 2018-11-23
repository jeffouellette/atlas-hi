#include "CalibUtils.h"
#include "Params.h"

#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

namespace JetCalibration {

TFile* xCalibSystematicsFile = NULL;
TFile* dataOverMCFile = NULL;
TFile* purityFile = NULL;


/**
 * Returns a new TH3D equal to tight - lnt.
 */
void SubtractLooseNonTight (TH3D* signal, const TH3D* tight, const TH3D* lnt, const TH3D* tightCounts, const TH3D* lntCounts) {
  const int nx = tight->GetNbinsX ();
  const int ny = tight->GetNbinsY ();
  const int nz = tight->GetNbinsZ ();
  for (int ix = 1; ix <= nx; ix++) {
   for (int iy = 1; iy <= ny; iy++) {
    double ntight, ntighterr, nlnt, nlnterr;
    ntight = tightCounts->IntegralAndError (ix, ix, iy, iy, 1, tightCounts->GetNbinsZ (), ntighterr);
    nlnt = lntCounts->IntegralAndError (ix, ix, iy, iy, 1, lntCounts->GetNbinsZ (), nlnterr);

    //const double ntight = tightCounts->GetBinContent (ix, iy);
    //const double ntighterr = tightCounts->GetBinError (ix, iy);
    //const double nlnt = lntCounts->GetBinContent (ix, iy); 
    //const double nlnterr = lntCounts->GetBinError (ix, iy);

    double w = 0, werr = 0;
    if (nlnt != 0) {
     w = ntight / nlnt;
     werr += pow (w * nlnterr / nlnt, 2);
    } 
    else
     w = 0;

    if (ntight != 0)
     werr += pow (w * ntighterr / ntight, 2);

    werr = sqrt (werr);

    for (int iz = 1; iz <= nz; iz++) {
     const double t = tight->GetBinContent (ix, iy, iz);
     const double terr = tight->GetBinError (ix, iy, iz);
     const double l = lnt->GetBinContent (ix, iy, iz);
     const double lerr = lnt->GetBinError (ix, iy, iz);

     signal->SetBinContent (ix, iy, iz, t - w * l);
     signal->SetBinError (ix, iy, iz, sqrt (pow (terr, 2) + pow (werr * l, 2) + pow (w * lerr, 2)));
    }
   }
  }
  return;
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
 * Returns the photon purities as a function of photon pT and eta.
 * Requires purityFile to be defined and open, otherwise will return 1 (perfect purity).
 */
double GetPurity (const double ppt, const double peta) {
  TFile* file = purityFile;

  if (!file || !file->IsOpen ()) {
   cout << "Warning: In Utils.C: Photon purity file not open! Will return 1." << endl;
   return 1;
  }

  if ((1.37 < TMath::Abs (peta) && TMath::Abs (peta) < 1.52) || 2.37 < TMath::Abs (peta)) {
   return 1;
  }

  TH2D* purities = (TH2D*)purityFile->Get ("purities");
  return purities->GetBinContent (purities->GetXaxis ()->FindBin (ppt), purities->GetYaxis ()->FindBin (peta));
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

} // end namespace
