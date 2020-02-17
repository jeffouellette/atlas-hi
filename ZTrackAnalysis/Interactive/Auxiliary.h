#ifndef __Auxiliary_h__
#define __Auxiliary_h__

#include "PhysicsAnalysis.h"
#include "Systematic.h"
#include "Params.h"

#include <AtlasUtils.h>

#include "math.h"

using namespace std;
using namespace atlashi;

typedef TGraphAsymmErrors TGAE;

void PlotPullDist (const PhysicsAnalysis* a, Systematic* sys = nullptr, const bool useTrkPt = true) {
  TCanvas* c = new TCanvas ("c_pull", "", 800, 600);
  TH1D* h_pull = new TH1D ("h_pull", Form (";Pull = (Y^{#it{ee}} (%s) - Y^{#it{#mu#mu}} (%s)) / #sigma;Entries", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}"), 20, -5, 5);

  int nPoints = 0;
  for (short iCent = 0; iCent < numCentBins; iCent++) {
  //for (short iCent = 0; iCent < 1; iCent++) {
    for (short iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
      TH1D* h_ee = (useTrkPt ? a->h_trk_pt_ptz_sub[0][iPtZ][iCent] : a->h_trk_xhz_ptz_sub[0][iPtZ][iCent]);
      TH1D* h_mumu = (useTrkPt ? a->h_trk_pt_ptz_sub[1][iPtZ][iCent] : a->h_trk_xhz_ptz_sub[1][iPtZ][iCent]);
      TGAE* g_ee_sys = (!sys ? nullptr : sys->GetTGAE (useTrkPt ? sys->h_trk_pt_ptz_sub[0][iPtZ][iCent] : sys->h_trk_xhz_ptz_sub[0][iPtZ][iCent]));
      TGAE* g_mumu_sys = (!sys ? nullptr : sys->GetTGAE (useTrkPt ? sys->h_trk_pt_ptz_sub[1][iPtZ][iCent] : sys->h_trk_xhz_ptz_sub[1][iPtZ][iCent]));
      const int nBinsX = h_ee->GetNbinsX ();
      for (short iX = 1; iX <= nBinsX; iX++) {
        float variance = pow (h_ee->GetBinError (iX), 2) + pow (h_mumu->GetBinError (iX), 2);
        if (sys)
          variance += pow (g_ee_sys->GetErrorY (iX-1), 2) + pow (g_mumu_sys->GetErrorY (iX-1), 2);

        const float pull = (h_ee->GetBinContent (iX) - h_mumu->GetBinContent (iX)) / sqrt (variance);
        h_pull->Fill (pull);
        nPoints++;
      } // end loop over iX
    } // end loop over iPtZ
  } // end loop over iCent

  TH1D* h_pull_expected = new TH1D ("h_pull_expected", Form (";Pull = (Y^{#it{ee}} (%s) - Y^{#it{#mu#mu}} (%s)) / #sigma;Entries", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}"), 20, -5, 5);
  TF1* gaus = new TF1 ("gaus", "gaus(0)", -5, 5);
  gaus->SetParameter (0, 1);
  gaus->FixParameter (1, 0);
  gaus->FixParameter (2, 1);
  for (short i = 0; i < 2000; i++) {
    h_pull_expected->Fill (gaus->GetRandom ());
  }
  SaferDelete (gaus);
  h_pull_expected->Scale (nPoints/2000.);
  h_pull_expected->SetLineColorAlpha (kRed, 0.35);
  h_pull_expected->SetLineStyle (2);
  h_pull_expected->SetLineWidth (3);
 
  h_pull->GetYaxis ()->SetRangeUser (0, 10);
  h_pull->Draw ("hist");
  h_pull_expected->Draw ("hist same");

  TF1* fit = new TF1 ("fit", "gaus(0)", -5, 5);
  fit->SetParameter (0, 1);
  fit->SetParameter (1, 0);
  fit->SetParameter (2, 0);
  h_pull->Fit (fit, "RN0Q");

  const float mu = fit->GetParameter (1);
  const float muErr = fit->GetParError (1);
  const float sigma = fit->GetParameter (2);
  const float sigmaErr = fit->GetParError (2);
  myText (0.20, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.055);
  myText (0.20, 0.80, kBlack, "Observed", 0.045);
  myText (0.25, 0.75, kBlack, Form ("#mu = %.2f #pm %.2f", mu, muErr), 0.045);
  myText (0.25, 0.70, kBlack, Form ("#sigma = %.2f #pm %.2f", sigma, sigmaErr), 0.045);
  myText (0.25, 0.65, kBlack, Form ("#chi^{2}/dof = %.2f / %i", fit->GetChisquare (), fit->GetNDF ()), 0.045);
  myText (0.20, 0.60, kRed, "Expected", 0.045);
  myText (0.25, 0.55, kRed, "#mu = 0", 0.045);
  myText (0.25, 0.50, kRed, "#sigma = 1", 0.045);
}

#endif
