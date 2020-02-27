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
  TCanvas* c = new TCanvas (Form ("c_pull_%s", useTrkPt ? "pTch" : "xhZ"), "", 800, 600);
  TH1D* h_pull = new TH1D (Form ("h_pull_%s", useTrkPt ? "pTch" : "xhZ"), Form (";Pull = (#LTY^{#it{ee}} (%s)#GT - #LTY^{#it{#mu#mu}} (%s)#GT) / #sigma;Entries", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}"), 20, -5, 5);

  int nPoints = 0;
  for (short iCent = 1; iCent < numCentBins; iCent++) {
  //for (short iCent = 0; iCent < 1; iCent++) {
    for (short iPtZ = 3; iPtZ < nPtZBins; iPtZ++) {
      TH1D* h_ee = (useTrkPt ? a->h_trk_pt_ptz_sub : a->h_trk_xhz_ptz_sub)[0][iPtZ][iCent];
      TH1D* h_mumu = (useTrkPt ? a->h_trk_pt_ptz_sub : a->h_trk_xhz_ptz_sub)[1][iPtZ][iCent];
      TGAE* g_ee_sys = (sys ? sys->GetTGAE ((useTrkPt ? sys->h_trk_pt_ptz_sub : sys->h_trk_xhz_ptz_sub)[0][iPtZ][iCent]) : nullptr);
      TGAE* g_mumu_sys = (sys ? sys->GetTGAE ((useTrkPt ? sys->h_trk_pt_ptz_sub : sys->h_trk_xhz_ptz_sub)[1][iPtZ][iCent]) : nullptr);
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

  TH1D* h_pull_expected = new TH1D (Form ("h_pull_expected_%s", useTrkPt ? "pTch" : "xhZ"), Form (";Pull = (#LTY^{#it{ee}} (%s)#GT - #LTY^{#it{#mu#mu}} (%s)#GT) / #sigma;Entries", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}"), 20, -5, 5);
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
 
  h_pull->GetYaxis ()->SetRangeUser (0, 11);
  h_pull->Draw ("hist");
  h_pull_expected->Draw ("hist same");

  TF1* fit = new TF1 ("fit", "gaus(0)", -5, 5);
  fit->SetParameter (0, 1);
  fit->SetParameter (1, 0);
  fit->SetParameter (2, 0);
  h_pull->Fit (fit, "RN0QL");

  const float mu = fit->GetParameter (1);
  const float muErr = fit->GetParError (1);
  const float sigma = fit->GetParameter (2);
  const float sigmaErr = fit->GetParError (2);
  myText (0.20, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.055);
  //myText (0.20, 0.80, kBlack, "#it{pp}, all points #it{p}_{T}^{Z} > 30 GeV", 0.045);
  myText (0.20, 0.80, kBlack, "Pb+Pb, all points #it{p}_{T}^{Z} > 30 GeV", 0.045);
  //myText (0.20, 0.74, kBlack, "Poisson uncertainties", 0.045);
  myText (0.20, 0.74, kBlack, "Generalized uncertainties", 0.045);

  myText (0.65, 0.85, kBlack, "Observed", 0.045);
  myText (0.70, 0.80, kBlack, Form ("#mu = %.2f #pm %.2f", mu, muErr), 0.045);
  myText (0.70, 0.75, kBlack, Form ("#sigma = %.2f #pm %.2f", sigma, sigmaErr), 0.045);
  myText (0.70, 0.70, kBlack, Form ("#chi^{2}/ndf = %.2f / %i", fit->GetChisquare (), fit->GetNDF ()), 0.045);
  myText (0.65, 0.65, kRed, "Expected", 0.045);
  myText (0.70, 0.60, kRed, "#mu = 0", 0.045);
  myText (0.70, 0.55, kRed, "#sigma = 1", 0.045);

  //c->SaveAs (Form ("%s/ChannelCompare/pull_%s_pp.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
  c->SaveAs (Form ("%s/ChannelCompare/pull_%s_PbPb.pdf", plotPath.Data (), useTrkPt ? "pTch" : "xhZ"));
}

#endif
