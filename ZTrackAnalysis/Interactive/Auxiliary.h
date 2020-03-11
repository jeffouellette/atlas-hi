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
      TH1D* h_ee = (useTrkPt ? a->h_trk_pt_ptz_iaa : a->h_trk_xhz_ptz_iaa)[0][iPtZ][iCent];
      TH1D* h_mumu = (useTrkPt ? a->h_trk_pt_ptz_iaa : a->h_trk_xhz_ptz_iaa)[1][iPtZ][iCent];
      TGAE* g_ee_sys = (sys ? sys->GetTGAE ((useTrkPt ? sys->h_trk_pt_ptz_iaa : sys->h_trk_xhz_ptz_iaa)[0][iPtZ][iCent]) : nullptr);
      TGAE* g_mumu_sys = (sys ? sys->GetTGAE ((useTrkPt ? sys->h_trk_pt_ptz_iaa : sys->h_trk_xhz_ptz_iaa)[1][iPtZ][iCent]) : nullptr);
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




/**
 * Plots nominal vs. up/down variations on the muon energy scale, then fits to a line. (This is how the muon ES systematic is evaluated.)
 */
void DoMuonESSystStudy (PhysicsAnalysis* a_nom, PhysicsAnalysis* a_up, PhysicsAnalysis* a_down, const short iPtZ = nPtZBins-1, const bool useTrkPt = true) {
  const char* canvasName = Form ("c_muonES_syst_%s_iPtZ%i", useTrkPt ? "pTch" : "xhZ", iPtZ);
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1600, 800);
    gDirectory->Add (c);
    c->cd ();
  }

  const double dPadY = 0.5;
  const double uPadY = 1. - dPadY;
  const int axisTextSize = 23;

  TH1D* h_nom = nullptr, *h_up = nullptr, *h_down = nullptr, *h_ratio = nullptr;
  TGraphAsymmErrors* g = nullptr;

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    const char* uPadName = Form ("uPad_%i", iCent);
    const char* dPadName = Form ("dPad_%i", iCent);

    TPad* uPad = new TPad (uPadName, "", (1./(numCentBins))*(iCent), dPadY, (1./(numCentBins))*(iCent+1), 1);
    TPad* dPad = new TPad (dPadName, "", (1./(numCentBins))*(iCent), 0, (1./(numCentBins))*(iCent+1), dPadY);

    uPad->SetTopMargin (0.04);
    uPad->SetBottomMargin (0);
    uPad->SetLeftMargin (0.17);
    uPad->SetRightMargin (0.06);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    dPad->SetLeftMargin (0.17);
    dPad->SetRightMargin (0.06);
    uPad->Draw ();
    dPad->Draw ();


    uPad->cd ();
    uPad->SetLogx ();
    uPad->SetLogy ();

    h_nom = (TH1D*) (useTrkPt ? a_nom->h_trk_pt_ptz : a_nom->h_trk_xhz_ptz)[2][iPtZ][iCent]->Clone ("h_nom");
    h_up = (TH1D*) (useTrkPt ? a_up->h_trk_pt_ptz : a_up->h_trk_xhz_ptz)[2][iPtZ][iCent]->Clone ("h_up");
    h_down = (TH1D*) (useTrkPt ? a_down->h_trk_pt_ptz : a_down->h_trk_xhz_ptz)[2][iPtZ][iCent]->Clone ("h_down");

    const float min = fmin (fmin (h_nom->GetMinimum (0), h_up->GetMinimum (0)), h_down->GetMinimum (0));
    const float max = fmax (fmax (h_nom->GetMaximum (),  h_up->GetMaximum ()), h_down->GetMaximum (0));

    g = a_nom->GetTGAE (h_nom);

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);

    useTrkPt ? g->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
    g->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
    g->GetYaxis ()->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T}d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x}_{hZ}d#Delta#phi");

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->Draw ("AP");

    g = a_up->GetTGAE (h_up);
    deltaize (g, 0.95, true);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kRed+1);
    g->SetLineColor (kRed+1);

    g->Draw ("P");

    g = a_down->GetTGAE (h_down);
    deltaize (g, 1.05, true);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kAzure-1);
    g->SetLineColor (kAzure-1);

    g->Draw ("P");


    if (iCent == 0) {
      myText (0.40, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.040/uPadY);
      myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.036/uPadY);
    }
    else myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.036/uPadY);

    if (iCent == 1) myText (0.45, 0.85, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.036/uPadY);
    else if (iCent == 2) {
      if (iPtZ == nPtZBins-1) myText (0.50, 0.85, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.036/uPadY);
      else                    myText (0.40, 0.85, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.036/uPadY);
    }
    else if (iCent == 3) {
      myMarkerTextNoLine (0.50, 0.9, kBlack, kFullCircle, "Central values", 1.4, 0.032/uPadY);
      myMarkerTextNoLine (0.50, 0.82, kRed+1, kOpenSquare, "Muons Up", 1.4, 0.032/uPadY);
      myMarkerTextNoLine (0.50, 0.74, kAzure-1, kOpenSquare, "Muons Down", 1.4, 0.032/uPadY);
    }


    dPad->cd ();
    dPad->SetLogx ();

    h_ratio = (TH1D*) h_up->Clone ("h_ratio");
    h_ratio->Reset ();
    for (int ix = 1; ix <= h_ratio->GetNbinsX (); ix++) {
      const float yd = h_nom->GetBinContent (ix);
      const float yde = h_nom->GetBinError (ix);
      const float yn = h_up->GetBinContent (ix);
      const float yne = h_up->GetBinError (ix);
      h_ratio->SetBinContent (ix, yn/yd);
      h_ratio->SetBinError (ix, fabs (yn/yd) * sqrt (fabs (pow (yde/yd, 2) + pow (yne/yn, 2) - 2.*yne*yne/(yd*yn))));
    }

    g = make_graph (h_ratio);
    deltaize (g, 0.95, true);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kRed+1);
    g->SetLineColor (kRed+1);

    useTrkPt ? g->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
    g->GetYaxis ()->SetRangeUser (0.95, 1.05);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
    g->GetYaxis ()->SetTitle ("Variation / Nominal");

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->GetYaxis ()->CenterTitle ();

    g->Draw ("AP");

    //TF1* fit1 = new TF1 ("fit1", "[0]", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    TF1* fit1 = new TF1 ("fit1", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    fit1->SetParameter (0, 1);
    fit1->SetParameter (1, 0);
    fit1->SetParameter (2, 0);
    h_ratio->Fit (fit1, "RN0Q");
    SaferDelete (h_ratio);

    fit1->SetLineColor (kRed+1);
    fit1->SetLineStyle (2);
    fit1->SetLineWidth (2);
    fit1->Draw ("same");

    //TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    //inv_fit1->SetParameter (0, 1/fit1->GetParameter (0));
    //inv_fit1->SetParameter (1, -fit1->GetParameter (1)/fit1->GetParameter (0));
    inv_fit1->SetParameter (0, 2-fit1->GetParameter (0));
    inv_fit1->SetParameter (1, -fit1->GetParameter (1));
    inv_fit1->SetParameter (2, -fit1->GetParameter (2));

    inv_fit1->SetLineColor (kRed+1);
    inv_fit1->SetLineStyle (2);
    inv_fit1->SetLineWidth (2);
    inv_fit1->Draw ("same");

    cout << "chi2/ndf = " << fit1->GetChisquare () << " / " << fit1->GetNDF () << endl;


    h_ratio = (TH1D*) h_down->Clone ("h_ratio");
    h_ratio->Reset ();
    for (int ix = 1; ix <= h_ratio->GetNbinsX (); ix++) {
      const float yd = h_nom->GetBinContent (ix);
      const float yde = h_nom->GetBinError (ix);
      const float yn = h_down->GetBinContent (ix);
      const float yne = h_down->GetBinError (ix);
      h_ratio->SetBinContent (ix, yn/yd);
      h_ratio->SetBinError (ix, fabs (yn/yd) * sqrt (fabs (pow (yde/yd, 2) + pow (yne/yn, 2) - 2.*yne*yne/(yd*yn))));
    }

    g = make_graph (h_ratio);
    deltaize (g, 1.05, true);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kAzure-1);
    g->SetLineColor (kAzure-1);

    g->Draw ("P");

    //TF1* fit2 = new TF1 ("fit2", "[0]", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    TF1* fit2 = new TF1 ("fit2", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    fit2->SetParameter (0, 1);
    fit2->SetParameter (1, 0);
    fit2->SetParameter (2, 0);
    h_ratio->Fit (fit2, "RN0Q");
    SaferDelete (h_ratio);

    fit2->SetLineColor (kAzure-1);
    fit2->SetLineStyle (2);
    fit2->SetLineWidth (2);
    fit2->Draw ("same");

    cout << "chi2/ndf = " << fit2->GetChisquare () << " / " << fit2->GetNDF () << endl;

    //TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    //inv_fit2->SetParameter (0, 1./fit2->GetParameter (0));
    //inv_fit2->SetParameter (1, -fit2->GetParameter (1)/fit2->GetParameter (0));
    inv_fit2->SetParameter (0, 2-fit2->GetParameter (0));
    inv_fit2->SetParameter (1, -fit2->GetParameter (1));
    inv_fit2->SetParameter (2, -fit2->GetParameter (2));

    inv_fit2->SetLineColor (kAzure-1);
    inv_fit2->SetLineStyle (2);
    inv_fit2->SetLineWidth (2);
    inv_fit2->Draw ("same");

    TLine* l = new TLine (useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], 1, useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]], 1);
    l->SetLineColor (kBlack);
    l->SetLineStyle (2);
    l->Draw ("same");

    SaferDelete (h_nom);
    SaferDelete (h_up);
    SaferDelete (h_down);
  }

  c->SaveAs (Form ("../Plots/LeptonESSystStudy/muonES_%s_iPtZ%i.pdf", useTrkPt ? "pTch" : "xhZ", iPtZ));
}




/**
 * Plots nominal vs. up/down variations on the electron energy scale, then fits to a line. (This is how the electron ES systematic is evaluated.)
 */
void DoElectronESSystStudy (PhysicsAnalysis* a_nom, PhysicsAnalysis* a_up, PhysicsAnalysis* a_down, const short iPtZ = nPtZBins-1, const bool useTrkPt = true) {
  const char* canvasName = Form ("c_electronES_syst_%s_iPtZ%i", useTrkPt ? "pTch" : "xhZ", iPtZ);
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1600, 800);
    gDirectory->Add (c);
    c->cd ();
  }

  const double dPadY = 0.5;
  const double uPadY = 1. - dPadY;
  const int axisTextSize = 23;

  TH1D* h_nom = nullptr, *h_up = nullptr, *h_down = nullptr, *h_ratio = nullptr;
  TGraphAsymmErrors* g = nullptr;

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    const char* uPadName = Form ("uPad_%i", iCent);
    const char* dPadName = Form ("dPad_%i", iCent);

    TPad* uPad = new TPad (uPadName, "", (1./(numCentBins))*(iCent), dPadY, (1./(numCentBins))*(iCent+1), 1);
    TPad* dPad = new TPad (dPadName, "", (1./(numCentBins))*(iCent), 0, (1./(numCentBins))*(iCent+1), dPadY);

    uPad->SetTopMargin (0.04);
    uPad->SetBottomMargin (0);
    uPad->SetLeftMargin (0.17);
    uPad->SetRightMargin (0.06);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    dPad->SetLeftMargin (0.17);
    dPad->SetRightMargin (0.06);
    uPad->Draw ();
    dPad->Draw ();


    uPad->cd ();
    uPad->SetLogx ();
    uPad->SetLogy ();

    h_nom = (TH1D*) (useTrkPt ? a_nom->h_trk_pt_ptz : a_nom->h_trk_xhz_ptz)[2][iPtZ][iCent]->Clone ("h_nom");
    h_up = (TH1D*) (useTrkPt ? a_up->h_trk_pt_ptz : a_up->h_trk_xhz_ptz)[2][iPtZ][iCent]->Clone ("h_up");
    h_down = (TH1D*) (useTrkPt ? a_down->h_trk_pt_ptz : a_down->h_trk_xhz_ptz)[2][iPtZ][iCent]->Clone ("h_down");

    const float min = fmin (fmin (h_nom->GetMinimum (0), h_up->GetMinimum (0)), h_down->GetMinimum (0));
    const float max = fmax (fmax (h_nom->GetMaximum (),  h_up->GetMaximum ()), h_down->GetMaximum (0));

    g = a_nom->GetTGAE (h_nom);

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);

    useTrkPt ? g->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
    g->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
    g->GetYaxis ()->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T}d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x}_{hZ}d#Delta#phi");

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->Draw ("AP");


    g = a_up->GetTGAE (h_up);
    deltaize (g, 0.95, true);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kRed+1);
    g->SetLineColor (kRed+1);

    g->Draw ("P");

    g = a_down->GetTGAE (h_down);
    deltaize (g, 1.05, true);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kAzure-1);
    g->SetLineColor (kAzure-1);

    g->Draw ("P");


    if (iCent == 0) {
      myText (0.40, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.040/uPadY);
      myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.036/uPadY);
    }
    else myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.036/uPadY);

    if (iCent == 1) myText (0.45, 0.85, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.036/uPadY);
    else if (iCent == 2) {
      if (iPtZ == nPtZBins-1) myText (0.50, 0.85, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.036/uPadY);
      else                    myText (0.40, 0.85, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.036/uPadY);
    }
    else if (iCent == 3) {
      myMarkerTextNoLine (0.50, 0.9, kBlack, kFullCircle, "Central values", 1.4, 0.032/uPadY);
      myMarkerTextNoLine (0.50, 0.82, kRed+1, kOpenSquare, "Electrons Up", 1.4, 0.032/uPadY);
      myMarkerTextNoLine (0.50, 0.74, kAzure-1, kOpenSquare, "Electrons Down", 1.4, 0.032/uPadY);
    }


    dPad->cd ();
    dPad->SetLogx ();

    h_ratio = (TH1D*) h_up->Clone ("h_ratio");
    h_ratio->Reset ();
    for (int ix = 1; ix <= h_ratio->GetNbinsX (); ix++) {
      const float yd = h_nom->GetBinContent (ix);
      const float yde = h_nom->GetBinError (ix);
      const float yn = h_up->GetBinContent (ix);
      const float yne = h_up->GetBinError (ix);
      h_ratio->SetBinContent (ix, yn/yd);
      h_ratio->SetBinError (ix, fabs (yn/yd) * sqrt (fabs (pow (yde/yd, 2) + pow (yne/yn, 2) - 2.*yne*yne/(yd*yn))));
    }

    g = make_graph (h_ratio);
    deltaize (g, 0.95, true);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kRed+1);
    g->SetLineColor (kRed+1);

    useTrkPt ? g->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]) : g->GetXaxis ()->SetLimits (xhZBins[iPtZ][0], xhZBins[iPtZ][nXhZBins[iPtZ]]);
    g->GetYaxis ()->SetRangeUser (0.95, 1.05);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
    g->GetYaxis ()->SetTitle ("Variation / Nominal");

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->GetYaxis ()->CenterTitle ();

    g->Draw ("AP");

    //TF1* fit1 = new TF1 ("fit1", "[0]", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    TF1* fit1 = new TF1 ("fit1", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    fit1->SetParameter (0, 1);
    fit1->SetParameter (1, 0);
    h_ratio->Fit (fit1, "RN0Q");
    SaferDelete (h_ratio);

    fit1->SetLineColor (kRed+1);
    fit1->SetLineStyle (2);
    fit1->SetLineWidth (2);
    fit1->Draw ("same");

    //TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    //inv_fit1->SetParameter (0, 1/fit1->GetParameter (0));
    //inv_fit1->SetParameter (1, -fit1->GetParameter (1)/fit1->GetParameter (0));
    inv_fit1->SetParameter (2, 2-fit1->GetParameter (0));
    inv_fit1->SetParameter (2, -fit1->GetParameter (1));
    inv_fit1->SetParameter (2, -fit1->GetParameter (2));

    inv_fit1->SetLineColor (kRed+1);
    inv_fit1->SetLineStyle (2);
    inv_fit1->SetLineWidth (2);
    inv_fit1->Draw ("same");

    cout << "chi2/ndf = " << fit1->GetChisquare () << " / " << fit1->GetNDF () << endl;


    h_ratio = (TH1D*) h_down->Clone ("h_ratio");
    h_ratio->Reset ();
    for (int ix = 1; ix <= h_ratio->GetNbinsX (); ix++) {
      const float yd = h_nom->GetBinContent (ix);
      const float yde = h_nom->GetBinError (ix);
      const float yn = h_down->GetBinContent (ix);
      const float yne = h_down->GetBinError (ix);
      h_ratio->SetBinContent (ix, yn/yd);
      h_ratio->SetBinError (ix, fabs (yn/yd) * sqrt (fabs (pow (yde/yd, 2) + pow (yne/yn, 2) - 2.*yne*yne/(yd*yn))));
    }

    g = make_graph (h_ratio);
    deltaize (g, 1.05, true);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kAzure-1);
    g->SetLineColor (kAzure-1);

    g->Draw ("P");

    //TF1* fit2 = new TF1 ("fit2", "[0]", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    TF1* fit2 = new TF1 ("fit2", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    fit2->SetParameter (0, 1);
    fit2->SetParameter (1, 0);
    fit2->SetParameter (2, 0);
    h_ratio->Fit (fit2, "RN0Q");
    SaferDelete (h_ratio);

    fit2->SetLineColor (kAzure-1);
    fit2->SetLineStyle (2);
    fit2->SetLineWidth (2);
    fit2->Draw ("same");

    cout << "chi2/ndf = " << fit2->GetChisquare () << " / " << fit2->GetNDF () << endl;

    //TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    //inv_fit2->SetParameter (0, 1./fit2->GetParameter (0));
    //inv_fit2->SetParameter (1, fit2->GetParameter (1)/fit2->GetParameter (0));
    inv_fit2->SetParameter (2, 2-fit2->GetParameter (0));
    inv_fit2->SetParameter (2, -fit2->GetParameter (1));
    inv_fit2->SetParameter (2, -fit2->GetParameter (2));

    inv_fit2->SetLineColor (kAzure-1);
    inv_fit2->SetLineStyle (2);
    inv_fit2->SetLineWidth (2);
    inv_fit2->Draw ("same");

    TLine* l = new TLine (useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], 1, useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]], 1);
    l->SetLineColor (kBlack);
    l->SetLineStyle (2);
    l->Draw ("same");

    SaferDelete (h_nom);
    SaferDelete (h_up);
    SaferDelete (h_down);
  }

  c->SaveAs (Form ("../Plots/LeptonESSystStudy/electronES_%s_iPtZ%i.pdf", useTrkPt ? "pTch" : "xhZ", iPtZ));
}




/**
 * Plots nominal (HILoose) vs. HITight variation on the track selection, then fits to a constant. (This is how the track quality WP systematic is evaluated.)
 */
void DoTrackWPSystStudy (PhysicsAnalysis* a_hiloose, PhysicsAnalysis* a_hitight, const short iPtZ = nPtZBins-1) {
  const char* canvasName = Form ("c_trackWPsyst_pTch_iPtZ%i", iPtZ);
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1600, 800);
    gDirectory->Add (c);
    c->cd ();
  }

  const double dPadY = 0.5;
  const double uPadY = 1. - dPadY;
  const int axisTextSize = 23;

  TH1D* h_hiloose = nullptr, *h_hitight = nullptr, *h_ratio = nullptr, *eff_hiloose = nullptr, *eff_hitight = nullptr, *eff_ratio = nullptr, *pur_ratio = nullptr;
  TGraphAsymmErrors* g = nullptr;

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    const char* uPadName = Form ("uPad_%i", iCent);
    const char* dPadName = Form ("dPad_%i", iCent);

    TPad* uPad = new TPad (uPadName, "", (1./(numCentBins))*(iCent), dPadY, (1./(numCentBins))*(iCent+1), 1);
    TPad* dPad = new TPad (dPadName, "", (1./(numCentBins))*(iCent), 0, (1./(numCentBins))*(iCent+1), dPadY);

    uPad->SetTopMargin (0.04);
    uPad->SetBottomMargin (0);
    uPad->SetLeftMargin (0.17);
    uPad->SetRightMargin (0.06);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.25);
    dPad->SetLeftMargin (0.17);
    dPad->SetRightMargin (0.06);
    uPad->Draw ();
    dPad->Draw ();


    uPad->cd ();
    uPad->SetLogx ();
    uPad->SetLogy ();

    h_hiloose = (TH1D*) a_hiloose->h_trk_pt_dphi_raw[2][iPtZ][1][iCent]->Clone ("h_hiloose");
    h_hitight = (TH1D*) a_hitight->h_trk_pt_dphi_raw[2][iPtZ][1][iCent]->Clone ("h_hitight");
    for (int idPhi = 2; idPhi < numPhiBins; idPhi++) {
      h_hiloose->Add (a_hiloose->h_trk_pt_dphi_raw[2][iPtZ][idPhi][iCent]);
      h_hitight->Add (a_hitight->h_trk_pt_dphi_raw[2][iPtZ][idPhi][iCent]);
    }
    for (int iX = 1; iX <= h_hiloose->GetNbinsX (); iX++)
      h_hiloose->SetBinError (iX, sqrt (h_hiloose->GetBinContent (iX)));
    for (int iX = 1; iX <= h_hitight->GetNbinsX (); iX++)
      h_hitight->SetBinError (iX, sqrt (h_hitight->GetBinContent (iX)));

    const float min = fmin (h_hiloose->GetMinimum (0), h_hitight->GetMinimum (0));
    const float max = fmax (h_hiloose->GetMaximum (),  h_hitight->GetMaximum ());

    g = a_hiloose->GetTGAE (h_hiloose);

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kBlack);
    g->SetLineColor (kBlack);

    g->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
    g->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    g->GetYaxis ()->SetTitle ("N_{ch}^{total}");

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->Draw ("AP");


    g = a_hitight->GetTGAE (h_hitight);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kRed+1);
    g->SetLineColor (kRed+1);

    g->Draw ("P");

    if (iCent == 0) {
      myText (0.40, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.040/uPadY);
      myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.036/uPadY);
      //if (iSpc == 0) myText (0.44, 0.80, kBlack, "#it{Z} #rightarrow e^{+}e^{-} events", 0.036/uPadY);
      //if (iSpc == 1) myText (0.44, 0.80, kBlack, "#it{Z} #rightarrow #mu^{+}#mu^{-} events", 0.036/uPadY);
    }
    else myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.036/uPadY);

    if (iCent == 1) myText (0.45, 0.85, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.036/uPadY);
    else if (iCent == 2) {
      if (iPtZ == nPtZBins-1) myText (0.50, 0.85, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.036/uPadY);
      else                    myText (0.40, 0.85, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.036/uPadY);
    }
    else if (iCent == 3) {
      myMarkerTextNoLine (0.55, 0.9, kBlack, kFullCircle, "HILoose", 1.4, 0.036/uPadY);
      myMarkerTextNoLine (0.55, 0.82, kRed+1, kOpenSquare, "HITight", 1.4, 0.036/uPadY);
    }


    dPad->cd ();
    dPad->SetLogx ();

    h_ratio = (TH1D*) h_hitight->Clone ("h_ratio");
    h_ratio->Reset ();
    //h_ratio->Divide (h_hiloose);
    for (int ix = 1; ix <= h_ratio->GetNbinsX (); ix++) {
      const float y1 = h_hiloose->GetBinContent (ix);
      const float y1e = h_hiloose->GetBinError (ix);
      const float y2 = h_hitight->GetBinContent (ix);
      const float y2e = h_hitight->GetBinError (ix);
      h_ratio->SetBinContent (ix, y1/y2);
      h_ratio->SetBinError (ix, fabs(y1/y2)*sqrt (fabs (pow (y1e/y1,2) + pow (y2e/y2,2) - 2.*y1e*y1e/(y1*y2))));
    }
    for (int ix = 1; ix <= h_ratio->GetNbinsX (); ix++) {
      const float passes = h_hitight->GetBinContent (ix);
      const float trials = h_hiloose->GetBinContent (ix);
      h_ratio->SetBinContent (ix, (passes+1) / (trials+2));
      h_ratio->SetBinError (ix, sqrt ((passes+1)*(passes+2)/((trials+2)*(trials+3)) - pow (passes+1, 2) / pow (trials+2, 2)));

      //h_ratio->SetBinError (ix, sqrt ((h_ratio->GetBinContent (ix)) * (1-h_ratio->GetBinContent (ix)) / h_hiloose->GetBinContent (ix)));
    }

    a_hiloose->LoadTrackingEfficiencies (true);
    a_hitight->LoadTrackingEfficiencies (true);
    a_hiloose->LoadTrackingPurities (true);
    a_hitight->LoadTrackingPurities (true);

    eff_hiloose = (TH1D*) a_hiloose->h2_num_trk_effs[iCent]->ProjectionY ("1");
    eff_ratio = (TH1D*) a_hiloose->h2_den_trk_effs[iCent]->ProjectionY ("2");
    eff_hiloose->Divide (eff_ratio);
    SaferDelete (eff_ratio);
    eff_hitight = (TH1D*) a_hitight->h2_num_trk_effs[iCent]->ProjectionY ("3");
    eff_ratio = (TH1D*) a_hitight->h2_den_trk_effs[iCent]->ProjectionY ("4");
    eff_hitight->Divide (eff_ratio);
    SaferDelete (eff_ratio);

    eff_ratio = (TH1D*) eff_hitight->Clone ("eff_ratio");
    eff_ratio->Divide (eff_hiloose);
    SaferDelete (eff_hiloose);
    SaferDelete (eff_hitight);

    eff_hiloose = (TH1D*) a_hiloose->h2_num_trk_purs[iCent]->ProjectionY ("1");
    eff_hiloose->Divide ((TH1D*) a_hiloose->h2_den_trk_purs[iCent]->ProjectionY ("2"));
    eff_hitight = (TH1D*) a_hitight->h2_num_trk_purs[iCent]->ProjectionY ("3");
    eff_hitight->Divide ((TH1D*) a_hitight->h2_den_trk_purs[iCent]->ProjectionY ("4"));

    pur_ratio = (TH1D*) eff_hitight->Clone ("pur_ratio");
    pur_ratio->Divide (eff_hiloose);
    SaferDelete (eff_hiloose);
    SaferDelete (eff_hitight);

    const int bin1 = eff_ratio->FindBin (h_ratio->GetBinCenter (1));
    const int bin2 = pur_ratio->FindBin (h_ratio->GetBinCenter (1));
    for (int ix = 1; ix <= h_ratio->GetNbinsX (); ix++) {
      h_ratio->SetBinContent (ix, h_ratio->GetBinContent (ix) * pur_ratio->GetBinContent (ix+bin2-1) / eff_ratio->GetBinContent (ix+bin1-1));
      h_ratio->SetBinError (ix, h_ratio->GetBinError (ix) * pur_ratio->GetBinContent (ix+bin2-1) / eff_ratio->GetBinContent (ix+bin1-1));
    }


    g = make_graph (h_ratio);
    SaferDelete (eff_ratio);
    SaferDelete (pur_ratio);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (kRed+1);
    g->SetLineColor (kRed+1);

    g->GetXaxis ()->SetLimits (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
    g->GetYaxis ()->SetRangeUser (0.90, 1.10);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    g->GetYaxis ()->SetTitle ("Performance Corr. Ratio");

    g->GetXaxis ()->SetTitleFont (43);
    g->GetXaxis ()->SetTitleSize (axisTextSize);
    g->GetXaxis ()->SetLabelFont (43);
    g->GetXaxis ()->SetLabelSize (axisTextSize);

    g->GetYaxis ()->SetTitleFont (43);
    g->GetYaxis ()->SetTitleSize (axisTextSize);
    g->GetYaxis ()->SetLabelFont (43);
    g->GetYaxis ()->SetLabelSize (axisTextSize);

    g->GetXaxis ()->SetTitleOffset (2.6 * g->GetXaxis ()->GetTitleOffset ());
    g->GetYaxis ()->SetTitleOffset (1.8 * g->GetYaxis ()->GetTitleOffset ());

    g->GetYaxis ()->CenterTitle ();

    g->Draw ("AP");

    TF1* fit1 = new TF1 ("fit1", "[0]", pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
    //TF1* fit1 = new TF1 ("fit1", "[0]+[1]*log(x)", pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
    fit1->SetParameter (0, 1);
    //fit1->SetParameter (1, 0);
    h_ratio->Fit (fit1, "RQN0");
    SaferDelete (h_ratio);

    fit1->SetLineColor (kRed+1);
    fit1->SetLineStyle (2);
    fit1->SetLineWidth (2);
    fit1->Draw ("same");

    TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]", pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
    //TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]-[1]*log(x)", pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
    inv_fit1->SetParameter (0, 1/fit1->GetParameter (0));
    //inv_fit1->SetParameter (1, fit1->GetParameter (1)/fit1->GetParameter (0));

    inv_fit1->SetLineColor (kRed+1);
    inv_fit1->SetLineStyle (2);
    inv_fit1->SetLineWidth (2);
    inv_fit1->Draw ("same");

    cout << "chi2/ndf = " << fit1->GetChisquare () << " / " << fit1->GetNDF () << endl;

    TLine* l = new TLine (pTchBins[iPtZ][0], 1, pTchBins[iPtZ][nPtchBins[iPtZ]], 1);
    l->SetLineColor (kBlack);
    l->SetLineStyle (2);
    l->Draw ("same");

    myText (0.22, 0.3, kBlack, Form ("Avg. = %s", FormatMeasurement (fit1->GetParameter (0), fit1->GetParError (0), 2)), 0.03/dPadY);

    SaferDelete (h_hiloose);
    SaferDelete (h_hitight);
  }
  c->SaveAs (Form ("../Plots/TrackWPSystStudy/trackWP_pTch_iPtZ%i.pdf", iPtZ));
}

#endif
