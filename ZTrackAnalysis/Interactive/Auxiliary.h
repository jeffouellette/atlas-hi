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

/**
 * Plots pull distributions between electron and muon channels. Can also provide a systematic uncertainty and add errors in quadrature.
 */
void PlotPullDist (const PhysicsAnalysis* data, const bool useTrkPt = true, const short pPtZ = -1, Systematic* sys = nullptr) {
  const bool doSpecificBin = (pPtZ >= 0 && pPtZ < nPtZBins);

  TCanvas* c = new TCanvas (Form ("c_pull_%s%s", useTrkPt ? "pTch" : "xhZ", doSpecificBin ? Form ("iPtZ%ipPtZ", pPtZ) : ""), "", 800, 600);
  TH1D* h_pull = new TH1D (Form ("h_pull_%s%s", useTrkPt ? "pTch" : "xhZ", doSpecificBin ? Form ("iPtZ%ipPtZ", pPtZ) : ""), Form (";Pull = (#LTY^{#it{ee}} (%s)#GT - #LTY^{#it{#mu#mu}} (%s)#GT) / #sigma;Entries", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}"), 21, -5, 5);

  const short iPtZLo = (doSpecificBin ? pPtZ : 2);
  const short iPtZHi = (doSpecificBin ? pPtZ+1 : 5);

  int nPoints = 0;
  for (short iCent = 1; iCent < numCentBins; iCent++) {
  //for (short iCent = 0; iCent < 1; iCent++) {
    for (short iPtZ = iPtZLo; iPtZ < iPtZHi; iPtZ++) {
      TH1D* h_ee = (useTrkPt ? data->h_trk_pt_ptz_sub : data->h_trk_xhz_ptz_sub)[0][iPtZ][iCent];
      TH1D* h_mumu = (useTrkPt ? data->h_trk_pt_ptz_sub : data->h_trk_xhz_ptz_sub)[1][iPtZ][iCent];
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

  TH1D* h_pull_expected = new TH1D (Form ("h_pull_expected_%s", useTrkPt ? "pTch" : "xhZ"), Form (";Pull = (#LTY^{#it{ee}} (%s)#GT - #LTY^{#it{#mu#mu}} (%s)#GT) / #sigma;Entries", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}", useTrkPt ? "#it{p}_{T}^{ch}" : "#it{x}_{hZ}"), 21, -5, 5);
  TF1* gaus = new TF1 ("gaus", "gaus(0)", -5, 5);
  gaus->SetParameter (0, 1);
  gaus->FixParameter (1, 0);
  gaus->FixParameter (2, 1);
  for (short i = 0; i < 2000; i++) {
    h_pull_expected->Fill (gaus->GetRandom ());
  }
  SaferDelete (&gaus);
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
  if (!doSpecificBin)           myText (0.20, 0.80, kBlack, "Pb+Pb, all points #it{p}_{T}^{Z} > 30 GeV", 0.045);
  else if (pPtZ != nPtZBins-1)  myText (0.20, 0.80, kBlack, Form ("Pb+Pb, all points %g < #it{p}_{T}^{Z} < %g GeV", zPtBins[pPtZ], zPtBins[pPtZ+1]), 0.045);
  else                          myText (0.20, 0.80, kBlack, Form ("Pb+Pb, all points #it{p}_{T}^{Z} > %g GeV", zPtBins[pPtZ]), 0.045);
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
    SaferDelete (&h_ratio);

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

    TF1* fit2 = new TF1 ("fit2", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
    fit2->SetParameter (0, 1);
    fit2->SetParameter (1, 0);
    fit2->SetParameter (2, 0);
    h_ratio->Fit (fit2, "RN0Q");
    SaferDelete (&h_ratio);

    fit2->SetLineColor (kAzure-1);
    fit2->SetLineStyle (2);
    fit2->SetLineWidth (2);
    fit2->Draw ("same");

    cout << "chi2/ndf = " << fit2->GetChisquare () << " / " << fit2->GetNDF () << endl;

    TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
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

    SaferDelete (&h_nom);
    SaferDelete (&h_up);
    SaferDelete (&h_down);
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
    SaferDelete (&h_ratio);

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
    SaferDelete (&h_ratio);

    fit2->SetLineColor (kAzure-1);
    fit2->SetLineStyle (2);
    fit2->SetLineWidth (2);
    fit2->Draw ("same");

    cout << "chi2/ndf = " << fit2->GetChisquare () << " / " << fit2->GetNDF () << endl;

    TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[iPtZ][0] : xhZBins[iPtZ][0], useTrkPt ? pTchBins[iPtZ][nPtchBins[iPtZ]] : xhZBins[iPtZ][nXhZBins[iPtZ]]);
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

    SaferDelete (&h_nom);
    SaferDelete (&h_up);
    SaferDelete (&h_down);
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
    SaferDelete (&eff_ratio);
    eff_hitight = (TH1D*) a_hitight->h2_num_trk_effs[iCent]->ProjectionY ("3");
    eff_ratio = (TH1D*) a_hitight->h2_den_trk_effs[iCent]->ProjectionY ("4");
    eff_hitight->Divide (eff_ratio);
    SaferDelete (&eff_ratio);

    eff_ratio = (TH1D*) eff_hitight->Clone ("eff_ratio");
    eff_ratio->Divide (eff_hiloose);
    SaferDelete (&eff_hiloose);
    SaferDelete (&eff_hitight);

    eff_hiloose = (TH1D*) a_hiloose->h2_num_trk_purs[iCent]->ProjectionY ("1");
    eff_hiloose->Divide ((TH1D*) a_hiloose->h2_den_trk_purs[iCent]->ProjectionY ("2"));
    eff_hitight = (TH1D*) a_hitight->h2_num_trk_purs[iCent]->ProjectionY ("3");
    eff_hitight->Divide ((TH1D*) a_hitight->h2_den_trk_purs[iCent]->ProjectionY ("4"));

    pur_ratio = (TH1D*) eff_hitight->Clone ("pur_ratio");
    pur_ratio->Divide (eff_hiloose);
    SaferDelete (&eff_hiloose);
    SaferDelete (&eff_hitight);

    const int bin1 = eff_ratio->FindBin (h_ratio->GetBinCenter (1));
    const int bin2 = pur_ratio->FindBin (h_ratio->GetBinCenter (1));
    for (int ix = 1; ix <= h_ratio->GetNbinsX (); ix++) {
      h_ratio->SetBinContent (ix, h_ratio->GetBinContent (ix) * pur_ratio->GetBinContent (ix+bin2-1) / eff_ratio->GetBinContent (ix+bin1-1));
      h_ratio->SetBinError (ix, h_ratio->GetBinError (ix) * pur_ratio->GetBinContent (ix+bin2-1) / eff_ratio->GetBinContent (ix+bin1-1));
    }


    g = make_graph (h_ratio);
    SaferDelete (&eff_ratio);
    SaferDelete (&pur_ratio);

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
    SaferDelete (&h_ratio);

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

    SaferDelete (&h_hiloose);
    SaferDelete (&h_hitight);
  }
  c->SaveAs (Form ("../Plots/TrackWPSystStudy/trackWP_pTch_iPtZ%i.pdf", iPtZ));
}




/**
 * Does low pT^Z channel comparison systematic study.
 */
void DoLowPtZSystStudy (PhysicsAnalysis* data, PhysicsAnalysis* bkg, const bool useTrkPt = true) {
  const char* canvasName = Form ("c_lowptz_syst_%s", useTrkPt ? "pTch" : "xhZ");
  const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
  TCanvas* c = nullptr;
  if (canvasExists)
    c = dynamic_cast <TCanvas*> (gDirectory->Get (canvasName));
  else {
    c = new TCanvas (canvasName, "", 1600, 800);
    gDirectory->Add (c);
  }
  c->cd ();

  const double dPadY = 0.5;
  const double uPadY = 1. - dPadY;
  const int axisTextSize = 23;

  TH1D* h_ee = nullptr, *h_mumu = nullptr, *h_comb = nullptr;
  TH1D* h_ee_bkg = nullptr, *h_mumu_bkg = nullptr, *h_comb_bkg = nullptr;
  TH1D* h_ee_diff = nullptr, *h_mumu_diff = nullptr;
  TGraphAsymmErrors* g = nullptr;

  TPad* ulPad = new TPad (Form ("ulPad_lowptz_syst_%s", useTrkPt ? "pTch" : "xhZ"), "", 0, dPadY, 0.25, 1);
  TPad* dlPad = new TPad (Form ("dlPad_lowptz_syst_%s", useTrkPt ? "pTch" : "xhZ"), "", 0, 0, 0.25, dPadY);
  TPad* uclPad = new TPad (Form ("uclPad_lowptz_syst_%s", useTrkPt ? "pTch" : "xhZ"), "", 0.25, dPadY, 0.5, 1);
  TPad* dclPad = new TPad (Form ("dclPad_lowptz_syst_%s", useTrkPt ? "pTch" : "xhZ"), "", 0.25, 0, 0.5, dPadY);
  TPad* ucrPad = new TPad (Form ("ucrPad_lowptz_syst_%s", useTrkPt ? "pTch" : "xhZ"), "", 0.5, dPadY, 0.75, 1);
  TPad* dcrPad = new TPad (Form ("dcrPad_lowptz_syst_%s", useTrkPt ? "pTch" : "xhZ"), "", 0.5, 0, 0.75, dPadY);
  TPad* urPad = new TPad (Form ("urPad_lowptz_syst_%s", useTrkPt ? "pTch" : "xhZ"), "", 0.75, dPadY, 1, 1);
  TPad* drPad = new TPad (Form ("drPad_lowptz_syst_%s", useTrkPt ? "pTch" : "xhZ"), "", 0.75, 0, 1, dPadY);

  ulPad->SetTopMargin (0.04);
  ulPad->SetBottomMargin (0.02);
  ulPad->SetLeftMargin (0.17);
  ulPad->SetRightMargin (0.06);
  dlPad->SetTopMargin (0.02);
  dlPad->SetBottomMargin (0.25);
  dlPad->SetLeftMargin (0.17);
  dlPad->SetRightMargin (0.06);
  uclPad->SetTopMargin (0.04);
  uclPad->SetBottomMargin (0.02);
  uclPad->SetLeftMargin (0.17);
  uclPad->SetRightMargin (0.06);
  dclPad->SetTopMargin (0.02);
  dclPad->SetBottomMargin (0.25);
  dclPad->SetLeftMargin (0.17);
  dclPad->SetRightMargin (0.06);
  ucrPad->SetTopMargin (0.04);
  ucrPad->SetBottomMargin (0.02);
  ucrPad->SetLeftMargin (0.17);
  ucrPad->SetRightMargin (0.06);
  dcrPad->SetTopMargin (0.02);
  dcrPad->SetBottomMargin (0.25);
  dcrPad->SetLeftMargin (0.17);
  dcrPad->SetRightMargin (0.06);
  urPad->SetTopMargin (0.04);
  urPad->SetBottomMargin (0.02);
  urPad->SetLeftMargin (0.17);
  urPad->SetRightMargin (0.06);
  drPad->SetTopMargin (0.02);
  drPad->SetBottomMargin (0.25);
  drPad->SetLeftMargin (0.17);
  drPad->SetRightMargin (0.06);

  ulPad->Draw ();
  dlPad->Draw ();
  uclPad->Draw ();
  dclPad->Draw ();
  ucrPad->Draw ();
  dcrPad->Draw ();
  urPad->Draw ();
  drPad->Draw ();

  TPad* uPads[4] = {ulPad, uclPad, ucrPad, urPad};
  TPad* dPads[4] = {dlPad, dclPad, dcrPad, drPad};

  TGAE** g_ee_arr = Get1DArray <TGAE*> (numCentBins);
  TGAE** g_mumu_arr = Get1DArray <TGAE*> (numCentBins);
  for (int iCent = 0; iCent < numCentBins; iCent++) {
    g_ee_arr[iCent] = new TGAE ();
    g_ee_arr[iCent]->SetName (Form ("g_ee_iCent%i", iCent));
    g_mumu_arr[iCent] = new TGAE ();
    g_mumu_arr[iCent]->SetName (Form ("g_mumu_iCent%i", iCent));
  }

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    uPads[iCent]->cd ();
    gPad->SetLogx ();
    gPad->SetLogy ();

    h_ee = (TH1D*) (useTrkPt ? data->h_trk_pt_ptz_sub : data->h_trk_xhz_ptz_sub)[0][2][iCent]->Clone ("h_ee");
    h_mumu = (TH1D*) (useTrkPt ? data->h_trk_pt_ptz_sub : data->h_trk_xhz_ptz_sub)[1][2][iCent]->Clone ("h_mumu");
    h_comb = (TH1D*) (useTrkPt ? data->h_trk_pt_ptz_sub : data->h_trk_xhz_ptz_sub)[2][2][iCent]->Clone ("h_comb");
    h_ee_bkg = (useTrkPt ? bkg->h_trk_pt_ptz : bkg->h_trk_xhz_ptz)[0][2][iCent];
    h_mumu_bkg = (useTrkPt ? bkg->h_trk_pt_ptz : bkg->h_trk_xhz_ptz)[1][2][iCent];
    h_comb_bkg = (useTrkPt ? bkg->h_trk_pt_ptz : bkg->h_trk_xhz_ptz)[2][2][iCent];

    if (!canvasExists) {
      const float min = (useTrkPt ? 1e-3 : 1e-2);
      const float max = (useTrkPt ? 2e1 : 8e2);

      TH1D* htemp = (TH1D*) h_ee->Clone ("htemp");
      htemp->Reset ();

      useTrkPt ? htemp->GetXaxis ()->SetLimits (pTchBins[2][0], pTchBins[2][nPtchBins[2]]) : htemp->GetXaxis ()->SetLimits (xhZBins[2][0], xhZBins[2][nXhZBins[2]]);
      htemp->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

      htemp->GetXaxis ()->SetMoreLogLabels ();

      htemp->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
      htemp->GetYaxis ()->SetTitle (useTrkPt ? "d^{2}Y / d#it{p}_{T}d#Delta#phi [GeV^{-1}]" : "d^{2}Y / d#it{x}_{hZ}d#Delta#phi");

      htemp->GetXaxis ()->SetTitleFont (43);
      htemp->GetXaxis ()->SetTitleSize (axisTextSize);
      htemp->GetXaxis ()->SetLabelFont (43);
      htemp->GetXaxis ()->SetLabelSize (axisTextSize);

      htemp->GetYaxis ()->SetTitleFont (43);
      htemp->GetYaxis ()->SetTitleSize (axisTextSize);
      htemp->GetYaxis ()->SetLabelFont (43);
      htemp->GetYaxis ()->SetLabelSize (axisTextSize);

      htemp->GetXaxis ()->SetTitleOffset (2.6 * htemp->GetXaxis ()->GetTitleOffset ());
      htemp->GetYaxis ()->SetTitleOffset (1.8 * htemp->GetYaxis ()->GetTitleOffset ());

      htemp->DrawCopy ("hist");
      SaferDelete (&htemp);
    }

    g = data->GetTGAE (h_ee);
    RecenterGraph (g);
    ResetXErrors (g);
    deltaize (g, 1.02, true);
    ResetXErrors (g);

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (colors[1]);
    g->SetLineColor (colors[1]);

    g->Draw ("P");

    g = data->GetTGAE (h_mumu);
    RecenterGraph (g);
    ResetXErrors (g);
    deltaize (g, 0.98, true);
    ResetXErrors (g);

    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (colors[2]);
    g->SetLineColor (colors[2]);

    g->Draw ("P");

    if (iCent == 0) {
      myText (0.36, 0.85, kBlack, "#bf{#it{ATLAS}} Internal", 0.040/uPadY);
      myText (0.24, 0.14, kBlack, "15 < #it{p}_{T}^{Z} < 30 GeV", 0.032/uPadY);
      myText (0.24, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.032/uPadY);
    }
    else
      myText (0.24, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.032/uPadY);

    if (iCent == 3) {
      myMarkerText (0.70, 0.85, colors[1], kFullCircle, "#it{ee}", 1.6, 0.032/uPadY);
      myMarkerText (0.70, 0.78, colors[2], kOpenCircle, "#it{#mu#mu}", 1.6, 0.032/uPadY);
    }


    dPads[iCent]->cd ();
    gPad->SetLogx ();

    if (!canvasExists) {
      TH1D* htemp = (TH1D*) h_ee->Clone ("htemp");
      htemp->Reset ();

      useTrkPt ? htemp->GetXaxis ()->SetLimits (pTchBins[2][0], pTchBins[2][nPtchBins[2]]) : htemp->GetXaxis ()->SetLimits (xhZBins[2][0], xhZBins[2][nXhZBins[2]]);
      if (iCent == 0) htemp->GetYaxis ()->SetRangeUser (-0.2, 0.2);
      else            htemp->GetYaxis ()->SetRangeUser (-1.5, 1.5);

      htemp->GetXaxis ()->SetMoreLogLabels ();

      htemp->GetXaxis ()->SetTitle (useTrkPt ? "#it{p}_{T}^{ ch} [GeV]" : "#it{x}_{hZ}");
      htemp->GetYaxis ()->SetTitle ("(Y_{#it{ll}} #minus Y_{comb.}) / Y_{comb.}");
      //htemp->GetYaxis ()->SetTitle (useTrkPt ? "(d^{2}Y_{#it{#mu#mu}}/d#it{p}_{T}d#Delta#phi #minus d^{2}Y_{#it{ee}}/d#it{p}_{T}d#Delta#phi)" : "(d^{2}Y_{#it{#mu#mu}}/d#it{x}d#Delta#phi #minus d^{2}Y_{#it{ee}}/d#it{x}d#Delta#phi)");
      //htemp->GetYaxis ()->SetTitle ("Muons #minus Electrons");

      htemp->GetXaxis ()->SetTitleFont (43);
      htemp->GetXaxis ()->SetTitleSize (axisTextSize);
      htemp->GetXaxis ()->SetLabelFont (43);
      htemp->GetXaxis ()->SetLabelSize (axisTextSize);

      htemp->GetYaxis ()->SetTitleFont (43);
      htemp->GetYaxis ()->SetTitleSize (axisTextSize);
      htemp->GetYaxis ()->SetLabelFont (43);
      htemp->GetYaxis ()->SetLabelSize (axisTextSize);

      htemp->GetXaxis ()->SetTitleOffset (2.6 * htemp->GetXaxis ()->GetTitleOffset ());
      htemp->GetYaxis ()->SetTitleOffset (1.8 * htemp->GetYaxis ()->GetTitleOffset ());

      htemp->GetYaxis ()->CenterTitle ();

      htemp->DrawCopy ("hist");
      SaferDelete (&htemp);

      //TLine* l = new TLine ();
      //l->SetLineColor (kBlack);
      //l->SetLineStyle (2);
      //l->DrawLine (useTrkPt ? pTchBins[2][0] : xhZBins[2][0], 1, useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]], 1);
      //l->DrawLine (useTrkPt ? pTchBins[2][0] : xhZBins[2][0], -1, useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]], -1);
    }

    TGAE* g_ee = g_ee_arr[iCent];
    TGAE* g_mumu = g_mumu_arr[iCent];
    h_ee_diff = (TH1D*) h_ee->Clone ("h_ee_diff");
    h_ee_diff->Reset ();
    h_mumu_diff = (TH1D*) h_mumu->Clone ("h_mumu_diff");
    h_mumu_diff->Reset ();

    const double wee = data->h_z_counts[0][2][iCent]->GetBinContent (2) / data->h_z_counts[2][2][iCent]->GetBinContent (2);
    const double wmm = data->h_z_counts[1][2][iCent]->GetBinContent (2) / data->h_z_counts[2][2][iCent]->GetBinContent (2);
    for (int ix = 1; ix <= h_ee_diff->GetNbinsX (); ix++) {
      const double yee = h_ee->GetBinContent (ix);
      const double yee_err = sqrt (pow (h_ee->GetBinError (ix), 2) + pow (h_ee_bkg->GetBinError (ix), 2));

      const double ymm = h_mumu->GetBinContent (ix);
      const double ymm_err = sqrt (pow (h_mumu->GetBinError (ix), 2) + pow (h_mumu_bkg->GetBinError (ix), 2));

      const double ycomb = h_comb->GetBinContent (ix);

      assert (yee != 0 && ymm != 0);

      const double fee = (yee-ycomb) / ycomb;
      const double fee_err = sqrt (pow (wmm*ymm*yee_err / pow (ycomb, 2), 2) + pow (-wmm*yee*ymm_err / pow (ycomb, 2), 2));

      h_ee_diff->SetBinContent (ix, fee);
      h_ee_diff->SetBinError (ix, fee_err);
      g_ee->SetPoint (g_ee->GetN (), h_ee_diff->GetBinCenter (ix), fee);
      g_ee->SetPointEYhigh (g_ee->GetN () - 1, fee_err);
      g_ee->SetPointEYlow (g_ee->GetN () - 1, fee_err);
      g_ee->SetPointEXhigh (g_ee->GetN () - 1, 0.5 * h_ee_diff->GetBinWidth (ix));
      g_ee->SetPointEXlow (g_ee->GetN () - 1, 0.5 * h_ee_diff->GetBinWidth (ix));

      const double fmm = (ymm-ycomb) / ycomb;
      const double fmm_err = sqrt (pow (wee*yee*ymm_err / pow (ycomb, 2), 2) + pow (-wee*ymm*yee_err / pow (ycomb, 2), 2));

      h_mumu_diff->SetBinContent (ix, fmm);
      h_mumu_diff->SetBinError (ix, fmm_err);
      g_mumu->SetPoint (g_mumu->GetN (), h_mumu_diff->GetBinCenter (ix), fmm);
      g_mumu->SetPointEYhigh (g_mumu->GetN () - 1, fmm_err);
      g_mumu->SetPointEYlow (g_mumu->GetN () - 1, fmm_err);
      g_mumu->SetPointEXhigh (g_mumu->GetN () - 1, 0.5 * h_mumu_diff->GetBinWidth (ix));
      g_mumu->SetPointEXlow (g_mumu->GetN () - 1, 0.5 * h_mumu_diff->GetBinWidth (ix));
    }

    g = make_graph (h_ee_diff);
    RecenterGraph (g);
    ResetXErrors (g);
    deltaize (g, 1.02, true);
    ResetXErrors (g);
    //if (iCent != 0) deltaize (g, 1+0.05*(iCent-2), true);

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (colors[1]);
    g->SetLineColor (colors[1]);

    g->Draw ("P");

    g = make_graph (h_mumu_diff);
    RecenterGraph (g);
    ResetXErrors (g);
    deltaize (g, 0.98, true);
    ResetXErrors (g);
    //if (iCent != 0) deltaize (g, 1+0.05*(iCent-2), true);

    g->SetMarkerStyle (kOpenCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (2);
    g->SetMarkerColor (colors[2]);
    g->SetLineColor (colors[2]);

    g->Draw ("P");

    SaferDelete (&h_ee);
    SaferDelete (&h_mumu);
    SaferDelete (&h_comb);
    SaferDelete (&h_ee_diff);
    SaferDelete (&h_mumu_diff);
  } // end loop over iCent

  for (int iCent = 0; iCent < numCentBins; iCent++) {

    dPads[iCent]->cd ();

    const bool doPP = (iCent == 0);

    g = g_ee_arr[iCent];

    //TF1* fit1 = new TF1 ("fit1", "[0]", useTrkPt ? pTchBins[2][0] : xhZBins[2][0], useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]]);
    TF1* fit1 = nullptr;
    if (doPP) {
      fit1 = new TF1 ("fit1", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[2][0] : xhZBins[2][0], useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]]);
      fit1->SetParameter (0, 0);
      fit1->SetParameter (1, 0);
      fit1->SetParameter (2, 0);
    }
    else {
      fit1 = new TF1 ("fit1", useTrkPt ? "[0]" : "[0]+[1]*log(x)", useTrkPt ? pTchBins[2][0] : xhZBins[2][0], useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]]);
      fit1->SetParameter (0, 0);
      if (!useTrkPt) fit1->SetParameter (1, 0);
    }
    g->Fit (fit1, "RN0Q");

    fit1->SetLineColor (colors[1]);
    fit1->SetLineStyle (2);
    fit1->SetLineWidth (2);
    fit1->Draw ("same");

    TF1* inv_fit1 = nullptr;
    if (doPP) {
      inv_fit1 = new TF1 ("inv_fit1", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[2][0] : xhZBins[2][0], useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]]);
      inv_fit1->SetParameter (0, -fit1->GetParameter (0));
      inv_fit1->SetParameter (1, -fit1->GetParameter (1));
      inv_fit1->SetParameter (2, -fit1->GetParameter (2));
    }
    else {
      inv_fit1 = new TF1 ("inv_fit1", useTrkPt ? "[0]" : "[0]+[1]*log(x)", useTrkPt ? pTchBins[2][0] : xhZBins[2][0], useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]]);
      inv_fit1->SetParameter (0, -fit1->GetParameter (0));
      if (!useTrkPt) inv_fit1->SetParameter (1, -fit1->GetParameter (1));
    }

    inv_fit1->SetLineColor (colors[1]);
    inv_fit1->SetLineStyle (2);
    inv_fit1->SetLineWidth (2);
    inv_fit1->Draw ("same");

    myText (0.52, 0.90, colors[1], Form ("#chi^{2}/ndf = %.2f/%i", fit1->GetChisquare (), fit1->GetNDF ()), 0.032/uPadY);

    cout << "chi2/ndf = " << fit1->GetChisquare () << " / " << fit1->GetNDF () << endl;

    g = g_mumu_arr[iCent];

    //TF1* fit2 = new TF1 ("fit2", "[0]", useTrkPt ? pTchBins[2][0] : xhZBins[2][0], useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]]);
    TF1* fit2 = nullptr;
    if (doPP) {
      fit2 = new TF1 ("fit2", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[2][0] : xhZBins[2][0], useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]]);
      fit2->SetParameter (0, 0);
      fit2->SetParameter (1, 0);
      fit2->SetParameter (2, 0);
    }
    else {
      fit2 = new TF1 ("fit2", useTrkPt ? "[0]" : "[0]+[1]*log(x)", useTrkPt ? pTchBins[2][0] : xhZBins[2][0], useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]]);
      fit2->SetParameter (0, 0);
      if (!useTrkPt) fit2->SetParameter (1, 0);
    }
    g->Fit (fit2, "RN0Q");

    fit2->SetLineColor (colors[2]);
    fit2->SetLineStyle (2);
    fit2->SetLineWidth (2);
    fit2->Draw ("same");

    TF1* inv_fit2 = nullptr;
    if (doPP) {
      inv_fit2 = new TF1 ("inv_fit2", "[0]+[1]*log(x)+[2]*log(x)*log(x)", useTrkPt ? pTchBins[2][0] : xhZBins[2][0], useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]]);
      inv_fit2->SetParameter (0, -fit2->GetParameter (0));
      inv_fit2->SetParameter (1, -fit2->GetParameter (1));
      inv_fit2->SetParameter (2, -fit2->GetParameter (2));
    }
    else {
      inv_fit2 = new TF1 ("inv_fit2", useTrkPt ? "[0]" : "[0]+[1]*log(x)", useTrkPt ? pTchBins[2][0] : xhZBins[2][0], useTrkPt ? pTchBins[2][nPtchBins[2]] : xhZBins[2][nXhZBins[2]]);
      inv_fit2->SetParameter (0, -fit2->GetParameter (0));
      if (!useTrkPt) inv_fit2->SetParameter (1, -fit2->GetParameter (1));
    }

    inv_fit2->SetLineColor (colors[2]);
    inv_fit2->SetLineStyle (2);
    inv_fit2->SetLineWidth (2);
    inv_fit2->Draw ("same");

    myText (0.52, 0.82, colors[2], Form ("#chi^{2}/ndf = %.2f/%i", fit2->GetChisquare (), fit2->GetNDF ()), 0.032/uPadY);
  } // end loop over iCent

  c->SaveAs (Form ("../Plots/LowPtZSystStudy/lowPtZ_%s.pdf", useTrkPt ? "pTch" : "xhZ"));
}

#endif
