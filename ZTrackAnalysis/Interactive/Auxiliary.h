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




void DoMuonESSystStudy (PhysicsAnalysis* a1, PhysicsAnalysis* a2, PhysicsAnalysis* a3, const short iPtZ = nPtZBins-1) {
  TCanvas* c = new TCanvas ("c", "", 1600, 800);
  const double dPadY = 0.5;
  const double uPadY = 1. - dPadY;
  const int axisTextSize = 23;

  TH1D* h1 = nullptr, *h2 = nullptr, *h3 = nullptr, *hrat = nullptr;
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

    h1 = a1->h_trk_pt_ptz[2][iPtZ][iCent];
    h2 = a2->h_trk_pt_ptz[2][iPtZ][iCent];
    h3 = a3->h_trk_pt_ptz[2][iPtZ][iCent];

    const float min = fmin (fmin (h1->GetMinimum (0), h2->GetMinimum (0)), h3->GetMinimum (0));
    const float max = fmax (fmax (h1->GetMaximum (),  h2->GetMaximum ()), h3->GetMaximum (0));

    g = a1->GetTGAE (h1);

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[2]);
    g->SetLineColor (colors[2]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtchBins[0], allPtchBins[maxNPtchBins]);
    g->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    //g->GetYaxis ()->SetTitle ("N_{ch}^{total}");
    g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T}d#Delta#phi [GeV^{-1}]");
    //g->GetYaxis ()->SetTitle ("Y / Y_{bkg}");

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


    g = a2->GetTGAE (h2);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[1]);
    g->SetLineColor (colors[1]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtchBins[0], allPtchBins[maxNPtchBins]);
    g->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->Draw ("P");

    g = a3->GetTGAE (h3);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[3]);
    g->SetLineColor (colors[3]);

    g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    //g->GetXaxis ()->SetLimits (allPtchBins[0], allPtchBins[maxNPtchBins]);
    g->GetYaxis ()->SetRangeUser (min, max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->Draw ("P");


    if (iCent == 0) {
      myText (0.44, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.040/uPadY);
      myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.036/uPadY);
    }
    else
      myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.036/uPadY);

    if (iCent == 1)
      myText (0.50, 0.88, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.036/uPadY);
    else if (iCent == 2) {
      if (iPtZ == nPtZBins-1)
        myText (0.50, 0.88, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.036/uPadY);
      else
        myText (0.50, 0.88, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.036/uPadY);
    }
    else if (iCent == 3) {
      myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "Central values", 1.4, 0.04/uPadY);
      myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "Muons Up", 1.4, 0.04/uPadY);
      myMarkerTextNoLine (0.50, 0.74, colors[3], kOpenSquare, "Muons Down", 1.4, 0.04/uPadY);
    }


    dPad->cd ();
    dPad->SetLogx ();

    hrat = (TH1D*) h2->Clone ("hrat");
    hrat->Reset ();
    for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
      const float y1 = h1->GetBinContent (ix);
      const float y1e = h1->GetBinError (ix);
      const float y2 = h2->GetBinContent (ix);
      const float y2e = h2->GetBinError (ix);
      hrat->SetBinContent (ix, y2/y1);
      hrat->SetBinError (ix, fabs (y2/y1) * sqrt (fabs (pow (y1e/y1, 2) + pow (y2e/y2, 2) - 2.*y2e*y2e/(y1*y2))));
    }

    g = make_graph (hrat);
    delete hrat;

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[1]);
    g->SetLineColor (colors[1]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtchBins[0], allPtchBins[maxNPtchBins]);
    g->GetYaxis ()->SetRangeUser (0.95, 1.05);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    //g->GetYaxis ()->SetTitle ("Z#rightarrow#mu#mu / Z#rightarrowee");
    g->GetYaxis ()->SetTitle ("Variation / Nominal");
    //g->GetYaxis ()->SetTitle ("Eff. Corrected / Uncorrected");
    //g->GetYaxis ()->SetTitle ("New ES / Old ES");
    //g->GetYaxis ()->SetTitle ("Unfolded / not unfolded");

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

    //TF1* fit1 = new TF1 ("fit1", "[0]", allPtchBins[0], allPtchBins[maxNPtchBins]);
    //TF1* fit1 = new TF1 ("fit1", "[0]+[1]*log(x)", allXHZBins[0], allXHZBins[maxNXHZBins]);
    TF1* fit1 = new TF1 ("fit1", "[0]+[1]*log(x)", allPtchBins[0], allPtchBins[maxNPtchBins]);
    fit1->SetParameter (0, 1);
    fit1->SetParameter (1, 0);
    g->Fit (fit1, "RQN0");

    fit1->SetLineColor (colors[1]);
    fit1->SetLineStyle (2);
    fit1->SetLineWidth (1);
    fit1->Draw ("same");

    //TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]", allPtchBins[0], allPtchBins[maxNPtchBins]);
    //TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]-[1]*log(x)", allXHZBins[0], allXHZBins[maxNXHZBins]);
    TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]-[1]*log(x)", allPtchBins[0], allPtchBins[maxNPtchBins]);
    inv_fit1->SetParameter (0, 1/fit1->GetParameter (0));
    inv_fit1->SetParameter (1, fit1->GetParameter (1)/fit1->GetParameter (0));

    inv_fit1->SetLineColor (colors[1]);
    inv_fit1->SetLineStyle (2);
    inv_fit1->SetLineWidth (1);
    inv_fit1->Draw ("same");

    cout << "chi2/ndf = " << fit1->GetChisquare () << " / " << fit1->GetNDF () << endl;


    hrat = (TH1D*) h3->Clone ("hrat");
    hrat->Reset ();
    for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
      const float y1 = h1->GetBinContent (ix);
      const float y1e = h1->GetBinError (ix);
      const float y3 = h3->GetBinContent (ix);
      const float y3e = h3->GetBinError (ix);
      hrat->SetBinContent (ix, y3/y1);
      //hrat->SetBinError (ix, y3e/y1);
      hrat->SetBinError (ix, fabs (y3/y1) * sqrt (fabs (pow (y1e/y1, 2) + pow (y3e/y3, 2) - 2.*y3e*y3e/(y1*y3))));
    }

    g = make_graph (hrat);
    delete hrat;

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[3]);
    g->SetLineColor (colors[3]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtchBins[0], allPtchBins[maxNPtchBins]);
    g->GetYaxis ()->SetRangeUser (0.95, 1.05);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->Draw ("P");

    //TF1* fit2 = new TF1 ("fit2", "[0]+[1]*log(x)", allXHZBins[0], allXHZBins[maxNXHZBins]);
    TF1* fit2 = new TF1 ("fit2", "[0]+[1]*log(x)", allPtchBins[0], allPtchBins[maxNPtchBins]);
    fit2->SetParameter (0, 1);
    fit2->SetParameter (1, 0);
    g->Fit (fit2, "RQN0");

    fit2->SetLineColor (colors[3]);
    fit2->SetLineStyle (2);
    fit2->SetLineWidth (1);
    fit2->Draw ("same");

    cout << "chi2/ndf = " << fit2->GetChisquare () << " / " << fit2->GetNDF () << endl;

    //TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]-[1]*log(x)", allXHZBins[0], allXHZBins[maxNPtchBins]);
    TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]-[1]*log(x)", allPtchBins[0], allPtchBins[maxNPtchBins]);
    inv_fit2->SetParameter (0, 1./fit2->GetParameter (0));
    inv_fit2->SetParameter (1, fit2->GetParameter (1)/fit2->GetParameter (0));

    inv_fit2->SetLineColor (colors[3]);
    inv_fit2->SetLineStyle (2);
    inv_fit2->SetLineWidth (1);
    inv_fit2->Draw ("same");

    TLine* l = new TLine (allPtchBins[0], 1, allPtchBins[maxNPtchBins], 1);
    l->SetLineColor (kBlack);
    l->SetLineStyle (2);
    l->Draw ("same");

    //myText (0.22, 0.3, kBlack, Form ("Avg. = %s", FormatMeasurement (fit1->GetParameter (0), fit1->GetParError (0), 2)), 0.03/dPadY);

    delete h1, h2;
  }

}



/*
void DoTrackWPSystStudy (const short iPtZ = nPtZBins-1) {
  TCanvas* c = new TCanvas ("c", "", 1600, 800);
  const double dPadY = 0.5;
  const double uPadY = 1. - dPadY;
  const int axisTextSize = 23;

  TH1D* h1 = nullptr, *h2 = nullptr, *h3 = nullptr, *hrat = nullptr, *eff1 = nullptr, *eff2 = nullptr, *effrat = nullptr, *purrat = nullptr;
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

    PhysicsAnalysis* a1, *a2 = nullptr, *a3 = nullptr;

    a1 = data18;
    a2 = data_muonPtUp;
    a3 = data_muonPtDown;

    h1 = a1->h_trk_pt_ptz[2][iPtZ][iCent];
    h2 = a2->h_trk_pt_ptz[2][iPtZ][iCent];
    h3 = a3->h_trk_pt_ptz[2][iPtZ][iCent];

    //float i1 = 0, i2 = 0;
    //for (int ix = 1; ix <= h1->GetNbinsX (); ix++) {
    //  cout << h2->GetBinContent (ix) << endl;
    //  i1 += h1->GetBinContent (ix) * h1->GetBinCenter (ix) * h1->GetBinWidth (ix);
    //}
    //for (int ix = 1; ix <= h1->GetNbinsX (); ix++) {
    //  cout << h2->GetBinContent (ix) << endl;
    //  i2 += h2->GetBinContent (ix) * h2->GetBinCenter (ix) * h2->GetBinWidth (ix);
    //}

    //i1 *= pi/4;
    //i2 *= pi/4;

    //cout << "iPtZ = " << iPtZ << ", iCent = " << iCent << endl;
    //cout << "i1 = " << i1 << endl;
    //cout << "i2 = " << i2 << endl;

    //h3 = a3->h_trk_xhz_ptz[iSpc][iPtZ][iCent];
    //h1 = (TH1D*) a1->h_z_trk_raw_pt[iSpc][iPtZ][1][iCent]->Clone ("h1");
    //h1->Add (a1->h_z_trk_raw_pt[iSpc][iPtZ][2][iCent]);
    //h2 = (TH1D*) a2->h_z_trk_raw_pt[iSpc][iPtZ][1][iCent]->Clone ("h2");
    //h2->Add (a2->h_z_trk_raw_pt[iSpc][iPtZ][2][iCent]);
    //h3 = (TH1D*) a3->h_z_trk_raw_pt[iSpc][iPtZ][1][iCent]->Clone ("h2");
    //h3->Add (a3->h_z_trk_raw_pt[iSpc][iPtZ][2][iCent]);
    //const float min = fmin (h1->GetMinimum (0), h2->GetMinimum (0));
    //const float max = fmax (h1->GetMaximum (),  h2->GetMaximum ());
    const float min = fmin (fmin (h1->GetMinimum (0), h2->GetMinimum (0)), h3->GetMinimum (0));
    const float max = fmax (fmax (h1->GetMaximum (),  h2->GetMaximum ()), h3->GetMaximum (0));

    g = a1->GetTGAE (h1);

    g->SetMarkerStyle (kFullCircle);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[2]);
    g->SetLineColor (colors[2]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtchBins[0], allPtchBins[maxNPtchBins]);
    g->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    //g->GetYaxis ()->SetTitle ("N_{ch}^{total}");
    g->GetYaxis ()->SetTitle ("d^{2}Y / d#it{p}_{T}d#Delta#phi [GeV^{-1}]");
    //g->GetYaxis ()->SetTitle ("Y / Y_{bkg}");

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


    g = a2->GetTGAE (h2);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[1]);
    g->SetLineColor (colors[1]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtchBins[0], allPtchBins[maxNPtchBins]);
    g->GetYaxis ()->SetRangeUser (0.5*min, 2*max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->Draw ("P");

    g = a3->GetTGAE (h3);

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[3]);
    g->SetLineColor (colors[3]);

    g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    //g->GetXaxis ()->SetLimits (allPtchBins[0], allPtchBins[maxNPtchBins]);
    g->GetYaxis ()->SetRangeUser (min, max);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->Draw ("P");


    if (iCent == 0) {
      myText (0.44, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.040/uPadY);
      myText (0.22, 0.06, kBlack, "#it{pp}, 5.02 TeV", 0.036/uPadY);
      if (iSpc == 0)
        myText (0.44, 0.80, kBlack, "Z #rightarrow e^{+}e^{-} events", 0.036/uPadY);
      if (iSpc == 1)
        myText (0.44, 0.80, kBlack, "Z #rightarrow #mu^{+}#mu^{-} events", 0.036/uPadY);
    }
    else
      myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.036/uPadY);

    if (iCent == 1)
      myText (0.50, 0.88, kBlack, "3#pi/4 < |#Delta#phi| < #pi", 0.036/uPadY);
    else if (iCent == 2) {
      if (iPtZ == nPtZBins-1)
        myText (0.50, 0.88, kBlack, Form ("#it{p}_{T}^{Z} > %g GeV", zPtBins[iPtZ]), 0.036/uPadY);
      else
        myText (0.50, 0.88, kBlack, Form ("%g < #it{p}_{T}^{Z} < %g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.036/uPadY);
    }
    else if (iCent == 3) {
      //myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "Z#rightarrowee", 1.4, 0.036/uPadY);
      //myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "Z#rightarrow#mu#mu", 1.4, 0.036/uPadY);
      //myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "No correction", 1.4, 0.036/uPadY);
      //myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "Eff. corrected", 1.4, 0.036/uPadY);
      //myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "Electrons", 1.4, 0.036/uPadY);
      //myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "Muons", 1.4, 0.04/uPadY);
      myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "Central values", 1.4, 0.04/uPadY);
      myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "Muons Up", 1.4, 0.04/uPadY);
      myMarkerTextNoLine (0.50, 0.74, colors[3], kOpenSquare, "Muons Down", 1.4, 0.04/uPadY);
      //myMarkerTextNoLine (0.50, 0.9, colors[2], kFullCircle, "#varepsilon_{Z} = 1", 1.4, 0.04/uPadY);
      //myMarkerTextNoLine (0.50, 0.82, colors[1], kOpenSquare, "#varepsilon_{Z} = #varepsilon_{trig} #times #varepsilon_{reco}", 1.4, 0.04/uPadY);
    }


    dPad->cd ();
    dPad->SetLogx ();

    hrat = (TH1D*) h2->Clone ("hrat");
    hrat->Reset ();
    //hrat->Divide (h1);
    for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
      const float y1 = h1->GetBinContent (ix);
      const float y1e = h1->GetBinError (ix);
      const float y2 = h2->GetBinContent (ix);
      const float y2e = h2->GetBinError (ix);
      hrat->SetBinContent (ix, y1/y2);
      hrat->SetBinError (ix, fabs(y1/y2)*sqrt (fabs (pow (y1e/y1,2) + pow (y2e/y2,2) - 2.*y1e*y1e/(y1*y2))));
    }
    //for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
    //  const float passes = h2->GetBinContent (ix);
    //  const float trials = h1->GetBinContent (ix);
    //  hrat->SetBinContent (ix, (passes+1) / (trials+2));
    //  hrat->SetBinError (ix, sqrt ((passes+1)*(passes+2)/((trials+2)*(trials+3)) - pow (passes+1, 2) / pow (trials+2, 2)));

    //  //hrat->SetBinError (ix, sqrt ((hrat->GetBinContent (ix)) * (1-hrat->GetBinContent (ix)) / h1->GetBinContent (ix)));
    //}

    //a1->LoadTrackingEfficiencies (true);
    //a2->LoadTrackingEfficiencies (true);
    //a1->LoadTrackingPurities (true);
    //a2->LoadTrackingPurities (true);

    //eff1 = (TH1D*) a1->h2_num_trk_effs[iCent]->ProjectionY ("1");
    //effrat = (TH1D*) a1->h2_den_trk_effs[iCent]->ProjectionY ("2");
    //eff1->Divide (effrat);
    //delete effrat;
    //eff2 = (TH1D*) a2->h2_num_trk_effs[iCent]->ProjectionY ("3");
    //effrat = (TH1D*) a2->h2_den_trk_effs[iCent]->ProjectionY ("4");
    //eff2->Divide (effrat);
    //delete effrat;

    //effrat = (TH1D*) eff2->Clone ("effrat");
    //effrat->Divide (eff1);
    //delete eff1, eff2;

    //eff1 = (TH1D*) a1->h2_num_trk_purs[iCent]->ProjectionY ("1");
    //eff1->Divide ((TH1D*) a1->h2_den_trk_purs[iCent]->ProjectionY ("2"));
    //eff2 = (TH1D*) a2->h2_num_trk_purs[iCent]->ProjectionY ("3");
    //eff2->Divide ((TH1D*) a2->h2_den_trk_purs[iCent]->ProjectionY ("4"));

    //purrat = (TH1D*) eff2->Clone ("purrat");
    //purrat->Divide (eff1);
    //delete eff1, eff2;

    //const int bin1 = effrat->FindBin (hrat->GetBinCenter (1));
    //const int bin2 = purrat->FindBin (hrat->GetBinCenter (1));
    //for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
    //  hrat->SetBinContent (ix, hrat->GetBinContent (ix) * purrat->GetBinContent (ix+bin2-1) / effrat->GetBinContent (ix+bin1-1));
    //  hrat->SetBinError (ix, hrat->GetBinError (ix) * purrat->GetBinContent (ix+bin2-1) / effrat->GetBinContent (ix+bin1-1));
    //}


    g = make_graph (hrat);
    delete hrat, effrat;

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[1]);
    g->SetLineColor (colors[1]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtchBins[0], allPtchBins[maxNPtchBins]);
    g->GetYaxis ()->SetRangeUser (0.95, 1.05);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
    //g->GetYaxis ()->SetTitle ("Z#rightarrow#mu#mu / Z#rightarrowee");
    g->GetYaxis ()->SetTitle ("Variation / Nominal");
    //g->GetYaxis ()->SetTitle ("Eff. Corrected / Uncorrected");
    //g->GetYaxis ()->SetTitle ("New ES / Old ES");
    //g->GetYaxis ()->SetTitle ("Unfolded / not unfolded");

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

    //TF1* fit1 = new TF1 ("fit1", "[0]", allPtchBins[0], allPtchBins[maxNPtchBins]);
    //TF1* fit1 = new TF1 ("fit1", "[0]+[1]*log(x)", allXHZBins[0], allXHZBins[maxNXHZBins]);
    TF1* fit1 = new TF1 ("fit1", "[0]+[1]*log(x)", allPtchBins[0], allPtchBins[maxNPtchBins]);
    fit1->SetParameter (0, 1);
    fit1->SetParameter (1, 0);
    g->Fit (fit1, "RQN0");

    fit1->SetLineColor (colors[1]);
    fit1->SetLineStyle (2);
    fit1->SetLineWidth (1);
    fit1->Draw ("same");

    //TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]", allPtchBins[0], allPtchBins[maxNPtchBins]);
    //TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]-[1]*log(x)", allXHZBins[0], allXHZBins[maxNXHZBins]);
    TF1* inv_fit1 = new TF1 ("inv_fit1", "[0]-[1]*log(x)", allPtchBins[0], allPtchBins[maxNPtchBins]);
    inv_fit1->SetParameter (0, 1/fit1->GetParameter (0));
    inv_fit1->SetParameter (1, fit1->GetParameter (1)/fit1->GetParameter (0));

    inv_fit1->SetLineColor (colors[1]);
    inv_fit1->SetLineStyle (2);
    inv_fit1->SetLineWidth (1);
    inv_fit1->Draw ("same");

    cout << "chi2/ndf = " << fit1->GetChisquare () << " / " << fit1->GetNDF () << endl;


    hrat = (TH1D*) h3->Clone ("hrat");
    hrat->Reset ();
    for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
      const float y1 = h1->GetBinContent (ix);
      const float y1e = h1->GetBinError (ix);
      const float y3 = h3->GetBinContent (ix);
      const float y3e = h3->GetBinError (ix);
      hrat->SetBinContent (ix, y1/y3);
      hrat->SetBinError (ix, (y1/y3)*sqrt (fabs (pow (y1e/y1,2) + pow (y3e/y3,2) - 2.*y1e*y1e/(y1*y3))));
    }

    //for (int ix = 1; ix <= hrat->GetNbinsX (); ix++) {
    //  hrat->SetBinError (ix, sqrt ((hrat->GetBinContent (ix)) * (1-hrat->GetBinContent (ix)) / h1->GetBinContent (ix)));
    //}

    //a1->LoadTrackingEfficiencies (true);
    //a2->LoadTrackingEfficiencies (true);
    //a1->LoadTrackingPurities (true);
    //a2->LoadTrackingPurities (true);

    //eff1 = (TH1D*) a1->h2_num_trk_effs[iCent]->ProjectionY ("1");
    //effrat = (TH1D*) a1->h2_den_trk_effs[iCent]->ProjectionY ("2");
    //eff1->Divide (effrat);
    //delete effrat;
    //eff2 = (TH1D*) a2->h2_num_trk_effs[iCent]->ProjectionY ("3");
    //effrat = (TH1D*) a2->h2_den_trk_effs[iCent]->ProjectionY ("4");
    //eff2->Divide (effrat);
    //delete effrat;

    //effrat = (TH1D*) eff2->Clone ("effrat");
    //effrat->Divide (eff1);
    //delete eff1, eff2;

    //eff1 = (TH1D*) a1->h2_num_trk_purs[iCent]->ProjectionY ("1");
    //eff1->Divide ((TH1D*) a1->h2_den_trk_purs[iCent]->ProjectionY ("2"));
    //eff2 = (TH1D*) a2->h2_num_trk_purs[iCent]->ProjectionY ("3");
    //eff2->Divide ((TH1D*) a2->h2_den_trk_purs[iCent]->ProjectionY ("4"));

    //purrat = (TH1D*) eff2->Clone ("purrat");
    //purrat->Divide (eff1);
    //delete eff1, eff2;

    //const int bin1 = effrat->FindBin (hrat->GetBinCenter (1));
    //const int bin2 = purrat->FindBin (hrat->GetBinCenter (1));
    //for (int iy = 1; iy <= hrat->GetNbinsX (); iy++)
    //  hrat->SetBinContent (iy, hrat->GetBinContent (iy) * purrat->GetBinContent (iy+bin2-1));// / effrat->GetBinContent (iy+bin1-1));

    g = make_graph (hrat);
    delete hrat;//, effrat;

    g->SetMarkerStyle (kOpenSquare);
    g->SetMarkerSize (1);
    g->SetLineWidth (1);
    g->SetMarkerColor (colors[3]);
    g->SetLineColor (colors[3]);

    //g->GetXaxis ()->SetLimits (allXHZBins[0], allXHZBins[maxNXHZBins]);
    g->GetXaxis ()->SetLimits (allPtchBins[0], allPtchBins[maxNPtchBins]);
    g->GetYaxis ()->SetRangeUser (0.95, 1.05);

    g->GetXaxis ()->SetMoreLogLabels ();

    g->Draw ("P");

    //TF1* fit2 = new TF1 ("fit2", "[0]+[1]*log(x)", allXHZBins[0], allXHZBins[maxNXHZBins]);
    TF1* fit2 = new TF1 ("fit2", "[0]+[1]*log(x)", allPtchBins[0], allPtchBins[maxNPtchBins]);
    fit2->SetParameter (0, 1);
    fit2->SetParameter (1, 0);
    g->Fit (fit2, "RQN0");

    fit2->SetLineColor (colors[3]);
    fit2->SetLineStyle (2);
    fit2->SetLineWidth (1);
    fit2->Draw ("same");

    cout << "chi2/ndf = " << fit2->GetChisquare () << " / " << fit2->GetNDF () << endl;

    //TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]-[1]*log(x)", allXHZBins[0], allXHZBins[maxNPtchBins]);
    TF1* inv_fit2 = new TF1 ("inv_fit2", "[0]-[1]*log(x)", allPtchBins[0], allPtchBins[maxNPtchBins]);
    inv_fit2->SetParameter (0, 1./fit2->GetParameter (0));
    inv_fit2->SetParameter (1, fit2->GetParameter (1)/fit2->GetParameter (0));

    inv_fit2->SetLineColor (colors[3]);
    inv_fit2->SetLineStyle (2);
    inv_fit2->SetLineWidth (1);
    inv_fit2->Draw ("same");

    TLine* l = new TLine (allPtchBins[0], 1, allPtchBins[maxNPtchBins], 1);
    l->SetLineColor (kBlack);
    l->SetLineStyle (2);
    l->Draw ("same");

    //myText (0.22, 0.3, kBlack, Form ("Avg. = %s", FormatMeasurement (fit1->GetParameter (0), fit1->GetParError (0), 2)), 0.03/dPadY);

    delete h1, h2;
  }

}
*/

#endif
