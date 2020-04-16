#ifndef __PlotJewel_C__
#define __PlotJewel_C__

#include "../Params.h"

#include "ArrayTemplates.h"

#include <AtlasUtils.h>

#include <TFile.h>
#include <TH2.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TColor.h>
#include <TLatex.h>

#include <iostream>
#include <string>

const Color_t jewelColor  = tcolor->GetColor (255, 170, 50);
const double  jewelAlpha  = 0.5;

void PlotJewel () {

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/Jewel/hists.root", rootPath.Data ()), "read");

  TGAE** g_jewel_xhz = Get1DArray <TGAE*> (nPtZBins);
  TGAE** g_jewel_pth = Get1DArray <TGAE*> (nPtZBins);

  TFile* jewelFile = new TFile ("../../rootFiles/Jewel/hists.root", "read");
  for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h_num = (TH1D*) jewelFile->Get (Form ("h_z_trk_pth_medium_iSpc2_iPtZ%i", iPtZ));
    TH1D* h_den = (TH1D*) jewelFile->Get (Form ("h_z_trk_pth_vacuum_iSpc2_iPtZ%i", iPtZ));
    
    h_num->Divide (h_den);
    g_jewel_pth[iPtZ] = make_graph (h_num);
    {
      TGAE* g = g_jewel_pth[iPtZ];
      double x, y;
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (g->GetErrorYhigh (i) < 0.05) g->SetPointEYhigh (i, 0.05); 
        if (g->GetErrorYlow (i) < 0.05) g->SetPointEYlow (i, 0.05); 
      }
    }
    
    h_num = (TH1D*) jewelFile->Get (Form ("h_z_trk_xhz_medium_iSpc2_iPtZ%i", iPtZ));
    h_den = (TH1D*) jewelFile->Get (Form ("h_z_trk_xhz_vacuum_iSpc2_iPtZ%i", iPtZ));
    
    h_num->Divide (h_den);
    g_jewel_xhz[iPtZ] = make_graph (h_num);

    {
      TGAE* g = g_jewel_xhz[iPtZ];
      double x, y;
      for (int i = 0; i < g->GetN (); i++) {
        g->GetPoint (i, x, y);
        if (g->GetErrorYhigh (i) < 0.05) g->SetPointEYhigh (i, 0.05); 
        if (g->GetErrorYlow (i) < 0.05) g->SetPointEYlow (i, 0.05); 
      }
    }
    
  }
  jewelFile->Close ();

  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();

  TCanvas* c = new TCanvas ("c", "", 1200, 800);
  c->Divide (3, 2);

  for (int iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
    c->cd (iPtZ-2+1);
    gPad->SetLogx ();

    {
      TH1D* h = new TH1D ("", "", nPtchBins[iPtZ], pTchBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{ch} [GeV]");
      xax->SetRangeUser (pTchBins[iPtZ][0], pTchBins[iPtZ][nPtchBins[iPtZ]]);
      //xax->SetTitleFont (43);
      //xax->SetTitleSize (36);
      //xax->SetTitleOffset (2.5);
      //xax->SetTitleOffset (1.2);
      //xax->SetLabelSize (0);
      xax->SetMoreLogLabels ();

      yax->SetTitle ("#it{I}_{AA} (#it{p}_{T}^{ch})");
      //yax->SetRangeUser (0.05, 10);
      double ymin, ymax;
      if (iPtZ == 3) { ymin = 0; ymax = 4.0; }
      if (iPtZ == 4) { ymin = 0; ymax = 5.0; }
      else           { ymin = 0; ymax = 4.0; }
      yax->SetRangeUser (ymin, ymax);
      //yax->SetTitleFont (43);
      //yax->SetTitleSize (36);
      //yax->SetTitleOffset (2.6);
      //yax->SetLabelFont (43);
      //yax->SetLabelSize (36);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      //const double yoff = ymin - 0.05 * (ymax-ymin);
      //tl->DrawLatex (1,  yoff, "1");
      //tl->DrawLatex (2,  yoff, "2");
      //tl->DrawLatex (3,  yoff, "3");
      //tl->DrawLatex (4,  yoff, "4");
      //tl->DrawLatex (5,  yoff, "5");
      //tl->DrawLatex (6,  yoff, "6");
      //tl->DrawLatex (7,  yoff, "7");
      //tl->DrawLatex (10, yoff, "10");
      //tl->DrawLatex (20, yoff, "20");
      //tl->DrawLatex (30, yoff, "30");
      //tl->DrawLatex (40, yoff, "40");
      //tl->DrawLatex (60, yoff, "60");

      if (iPtZ < nPtZBins-1) myText (0.60, 0.84, kBlack, Form ("#it{p}_{T}^{Z} = %g-%g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
      else                   myText (0.60, 0.84, kBlack, Form ("#it{p}_{T}^{Z} = %g+ GeV", zPtBins[iPtZ]), 0.05);

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (pTchBins[iPtZ][0], 1, pTchBins[iPtZ][nPtchBins[iPtZ]], 1);
    }

    TGAE* g = g_jewel_pth[iPtZ];

    g->SetFillColorAlpha (jewelColor, jewelAlpha);
    g->Draw ("3");

    c->cd (iPtZ-2+4);
    gPad->SetLogx ();

    {
      TH1D* h = new TH1D ("", "", nXhZBins[iPtZ], xhZBins[iPtZ]);

      TAxis* xax = h->GetXaxis ();
      TAxis* yax = h->GetYaxis ();

      xax->SetTitle ("#it{x}_{hZ}");
      xax->SetRangeUser (xhZBins[iPtZ][0], 1.);
      //xax->SetTitleFont (43);
      //xax->SetTitleSize (36);
      //xax->SetTitleOffset (2.5);
      //xax->SetTitleOffset (1.0);
      //xax->SetLabelSize (0);
      xax->SetMoreLogLabels ();

      yax->SetTitle ("#it{I}_{AA} (#it{x}_{hZ})");
      //yax->SetRangeUser (0.05, 10);
      double ymin, ymax;
      if (iPtZ == 3) { ymin = -0.1; ymax = 4.0; }
      if (iPtZ == 4) { ymin = -0.1; ymax = 5.0; }
      else           { ymin = 0;    ymax = 4.0; }
      yax->SetRangeUser (ymin, ymax);
      //yax->SetTitleFont (43);
      //yax->SetTitleSize (36);
      //yax->SetTitleOffset (2.6);
      //yax->SetLabelFont (43);
      //yax->SetLabelSize (36);

      h->SetLineWidth (0);

      h->DrawCopy ("");
      SaferDelete (&h);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
      //const double yoff = 0.034;

      //const double yoff = ymin - 0.05 * (ymax-ymin);
      //if (iPtZ > 2) {
      //  tl->DrawLatex (5e-2,  yoff, "5#times10^{-2}");
      //  if (iPtZ > 3) tl->DrawLatex (2e-2,  yoff, "2#times10^{-2}");
      //}
      //tl->DrawLatex (1e-1,  yoff, "10^{-1}");
      //tl->DrawLatex (2e-1,  yoff, "2#times10^{-1}");
      //tl->DrawLatex (5e-1,  yoff, "5#times10^{-1}");
      //tl->DrawLatex (1,     yoff, "1");

      if (iPtZ < nPtZBins-1) myText (0.60, 0.84, kBlack, Form ("#it{p}_{T}^{Z} = %g-%g GeV", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.05);
      else                   myText (0.60, 0.84, kBlack, Form ("#it{p}_{T}^{Z} = %g+ GeV", zPtBins[iPtZ]), 0.05);

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine (xhZBins[iPtZ][0], 1, xhZBins[iPtZ][nXhZBins[iPtZ]], 1);
    }
    g = g_jewel_xhz[iPtZ];

    g->SetFillColorAlpha (jewelColor, jewelAlpha);
    g->Draw ("3");
  }
  c->SaveAs ("jewel_allIAAs.pdf");

}


#endif
