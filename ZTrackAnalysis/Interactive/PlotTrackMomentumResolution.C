#ifndef __PlotTrackMomentumResolution_C__
#define __PlotTrackMomentumResolution_C__

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLine.h>
#include <TCanvas.h>

#include <AtlasUtils.h>
#include <Utilities.h>

#include "Params.h"

using namespace atlashi;

typedef TGraphAsymmErrors TGAE;

const int numFinerEtachBins = 40;
const double* finerEtachBins = linspace (-2.5, 2.5, numFinerEtachBins);

const double etachBins[6] = {0, 0.5, 1.0, 1.5, 2.0, 2.5};
const int numEtachBins = sizeof (etachBins) / sizeof (etachBins[0]) - 1;

const double pTchBinsPP[30] = {0.5, 0.7, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 60};
const int numPtchBinsPP = sizeof (pTchBinsPP) / sizeof (pTchBinsPP[0]) - 1;
const double pTchBinsPbPb[28] = {0.5, 0.7, 1, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10, 15, 20, 30, 40, 60};
const int numPtchBinsPbPb = sizeof (pTchBinsPbPb) / sizeof (pTchBinsPbPb[0]) - 1;

const bool isPbPb = false;

TH1D**** h_tms = nullptr;
TH2D** h2_avg_tms = nullptr;
TH2D** h2_avg_tmr = nullptr;
TH1D*** h_avg_tms = nullptr;
TH1D*** h_avg_tmr = nullptr;


void PlotTrackMomentumResolution () {

  TFile* inFile = new TFile ("../rootFiles/TrackingMomentumResolution/Nominal/trackingMomentumResolution.root", "read");


  h_tms = new TH1D***[numCentBins];
  h2_avg_tms = new TH2D*[numCentBins];
  h2_avg_tmr = new TH2D*[numCentBins];
  h_avg_tms = new TH1D**[numCentBins];
  h_avg_tmr = new TH1D**[numCentBins];

  for (int iCent = 0; iCent < 1; iCent++) {
    h2_avg_tms[iCent] = new TH2D (Form ("h2_avg_tms_iCent%i", iCent), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};TMS [%]", (isPbPb ? numPtchBinsPP : numPtchBinsPbPb), (isPbPb ? pTchBinsPP : pTchBinsPbPb), numFinerEtachBins, finerEtachBins);
    h2_avg_tmr[iCent] = new TH2D (Form ("h2_avg_tmr_iCent%i", iCent), ";#it{p}_{T}^{truth} [GeV];#eta_{truth};TMR [%]", (isPbPb ? numPtchBinsPP : numPtchBinsPbPb), (isPbPb ? pTchBinsPP : pTchBinsPbPb), numFinerEtachBins, finerEtachBins);

    h_avg_tms[iCent] = new TH1D*[numEtachBins+1];
    h_avg_tmr[iCent] = new TH1D*[numEtachBins+1];
    for (int iEta = 0; iEta <= numEtachBins; iEta++) {
      h_avg_tms[iCent][iEta] = new TH1D (Form ("h_avg_tms_iCent%i_iEta%i", iCent, iEta), ";#it{p}_{T}^{truth} [GeV];TMS [%]", (isPbPb ? numPtchBinsPP : numPtchBinsPbPb), (isPbPb ? pTchBinsPP : pTchBinsPbPb));
      h_avg_tmr[iCent][iEta] = new TH1D (Form ("h_avg_tmr_iCent%i_iEta%i", iCent, iEta), ";#it{p}_{T}^{truth} [GeV];TMR [%]", (isPbPb ? numPtchBinsPP : numPtchBinsPbPb), (isPbPb ? pTchBinsPP : pTchBinsPbPb));
    } // end loop over iEta


    h_tms[iCent] = new TH1D**[(isPbPb ? numPtchBinsPP : numPtchBinsPbPb)];
    for (int iPtch = 0; iPtch < (isPbPb ? numPtchBinsPP : numPtchBinsPbPb); iPtch++) {
      h_tms[iCent][iPtch] = new TH1D*[numFinerEtachBins];
      for (int iEta = 0; iEta < numFinerEtachBins; iEta++) {
        h_tms[iCent][iPtch][iEta] = (TH1D*) inFile->Get (Form ("h_tms_iCent%i_iPtch%i_iEta%i", iCent, iPtch, iEta));
      } // end loop over iEta
    } // end loop over iPtch

  } // end loop over iCent



  for (int iCent = 0; iCent < 1; iCent++) {

    for (int iPtch = 0; iPtch < (isPbPb ? numPtchBinsPP : numPtchBinsPbPb); iPtch++) {

      TH1D** h_tms_integratedEta = new TH1D*[numEtachBins+1];
      for (int iEta = 0; iEta <= numEtachBins; iEta++) {
        h_tms_integratedEta[iEta] = new TH1D (Form ("h_tms_integratedEta_iEta%i", iEta), "#it{p}_{T}^{reco} / #it{p}_{T}^{truth}", 200, 0.5, 1.5);
        h_tms_integratedEta[iEta]->Sumw2 ();
      }

      for (int iFinerEta = 0; iFinerEta < numFinerEtachBins; iFinerEta++) {

        // first add to eta-integrated plot
        const float binCenter = 0.5 * fabs (finerEtachBins[iFinerEta] + finerEtachBins[iFinerEta+1]);
        int iEta = 0;
        while (iEta < numEtachBins) {
          if (etachBins[iEta] < binCenter && binCenter < etachBins[iEta+1])
            break;
          iEta++;
        }

        h_tms_integratedEta[iEta]->Add (h_tms[iCent][iPtch][iFinerEta]);
        h_tms_integratedEta[numEtachBins]->Add (h_tms[iCent][iPtch][iFinerEta]);

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.5, 1.5);
        fit->SetParameter (0, h_tms[iCent][iPtch][iFinerEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_tms[iCent][iPtch][iFinerEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", mean-4*sigma, mean+4*sigma);

        fit->SetParameter (0, h_tms[iCent][iPtch][iFinerEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_tms[iCent][iPtch][iFinerEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double tms = mean;
        const double tms_err = mean_err;
        
        const double tmr = sigma / mean;
        const double tmr_err = fabs (tmr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h2_avg_tms[iCent]->SetBinContent (iPtch+1, iFinerEta+1, tms * 100);
        h2_avg_tms[iCent]->SetBinError (iPtch+1, iFinerEta+1, tms_err * 100);
        h2_avg_tmr[iCent]->SetBinContent (iPtch+1, iFinerEta+1, tmr * 100);
        h2_avg_tmr[iCent]->SetBinError (iPtch+1, iFinerEta+1, tmr_err * 100);
      }

      for (int iEta = 0; iEta <= numEtachBins; iEta++) {

        TF1* fit = new TF1 ("fit", "gaus(0)", 0.5, 1.5);
        fit->SetParameter (0, h_tms_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, 1);
        fit->SetParameter (2, 1);

        h_tms_integratedEta[iEta]->Fit (fit, "RN0Q");

        double mean = fit->GetParameter (1);
        double sigma = fit->GetParameter (2);

        delete fit;
        fit = new TF1 ("fit", "gaus(0)", mean-4*sigma, mean+4*sigma);

        fit->SetParameter (0, h_tms_integratedEta[iEta]->Integral ());
        fit->SetParameter (1, mean);
        fit->SetParameter (2, sigma);

        h_tms_integratedEta[iEta]->Fit (fit, "RN0Q");

        mean = fit->GetParameter (1);
        sigma = fit->GetParameter (2);
        double mean_err = fit->GetParError (1);
        double sigma_err = fit->GetParError (2);

        delete fit;

        const double tms = mean;
        const double tms_err = mean_err;

        const double tmr = sigma / mean;
        const double tmr_err = fabs (tmr) * sqrt (pow (mean_err/mean, 2) + pow (sigma_err/sigma, 2));

        h_avg_tms[iCent][iEta]->SetBinContent (iPtch+1, tms * 100);
        h_avg_tms[iCent][iEta]->SetBinError (iPtch+1, tms_err * 100);
        h_avg_tmr[iCent][iEta]->SetBinContent (iPtch+1, tmr * 100);
        h_avg_tmr[iCent][iEta]->SetBinError (iPtch+1, tmr_err * 100);

        delete h_tms_integratedEta[iEta];
      } // end loop over iEta
  
      delete[] h_tms_integratedEta;
    } // end loop over iPtch
  } // end loop over iCent


  TLatex* tl = new TLatex ();
  TLine* l = new TLine ();


  {
    TCanvas* c1 = new TCanvas ("c1", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c1->SetLeftMargin (lMargin);
    c1->SetRightMargin (rMargin);
    c1->SetBottomMargin (bMargin);
    c1->SetTopMargin (tMargin);

    {
      gPad->SetLogx ();

      TH1D* htemp = new TH1D ("htemp", "", 1, 0.5, 60);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("TMS [%]");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 95;
      const double ymax = 105;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");

      l->SetLineStyle (2);
      l->SetLineWidth (2);
      //l->SetLineColor (kPink-8);
      l->DrawLine ((isPbPb ? pTchBinsPP : pTchBinsPbPb)[0], 100, (isPbPb ? pTchBinsPP : pTchBinsPbPb)[(isPbPb ? numPtchBinsPP : numPtchBinsPbPb)], 100);
    }
     

    for (int iEta = 0; iEta < numEtachBins; iEta++) {
      TGAE* g = make_graph (h_avg_tms[0][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_avg_tms[0][numEtachBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.30, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV");
    tl->DrawLatexNDC (0.30, 0.800, "Powheg + Pythia #it{Z} #rightarrow #it{ll}");
    for (int iEta = 0; iEta < numEtachBins; iEta++)
      myMarkerTextNoLine (0.3, 0.39-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.04);
    myMarkerTextNoLine (0.3, 0.44, colors[0], kFullCircle, "0 < |#eta| < 2.5", 1.2, 0.04);

    c1->SaveAs ("../Plots/TrackMomentumScale.pdf");
  }



  {
    TCanvas* c2 = new TCanvas ("c2", "", 800, 800);
    const double lMargin = 0.15;
    const double rMargin = 0.04;
    const double bMargin = 0.15;
    const double tMargin = 0.04;

    c2->SetLeftMargin (lMargin);
    c2->SetRightMargin (rMargin);
    c2->SetBottomMargin (bMargin);
    c2->SetTopMargin (tMargin);

    {
      gPad->SetLogx ();

      TH1D* htemp = new TH1D ("htemp", "", 1, 0.5, 60);
  
      TAxis* xax = htemp->GetXaxis ();
      TAxis* yax = htemp->GetYaxis ();

      xax->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      xax->SetTitleOffset (0.9 * xax->GetTitleOffset ());
      xax->SetLabelSize (0);

      yax->SetTitle ("TMR [%]");
      yax->SetLabelOffset (1.8 * yax->GetLabelOffset ());
      const double ymin = 0;
      const double ymax = 20;
      yax->SetRangeUser (ymin, ymax);

      htemp->SetLineWidth (0);

      htemp->DrawCopy ("");
      SaferDelete (&htemp);

      tl->SetTextFont (43);
      tl->SetTextSize (36);
      tl->SetTextAlign (21);
      
      const double yoff = ymin - 0.05 * (ymax-ymin) / (1.-tMargin-bMargin);
      tl->DrawLatex (1,  yoff, "1");
      tl->DrawLatex (2,  yoff, "2");
      tl->DrawLatex (3,  yoff, "3");
      tl->DrawLatex (4,  yoff, "4");
      tl->DrawLatex (5,  yoff, "5");
      tl->DrawLatex (6,  yoff, "6");
      tl->DrawLatex (7,  yoff, "7");
      tl->DrawLatex (10, yoff, "10");
      tl->DrawLatex (20, yoff, "20");
      tl->DrawLatex (30, yoff, "30");
      tl->DrawLatex (40, yoff, "40");
      tl->DrawLatex (60, yoff, "60");
    }
     

    for (int iEta = 0; iEta < numEtachBins; iEta++) {
      TGAE* g = make_graph (h_avg_tmr[0][iEta]);

      g->SetMarkerColor (colors[iEta+1]);
      g->SetLineColor (colors[iEta+1]);
      g->SetMarkerStyle (kOpenCircle);
     
      g->Draw ("P"); 
    }
    {
      TGAE* g = make_graph (h_avg_tmr[0][numEtachBins]);

      g->SetMarkerColor (colors[0]);
      g->SetLineColor (colors[0]);
      g->SetMarkerStyle (kFullCircle);
     
      g->Draw ("P"); 
    }

    tl->SetTextColor (kBlack);
    tl->SetTextAlign (11);

    tl->SetTextSize (32);
    tl->DrawLatexNDC (0.30, 0.890, "#bf{#it{ATLAS}} Simulation Internal");
    tl->SetTextSize (30);
    tl->DrawLatexNDC (0.30, 0.845, "#it{pp}, #sqrt{s} = 5.02 TeV");
    tl->DrawLatexNDC (0.30, 0.800, "Powheg + Pythia #it{Z} #rightarrow #it{ll}");

    for (int iEta = 0; iEta < numEtachBins; iEta++)
      myMarkerTextNoLine (0.5, 0.69-0.04*iEta, colors[iEta+1], kOpenCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.04);
    myMarkerTextNoLine (0.5, 0.74, colors[0], kFullCircle, "All #eta", 1.2, 0.04);

    c2->SaveAs ("../Plots/TrackMomentumResolution.pdf");
  }


}

#endif
