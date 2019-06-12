#ifndef __Systematics_h__
#define __Systematics_h__

#include "Params.h"
#include "Analysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class Systematics : public Analysis {

  private:
  vector<Analysis*> variations;

  void NullifyErrors ();

  public:
  Systematics (Analysis* nom) : Analysis (){
    name = "systematics";
    directory = "Systematics/";
    SetupDirectories (directory, "ZTrackAnalysis/");

    CopyAnalysis (nom, true);
  }

  void AddVariation (Analysis* a);
  void CombineErrors ();
  void AddSystematic (Systematics* s);

};


////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: AddVariation (Analysis* a) {
  variations.push_back (a);
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Sets the errors in this systematic as a combination of all added variations.
// Takes the maximum error for each point.
// Intended for combining up & down variations, but expandable for additional categories
// (e.g. track quality criteria)
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: CombineErrors () {

  NullifyErrors ();

  for (Analysis* a : variations) {

    a->SubtractBackground ();
    a->CalculateIAA ();
    a->CalculateICP ();

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

        // Hadron yield systematics, signal & signal+bkg levels
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 0; iCent < numCentBins; iCent++) {
            for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
              TH1D* sys = h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              TH1D* var = a->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              if (sys && var) CalcSystematics (sys, var);

              
              sys = h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              var = a->h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              if (sys && var) CalcSystematics (sys, var);

              sys = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              var = a->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent];
              if (sys && var) CalcSystematics (sys, var);
            } // end loop over xztrk bins
          } // end loop over cents
        } // end loop over phi

        // IAA, ICP systematics
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          for (short iCent = 1; iCent < numCentBins; iCent++) {
            TH1D* sys = h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];
            TH1D* var = a->h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
          } // end loop over cents
          for (short iCent = 2; iCent < numCentBins; iCent++) {
            TH1D* sys = h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];
            TH1D* var = a->h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];
            if (sys && var) CalcSystematics (sys, var);
          } // end loop over cents
        } // end loop over phi

      } // end loop over pT^Z bins
    } // end loop over species

  } // end loop over variations
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Addition of independent systematics in quadrature; adds the errors of s to this systematic.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: NullifyErrors () {

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

      // Hadron yield systematics, signal & signal+bkg levels
      for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
            TH1D* master = h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master) ResetHistErrors (master);

            master = h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master) ResetHistErrors (master);

            master = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master) ResetHistErrors (master);
          } // end loop over xztrk bins
        } // end loop over cents
      } // end loop over phi

      // IAA, ICP systematics
      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];
          if (master) ResetHistErrors (master);
        } // end loop over cents
        for (short iCent = 2; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];
          if (master) ResetHistErrors (master);
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species

}




////////////////////////////////////////////////////////////////////////////////////////////////
// Addition of independent systematics in quadrature; adds the errors of s to this systematic.
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: AddSystematic (Systematics* s) {

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) { 

      // Hadron yield systematics, signal & signal+bkg levels
      for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 0; iCent < numCentBins; iCent++) {
          for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
            TH1D* master = h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            TH1D* sys = s->h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master && sys) AddSystematics (master, sys);

            master = h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            sys = s->h_z_trk_pt_sub[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master && sys) AddSystematics (master, sys);

            master = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            sys = s->h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            if (master && sys) AddSystematics (master, sys);
          } // end loop over xztrk bins
        } // end loop over cents
      } // end loop over phi

      // IAA, ICP systematics
      for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
        for (short iCent = 1; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];
          TH1D* sys = s->h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];
          if (master && sys) AddSystematics (master, sys);
        } // end loop over cents
        for (short iCent = 2; iCent < numCentBins; iCent++) {
          TH1D* master = h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];
          TH1D* sys = s->h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];
          if (master && sys) AddSystematics (master, sys);
        } // end loop over cents
      } // end loop over phi
    } // end loop over pT^Z bins
  } // end loop over species

}



/*
////////////////////////////////////////////////////////////////////////////////////////////////
// Plot pTtrk binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: PlotTrkYield (const short pSpc, const short pPtZ) {
  if (!backgroundSubtracted)
    SubtractBackground ();

  const double padRatio = 0.9; // ratio of size of upper pad to middle & lower pads. Used to scale plots and font sizes equally.
  const double dPadY = padRatio / (2*padRatio+1.0);
  const double mPadY = padRatio / (2*padRatio+1.0);
  const double uPadY = 1.0 - mPadY - dPadY;
  const int axisTextSize = 23;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
      continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("TrkYieldCanvas_%s_iPtZ%i", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 400*numCentBins, 1000);
        gDirectory->Add (c);
      }

      for (short iCent = 0; iCent < numCentBins; iCent++) {
        c->cd ();

        const char* topPadName = Form ("p_top_%s_iPtZ%i_iCent%i", spc, iPtZ, iCent);
        const char* middlePadName = Form ("p_middle_%s_iPtZ%i_iCent%i", spc, iPtZ, iCent);
        const char* bottomPadName = Form ("p_bottom_%s_iPtZ%i_iCent%i", spc, iPtZ, iCent);

        TPad* topPad = nullptr, *middlePad = nullptr, *bottomPad = nullptr;
        if (!canvasExists) {
          topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY+mPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
          middlePad = new TPad (middlePadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, dPadY+mPadY);
          bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);

          gDirectory->Add (topPad);
          gDirectory->Add (middlePad);
          gDirectory->Add (bottomPad);

          topPad->SetTopMargin (0.04);
          topPad->SetBottomMargin (0);
          topPad->SetLeftMargin (0.17);
          topPad->SetRightMargin (0.06);
          middlePad->SetTopMargin (0);
          middlePad->SetBottomMargin (0);
          middlePad->SetLeftMargin (0.17);
          middlePad->SetRightMargin (0.06);
          bottomPad->SetTopMargin (0);
          bottomPad->SetBottomMargin (0.20);
          bottomPad->SetLeftMargin (0.17);
          bottomPad->SetRightMargin (0.06);
          topPad->Draw ();
          middlePad->Draw ();
          bottomPad->Draw ();
        }
        else {
          topPad = dynamic_cast<TPad*> (gDirectory->Get (topPadName));
          middlePad = dynamic_cast<TPad*> (gDirectory->Get (middlePadName));
          bottomPad = dynamic_cast<TPad*> (gDirectory->Get (bottomPadName));
        }

        topPad->cd ();
        GetDrawnObjects ();
        bool plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
        gPad->SetLogx ();
        gPad->SetLogy ();

        double min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          TH1D* h = h_z_trk_pt[iSpc][iPtZ][0][iPhi][iCent];
          if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
          if (h->GetMaximum () > max)  max = h->GetMaximum ();
        } // end loop over phi
        min = (min > 0 ? 0.5*min : 0.1);
        max = (max > 0 ? 2*max : 1);
        SetMinAndMax (min, max);

        for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
          TH1D* h = h_z_trk_pt[iSpc][iPtZ][0][iPhi][iCent];

          h->SetFillColorAlpha (fillColors[iPhi], 0.8);
          h->SetMarkerSize (0);
          h->SetLineColor (kBlack);
          //h->SetLineWidth (0);

          h->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
          h->GetYaxis ()->SetTitle ("Per-trigger yield Y (#it{p}_{T})");

          h->GetXaxis ()->SetTitleFont (43);
          h->GetXaxis ()->SetTitleSize (axisTextSize);
          h->GetXaxis ()->SetLabelFont (43);
          h->GetXaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleFont (43);
          h->GetYaxis ()->SetTitleSize (axisTextSize);
          h->GetYaxis ()->SetLabelFont (43);
          h->GetYaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

          h->Draw (plotNewAxes && iPhi == numPhiBins-1 ? "e2" : "e2 same");
        }

        if (!canvasExists)
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++)
            LabelTrkYield (iCent, iPhi);

        if (!plotSignal)
          continue;

        middlePad->cd ();
        GetDrawnObjects ();
        plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
        gPad->SetLogx ();
        gPad->SetLogy ();

        min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          TH1D* h = h_z_trk_pt_sub[iSpc][iPtZ][0][iPhi][iCent];
          if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
          if (h->GetMaximum () > max) max = h->GetMaximum ();
        } // end loop over phi
        float delta = log10 (max) - log10 (min);
        min = pow (10, log10 (min) - 0.1*delta);
        max = pow (10, log10 (max) + 0.1*delta);
        SetMinAndMax (min, max);

        for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
          TH1D* h = h_z_trk_pt_sub[iSpc][iPtZ][0][iPhi][iCent];

          h->SetFillColorAlpha (fillColors[iPhi], 0.8);
          h->SetLineColor (kBlack);
          h->SetMarkerSize (0);
          //h->SetLineWidth (0);

          h->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
          h->GetYaxis ()->SetTitle ("Signal Yield");

          h->GetXaxis ()->SetTitleFont (43);
          h->GetXaxis ()->SetTitleSize (axisTextSize);
          h->GetXaxis ()->SetLabelFont (43);
          h->GetXaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleFont (43);
          h->GetYaxis ()->SetTitleSize (axisTextSize);
          h->GetYaxis ()->SetLabelFont (43);
          h->GetYaxis ()->SetLabelSize (axisTextSize);

          h->GetXaxis ()->SetTitleOffset (1.5 * h->GetXaxis ()->GetTitleOffset ());
          h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

          h->Draw (plotNewAxes && iPhi == numPhiBins-1 ? "e2" : "e2 same");
        } // end loop over phi
        

        bottomPad->cd ();
        GetDrawnObjects ();
        plotNewAxes = (drawnHists.size () == 0 && drawnGraphs.size () == 0);
        gPad->SetLogx ();
        gPad->SetLogy ();

        min = 1e30, max = 0;
        GetMinAndMax (min, max, true);
        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          TH1D* h = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][0][iPhi][iCent];
          if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
          if (h->GetMaximum () > max) max = h->GetMaximum ();
        } // end loop over phi
        //delta = max - min;
        //min = min - 0.3*delta;
        //max = max + 0.3*delta;
        delta = log10 (max) - log10 (min);
        min = pow (10, log10 (min) - 0.1*delta);
        max = pow (10, log10 (max) + 0.1*delta);
        SetMinAndMax (min, max);

        for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
          TH1D* h = h_z_trk_pt_sig_to_bkg[iSpc][iPtZ][0][iPhi][iCent];

          h->SetFillColorAlpha (fillColors[iPhi], 0.8);
          h->SetLineColor (kBlack);
          h->SetMarkerSize (0);
          //h->SetLineWidth (0);

          h->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
          h->GetYaxis ()->SetRangeUser (min, max);

          h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
          h->GetYaxis ()->SetTitle ("Signal / Bkg.");

          h->GetXaxis ()->SetTitleFont (43);
          h->GetXaxis ()->SetTitleSize (axisTextSize);
          h->GetXaxis ()->SetLabelFont (43);
          h->GetXaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleFont (43);
          h->GetYaxis ()->SetTitleSize (axisTextSize);
          h->GetYaxis ()->SetLabelFont (43);
          h->GetYaxis ()->SetLabelSize (axisTextSize);

          h->GetXaxis ()->SetTitleOffset (1.5 * h->GetXaxis ()->GetTitleOffset ());
          h->GetYaxis ()->SetTitleOffset (2.2 * h->GetYaxis ()->GetTitleOffset ());

          h->Draw (plotNewAxes && iPhi == numPhiBins-1 ? "e2" : "e2 same");
        } // end loop over phi
        
      } // end loop over cents
      
      c->SaveAs (Form ("%s/pTTrk_dists_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: PlotIAARatios (const short pSpc, const short pPtZ) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!iaaCalculated)
    CalculateIAA ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_z_trk_pt_iaa_%s_iPtZ%i", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numCentBins-1), 500);
        c->Divide (numCentBins-1, 1);
        gDirectory->Add (c);
      }

      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);
        gPad->SetLogx ();

        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          TH1D* h = h_z_trk_pt_iaa[iSpc][iPtZ][0][iPhi][iCent];

          h->SetFillColorAlpha (fillColors[iPhi], 0.5);
          h->SetMarkerSize (0);
          h->SetLineColor (kBlack);
          h->SetLineWidth (0);

          h->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
          h->GetYaxis ()->SetRangeUser (0, 1.6);

          h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
          h->GetYaxis ()->SetTitle ("I_{AA}");

          h->GetXaxis ()->SetTitleFont (43);
          h->GetXaxis ()->SetTitleSize (axisTextSize);
          h->GetXaxis ()->SetLabelFont (43);
          h->GetXaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleFont (43);
          h->GetYaxis ()->SetTitleSize (axisTextSize);
          h->GetYaxis ()->SetLabelFont (43);
          h->GetYaxis ()->SetLabelSize (axisTextSize);

          h->GetXaxis ()->SetTitleOffset (0.8 * h->GetXaxis ()->GetTitleOffset ());
          h->GetYaxis ()->SetTitleOffset (0.9 * h->GetYaxis ()->GetTitleOffset ());

          h->Draw (!canvasExists && iPhi == 1 ? "e2" : "e2 same");

          LabelIAARatios (iCent, iPhi);
        } // end loop over phi
      } // end loop over cents
      c->cd (1);

      for (short iCent = 1; iCent < numCentBins; iCent++) {
        c->cd (iCent);
        myText (0.22, 0.24, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
      } // end loop over cents

      c->SaveAs (Form ("%s/iaa_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between central and peripheral Pb+Pb
////////////////////////////////////////////////////////////////////////////////////////////////
void Systematics :: PlotICPRatios (const short pSpc, const short pPtZ) {
  if (!backgroundSubtracted)
    SubtractBackground ();
  if (!icpCalculated)
    CalculateICP ();

  const int axisTextSize = 28;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    if (pSpc != -1 && iSpc != pSpc)
       continue; // allows user to define which plots should be made
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      if (pPtZ != -1 && iPtZ != pPtZ)
        continue; // allows user to define which plots should be made

      const char* canvasName = Form ("c_z_trk_pt_icp_%s_iPtZ%i", spc, iPtZ);
      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
      TCanvas* c = nullptr;
      if (canvasExists)
        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
      else {
        c = new TCanvas (canvasName, "", 500*(numCentBins-2), 500);
        c->Divide (numCentBins-2, 1);
        gDirectory->Add (c);
      }

      for (short iCent = 2; iCent < numCentBins; iCent++) {
        c->cd (iCent-1);
        gPad->SetLogx ();

        for (int iPhi = 1; iPhi < numPhiBins; iPhi++) {
          TH1D* h = h_z_trk_pt_icp[iSpc][iPtZ][0][iPhi][iCent];

          h->SetFillColorAlpha (fillColors[iPhi], 0.5);
          h->SetMarkerSize (0);
          h->SetLineColor (kBlack);
          h->SetLineWidth (0);

          h->GetXaxis ()->SetLimits (ptTrkBins[0], ptTrkBins[nPtTrkBins]);
          h->GetYaxis ()->SetRangeUser (0, 1.6);

          h->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
          h->GetYaxis ()->SetTitle ("I_{CP}");

          h->GetXaxis ()->SetTitleFont (43);
          h->GetXaxis ()->SetTitleSize (axisTextSize);
          h->GetXaxis ()->SetLabelFont (43);
          h->GetXaxis ()->SetLabelSize (axisTextSize);

          h->GetYaxis ()->SetTitleFont (43);
          h->GetYaxis ()->SetTitleSize (axisTextSize);
          h->GetYaxis ()->SetLabelFont (43);
          h->GetYaxis ()->SetLabelSize (axisTextSize);

          h->GetXaxis ()->SetTitleOffset (0.8 * h->GetXaxis ()->GetTitleOffset ());
          h->GetYaxis ()->SetTitleOffset (0.9 * h->GetYaxis ()->GetTitleOffset ());

          h->Draw (!canvasExists && iPhi == 1 ? "e2" : "e2 same");

          LabelICPRatios (iCent, iPhi);
        } // end loop over phi
      } // end loop over cents

      c->cd (1);
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        c->cd (iCent-1);
        myText (0.22, 0.24, kBlack, Form ("%i-%i%% / %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1], (int)centCuts[1], (int)centCuts[0]), 0.06);
      } // end loop over cents

      c->SaveAs (Form ("%s/icp_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
    } // end loop over pT^Z bins
  } // end loop over species
}
*/



#endif
