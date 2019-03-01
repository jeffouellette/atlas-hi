#ifndef __Signal_h__
#define __Signal_h__

#include "Params.h"
#include "Analysis.h"
#include "ZTrackUtilities.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class Signal : public Analysis {

  private:
  Analysis* obs;
  Analysis* bkg;

  public:
  //TH1D***** ZTracksSubYields = nullptr;
  //TH1D***** ZTracksIAARatios = nullptr;
  //TH1D***** ZTracksICPRatios = nullptr;

  Signal (Analysis* _obs, Analysis* _bkg);

  void SubtractBackground ();

  //void CalculateIAA ();
  //void CalculateICP ();

  //void PlotTrkYield (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  //void PlotIAARatios (const short pSpc = 2, const short pPtZ = nPtZBins-1);
  //void PlotICPRatios (const short pSpc = 2, const short pPtZ = nPtZBins-1);
};


////////////////////////////////////////////////////////////////////////////////////////////////
// Construct a signal object
////////////////////////////////////////////////////////////////////////////////////////////////
Signal::Signal (Analysis* _obs, Analysis* _bkg) : Analysis (_obs) {

  SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");
  assert (!_obs);
  assert (!_bkg);

  obs = _obs;
  bkg = _bkg;

  ZTracksSubYields = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  ZTracksIAARatios = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
  ZTracksICPRatios = Get4DArray <TH1D*> (3, nPtZBins, numPhiBins, numCentBins); // iSpc, iPtZ, iPhi, iCent
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Subtracts a background from the observed track yields
////////////////////////////////////////////////////////////////////////////////////////////////
void Signal::SubtractBackground () {
  if (backgroundSubtracted)
    return;
  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          TH1D* subYield = new TH1D (Form ("subYieldPt_%s_iPtZ%i_iPhi%i_iCent%i_signal_%s_%s", spc, iPtZ, iPhi, iCent, obs->Name (), bkg->Name ()), "", nPtTrkBins, ptTrkBins);
          subYield->Sumw2 ();

          subYield->Add (obs->ZTracksPt[iSpc][iPtZ][iPhi][iCent]);
          subYield->Add (bkg->ZTracksPt[iSpc][iPtZ][iPhi][iCent], -1);

          ZTracksSubYields[iSpc][iPtZ][iPhi][iCent] = subYield;
        } // end loop over phi
      } // end loop over centralities
    } // end loop over pT^Z bins
  } // end loop over species
  backgroundSubtracted = true;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//// Plot pTtrk binned in dPhi
//////////////////////////////////////////////////////////////////////////////////////////////////
//void Signal::PlotTrkYield (const short pSpc, const short pPtZ) {
//  const double padRatio = 1.1; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
//  const double dPadY = 1.0/ (padRatio+1.0);
//  const double uPadY = 1.0 - dPadY;
//
//  for (short iSpc = 0; iSpc < 3; iSpc++) {
//    if (pSpc != -1 && iSpc != pSpc)
//       continue; // allows user to define which plots should be made
//    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//      if (pPtZ != -1 && iPtZ != pPtZ)
//        continue; // allows user to define which plots should be made
//
//      const char* canvasName = Form ("TrkYieldCanvas_%s_iPtZ%i", spc, iPtZ);
//      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
//      TCanvas* c = nullptr;
//      if (canvasExists)
//        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
//      else {
//        c = new TCanvas (canvasName, "", 860*numCentBins, 1000);
//        gDirectory->Add (c);
//      }
//
//      for (short iCent = 0; iCent < numCentBins; iCent++) {
//        c->cd ();
//
//        const char* topPadName = Form ("topPad_%s_iPtZ%i_iCent%i", spc, iPtZ, iCent);
//        const char* bottomPadName = Form ("bottomPad_%s_iPtZ%i_iCent%i", spc, iPtZ, iCent);
//
//        TPad* topPad = nullptr, *bottomPad = nullptr;
//        if (!canvasExists) {
//          topPad = new TPad (topPadName, "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
//          bottomPad = new TPad (bottomPadName, "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);
//
//          gDirectory->Add (topPad);
//          gDirectory->Add (bottomPad);
//
//          topPad->SetTopMargin (0.04);
//          topPad->SetBottomMargin (0);
//          topPad->SetLeftMargin (0.17);
//          topPad->SetRightMargin (0.06);
//          bottomPad->SetTopMargin (0);
//          bottomPad->SetBottomMargin (0.20);
//          bottomPad->SetLeftMargin (0.17);
//          bottomPad->SetRightMargin (0.06);
//          topPad->Draw ();
//          bottomPad->Draw ();
//        }
//        else {
//          topPad = dynamic_cast<TPad*> (gDirectory->Get (topPadName));
//          bottomPad = dynamic_cast<TPad*> (gDirectory->Get (bottomPadName));
//        } 
//
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//
//          TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iPhi][iCent];
//          TH1D* countsHist = ZCounts[iSpc][iPtZ][iCent];
//          double yieldNormFactor = countsHist->GetBinContent (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);
//          double yieldNormFactorError = countsHist->GetBinError (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);
//
//          //RescaleWithError (thisHist, yieldNormFactor, yieldNormFactorError);
//          if (yieldNormFactor > 0)
//            thisHist->Scale (1. / yieldNormFactor);
//        } // end loop over phi
//
//        topPad->cd ();
//        GetDrawnObjects ();
//        gPad->SetLogx ();
//        gPad->SetLogy ();
//
//        double min = 1e30, max = 0;
//        GetMinAndMax (min, max, true);
//        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
//          TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iPhi][iCent];
//          if (thisHist->GetMinimum (0) < min) min = thisHist->GetMinimum (0);
//          if (thisHist->GetMaximum () > max)  max = thisHist->GetMaximum ();
//        } // end loop over phi
//        min *= 0.5;
//        max *= 2;
//        SetMinAndMax (min, max);
//
//        if (plotFill) {
//          for (int iPhi = numPhiBins-1; iPhi >= 0; iPhi--) {
//            TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iPhi][iCent];
//
//            thisHist->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
//            thisHist->SetMarkerSize (0);
//            thisHist->SetLineColor (kBlack);
//            thisHist->SetLineWidth (0);
//            thisHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
//            thisHist->GetYaxis ()->SetTitle ("Per-trigger yield Y (#it{p}_{T})");
//            thisHist->GetXaxis ()->SetTitleSize (0.03 / uPadY);
//            thisHist->GetYaxis ()->SetTitleSize (0.03 / uPadY);
//            thisHist->GetXaxis ()->SetLabelSize (0.03 / uPadY);
//            thisHist->GetYaxis ()->SetLabelSize (0.03 / uPadY);
//            thisHist->GetYaxis ()->SetTitleOffset (2.6 * uPadY);
//
//            thisHist->GetYaxis ()->SetRangeUser (min, max);
//            thisHist->DrawCopy (!canvasExists && iPhi == numPhiBins-1 ? "bar" : "bar same");
//            thisHist->SetLineWidth (1);
//            thisHist->Draw ("hist same");
//          }
//          gPad->RedrawAxis ();
//        }
//        else {
//          for (int iPhi = 0; iPhi < numPhiBins-1; iPhi++) {
//            const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);
//            TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iPhi][iCent];
//            
//            TGraphAsymmErrors* thisGraph = make_graph (thisHist);
//            RecenterGraph (thisGraph);
//            ResetXErrors (thisGraph);
//
//            thisGraph->SetMarkerStyle (markerStyle);
//            thisGraph->SetMarkerColor (colors[iPhi]);
//            thisGraph->SetLineColor (colors[iPhi]);
//            thisGraph->SetMarkerSize (1);
//            thisGraph->SetLineWidth (2);
//  
//            thisGraph->GetYaxis ()->SetRangeUser (0.5*min, 2*max);
//            
//            thisGraph->Draw (!canvasExists && iPhi == numPhiBins-1 ? "AP" : "P");
//          } // end loop over phi
//        }
//
//        if (!canvasExists)
//          for (int iPhi = 0; iPhi < numPhiBins; iPhi++)
//            LabelTrkYield (iCent, iPhi);
//
//        bottomPad->cd ();
//        GetDrawnObjects ();
//        gPad->SetLogx ();
//        gPad->SetLogy ();
//
//        min = 1e30, max = 0;
//        GetMinAndMax (min, max, true);
//        for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
//          TH1D* thisHist = ZTracksSubYields[iSpc][iPtZ][iPhi][iCent];
//          if (thisHist->GetMinimum (0) < min) min = thisHist->GetMinimum (0);
//          if (thisHist->GetMaximum () > max) max = thisHist->GetMaximum ();
//        } // end loop over phi
//        float delta = log10 (max) - log10 (min);
//        min = pow (10, log10 (min) - 0.1*delta);
//        max = pow (10, log10 (max) + 0.1*delta);
//        SetMinAndMax (min, max);
//
//        if (plotFill) {
//          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
//            TH1D* subYield = ZTracksSubYields[iSpc][iPtZ][iPhi][iCent];
//
//            subYield->GetYaxis ()->SetRangeUser (min, max);
//
//            subYield->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
//            subYield->SetLineColor (kBlack);
//            subYield->SetMarkerSize (0);
//            subYield->SetLineWidth (0);
//            subYield->SetMarkerStyle (kFullCircle);
//            subYield->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
//            subYield->GetYaxis ()->SetTitle ("Signal Yield");
//            subYield->GetXaxis ()->SetTitleSize (0.03 / uPadY);
//            subYield->GetYaxis ()->SetTitleSize (0.03 / uPadY);
//            subYield->GetXaxis ()->SetLabelSize (0.03 / uPadY);
//            subYield->GetYaxis ()->SetLabelSize (0.03 / uPadY);
//            subYield->GetXaxis ()->SetTitleOffset (1.5);
//            subYield->GetYaxis ()->SetTitleOffset (2.8 * uPadY);
//
//            subYield->DrawCopy (!canvasExists && iPhi == numPhiBins-1 ? "bar" : "bar same");
//            subYield->SetLineWidth (1);
//            subYield->Draw ("hist same");
//          } // end loop over phi
//          gPad->RedrawAxis ();
//        }
//        else {
//          for (int iPhi = numPhiBins-1; iPhi >= 1; iPhi--) {
//            const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);
//            TGraphAsymmErrors* subYield = make_graph (ZTracksSubYields[iSpc][iPtZ][iPhi][iCent]);
//            RecenterGraph (subYield);
//            ResetXErrors (subYield);
//
//            float delta = log10 (max) - log10 (min);
//            subYield->GetYaxis ()->SetRangeUser (pow (10, log10 (min) - 0.1*delta), pow (10, log10 (max) + 0.1*delta));
//            subYield->SetMarkerStyle (markerStyle);
//            subYield->SetLineColor (colors[iPhi]);
//            subYield->SetMarkerColor (colors[iPhi]);
//            subYield->SetMarkerSize (1);
//            subYield->SetLineWidth (2);
//            subYield->SetMarkerStyle (kFullCircle);
//
//            subYield->Draw (!canvasExists && iPhi == numPhiBins-1 ? "AP" : "P");
//          } // end loop over phi
//        }
//      } // end loop over cents
//      
//      c->SaveAs (Form ("%s/pTTrk_dists_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
//    } // end loop over pT^Z bins
//  } // end loop over species
//}


//////////////////////////////////////////////////////////////////////////////////////////////////
//// Calculates subtracted yield ratios between Pb+Pb and pp
//////////////////////////////////////////////////////////////////////////////////////////////////
//void Signal::CalculateIAA () {
//  for (short iSpc = 0; iSpc < 3; iSpc++) {
//    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//      for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
//        TH1D* ppHist = ZTracksSubYields[iSpc][iPtZ][iPhi][0];
//
//        for (short iCent = 1; iCent < numCentBins; iCent++) {
//          if (!ZTracksIAARatios[iSpc][iPtZ][iPhi][iCent]) {
//            TH1D* PbPbHist = (TH1D*)(ZTracksSubYields[iSpc][iPtZ][iPhi][iCent]->Clone ());
//            PbPbHist->Divide (ppHist);
//            ZTracksIAARatios[iSpc][iPtZ][iPhi][iCent] = PbPbHist;
//          } 
//        } // end loop over cents
//      } // end loop over phi
//    } // end loop over pT^Z bins
//  } // end loop over species
//}


//////////////////////////////////////////////////////////////////////////////////////////////////
//// Plots subtracted yield ratios between Pb+Pb and pp
//////////////////////////////////////////////////////////////////////////////////////////////////
//void Signal::PlotIAARatios (const short pSpc, const short pPtZ) {
//  for (short iSpc = 0; iSpc < 3; iSpc++) {
//    if (pSpc != -1 && iSpc != pSpc)
//       continue; // allows user to define which plots should be made
//    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//      if (pPtZ != -1 && iPtZ != pPtZ)
//        continue; // allows user to define which plots should be made
//
//      const char* canvasName = Form ("IAACanvas_%s_iPtZ%i", spc, iPtZ);
//      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
//      TCanvas* c = nullptr;
//      if (canvasExists)
//        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
//      else {
//        c = new TCanvas (canvasName, "", 860*numCentBins, 1000);
//        c->Divide (numCentBins-1, 1);
//        gDirectory->Add (c);
//      }
//
//      for (short iCent = 1; iCent < numCentBins; iCent++) {
//        c->cd (iCent);
//        gPad->SetLogx ();
//
//        if (plotFill) {
//          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
//            TH1D* ratioHist = ZTracksIAARatios[iSpc][iPtZ][iPhi][iCent];
//
//            ratioHist->GetYaxis ()->SetRangeUser (0, 3);
//            ratioHist->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
//            ratioHist->SetMarkerSize (0);
//            ratioHist->SetLineColor (kBlack);
//            ratioHist->SetLineWidth (0);
//            ratioHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
//            ratioHist->GetYaxis ()->SetTitle ("I_{AA}");
//            ratioHist->GetXaxis ()->SetTitleSize (0.06);
//            ratioHist->GetYaxis ()->SetTitleSize (0.06);
//            ratioHist->GetXaxis ()->SetLabelSize (0.06);
//            ratioHist->GetYaxis ()->SetLabelSize (0.06);
//            ratioHist->GetXaxis ()->SetTitleOffset (1.2);
//            ratioHist->GetYaxis ()->SetTitleOffset (1.2);
//
//            ratioHist->DrawCopy (!canvasExists && iPhi == 1 ? "bar" : "bar same");
//            ratioHist->SetLineWidth (1);
//            ratioHist->Draw ("hist same");
//          } // end loop over phi
//          gPad->RedrawAxis ();
//        }
//        else {
//          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
//            TGraphAsymmErrors* ratioGraph = make_graph (ZTracksIAARatios[iSpc][iPtZ][iPhi][iCent]);
//            RecenterGraph (ratioGraph);
//            ResetXErrors (ratioGraph);
//            deltaize (ratioGraph, 1+(0.5*(numPhiBins-1)-iPhi)*0.04, true); // 2.5 = 0.5*(numPhiBins-1)
//
//            ratioGraph->GetYaxis ()->SetRangeUser (0, 3);
//            ratioGraph->SetLineColor (colors[iPhi]);
//            ratioGraph->SetMarkerColor (colors[iPhi]);
//            ratioGraph->SetMarkerSize (1.2);
//            ratioGraph->SetLineWidth (2);
//            ratioGraph->SetMarkerStyle (kFullCircle);
//
//            ratioGraph->Draw (!canvasExists && iPhi == 1 ? "AP" : "P");
//
//            if (iCent == 1)
//              myText (0.49, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
//            else if (iCent == numCentBins-1) {
//              if (iPhi == 1) {
//                myText (0.44, 0.91, kBlack, "Data", 0.05);
//                myText (0.577, 0.91, kBlack, "MC", 0.05);
//                myText (0.685, 0.91, kBlack, "#Delta#phi", 0.05);
//              }
//              TVirtualPad* cPad = gPad; // store current pad
//              const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
//              const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
//              TBox* b = TBoxNDC (0.61-0.024, 0.91-0.06*iPhi-0.016, 0.61+0.024, 0.91-0.06*iPhi+0.016);
//              b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
//              b->Draw ("l");
//              myMarkerText (0.512, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
//              myText (0.68, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.05);
//              cPad->cd ();
//            }
//          } // end loop over phi
//        }
//      } // end loop over cents
//      c->cd (1);
//
//      for (short iCent = 1; iCent < numCentBins; iCent++) {
//        c->cd (iCent);
//        myText (0.22, 0.24, kBlack, Form ("%i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.06);
//      } // end loop over cents
//
//      c->SaveAs (Form ("%s/iaa_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
//    } // end loop over pT^Z bins
//  } // end loop over species
//}


//////////////////////////////////////////////////////////////////////////////////////////////////
//// Calculates subtracted yield ratios between central and peripheral Pb+Pb
//////////////////////////////////////////////////////////////////////////////////////////////////
//void Signal::CalculateICP () {
//  for (short iSpc = 0; iSpc < 3; iSpc++) {
//    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//      for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
//        TH1D* periphHist = ZTracksSubYields[iSpc][iPtZ][iPhi][1];
//
//        for (short iCent = 2; iCent < numCentBins; iCent++) {
//          if (!ZTracksICPRatios[iSpc][iPtZ][iPhi][iCent]) {
//            TH1D* centHist = (TH1D*)(ZTracksSubYields[iSpc][iPtZ][iPhi][iCent]->Clone ());
//            centHist->Divide (periphHist);
//            ZTracksICPRatios[iSpc][iPtZ][iPhi][iCent] = centHist;
//          } 
//        } // end loop over cents
//      } // end loop over phi
//    } // end loop over pT^Z bins
//  } // end loop over species
//  return;
//}


//////////////////////////////////////////////////////////////////////////////////////////////////
//// Plots subtracted yield ratios between central and peripheral Pb+Pb
//////////////////////////////////////////////////////////////////////////////////////////////////
//void Signal::PlotICPRatios (const short pSpc, const short pPtZ) {
//  for (short iSpc = 0; iSpc < 3; iSpc++) {
//    if (pSpc != -1 && iSpc != pSpc)
//       continue; // allows user to define which plots should be made
//    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
//    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
//      if (pPtZ != -1 && iPtZ != pPtZ)
//        continue; // allows user to define which plots should be made
//
//      const char* canvasName = Form ("ICPCanvas_%s_iPtZ%i", spc, iPtZ);
//      const bool canvasExists = (gDirectory->Get (canvasName) != nullptr);
//      TCanvas* c = nullptr;
//      if (canvasExists)
//        c = dynamic_cast<TCanvas*>(gDirectory->Get (canvasName));
//      else {
//        c = new TCanvas (canvasName, "", 600*(numCentBins-1), 500);
//        c->Divide (numCentBins-2, 1);
//        gDirectory->Add (c);
//      }
//
//      for (short iCent = 2; iCent < numCentBins; iCent++) {
//        c->cd (iCent-1);
//        gPad->SetLogx ();
//
//        if (plotFill) {
//          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
//            TH1D* ratioHist = ZTracksICPRatios[iSpc][iPtZ][iPhi][iCent];
//
//            ratioHist->GetYaxis ()->SetRangeUser (0, 3);
//            ratioHist->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
//            ratioHist->SetMarkerSize (0);
//            ratioHist->SetLineColor (kBlack);
//            ratioHist->SetLineWidth (0);
//            ratioHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
//            ratioHist->GetYaxis ()->SetTitle ("I_{CP}");
//            ratioHist->GetXaxis ()->SetTitleSize (0.06);
//            ratioHist->GetYaxis ()->SetTitleSize (0.06);
//            ratioHist->GetXaxis ()->SetLabelSize (0.06);
//            ratioHist->GetYaxis ()->SetLabelSize (0.06);
//            ratioHist->GetXaxis ()->SetTitleOffset (1.2);
//            ratioHist->GetYaxis ()->SetTitleOffset (1.2);
//
//            ratioHist->DrawCopy (!canvasExists && iPhi == 1 : "bar" : "bar same");
//            ratioHist->SetLineWidth (1);
//            ratioHist->Draw ("hist same");
//          } // end loop over phi
//          gPad->RedrawAxis ();
//        }
//        else {
//          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
//            c->cd (iCent-1);
//            
//            TGraphAsymmErrors* ratioGraph = make_graph (ZTracksICPRatios[iSpc][iPtZ][iPhi][iCent]);
//            RecenterGraph (ratioGraph);
//            ResetXErrors (ratioGraph);
//            deltaize (ratioGraph, 1+(0.5*(numPhiBins-1)-iPhi)*0.04, true); // 2.5 = 0.5*(numPhiBins-1)
//
//            //ratioGraph->GetYaxis ()->SetRangeUser (min - 0.2*(max-min), max + 0.2*(max-min));
//            ratioGraph->GetYaxis ()->SetRangeUser (0, 3);
//            ratioGraph->SetLineColor (colors[iPhi]);
//            ratioGraph->SetMarkerColor (colors[iPhi]);
//            ratioGraph->SetMarkerSize (1.2);
//            ratioGraph->SetLineWidth (2);
//            ratioGraph->SetMarkerStyle (kFullCircle);
//
//            ratioGraph->Draw (!canvasExists && iPhi == 1 ? "AP" : "P");
//
//            if (iCent == 2)
//              myText (0.59, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
//            else if (iCent == numCentBins-1) {
//              if (iPhi == 1) {
//                myText (0.555, 0.91, kBlack, "Data", 0.05);
//                myText (0.682, 0.91, kBlack, "MC", 0.05);
//                myText (0.785, 0.91, kBlack, "#Delta#phi", 0.05);
//              }
//              TVirtualPad* cPad = gPad; // store current pad
//              const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
//              const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
//              TBox* b = TBoxNDC (0.71-0.024, 0.91-0.06*iPhi-0.016, 0.71+0.024, 0.91-0.06*iPhi+0.016);
//              b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
//              b->Draw ("l");
//              myMarkerText (0.612, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
//              myText (0.78, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.05);
//              //myText (0.56, 0.90-0.07*(iPhi-1), colors[iPhi], Form ("(%s, %s)", lo, hi), 0.06);
//              cPad->cd ();
//            }
//          } // end loop over phi
//        }
//      } // end loop over cents
//
//      c->cd (1);
//      for (short iCent = 2; iCent < numCentBins; iCent++) {
//        c->cd (iCent-1);
//        myText (0.22, 0.24, kBlack, Form ("%i-%i%% / %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1], (int)centCuts[1], (int)centCuts[0]), 0.06);
//      } // end loop over cents
//
//      c->SaveAs (Form ("%s/icp_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));
//    } // end loop over pT^Z bins
//  } // end loop over species
//}

#endif
