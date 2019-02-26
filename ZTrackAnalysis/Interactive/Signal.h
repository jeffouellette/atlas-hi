#ifndef __Signal_h__
#define __Signal_h__

#include "Params.h"
#include "Analysis.h"
#include "ZTrackUtilities.h"
#include "ZTrackAnalysis.h"
#include "MinbiasAnalysis.h"

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
  TH1D****** ZTracksSubYields = nullptr;
  TH1D****** ZTracksIAARatios = nullptr;
  TH1D****** ZTracksICPRatios = nullptr;

  Signal (Analysis* _obs, Analysis* _bkg) : Analysis (_obs) {
    SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");
    assert (!_obs);
    assert (!_bkg);

    obs = _obs;
    bkg = _bkg;

    ZTracksSubYields = Get5DArray <TH1D*> (3, nPtZBins, 2, numPhiBins, numCentBins); // iSpc, iPtZ, iData, iBkg, iPhi, iCent
    ZTracksIAARatios = Get5DArray <TH1D*> (3, nPtZBins, 2, numPhiBins, numCentBins); // iSpc, iPtZ, iData, iBkg, iPhi, iCent
    ZTracksICPRatios = Get5DArray <TH1D*> (3, nPtZBins, 2, numPhiBins, numCentBins); // iSpc, iPtZ, iData, iBkg, iPhi, iCent
  }

  void GenerateHistograms ();
  void PlotTrkYield ();
  void PlotIAARatios ();
  void PlotICPRatios ();
};


void Signal::GenerateHistograms () {

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        for (short iData = 0; iData < 2; iData++) {
          const char* data = iData == 0 ? "data" : "mc";
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            TH1D* thisHist = obs->ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent];
            double yieldNormFactor = obs->ZCounts[iSpc][iPtZ][iData][iCent]->GetBinContent (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);
            double yieldNormFactorError = obs->ZCounts[iSpc][iPtZ][iData][iCent]->GetBinError (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);
            //RescaleWithError (thisHist, yieldNormFactor, yieldNormFactorError);
            if (yieldNormFactor > 0)
              thisHist->Scale (1. / yieldNormFactor);

            thisHist = bkg->ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent];
            yieldNormFactor = bkg->ZCounts[iSpc][iPtZ][iData][iCent]->GetBinContent (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);
            yieldNormFactorError = bkg->ZCounts[iSpc][iPtZ][iData][iCent]->GetBinError (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);
            //RescaleWithError (thisHist, yieldNormFactor, yieldNormFactorError);
            if (yieldNormFactor > 0)
              thisHist->Scale (1. / yieldNormFactor);

            TH1D* subYield = new TH1D (Form ("subYieldPt_%s_%s_iPtZ%i_iPhi%i_iCent%i", spc, data, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
            subYield->Sumw2 ();

            subYield->Add (obs->ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent]);
            subYield->Add (bkg->ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent], -1);

            ZTracksSubYields[iSpc][iPtZ][iData][iPhi][iCent] = subYield;
          } // end loop over phi
        } // end loop over data types
      } // end loop over data types
    } // end loop over pT^Z bins
  } // end loop over species

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Plot pTtrk binned in dPhi
////////////////////////////////////////////////////////////////////////////////////////////////
void Signal::PlotTrkYield () {
  TCanvas* c = new TCanvas ("TrkYieldCanvas", "", 860*numCentBins, 1000);
  c->cd ();

  const double padRatio = 1.1; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = 1.0/ (padRatio+1.0);
  const double uPadY = 1.0 - dPadY;

  TPad** topPads = Get1DArray <TPad*> (numCentBins);
  TPad** bottomPads = Get1DArray <TPad*> (numCentBins);
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    c->cd ();

    TPad* topPad = new TPad ("topPad", "", 0+(1./numCentBins)*iCent, dPadY, (1./numCentBins)+(1./numCentBins)*iCent, 1);
    TPad* bottomPad = new TPad ("bottomPad", "", 0+(1./numCentBins)*iCent, 0, (1./numCentBins)+(1./numCentBins)*iCent, dPadY);
    topPads[iCent] = topPad;
    bottomPads[iCent] = bottomPad;

    topPad->SetTopMargin (0.04);
    topPad->SetBottomMargin (0);
    topPad->SetLeftMargin (0.17);
    topPad->SetRightMargin (0.06);
    bottomPad->SetTopMargin (0);
    bottomPad->SetBottomMargin (0.20);
    bottomPad->SetLeftMargin (0.17);
    bottomPad->SetRightMargin (0.06);
    topPad->Draw ();
    bottomPad->Draw ();
  }

  //double***** yieldNormFactor      = Get5DArray <double> (3, nPtZBins, 2, numCentBins, numPhiBins);
  //double***** yieldNormFactorError = Get5DArray <double> (3, nPtZBins, 2, numCentBins, numPhiBins);
  //for (short iSpc = 0; iSpc < 3; iSpc++) {
  //  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
  //    for (short iData = 0; iData < 2; iData++) {
  //      for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
  //        for (short iCent = 0; iCent < numCentBins; iCent++) {
  //          double f = 0, e = 0;
  //          TH1D* h = ZPtSpecs[iData][iCent][iSpc];
  //          if (iPhi < numPhiBins-1)
  //            f = h->IntegralAndError (h->FindBin (zPtBins[iPtZ]), h->FindBin (zPtBins[iPtZ+1])-1, e);
  //          else
  //            f = h->IntegralAndError (1, 3, e);
  //          yieldNormFactor[iSpc][iPtZ][iData][iCent][iPhi]      = f * (phiHighBins[iPhi] - phiLowBins[iPhi]);
  //          yieldNormFactorError[iSpc][iPtZ][iData][iCent][iPhi] = e * (phiHighBins[iPhi] - phiLowBins[iPhi]);
  //        } // end loop over cents
  //      } // end loop over phi
  //    } // end loop over data types
  //  } // end loop over pT^Z bins
  //} // end loop over species

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 0; iCent < numCentBins; iCent++) {

        TPad* topPad = topPads[iCent];
        TPad* bottomPad = bottomPads[iCent];

        //for (short iData = 0; iData < 2; iData++) {
        //  for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {

        //    TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent];

        //    //RescaleWithError (thisHist, yieldNormFactor[iSpc][iPtZ][iData][iCent][iPhi], yieldNormFactorError[iSpc][iPtZ][iData][iCent][iPhi]);
        //    if (yieldNormFactor[iSpc][iPtZ][iData][iCent][iPhi] > 0)
        //      thisHist->Scale (1. / yieldNormFactor[iSpc][iPtZ][iData][iCent][iPhi]);
        //  } // end loop over phi
        //} // end loop over data types
        

        double min = 1e30;
        double max = 0;
        for (short iData = 0; iData < 2; iData++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent];
            if (thisHist->GetMinimum (0) < min) min = thisHist->GetMinimum (0);
            if (thisHist->GetMaximum () > max)  max = thisHist->GetMaximum ();
          } // end loop over phi
        } // end loop over data types

        topPad->cd ();
        gPad->SetLogx ();
        gPad->SetLogy ();

        for (int iPhi = numPhiBins-2; iPhi >= 0; iPhi--) {
          TH1D* thisHist = ZTracksPt[iSpc][iPtZ][1][iPhi][iCent];

          thisHist->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
          thisHist->SetMarkerSize (0);
          thisHist->SetLineColor (kBlack);
          thisHist->SetLineWidth (0);
          thisHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
          thisHist->GetYaxis ()->SetTitle ("Per-trigger yield Y (#it{p}_{T})");
          thisHist->GetXaxis ()->SetTitleSize (0.03 / uPadY);
          thisHist->GetYaxis ()->SetTitleSize (0.03 / uPadY);
          thisHist->GetXaxis ()->SetLabelSize (0.03 / uPadY);
          thisHist->GetYaxis ()->SetLabelSize (0.03 / uPadY);
          thisHist->GetYaxis ()->SetTitleOffset (2.6 * uPadY);

          thisHist->GetYaxis ()->SetRangeUser (0.5*min, 2*max);
          if (iPhi == numPhiBins-2)
            thisHist->DrawCopy ("bar");
          else
            thisHist->DrawCopy ("bar same");
          thisHist->SetLineWidth (1);
          thisHist->Draw ("hist same");
        }
        gPad->RedrawAxis ();

        for (int iPhi = 0; iPhi < numPhiBins-1; iPhi++) {
          const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);
          TH1D* thisHist = ZTracksPt[iSpc][iPtZ][0][iPhi][iCent];
          
          TGraphAsymmErrors* thisGraph = make_graph (thisHist);
          RecenterGraph (thisGraph);
          ResetXErrors (thisGraph);

          thisGraph->SetMarkerStyle (markerStyle);
          thisGraph->SetMarkerColor (colors[iPhi]);
          thisGraph->SetLineColor (colors[iPhi]);
          thisGraph->SetMarkerSize (1);
          thisGraph->SetLineWidth (2);
  
          thisGraph->GetYaxis ()->SetRangeUser (0.5*min, 2*max);
          thisGraph->Draw ("P");

          if (iCent == 0)
            myText (0.485, 0.903, kBlack, "#bf{#it{ATLAS}} Internal", 0.034/uPadY);
          else if (iCent == numCentBins-1) {
            if (iPhi == 0) {
              myText (0.43, 0.91, kBlack, "Data", 0.03/uPadY);
              myText (0.574, 0.91, kBlack, "MC", 0.03/uPadY);
              myText (0.683, 0.91, kBlack, "#Delta#phi", 0.03/uPadY);
            }

            TVirtualPad* cPad = gPad; // store current pad
            const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
            const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
            TBox* b = TBoxNDC (0.61-0.024, 0.85-0.06*iPhi-0.016, 0.61+0.024, 0.85-0.06*iPhi+0.016);
            b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            b->Draw ("l");
            myMarkerText (0.512, 0.852-0.06*iPhi, colors[iPhi], markerStyle, "", 1.4, 0.03/uPadY);
            myText (0.68, 0.85-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.03/uPadY);
            cPad->cd ();
          }

        } // end loop over phi

        bottomPad->cd ();
        gPad->SetLogx ();
        gPad->SetLogy ();

        //for (short iData = 0; iData < 2; iData++) {
        //  const char* data = iData == 0 ? "data" : "mc";
        //  for (int iPhi = numPhiBins-2; iPhi >= 0; iPhi--) {
        //    TH1D* subYield = new TH1D (Form ("subYieldPt_%s_%s_iPtZ%i_iBkg%i_iPhi%i_iCent%i", spc, data, iPtZ, iBkg, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
        //    subYield->Sumw2 ();
        //    //subYield->Rebin (8);

        //    subYield->Add (ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent]);
        //    if (iBkg == 0)
        //      subYield->Add (ZTracksPt[iSpc][iPtZ][iData][0][iCent], -1);
        //    else
        //      subYield->Add (ZTracksPt[iSpc][iPtZ][iData][numPhiBins-1][iCent], -1);

        //    ZTracksSubYields[iSpc][iPtZ][iData][iBkg][iPhi][iCent] = subYield;
        //  } // end loop over phi
        //} // end loop over data types
        
        min = 1e30;
        max = 0;
        for (short iData = 0; iData < 2; iData++) {
          for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
            TH1D* thisHist = ZTracksSubYields[iSpc][iPtZ][iData][iPhi][iCent];
            if (thisHist->GetMinimum (0) < min) min = thisHist->GetMinimum (0);
            if (thisHist->GetMaximum () > max) max = thisHist->GetMaximum ();
          } // end loop over phi
        } // end loop over data types

        for (int iPhi = numPhiBins-2; iPhi >= 1; iPhi--) {
          TH1D* subYield = ZTracksSubYields[iSpc][iPtZ][1][iPhi][iCent];

          float delta = log10 (max) - log10 (min);
          subYield->GetYaxis ()->SetRangeUser (pow (10, log10 (min) - 0.1*delta), pow (10, log10 (max) + 0.1*delta));

          subYield->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
          subYield->SetLineColor (kBlack);
          subYield->SetMarkerSize (0);
          subYield->SetLineWidth (0);
          subYield->SetMarkerStyle (kFullCircle);
          subYield->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
          subYield->GetYaxis ()->SetTitle ("Signal Yield");
          subYield->GetXaxis ()->SetTitleSize (0.03 / uPadY);
          subYield->GetYaxis ()->SetTitleSize (0.03 / uPadY);
          subYield->GetXaxis ()->SetLabelSize (0.03 / uPadY);
          subYield->GetYaxis ()->SetLabelSize (0.03 / uPadY);
          subYield->GetXaxis ()->SetTitleOffset (1.5);
          subYield->GetYaxis ()->SetTitleOffset (2.8 * uPadY);

          if (iPhi == numPhiBins-2)
            subYield->DrawCopy ("bar");
          else
            subYield->DrawCopy ("bar same");
          subYield->SetLineWidth (1);
          subYield->Draw ("hist same");
        } // end loop over phi
        gPad->RedrawAxis ();

        for (int iPhi = numPhiBins-2; iPhi >= 1; iPhi--) {
          const Style_t markerStyle = (iPhi == 0 ? kFullSquare : kFullCircle);
          TGraphAsymmErrors* subYield = make_graph (ZTracksSubYields[iSpc][iPtZ][0][iPhi][iCent]);
          RecenterGraph (subYield);
          ResetXErrors (subYield);

          float delta = log10 (max) - log10 (min);
          subYield->GetYaxis ()->SetRangeUser (pow (10, log10 (min) - 0.1*delta), pow (10, log10 (max) + 0.1*delta));
          subYield->SetMarkerStyle (markerStyle);
          subYield->SetLineColor (colors[iPhi]);
          subYield->SetMarkerColor (colors[iPhi]);
          subYield->SetMarkerSize (1);
          subYield->SetLineWidth (2);
          subYield->SetMarkerStyle (kFullCircle);

          subYield->Draw ("P");
        } // end loop over phi
      } // end loop over cents
      

      topPads[0]->cd ();
      for (short iCent = 0; iCent < numCentBins; iCent++) {
        topPads[iCent]->cd ();
        if (iCent == 0)
          myText (0.22, 0.06, kBlack, "#it{pp}", 0.03 / uPadY);// or d#it{p}_{T}^{1,2} > 5 GeV", 0.03 / uPadY);
        else
          myText (0.22, 0.06, kBlack, Form ("Pb+Pb, %i-%i%%", (int)centCuts[iCent], (int)centCuts[iCent-1]), 0.03 / uPadY);
      }

      c->SaveAs (Form ("%s/pTTrk_dists_%s_iPtZ%i.pdf", plotPath.Data (), spc, iPtZ));

    } // end loop over pT^Z bins
  } // end loop over species
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plots subtracted yield ratios between Pb+Pb and pp
////////////////////////////////////////////////////////////////////////////////////////////////
void Signal::PlotIAARatios () {
  TCanvas* c = new TCanvas ("IAACanvas", "", 800*(numCentBins-1), 500);
  c->cd ();

  c->Divide (numCentBins-1, 1);

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iData = 0; iData < 2; iData++) {
        for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
          TH1D* ppHist = ZTracksSubYields[iSpc][iPtZ][0][iPhi][0];

          for (short iCent = 1; iCent < numCentBins; iCent++) {
            if (!ZTracksIAARatios[iSpc][iPtZ][iData][iPhi][iCent]) {
              TH1D* PbPbHist = (TH1D*)(ZTracksSubYields[iSpc][iPtZ][iData][iPhi][iCent]->Clone ());
              PbPbHist->Divide (ppHist);
              ZTracksIAARatios[iSpc][iPtZ][iData][iPhi][iCent] = PbPbHist;
            } 
          } // end loop over cents
        } // end loop over phi
      } // end loop over data types
    } // end loop over pT^Z bins
  } // end loop over species

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 1; iCent < numCentBins; iCent++) {
        for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
          c->cd (iCent);
          gPad->SetLogx ();
          
          TH1D* ratioHist = ZTracksIAARatios[iSpc][iPtZ][1][iPhi][iCent];

          //ratioHist->GetYaxis ()->SetRangeUser (min - 0.2*(max-min), max + 0.2*(max-min));
          ratioHist->GetYaxis ()->SetRangeUser (0, 3);
          ratioHist->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
          ratioHist->SetMarkerSize (0);
          ratioHist->SetLineColor (kBlack);
          ratioHist->SetLineWidth (0);
          ratioHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
          ratioHist->GetYaxis ()->SetTitle ("I_{AA}");
          ratioHist->GetXaxis ()->SetTitleSize (0.06);
          ratioHist->GetYaxis ()->SetTitleSize (0.06);
          ratioHist->GetXaxis ()->SetLabelSize (0.06);
          ratioHist->GetYaxis ()->SetLabelSize (0.06);
          ratioHist->GetXaxis ()->SetTitleOffset (1.2);
          ratioHist->GetYaxis ()->SetTitleOffset (1.2);

          if (iPhi == 1)
            ratioHist->DrawCopy ("bar");
          else
            ratioHist->DrawCopy ("bar same");
          ratioHist->SetLineWidth (1);
          ratioHist->Draw ("hist same");
        } // end loop over phi
        gPad->RedrawAxis ();

        for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
          c->cd (iCent);
          gPad->SetLogx ();
          //gPad->SetLogy ();
          
          TGraphAsymmErrors* ratioGraph = make_graph (ZTracksIAARatios[iSpc][iPtZ][0][iPhi][iCent]);
          RecenterGraph (ratioGraph);
          ResetXErrors (ratioGraph);
          deltaize (ratioGraph, 1+(0.5*(numPhiBins-2)+0.5-iPhi)*0.04, true); // 2.5 = 0.5*(numPhiBins-2)+0.5

          //ratioGraph->GetYaxis ()->SetRangeUser (min - 0.2*(max-min), max + 0.2*(max-min));
          ratioGraph->GetYaxis ()->SetRangeUser (0, 3);
          ratioGraph->SetLineColor (colors[iPhi]);
          ratioGraph->SetMarkerColor (colors[iPhi]);
          ratioGraph->SetMarkerSize (1.2);
          ratioGraph->SetLineWidth (2);
          ratioGraph->SetMarkerStyle (kFullCircle);

          ratioGraph->Draw ("P");

          if (iCent == 1)
            myText (0.49, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
          else if (iCent == numCentBins-1) {
            if (iPhi == 1) {
              myText (0.44, 0.91, kBlack, "Data", 0.05);
              myText (0.577, 0.91, kBlack, "MC", 0.05);
              myText (0.685, 0.91, kBlack, "#Delta#phi", 0.05);
            }
            TVirtualPad* cPad = gPad; // store current pad
            const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
            const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
            TBox* b = TBoxNDC (0.61-0.024, 0.91-0.06*iPhi-0.016, 0.61+0.024, 0.91-0.06*iPhi+0.016);
            b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            b->Draw ("l");
            myMarkerText (0.512, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
            myText (0.68, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.05);
            //myText (0.56, 0.90-0.07*(iPhi-1), colors[iPhi], Form ("(%s, %s)", lo, hi), 0.06);
            cPad->cd ();
          }
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
void Signal::PlotICPRatios () {
  TCanvas* c = new TCanvas ("ICPCanvas", "", 600*(numCentBins-2), 500);
  c->cd ();

  c->Divide (numCentBins-2, 1);

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iData = 0; iData < 2; iData++) {
        for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
          TH1D* periphHist = ZTracksSubYields[iSpc][iPtZ][iData][iPhi][1];

          for (short iCent = 2; iCent < numCentBins; iCent++) {
            if (!ZTracksICPRatios[iSpc][iPtZ][iData][iPhi][iCent]) {
              TH1D* centHist = (TH1D*)(ZTracksSubYields[iSpc][iPtZ][iData][iPhi][iCent]->Clone ());
              centHist->Divide (periphHist);
              ZTracksICPRatios[iSpc][iPtZ][iData][iPhi][iCent] = centHist;
            } 
          } // end loop over cents
        } // end loop over phi
      } // end loop over data types
    } // end loop over pT^Z bins
  } // end loop over species

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
    for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      for (short iCent = 2; iCent < numCentBins; iCent++) {
        for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
          c->cd (iCent-1);
          gPad->SetLogx ();
          
          TH1D* ratioHist = ZTracksICPRatios[iSpc][iPtZ][1][iPhi][iCent];

          //ratioHist->GetYaxis ()->SetRangeUser (min - 0.2*(max-min), max + 0.2*(max-min));
          ratioHist->GetYaxis ()->SetRangeUser (0, 3);
          ratioHist->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
          ratioHist->SetMarkerSize (0);
          ratioHist->SetLineColor (kBlack);
          ratioHist->SetLineWidth (0);
          ratioHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{ ch} [GeV]");
          ratioHist->GetYaxis ()->SetTitle ("I_{CP}");
          ratioHist->GetXaxis ()->SetTitleSize (0.06);
          ratioHist->GetYaxis ()->SetTitleSize (0.06);
          ratioHist->GetXaxis ()->SetLabelSize (0.06);
          ratioHist->GetYaxis ()->SetLabelSize (0.06);
          ratioHist->GetXaxis ()->SetTitleOffset (1.2);
          ratioHist->GetYaxis ()->SetTitleOffset (1.2);

          if (iPhi == 1)
            ratioHist->DrawCopy ("bar");
          else
            ratioHist->DrawCopy ("bar same");
          ratioHist->SetLineWidth (1);
          ratioHist->Draw ("hist same");
        } // end loop over phi
        gPad->RedrawAxis ();

        for (int iPhi = 1; iPhi < numPhiBins-1; iPhi++) {
          c->cd (iCent-1);
          
          TGraphAsymmErrors* ratioGraph = make_graph (ZTracksICPRatios[iSpc][iPtZ][0][iPhi][iCent]);
          RecenterGraph (ratioGraph);
          ResetXErrors (ratioGraph);
          deltaize (ratioGraph, 1+(0.5*(numPhiBins-2)+0.5-iPhi)*0.04, true); // 2.5 = 0.5*(numPhiBins-2)+0.5

          //ratioGraph->GetYaxis ()->SetRangeUser (min - 0.2*(max-min), max + 0.2*(max-min));
          ratioGraph->GetYaxis ()->SetRangeUser (0, 3);
          ratioGraph->SetLineColor (colors[iPhi]);
          ratioGraph->SetMarkerColor (colors[iPhi]);
          ratioGraph->SetMarkerSize (1.2);
          ratioGraph->SetLineWidth (2);
          ratioGraph->SetMarkerStyle (kFullCircle);

          ratioGraph->Draw ("P");

          if (iCent == 2)
            myText (0.59, 0.90, kBlack, "#bf{#it{ATLAS}}  Internal", 0.06);
          else if (iCent == numCentBins-1) {
            if (iPhi == 1) {
              myText (0.555, 0.91, kBlack, "Data", 0.05);
              myText (0.682, 0.91, kBlack, "MC", 0.05);
              myText (0.785, 0.91, kBlack, "#Delta#phi", 0.05);
            }
            TVirtualPad* cPad = gPad; // store current pad
            const char* lo = phiLowBins[iPhi] != 0 ? (phiLowBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiLowBins[iPhi]*8/pi)) : "0";
            const char* hi = phiHighBins[iPhi] != pi ? (phiHighBins[iPhi] == pi/2 ? "#pi/2" : Form ("%.0f#pi/8", phiHighBins[iPhi]*8/pi)) : "#pi";
            TBox* b = TBoxNDC (0.71-0.024, 0.91-0.06*iPhi-0.016, 0.71+0.024, 0.91-0.06*iPhi+0.016);
            b->SetFillColorAlpha (fillColors[iPhi], fillAlpha);
            b->Draw ("l");
            myMarkerText (0.612, 0.912-0.06*iPhi, colors[iPhi], kFullCircle, "", 1.4, 0.05);
            myText (0.78, 0.91-0.06*iPhi, kBlack, Form ("(%s, %s)", lo, hi), 0.05);
            //myText (0.56, 0.90-0.07*(iPhi-1), colors[iPhi], Form ("(%s, %s)", lo, hi), 0.06);
            cPad->cd ();
          }
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

#endif
