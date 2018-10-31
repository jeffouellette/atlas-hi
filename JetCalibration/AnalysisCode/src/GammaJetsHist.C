#include "GammaJetsHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TVirtualFitter.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

using namespace atlashi;

namespace JetCalibration {


void GammaJetsHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("GammaJets/", "JetCalibration/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TH3D****** gJetHists = Get5DArray <TH3D*> (2, 3, 3, 3, 3); // iYear, iPer, iData, iErr, iPCut
  TH2D****** gJetCounts = Get5DArray <TH2D*> (2, 3, 3, 3, 2); // iYear, iPer, iData, iPCut, iWgt

  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0)
    continue; // ignore old insitu factors if not desired

   const TString year = (iYear == 0 ? "2015" : "2016");

   for (short iPer = 0; iPer < 3; iPer++) {
    const TString period = (iPer == 0 ? "periodA" : (iPer == 1 ? "periodB" : "periodAB"));

    for (short iData = 0; iData < 3; iData++) { // iData is 0 for data, 1 for MC
     if (iData == 2 && skipSignalMC)
      continue; // ignore signal MC if not desired

     const TString dataType = (iData == 0 ? "data" : (iData == 1 ? "mc_overlay" : "mc_signalonly"));

     for (short iPCut = 0; iPCut < 3; iPCut++) {
      const TString pCut = (iPCut == 0 ? "tight" : (iPCut == 1 ? "lnt" : "signal")); // lnt = loose non-tight

      for (short iErr = 0; iErr < 3; iErr++) {
       if (iErr != 1 && iData != 0)
        continue; // ignore systematics for MC

       const TString error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

       const TString name = Form ("gJetPtRatio_%s_%s_%s_%s_%s", year.Data (), dataType.Data (), error.Data (), period.Data (), pCut.Data ());
       gJetHists[iYear][iPer][iData][iErr][iPCut] = (TH3D*)inFile->Get (name);
      }

      for (short iWgt = 0; iWgt < 2; iWgt++) {
       const TString weight = (iWgt == 0 ? "unweighted" : "weighted");

       const TString name = Form ("gJetCounts_%s_%s_%s_%s_%s", year.Data (), dataType.Data (), period.Data (), pCut.Data (), weight.Data ());
       gJetCounts[iYear][iPer][iData][iPCut][iWgt] = (TH2D*)inFile->Get (name);
      }
     }
    }
   }
  }


  TLine* glines[5] = {};
  TLine* getalines[5] = {};
  TLine* xlines[5] = {};
  for (short i = 0; i < 5; i++) {
   const float dg = 0.05;
   const float dx = 0.2;

   glines[i] = new TLine (pbins[0], 1.0-1*dg+dg*i, pbins[numpbins], 1.0-1*dg+dg*i);
   getalines[i] = new TLine (etabins[0], 1.0-1*dg+dg*i, etabins[numetabins], 1.0-1*dg+dg*i);
   xlines[i] = new TLine (xjrefbins[0], 1.0-2*dx+dx*i, xjrefbins[numxjrefbins], 1.0-2*dx+dx*i);

   if (1.0-1*dg+dg*i == 1) glines[i]->SetLineStyle (1);
   else glines[i]->SetLineStyle (3);
   if (1.0-1*dg+dg*i == 1) getalines[i]->SetLineStyle (1);
   else getalines[i]->SetLineStyle (3);
   if (1.0-2*dx+dx*i == 1) xlines[i]->SetLineStyle (1);
   else xlines[i]->SetLineStyle (3);
  }


  /**** Canvas definitions ****/
  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);
  const double padRatio = 1.5; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = 1.0/ (padRatio+1.0);
  const double uPadY = 1.0 - dPadY;
  TPad* topPad = new TPad ("topPad", "", 0, dPadY, 1, 1);
  TPad* bottomPad = new TPad ("bottomPad", "", 0, 0, 1, dPadY);
  topPad->SetBottomMargin (0);
  topPad->SetLeftMargin (-0.20);
  bottomPad->SetTopMargin (0);
  bottomPad->SetBottomMargin (0.30);
  bottomPad->SetLeftMargin (-0.20);
  topPad->Draw ();
  bottomPad->Draw ();


  /**** Define local histograms, graphs, etc. ****/
  TH1D *vJetHist = NULL, *vJetHist_mc_overlay = NULL, *vJetHist_mc_signal = NULL, *vJetHist_lo = NULL, *vJetHist_hi = NULL, *vJetHist_rat = NULL, *vJetHist_rat_lo = NULL, *vJetHist_rat_hi = NULL, *vJetRatioCI = NULL;
  TH2D *proj = NULL, *proj_mc_overlay = NULL, *proj_mc_signal = NULL, *proj_lo = NULL, *proj_hi = NULL;
  TGraphAsymmErrors *vJetGraph_sys = NULL, *vJetGraph_rat_sys = NULL;
  TF1* vJetRatioFit = NULL;
  char* plotName;

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
   const TString period = (iPer == 0 ? "Period A" : (iPer == 1 ? "Period B" : "Period A+B"));
   const TString perType = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));

   for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta

    int eta_lo = iEta+1;
    int eta_hi = iEta+1;
    if (iEta == numetabins) {
     eta_lo = eta_lo_comb;
     eta_hi = eta_hi_comb;
    }

    /**** Plots GammaJet info as a function of p_T^ref****/
    for (short iYear = 0; iYear < 2; iYear++) {
     if (skipOldInsitu && iYear == 0) continue;

     const Style_t dataStyle = (iYear == 0 ? 24 : 20);
     const Style_t mcOverlayStyle = 33;
     const Style_t mcSignalStyle = 34;

     topPad->cd ();
     topPad->SetLogx ();

     proj = Project2D ("", gJetHists[iYear][iPer][0][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist = GetProfileX ("vJetHist", proj, numpbins, pbins, false);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     vJetHist->SetYTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
     double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;
     vJetHist->SetAxisRange (middle - 0.35, middle + 0.35, "Y");
     vJetHist->SetMarkerStyle (dataStyle);
     vJetHist->SetMarkerColor (dataColor);
     vJetHist->SetLineColor (dataColor);
     vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", gJetHists[iYear][iPer][0][0][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numpbins, pbins, false);

     proj_hi = Project2D ("", gJetHists[iYear][iPer][0][2][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numpbins, pbins, false);

     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     vJetGraph_sys->SetFillColor (dataColor);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", gJetHists[iYear][iPer][1][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numpbins, pbins, false);

     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

     if (iPer != 1 && !skipSignalMC) {
      proj_mc_signal = Project2D ("", gJetHists[iYear][iPer][2][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
      vJetHist_mc_signal = GetProfileX ("vJetHist_mc_signal", proj_mc_signal, numpbins, pbins, false);

      vJetHist_mc_signal->SetMarkerStyle (mcSignalStyle);
      vJetHist_mc_signal->SetMarkerColor (mcSignalColor);
      vJetHist_mc_signal->SetLineColor (mcSignalColor);
     }

     //for (short iErr = 0; iErr < 3; iErr++) {
     // const TString error = (iErr == 0 ? "syslo" : (iErr == 1 ? "stat" : "syshi"));
     // TString periodStr = "periodA";
     // if (iPer == 1) periodStr = "periodB";
     // else if (iPer == 2) periodStr = "periodAB";
     // gJetHistDifference[iPer][iEta][iErr] = new TH1D (Form ("gJetPtRatio_diff%i_%s_%s", iEta, error.Data (), periodStr.Data ()), ";#it{p}_{T}^{ref} #left[GeV#right]", numetabins, etabins);
     // for (short iEta = 1; iEta <= numetabins; iEta++) {
     //  double dataVal, dataErr;
     //  switch (iErr) {
     //   case 0:
     //    dataVal = vJetHist_lo->GetBinContent (iEta);
     //    dataErr = vJetHist_lo->GetBinError (iEta);
     //    break;
     //   case 2:
     //    dataVal = vJetHist_hi->GetBinContent (iEta);
     //    dataErr = vJetHist_hi->GetBinError (iEta);
     //    break;
     //   default:
     //    dataVal = vJetHist->GetBinContent (iEta);
     //    dataErr = vJetHist->GetBinError (iEta);
     //  } 
     //  gJetHistDifference[iPer][iEta][iErr]->SetBinContent (iEta, dataVal - vJetHist_mc_overlay->GetBinContent (iEta));
     //  gJetHistDifference[iPer][iEta][iErr]->SetBinError (iEta, TMath::Sqrt (TMath::Power (dataErr,2) + TMath::Power (vJetHist_mc_overlay->GetBinError (iEta),2)));
     // }
     // gJetHistDifference[iPer][iEta][iErr]->Write ();
     //}
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     if (skipOldInsitu || iYear == 0) {
      vJetHist->DrawCopy ("e1 x0");
      vJetHist_mc_overlay->DrawCopy ("same e1 x0"); // insitu factors are not applied to MC
      if (iPer != 1 && !skipSignalMC) vJetHist_mc_signal->DrawCopy ("same e1 x0");
     }
     else {
      vJetHist->DrawCopy ("same e1 x0");
     }
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");
     //for (TLine* line : dplines) line->Draw ();

     if (skipOldInsitu || iYear == 0) {
      int countsData = 0, countsMC = 0, countsMC_sig = 0;
      if (exclusive && iEta == numetabins) {
       countsData = gJetCounts[iYear][iPer][0][0][0]->Integral () - gJetCounts[iYear][iPer][0][0][0]->Integral (1, gJetCounts[iYear][iPer][0][0][0]->GetNbinsX(), eta_lo, eta_hi);
       countsMC = gJetCounts[iYear][iPer][1][0][0]->Integral () - gJetCounts[iYear][iPer][1][0][0]->Integral (1, gJetCounts[iYear][iPer][1][0][0]->GetNbinsX(), eta_lo, eta_hi);
       if (!skipSignalMC)
        countsMC_sig = gJetCounts[iYear][iPer][2][0][0]->Integral () - gJetCounts[iYear][iPer][2][0][0]->Integral (1, gJetCounts[iYear][iPer][2][0][0]->GetNbinsX(), eta_lo, eta_hi);
      }
      else {
       countsData = gJetCounts[iYear][iPer][0][0][0]->Integral (1, gJetCounts[iYear][iPer][0][0][0]->GetNbinsX(), eta_lo, eta_hi);
       countsMC = gJetCounts[iYear][iPer][1][0][0]->Integral (1, gJetCounts[iYear][iPer][1][0][0]->GetNbinsX(), eta_lo, eta_hi);
       if (!skipSignalMC)
        countsMC_sig = gJetCounts[iYear][iPer][2][0][0]->Integral (1, gJetCounts[iYear][iPer][2][0][0]->GetNbinsX(), eta_lo, eta_hi);
      }

      myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countsData), 1.25, 0.04/uPadY);
      myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
      if (iPer != 1 && !skipSignalMC)
       myMarkerText (0.175, 0.74, mcSignalColor, kFullCross, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", countsMC_sig), 1.25, 0.04/uPadY);
      if (!skipOldInsitu) {
       myMarkerText (0.65, 0.66, dataColor, kFullCircle, "2016 Insitu Factors", 1.25, 0.04/uPadY);
       myMarkerText (0.65, 0.59, dataColor, kOpenCircle, "2015 Insitu Factors", 1.25, 0.04/uPadY);
      }
      myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      if (eta_lo != 1 || eta_hi != numetabins)
       if (!exclusive || iEta != numetabins)
        myText (0.155, 0.15, dataColor, Form ("#gamma + Jet, %g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04/uPadY);
       else
        myText (0.155, 0.15, dataColor, Form ("#gamma + Jet, %g < #left|#eta_{det}^{Jet}#right|", etabins[eta_hi]), 0.04/uPadY);
      else
       myText (0.155, 0.15, dataColor, "#gamma + Jet", 0.04/uPadY);
      myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);
     }

     bottomPad->cd ();
     bottomPad->SetLogx ();

     for (int iMC = 1; iMC < 3; iMC++) { // loops over overlay then signal MC
      if (iMC == 2 && (skipSignalMC || iPer == 1)) continue; // no signal only for period B!

      const TString mcType = (iMC == 1 ? "overlay":"signal");
      TH2D* proj_mc = (iMC == 1 ? proj_mc_overlay:proj_mc_signal);

      vJetHist_rat = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_%s_iEta%i", mcType.Data (), iEta)), proj, proj_mc, numpbins, pbins, false, "x");
      vJetHist_rat_lo = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_lo_%s_iEta%i", mcType.Data (), iEta)), proj_lo, proj_mc, numpbins, pbins, false, "x");
      vJetHist_rat_hi = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_hi_%s_iEta%i", mcType.Data (), iEta)), proj_hi, proj_mc, numpbins, pbins, false, "x");

      vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
      CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
      if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
      if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
      //vJetGraph_rat_sys->SetFillColor (iMC == 1 ? dataColor:mcSignalColor);
      vJetGraph_rat_sys->SetFillColor (dataColor);
      vJetGraph_rat_sys->SetFillStyle (3001);

      vJetHist_rat->SetXTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
      vJetHist_rat->SetYTitle ("Data / MC");
      vJetHist_rat->SetAxisRange (0.91, 1.09, "Y");
      vJetHist_rat->SetMarkerStyle (iMC == 1 ? dataStyle:mcSignalStyle);
      vJetHist_rat->SetMarkerColor (iMC == 1 ? dataColor:mcSignalColor);
      vJetHist_rat->SetLineColor (iMC == 1 ? dataColor:mcSignalColor);
      vJetHist_rat->GetYaxis ()->SetNdivisions (405);
      vJetHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
      vJetHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
      vJetHist_rat->GetXaxis ()->SetTitleOffset (1);
      vJetHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
      vJetHist_rat->GetYaxis ()->CenterTitle (true);
      vJetHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
      vJetHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
      vJetHist_rat->GetXaxis ()->SetTickLength (0.08);

      if ((skipOldInsitu || iYear == 0) && iMC == 1) vJetHist_rat->DrawCopy ("e1 x0");
      else vJetHist_rat->DrawCopy ("same e1 x0");
      ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
      for (TLine* line : glines) line->Draw ();
      //for (TLine* line : dplines_bottom) line->Draw ();

      if (vJetHist_rat->Integral () != 0.) {
       // add systematic errors
       for (int ix = 1; ix < vJetHist_rat->GetNbinsX (); ix++) {
        const double sys_err = max (vJetGraph_rat_sys->GetErrorYhigh (ix-1), vJetGraph_rat_sys->GetErrorYlow (ix-1));
        vJetHist_rat->SetBinError (ix, sqrt (pow (vJetHist_rat->GetBinError (ix), 2) + pow (sys_err, 2)));
       }

       const TString func = "[0] + [1]*((log(x)-[3])/[4]) + [2]*(2*((log(x)-[3])/[4])^2-1)";
       vJetRatioFit = new TF1 ("vJetPtDataMCRatioFit", func, pbins[0], pbins[numpbins]);
       vJetRatioFit->SetParameter (0, 1);
       vJetRatioFit->SetParameter (1, 0);
       vJetRatioFit->SetParameter (2, 0);
       vJetRatioFit->FixParameter (3, 0.5 * (log (pbins[numpbins]) + log (pbins[0])));
       vJetRatioFit->FixParameter (4, 0.5 * (log (pbins[numpbins]) - log (pbins[0])));
       vJetHist_rat->Fit (vJetRatioFit, "RN0Q");

       vJetRatioCI = new TH1D ("vJetPtDataMCRatioCI", "", numpbins, pbins);
       (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (vJetRatioCI);

       vJetRatioFit->SetLineColor (2);
       ( (TF1*)vJetRatioFit->Clone ())->Draw ("same");
       vJetRatioCI->SetMarkerStyle (kDot);
       vJetRatioCI->SetFillColorAlpha (2, 0.4);
       vJetRatioCI->DrawCopy ("e3 same");

       if (vJetRatioFit) { delete vJetRatioFit; vJetRatioFit = NULL; }
       if (vJetRatioCI) { delete vJetRatioCI; vJetRatioCI = NULL; }
      }

      if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
      if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
     }
     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetHist_mc_signal) { delete vJetHist_mc_signal; vJetHist_mc_signal = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_mc_signal) { delete proj_mc_signal; proj_mc_signal = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }
    } // end loop over insitu configurations

    if (iEta < numetabins) plotName = Form ("gamma_jet_iEta%i.pdf", iEta);
    else plotName = Form ("gamma_jet_iEta_combined.pdf");
    switch (iPer) {
     case 0:
      canvas->SaveAs (Form ("%s/PeriodA/%s", plotPath.Data (), plotName));
      break;
     case 1:
      canvas->SaveAs (Form ("%s/PeriodB/%s", plotPath.Data (), plotName));
      break;
     case 2:
      canvas->SaveAs (Form ("%s/PeriodAB/%s", plotPath.Data (), plotName));
      break;
    } // end switch
    //for (short iErr = 0; iErr < 3; iErr++)
    // if (gJetHistDifference[iPer][iEta][iErr])
    //  delete gJetHistDifference[iPer][iEta][iErr];

    if (!plot_xjref || iEta == numetabins) continue;

    /**** Plots xjref distributions, binned by ptref ****/
    for (short iP = 0; iP < numpbins; iP++) {
     const double pref_lo = pbins[iP];
     const double pref_hi =  pbins[iP+1];

     for (short iYear = 0; iYear < 2; iYear++) {
      if (skipOldInsitu && iYear == 0) continue;

      const Style_t dataStyle = (iYear == 0 ? 24 : 20);
      const Style_t mcOverlayStyle = (iYear == 0 ? 27 : 33);
      //const Style_t mcSignalStyle = (iYear == 0 ? 28 : 34);

      topPad->cd ();
      topPad->SetLogx (0);

      proj = Project2D ("", gJetHists[iYear][iPer][0][1][2], "x", "z", eta_lo, eta_hi);
      vJetHist = proj->ProjectionY ("vJetProjection", iP, iP);

      vJetHist->Rebin (rebinFactor);
      if (vJetHist->Integral () != 0) vJetHist->Scale (1./vJetHist->Integral ());
      vJetHist->SetXTitle ("#it{x}_{J}^{ref}");
      vJetHist->SetYTitle ("Counts / Total");
      vJetHist->SetMarkerStyle (dataStyle);
      vJetHist->SetMarkerColor (dataColor);
      vJetHist->SetLineColor (dataColor);
      //vJetHist->GetXaxis ()->SetRangeUser (0., 2.);
      vJetHist->GetYaxis ()->SetRangeUser (0., 0.6);//vJetHist->GetYaxis ()->GetXmax ());
      
      vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
      vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
      vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
      vJetHist->GetYaxis ()->SetTitleOffset (1.1*uPadY);

      proj_mc_overlay = Project2D ("", gJetHists[iYear][iPer][1][1][2], "x", "z", eta_lo, eta_hi);
      vJetHist_mc_overlay = proj_mc_overlay->ProjectionY ("vJetProjection_mc", iP, iP);

      vJetHist_mc_overlay->Rebin (rebinFactor);
      if (vJetHist_mc_overlay->Integral () != 0) vJetHist_mc_overlay->Scale (1./vJetHist_mc_overlay->Integral ()); 

      vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
      vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
      vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

      if (skipOldInsitu || iYear == 0) {
       vJetHist->DrawCopy ("e1 x0");
       vJetHist_mc_overlay->DrawCopy ("same e1 x0"); // insitu factors are not applied to MC
      }
      else {
       vJetHist->DrawCopy ("same e1 x0");
      }

      float mean, mean_err, mean_mc, mean_mc_err;
   
      if (useGaussian) {
       TF1* gaus_data = new TF1 ("gaus_data", "gaus (0)", 0, 4.0);
       vJetHist->Fit (gaus_data, "Q0R");
       TF1* gaus_mc = new TF1 ("gaus_mc", "gaus (0)", 0, 4.0);
       vJetHist_mc_overlay->Fit (gaus_mc, "Q0R");
       mean = gaus_data->GetParameter (1);
       mean_err = gaus_data->GetParError (1);
       mean_mc = gaus_mc->GetParameter (1);
       mean_mc_err = gaus_mc->GetParError (1);
       if (gaus_data) delete gaus_data;
       if (gaus_mc) delete gaus_mc;
      }
      else {
       mean = vJetHist->GetMean ();
       mean_err = vJetHist->GetMeanError ();
       mean_mc = vJetHist_mc_overlay->GetMean ();
       mean_mc_err = vJetHist_mc_overlay->GetMeanError ();
      }

      if (iYear == 0) {
       const int countData = gJetCounts[iYear][iPer][0][0][0]->GetBinContent (iP, iEta);
       const int countMC = gJetCounts[iYear][iPer][1][0][0]->GetBinContent (iP, iEta);

       myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countData), 1.25, 0.04/uPadY);
       myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV %s (%i events)", (runValidation ? "":"with #it{p}-Pb Overlay"), countMC), 1.25, 0.04/uPadY);

       myText (0.155, 0.73, dataColor, Form ("<#it{x}_{J}^{ref}>^{data} = %.2f #pm %.2f", mean, mean_err), 0.04/uPadY);
       myText (0.155, 0.64, dataColor, Form ("<#it{x}_{J}^{ref}>^{MC} = %.2f #pm %.2f", mean_mc, mean_mc_err), 0.04/uPadY);

       myText (0.155, 0.43, dataColor, "#gamma + Jet", 0.04/uPadY);
       myText (0.155, 0.34, dataColor, Form ("%g < #it{p}_{T}^{ref} < %g", pref_lo, pref_hi), 0.04/uPadY);
       myText (0.155, 0.25, dataColor, period.Data (), 0.04/uPadY);
       myText (0.155, 0.16, dataColor, Form ("%g < #eta_{det}^{Jet} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
       myText (0.155, 0.08, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      }

      bottomPad->cd ();
      bottomPad->SetLogx (false);
      vJetHist->Divide (vJetHist_mc_overlay);

      vJetHist->SetYTitle ("Data / MC");
      vJetHist->SetAxisRange (0.45, 1.65, "Y");
      vJetHist->GetYaxis ()->SetNdivisions (605);
      vJetHist->GetXaxis ()->SetTitleSize (0.04/dPadY);
      vJetHist->GetYaxis ()->SetTitleSize (0.04/dPadY);
      vJetHist->GetXaxis ()->SetTitleOffset (1);
      vJetHist->GetYaxis ()->SetTitleOffset (1.1*dPadY);
      vJetHist->GetYaxis ()->CenterTitle (true);
      vJetHist->GetXaxis ()->SetLabelSize (0.04/dPadY);
      vJetHist->GetYaxis ()->SetLabelSize (0.04/dPadY);
      vJetHist->GetXaxis ()->SetTickLength (0.08);

      if (skipOldInsitu || iYear == 0) vJetHist->DrawCopy ("e1 x0"); 
      else vJetHist->DrawCopy ("same e1 x0"); 
      for (TLine* line : xlines) line->Draw ();

      if (vJetHist) { delete vJetHist; vJetHist = NULL; }
      if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
      if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
      if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

      if (proj) { delete proj; proj = NULL; }
      if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
      if (proj_lo) { delete proj_lo; proj_lo = NULL; }
      if (proj_hi) { delete proj_hi; proj_hi = NULL; }
     } // end loop over insitu configurations

     plotName = Form ("pref_slices/gamma_jet_iEta%i_iP%i.pdf", iEta, iP);
     switch (iPer) {
      case 0:
       canvas->SaveAs (Form ("%s/PeriodA/%s", plotPath.Data (), plotName));
       break;
      case 1:
       canvas->SaveAs (Form ("%s/PeriodB/%s", plotPath.Data (), plotName));
       break;
      case 2:
       canvas->SaveAs (Form ("%s/PeriodAB/%s", plotPath.Data (), plotName));
       break;
     } // end switch
      
    } // end loop over pT bins
   } // end loop over etabins


   /**** Now loop over pT bins and plot response as function of eta^jet ****/
   for (short iP = 0; iP <= numpbins; iP++) {

    int p_lo = iP+1;
    int p_hi = iP+1;
    if (iP == numpbins) {
     p_lo = p_lo_comb;
     p_hi = p_hi_comb;
    }


    /**** Plots GammaJet info as a function of eta^jet****/
    for (short iYear = 0; iYear < 2; iYear++) {
     if (skipOldInsitu && iYear == 0) continue;

     const Style_t dataStyle = (iYear == 0 ? 24 : 20);
     const Style_t mcOverlayStyle = 33;
     const Style_t mcSignalStyle = 34;

     topPad->cd ();
     topPad->SetLogx (false);

     proj = Project2D ("", gJetHists[iYear][iPer][0][1][2], "y", "z", p_lo, p_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numetabins, etabins, false);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     vJetHist->SetYTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
     double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;
     vJetHist->SetAxisRange (middle - 0.35, middle + 0.35, "Y");
     vJetHist->SetMarkerStyle (dataStyle);
     vJetHist->SetMarkerColor (dataColor);
     vJetHist->SetLineColor (dataColor);
     vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", gJetHists[iYear][iPer][0][0][2], "y", "z", p_lo, p_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numetabins, etabins, false);

     proj_hi = Project2D ("", gJetHists[iYear][iPer][0][2][2], "y", "z", p_lo, p_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, false);

     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     vJetGraph_sys->SetFillColor (dataColor);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", gJetHists[iYear][iPer][1][1][2], "y", "z", p_lo, p_hi);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numetabins, etabins, false);

     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

     if (iPer != 1 && !skipSignalMC) {
      proj_mc_signal = Project2D ("", gJetHists[iYear][iPer][2][1][2], "y", "z", p_lo, p_hi);
      vJetHist_mc_signal = GetProfileX ("vJetHist_mc_signal", proj_mc_signal, numetabins, etabins, false);

      vJetHist_mc_signal->SetMarkerStyle (mcSignalStyle);
      vJetHist_mc_signal->SetMarkerColor (mcSignalColor);
      vJetHist_mc_signal->SetLineColor (mcSignalColor);
     }

     //for (short iErr = 0; iErr < 3; iErr++) {
     // const TString error = (iErr == 0 ? "syslo" : (iErr == 1 ? "stat" : "syshi"));
     // TString periodStr = "periodA";
     // if (iPer == 1) periodStr = "periodB";
     // else if (iPer == 2) periodStr = "periodAB";
     // gJetHistDifference[iPer][iP][iErr] = new TH1D (Form ("gJetPtRatio_diff%i_%s_%s", iP, error.Data (), periodStr.Data ()), ";#it{p}_{T}^{ref} #left[GeV#right]", numetabins, etabins);
     // for (short iP = 1; iP <= numpbins; iP++) {
     //  double dataVal, dataErr;
     //  switch (iErr) {
     //   case 0:
     //    dataVal = vJetHist_lo->GetBinContent (iP);
     //    dataErr = vJetHist_lo->GetBinError (iP);
     //    break;
     //   case 2:
     //    dataVal = vJetHist_hi->GetBinContent (iP);
     //    dataErr = vJetHist_hi->GetBinError (iP);
     //    break;
     //   default:
     //    dataVal = vJetHist->GetBinContent (iP);
     //    dataErr = vJetHist->GetBinError (iP);
     //  } 
     //  gJetHistDifference[iPer][iP][iErr]->SetBinContent (iP, dataVal - vJetHist_mc_overlay->GetBinContent (iP));
     //  gJetHistDifference[iPer][iP][iErr]->SetBinError (iP, TMath::Sqrt (TMath::Power (dataErr,2) + TMath::Power (vJetHist_mc_overlay->GetBinError (iP),2)));
     // }
     // gJetHistDifference[iPer][iP][iErr]->Write ();
     //}
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     if (skipOldInsitu || iYear == 0) {
      vJetHist->DrawCopy ("e1 x0");
      vJetHist_mc_overlay->DrawCopy ("same e1 x0"); // insitu factors are not applied to MC
      if (iPer != 1 && !skipSignalMC) vJetHist_mc_signal->DrawCopy ("same e1 x0");
     }
     else {
      vJetHist->DrawCopy ("same e1 x0");
     }
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");
     //for (TLine* line : dplines) line->Draw ();

     if (skipOldInsitu || iYear == 0) {
      int countsData = 0, countsMC = 0, countsMC_sig = 0;
      countsData = gJetCounts[iYear][iPer][0][0][0]->Integral (p_lo, p_hi, 1, gJetCounts[iYear][iPer][0][0][0]->GetNbinsY());
      countsMC = gJetCounts[iYear][iPer][1][0][0]->Integral (p_lo, p_hi, 1, gJetCounts[iYear][iPer][1][0][0]->GetNbinsY());
      if (!skipSignalMC)
       countsMC_sig = gJetCounts[iYear][iPer][2][0][0]->Integral (p_lo, p_hi, 1, gJetCounts[iYear][iPer][2][0][0]->GetNbinsY());

      myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countsData), 1.25, 0.04/uPadY);
      myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
      if (iPer != 1 && !skipSignalMC)
       myMarkerText (0.175, 0.74, mcSignalColor, kFullCross, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", countsMC_sig), 1.25, 0.04/uPadY);
      if (!skipOldInsitu) {
       myMarkerText (0.65, 0.66, dataColor, kFullCircle, "2016 Insitu Factors", 1.25, 0.04/uPadY);
       myMarkerText (0.65, 0.59, dataColor, kOpenCircle, "2015 Insitu Factors", 1.25, 0.04/uPadY);
      }
      myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      if (p_lo != 1 || p_hi != numpbins)
       myText (0.155, 0.15, dataColor, Form ("#gamma + Jet, %g < #it{p}_{T}^{#gamma} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
      else
       myText (0.155, 0.15, dataColor, "#gamma + Jet", 0.04/uPadY);
      myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);
     }

     bottomPad->cd ();
     bottomPad->SetLogx (false);

     for (int iMC = 1; iMC < 3; iMC++) { // loops over overlay then signal MC
      if (iMC == 2 && (skipSignalMC || iPer == 1)) continue; // no signal only for period B!

      const TString mcType = (iMC == 1 ? "overlay":"signal");
      TH2D* proj_mc = (iMC == 1 ? proj_mc_overlay:proj_mc_signal);

      vJetHist_rat = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_iP%i", iP)), proj, proj_mc, numetabins, etabins, false, "x");
      vJetHist_rat_lo = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_lo_iP%i", iP)), proj_lo, proj_mc, numetabins, etabins, false, "x");
      vJetHist_rat_hi = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_hi_iP%i", iP)), proj_hi, proj_mc, numetabins, etabins, false, "x");

      vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
      CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
      if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
      if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
      //vJetGraph_rat_sys->SetFillColor (iMC == 1 ? dataColor:mcSignalColor);
      vJetGraph_rat_sys->SetFillColor (dataColor);
      vJetGraph_rat_sys->SetFillStyle (3001);

      vJetHist_rat->SetXTitle ("Jet #eta");
      vJetHist_rat->SetYTitle ("Data / MC");
      vJetHist_rat->SetAxisRange (0.91, 1.09, "Y");
      vJetHist_rat->SetMarkerStyle (iMC == 1 ? dataStyle:mcSignalStyle);
      vJetHist_rat->SetMarkerColor (iMC == 1 ? dataColor:mcSignalColor);
      vJetHist_rat->SetLineColor (iMC == 1 ? dataColor:mcSignalColor);
      vJetHist_rat->GetYaxis ()->SetNdivisions (405);
      vJetHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
      vJetHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
      vJetHist_rat->GetXaxis ()->SetTitleOffset (1);
      vJetHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
      vJetHist_rat->GetYaxis ()->CenterTitle (true);
      vJetHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
      vJetHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
      vJetHist_rat->GetXaxis ()->SetTickLength (0.08);

      if ((skipOldInsitu || iYear == 0) && iMC == 1) vJetHist_rat->DrawCopy ("e1 x0");
      else vJetHist_rat->DrawCopy ("same e1 x0");
      ((TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
      for (TLine* line : getalines) line->Draw ();
      //for (TLine* line : dplines_bottom) line->Draw ();

      if (vJetHist_rat->Integral () != 0.) {
       // add systematic errors
       for (int ix = 1; ix < vJetHist_rat->GetNbinsX (); ix++) {
        const double sys_err = max (vJetGraph_rat_sys->GetErrorYhigh (ix-1), vJetGraph_rat_sys->GetErrorYlow (ix-1));
        vJetHist_rat->SetBinError (ix, sqrt (pow (vJetHist_rat->GetBinError (ix), 2) + pow (sys_err, 2)));
       }

       const TString func = "[0] + [1]*((x-[3])/[4]) + [2]*(2*((x-[3])/[4])^2-1)";
       TF1* vJetRatioFit = new TF1 ("vJetPtDataMCRatioFit", func, etabins[0], etabins[numetabins]);
       vJetRatioFit->SetParameter (0, 1);
       vJetRatioFit->SetParameter (1, 0);
       vJetRatioFit->SetParameter (2, 0);
       vJetRatioFit->FixParameter (3, 0.5 * (etabins[numetabins] + etabins[0]));
       vJetRatioFit->FixParameter (4, 0.5 * (etabins[numetabins] - etabins[0]));
       vJetHist_rat->Fit (vJetRatioFit, "RN0Q");

       vJetRatioCI = new TH1D ("vJetPtDataMCRatioCI", "", numetabins, etabins);
       (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (vJetRatioCI);

       vJetRatioFit->SetLineColor (2);
       ( (TF1*)vJetRatioFit->Clone ())->Draw ("same");
       vJetRatioCI->SetMarkerStyle (kDot);
       vJetRatioCI->SetFillColorAlpha (2, 0.4);
       vJetRatioCI->DrawCopy ("e3 same");

       if (vJetRatioFit) { delete vJetRatioFit; vJetRatioFit = NULL; }
       if (vJetRatioCI) { delete vJetRatioCI; vJetRatioCI = NULL; }
      }

      if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
      if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
     } // end loop over MC types
     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetHist_mc_signal) { delete vJetHist_mc_signal; vJetHist_mc_signal = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_mc_signal) { delete proj_mc_signal; proj_mc_signal = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     //for (short iErr = 0; iErr < 3; iErr++)
     // if (gJetHistDifference[iPer][iP][iErr])
     //  delete gJetHistDifference[iPer][iP][iErr];
    } // end loop over insitu configurations

    if (iP < numpbins) plotName = Form ("gamma_jet_iP%i.pdf", iP);
    else plotName = Form ("gamma_jet_iP_combined.pdf");
    switch (iPer) {
     case 0:
      canvas->SaveAs (Form ("%s/PeriodA/%s", plotPath.Data (), plotName));
      break;
     case 1:
      canvas->SaveAs (Form ("%s/PeriodB/%s", plotPath.Data (), plotName));
      break;
     case 2:
      canvas->SaveAs (Form ("%s/PeriodAB/%s", plotPath.Data (), plotName));
      break;
    } // end switch

   } // end loop over pT bins
  } // end loop over beam periods


//  /**** Plots systematic errors vs jet pt ****/
//  for (short iPer = 0; iPer < 2; iPer++) {
//   const TString period = (iPer == 0 ? "Period A":"Period B");
//
//   for (short iEta = 0; iEta < numetabins; iEta++) {
//    topPad->cd ();
//    TH2D* thisHist = gJetHistsSys[iPer][iEta][0][1];
//    TH1D* rmsHist = new TH1D (Form ("rms_iEta%i_%s", iEta, (iPer==0?"pPb":"Pbp")), "", numpbins, pbins);
//    for (short pbin = 0; pbin < numpbins; pbin++) {
//     float rms = 0;
//     float sumWeights = 0;
//     for (short sigbin = 0; sigbin < numSigmaBins; sigbin++) {
//      const float sig = thisHist->GetYaxis ()->GetBinCenter (sigbin+1);
//      const float weight = thisHist->GetBinContent (pbin+1, sigbin+1);
//      rms += pow (sig, 2) * weight;
//      sumWeights += weight;
//     }
//     if (sumWeights > 0) rms = sqrt (rms) / sqrt (sumWeights);
//     rmsHist->SetBinContent (pbin+1, rms);
//    }
//    topPad->SetLogz ();
//    thisHist->Draw ("col");
//    thisHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
//    thisHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
//    thisHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
//    thisHist->GetYaxis ()->SetTitleOffset (1.1*uPadY);
//    myText (0.72, 0.89, dataColor, period.Data (), 0.04/uPadY);
//    myText (0.72, 0.8, dataColor, Form ("%g < #eta_{det}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
//
//    bottomPad->cd ();
//    rmsHist->SetXTitle ("#it{p}_{T}^{jet} #left[GeV#right]");
//    //rmsHist->SetYTitle ("RMS (#sigma)/#it{p}_{T}^{jet}");
//    rmsHist->SetYTitle ("RMS");
//    rmsHist->SetAxisRange (0, 0.09, "Y");
//    rmsHist->GetYaxis ()->SetNdivisions (405);
//    rmsHist->GetXaxis ()->SetTitleSize (0.04/dPadY);
//    rmsHist->GetYaxis ()->SetTitleSize (0.04/dPadY);
//    rmsHist->GetXaxis ()->SetTitleOffset (1);
//    rmsHist->GetYaxis ()->SetTitleOffset (1.1*dPadY);
//    rmsHist->GetYaxis ()->CenterTitle (true);
//    rmsHist->GetXaxis ()->SetLabelSize (0.04/dPadY);
//    rmsHist->GetYaxis ()->SetLabelSize (0.04/dPadY);
//    rmsHist->GetXaxis ()->SetTickLength (0.08);
//    rmsHist->Draw ("hist");
//    canvas->SaveAs (Form ("%s/Period%s/jetSystematics_iEta%i.pdf", plotPath.Data (), (iPer==0 ? "A":"B"), iEta));
//    if (rmsHist) delete rmsHist;
//   }
//  }
//  outFile->Write ();
//  if (outFile) delete outFile;


//  /**** Plot 2d histograms ****/
//  TCanvas* th2canvas = new TCanvas ("th2canvas", "", 800, 600);
//  const double padRatio_th2 = 1; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
//  const double rPadX = 1.0/ (padRatio_th2+1.0);
//  const double lPadX = 1.0 - rPadX;
//  TPad* leftPad = new TPad ("leftPad", "", 0, 0, lPadX, 1);
//  TPad* rightPad = new TPad ("rightPad", "", lPadX, 0, 1, 1);
//  leftPad->SetRightMargin (0);
//  rightPad->SetLeftMargin (0);
//  leftPad->SetTopMargin (-0.05);
//  rightPad->SetTopMargin (-0.05);
//  leftPad->SetBottomMargin (-0.10);
//  rightPad->SetBottomMargin (-0.10);
//  leftPad->SetLeftMargin (-0.20);
//  rightPad->SetRightMargin (-0.10);
//  leftPad->Draw ();
//  rightPad->Draw ();
//  for (short iPer = 0; iPer < 3; iPer++) {
//   const TString period = (iPer == 0 ? "Period A": (iPer == 1 ? "Period B":"Period A+B"));
//
//   for (short iEta = 0; iEta < numetabins; iEta++) {
//    TH2D* dataHist = Project2D ("", gJetHists[iYear][iPer][0][1], "x", "z", iEta, iEta);
//    TH2D* mcHist = Project2D ("", gJetHists[iYear][iPer][1][1], "x", "z", iEta, iEta);
//
//    dataHist->Scale (1./dataHist->Integral ());
//    mcHist->Scale (1./mcHist->Integral ());
//
//    rightPad->cd ();
//    rightPad->SetLogx ();
//    rightPad->SetLogz ();
//    dataHist->GetXaxis ()->SetTitleSize (0.02/rPadX);
//    dataHist->GetYaxis ()->SetTitleSize (0.02/rPadX);
//    dataHist->GetXaxis ()->SetTitleOffset (1);
//    dataHist->GetYaxis ()->SetTitleOffset (1);
//    dataHist->GetXaxis ()->SetLabelSize (0.02/rPadX);
//    dataHist->GetYaxis ()->SetLabelSize (0.02/rPadX);
//    dataHist->Draw ("col");
//    myText (0.1, 0.15, dataColor, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", nGammaJet[iYear][iPer][0][iEta][numpbins]), 0.02/rPadX);
//    if (iPer != 2) myText (0.6, 0.85, dataColor, Form ("%g < #eta_{det}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.02/rPadX);
//    else myText (0.6, 0.85, dataColor, Form ("%g < #eta_{det}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.02/rPadX);
//    myText (0.6, 0.8, dataColor, period.Data (), 0.02/rPadX);
//
//    leftPad->cd ();
//    leftPad->SetLogx ();
//    leftPad->SetLogz ();
//    mcHist->GetXaxis ()->SetTitleSize (0.02/lPadX);
//    mcHist->GetYaxis ()->SetTitleSize (0.02/lPadX);
//    mcHist->GetXaxis ()->SetTitleOffset (1);
//    mcHist->GetYaxis ()->SetTitleOffset (1);
//    mcHist->GetXaxis ()->SetLabelSize (0.02/lPadX);
//    mcHist->GetYaxis ()->SetLabelSize (0.02/lPadX);
//    mcHist->Draw ("col");
//    myText (0.2, 0.15, dataColor, Form ("Pythia8 #it{pp} 8.16 TeV %s (%i events)", (runValidation ? "":"with #it{p}-Pb Overlay"), nGammaJet[iYear][iPer][1][iEta][numpbins]), 0.02/lPadX);
//
//    plotName = Form ("gamma_jet%i_th2.pdf", iEta);
//    switch (iPer) {
//     case 0:
//      th2canvas->SaveAs (Form ("%s/PeriodA/%s", plotPath.Data (), plotName));
//      break;
//     case 1:
//      th2canvas->SaveAs (Form ("%s/PeriodB/%s", plotPath.Data (), plotName));
//      break;
//     case 2:
//      th2canvas->SaveAs (Form ("%s/PeriodAB/%s", plotPath.Data (), plotName));
//      break;
//    }
//
//    if (dataHist) { delete dataHist; dataHist = NULL; }
//    if (mcHist) { delete mcHist; mcHist = NULL; }
//   }
//  }


  return;
}

} // end namespace
