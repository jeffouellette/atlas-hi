#include "ZeeJetsHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TVirtualFitter.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TVectorT.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

namespace JetCalibration {

void ZeeJetsHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("ZeeJets/", "JetCalibration/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TH3D***** zeeJetHists = Get4DArray <TH3D*> (2, 3, 2, 3); // iYear, iPer, iData, iErr
  TH2D**** zeeJetCounts = Get3DArray <TH2D*> (2, 3, 2);

  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0)
    continue; // ignore old insitu factors if not desired

   const TString year = (iYear == 0 ? "2015" : "2016");

   for (short iPer = 0; iPer < 3; iPer++) {
    const TString period = (iPer == 0 ? "periodA" : (iPer == 1 ? "periodB" : "periodAB"));

    for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
     const TString dataType = (iData == 0 ? "data" : "mc_overlay");

     for (short iErr = 0; iErr < 3; iErr++) {
      if (iErr != 1 && iData != 0)
       continue;

      const TString error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

      const TString name = Form ("zeeJetPtRatio_%s_%s_%s_%s", year.Data (), dataType.Data (), error.Data (), period.Data ());
      zeeJetHists[iYear][iPer][iData][iErr] = (TH3D*)inFile->Get (name);
     }

     const TString name = Form ("zeeJetCounts_%s_%s_%s", year.Data (), dataType.Data (), period.Data ());
     zeeJetCounts[iYear][iPer][iData] = (TH2D*)inFile->Get (name);
    }
   }
  }


  TLine* zlines[5] = {};
  TLine* zetalines[5] = {};
  for (short i = 0; i < 5; i++) {
   const float dz = 0.05;

   zlines[i] = new TLine (pbins[0], 1.0-2*dz+dz*i, pbins[numpbins], 1.0-2*dz+dz*i);
   zetalines[i] = new TLine (etabins[0], 1.0-2*dz+dz*i, etabins[numetabins], 1.0-2*dz+dz*i);

   if (1.0-2*dz+dz*i == 1) zlines[i]->SetLineStyle (1);
   else zlines[i]->SetLineStyle (3);
   if (1.0-2*dz+dz*i == 1) zetalines[i]->SetLineStyle (1);
   else zetalines[i]->SetLineStyle (3);
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
  TH1D *vJetHist = NULL, *vJetHist_mc_overlay = NULL, *vJetHist_lo = NULL, *vJetHist_hi = NULL, *vJetHist_rat = NULL, *vJetHist_rat_lo = NULL, *vJetHist_rat_hi = NULL, *vJetRatioCI = NULL;
  TH2D *proj = NULL, *proj_mc_overlay = NULL, *proj_lo = NULL, *proj_hi = NULL;
  TF1* vJetRatioFit = NULL;
  TGraphAsymmErrors *vJetGraph_sys = NULL, *vJetGraph_rat_sys = NULL;
  char* plotName;

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
   TString period = "Period A";
   if (iPer == 1) period = "Period B";
   else if (iPer == 2) period = "Period A+B";

   TString perType = "pA";
   if (iPer == 1) perType = "pB";
   else if (iPer == 2) perType = "pAB";

   for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta

    const int eta_lo = (iEta != numetabins ? iEta+1 : eta_lo_comb);
    const int eta_hi = (iEta != numetabins ? iEta+1 : eta_hi_comb);

    /**** Plots ZeeJet info ****/
    for (short iYear = 0; iYear < 2; iYear++) {
     if (skipOldInsitu && iYear == 0) continue;

     const Style_t dataStyle = (iYear == 0 ? 24 : 20);
     const Style_t mcOverlayStyle = 33;
     //const Style_t mcSignalStyle = 34;

     topPad->cd ();
     topPad->SetLogx ();

     proj = Project2D ("", zeeJetHists[1][iPer][0][1], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist = GetProfileX ("vJetHist", proj, numpbins, pbins, true);

     vJetGraph_sys = make_graph (vJetHist); // for plotting systematics
     vJetHist->SetYTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
     vJetHist->SetMarkerStyle (dataStyle);
     vJetHist->SetMarkerColor (kBlack);
     vJetHist->SetLineColor (kBlack);
     vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", zeeJetHists[1][iPer][0][0], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numpbins, pbins, true);

     proj_hi = Project2D ("", zeeJetHists[1][iPer][0][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numpbins, pbins, true);

     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }
     vJetGraph_sys->SetFillColor (kBlack);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", zeeJetHists[1][iPer][1][1], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numpbins, pbins, true);

     // set y axis range according to MC
     double middle = 0.05 * floor (20 * vJetHist_mc_overlay->Integral () / vJetHist_mc_overlay->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;
     vJetHist->SetAxisRange (middle - 0.35, middle + 0.35, "Y");

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
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     int countsData = 0, countsMC = 0;
     if (exclusive && iEta == numetabins) {
      countsData = zeeJetCounts[1][iPer][0]->Integral () - zeeJetCounts[1][iPer][0]->Integral (1, zeeJetCounts[1][iPer][0]->GetNbinsX(), eta_lo, eta_hi);
      countsMC = zeeJetCounts[1][iPer][1]->Integral () - zeeJetCounts[1][iPer][1]->Integral (1, zeeJetCounts[1][iPer][1]->GetNbinsX(), eta_lo, eta_hi);
     }
     else {
      countsData = zeeJetCounts[1][iPer][0]->Integral (1, zeeJetCounts[1][iPer][0]->GetNbinsX(), eta_lo, eta_hi);
      countsMC = zeeJetCounts[1][iPer][1]->Integral (1, zeeJetCounts[1][iPer][1]->GetNbinsX(), eta_lo, eta_hi);
     }

     myMarkerText (0.175, 0.88, kBlack, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countsData), 1.25, 0.04/uPadY);
     myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
     myText (0.155, 0.22, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
     if (eta_lo != 1 || eta_hi != numetabins)
      if (!(exclusive && iEta == numetabins))
       myText (0.155, 0.15, kBlack, Form ("Z (ee) + Jet, %g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04/uPadY);
      else
       myText (0.155, 0.15, kBlack, Form ("Z (ee) + Jet, %g < #left|#eta_{det}^{Jet}#right|", etabins[eta_hi]), 0.04/uPadY);
     else
      myText (0.155, 0.15, kBlack, "Z (ee) + Jet", 0.04/uPadY);
     myText (0.155, 0.08, kBlack, period.Data (), 0.04/uPadY);

     bottomPad->cd ();
     bottomPad->SetLogx ();

     vJetHist_rat = GetDataOverMC (TString (Form ("zeeJetPtDataMCRatio_iEta%i", iEta)), proj, proj_mc_overlay, numpbins, pbins, true, "x");
     vJetHist_rat_lo = GetDataOverMC (TString (Form ("zeeJetPtDataMCRatio_lo_iEta%i", iEta)), proj_lo, proj_mc_overlay, numpbins, pbins, true, "x");
     vJetHist_rat_hi = GetDataOverMC (TString (Form ("zeeJetPtDataMCRatio_hi_iEta%i", iEta)), proj_hi, proj_mc_overlay, numpbins, pbins, true, "x");

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     vJetGraph_rat_sys = make_graph (vJetHist_rat);
     CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
     if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
     vJetGraph_rat_sys->SetFillColor (kBlack);
     vJetGraph_rat_sys->SetFillStyle (3001);

     vJetHist_rat->SetXTitle ("#it{p}_{T}^{Z} #left[GeV#right]");
     vJetHist_rat->SetYTitle ("Data / MC");
     vJetHist_rat->SetAxisRange (0.85, 1.15, "Y");
     //vJetHist_rat->SetAxisRange (0.75, 1.35, "Y");
     vJetHist_rat->SetMarkerStyle (dataStyle);
     vJetHist_rat->SetMarkerColor (kBlack);
     vJetHist_rat->SetLineColor (kBlack);
     vJetHist_rat->GetYaxis ()->SetNdivisions (405);
     vJetHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTitleOffset (1);
     vJetHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
     vJetHist_rat->GetYaxis ()->CenterTitle (true);
     vJetHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTickLength (0.08);

     if (skipOldInsitu || iYear == 0) vJetHist_rat->DrawCopy ("e1 x0");
     else vJetHist_rat->DrawCopy ("same e1 x0");
     ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     for (TLine* line : zlines) line->Draw ();

     //if (vJetHist_rat->Integral () != 0.) {
     // // add systematic errors
     // for (int ix = 1; ix < vJetHist_rat->GetNbinsX (); ix++) {
     //  const double sys_err = max (vJetGraph_rat_sys->GetErrorYhigh (ix-1), vJetGraph_rat_sys->GetErrorYlow (ix-1));
     //  vJetHist_rat->SetBinError (ix, sqrt (pow (vJetHist_rat->GetBinError (ix), 2) + pow (sys_err, 2)));
     // }

     // const TString func = "[0] + [1]*((log(x)-[3])/[4]) + [2]*(2*((log(x)-[3])/[4])^2-1)";
     // vJetRatioFit = new TF1 ("vJetPtDataMCRatioFit", func, pbins[0], pbins[numpbins]);
     // vJetRatioFit->SetParameter (0, 1);
     // vJetRatioFit->SetParameter (1, 0);
     // vJetRatioFit->SetParameter (2, 0);
     // vJetRatioFit->FixParameter (3, 0.5 * (log (pbins[numpbins]) + log (pbins[0])));
     // vJetRatioFit->FixParameter (4, 0.5 * (log (pbins[numpbins]) - log (pbins[0])));
     // vJetHist_rat->Fit (vJetRatioFit, "RN0Q");

     // vJetRatioCI = new TH1D ("vJetPtDataMCRatioCI", "", numpbins, pbins);
     // (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (vJetRatioCI);

     // vJetRatioFit->SetLineColor (2);
     // ( (TF1*)vJetRatioFit->Clone ())->Draw ("same");
     // vJetRatioCI->SetMarkerStyle (kDot);
     // vJetRatioCI->SetFillColorAlpha (2, 0.4);
     // vJetRatioCI->DrawCopy ("e3 same");

     // if (vJetRatioFit) { delete vJetRatioFit; vJetRatioFit = NULL; }
     // if (vJetRatioCI) { delete vJetRatioCI; vJetRatioCI = NULL; }
     //}

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
    } // end loop of insitu configurations
     

    if (iEta < numetabins) plotName = Form ("z_ee_jet_iEta%i.pdf", iEta);
    else plotName = Form ("z_ee_jet_iEta_combined.pdf");

    canvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));
   } // end loop over etabins

   /**** Now loop over pT bins and plot response as function of eta^jet ****/
   for (short iP = 0; iP <= numpbins; iP++) {

    int p_lo = iP+1;
    int p_hi = iP+1;
    if (iP == numpbins) {
     p_lo = 1;
     p_hi = 7;//numpbins;
    }

    /**** Plots ZeeJet info ****/
    for (short iYear = 0; iYear < 2; iYear++) {
     if (skipOldInsitu && iYear == 0) continue;

     const Style_t dataStyle = (iYear == 0 ? 24 : 20);
     const Style_t mcOverlayStyle = 33;
     //const Style_t mcSignalStyle = 34;

     topPad->cd ();
     topPad->SetLogx (false);

     proj = Project2D ("", zeeJetHists[1][iPer][0][1], "y", "z", p_lo, p_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numetabins, etabins, true);

     vJetGraph_sys = make_graph (vJetHist); // for plotting systematics
     vJetHist->SetYTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
     vJetHist->SetMarkerStyle (dataStyle);
     vJetHist->SetMarkerColor (kBlack);
     vJetHist->SetLineColor (kBlack);
     vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", zeeJetHists[1][iPer][0][0], "y", "z", p_lo, p_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numetabins, etabins, true);

     proj_hi = Project2D ("", zeeJetHists[1][iPer][0][2], "y", "z", p_lo, p_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, true);

     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }
     vJetGraph_sys->SetFillColor (kBlack);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", zeeJetHists[1][iPer][1][1], "y", "z", p_lo, p_hi);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numetabins, etabins, true);

     // set y axis range according to MC
     double middle = 0.05 * floor (20 * vJetHist_mc_overlay->Integral () / vJetHist_mc_overlay->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;
     vJetHist->SetAxisRange (middle - 0.35, middle + 0.35, "Y");

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
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     const int countsData = zeeJetCounts[1][iPer][0]->Integral (p_lo, p_hi, 1, zeeJetCounts[1][iPer][0]->GetNbinsY());
     const int countsMC = zeeJetCounts[1][iPer][1]->Integral (p_lo, p_hi, 1, zeeJetCounts[1][iPer][1]->GetNbinsY());

     myMarkerText (0.175, 0.88, kBlack, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countsData), 1.25, 0.04/uPadY);
     myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
     myText (0.155, 0.22, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
     if (p_lo != 1 || p_hi != numpbins)
      myText (0.155, 0.15, kBlack, Form ("Z (ee) + Jet, %g < #it{p}_{T}^{Z} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
     else
      myText (0.155, 0.15, kBlack, "Z (ee) + Jet", 0.04/uPadY);
     myText (0.155, 0.08, kBlack, period.Data (), 0.04/uPadY);

     bottomPad->cd ();
     bottomPad->SetLogx (false);

     vJetHist_rat = GetDataOverMC (TString (Form ("zeeJetPtDataMCRatio_iP%i", iP)), proj, proj_mc_overlay, numetabins, etabins, true, "x");
     vJetHist_rat_lo = GetDataOverMC (TString (Form ("zeeJetPtDataMCRatio_lo_iP%i", iP)), proj_lo, proj_mc_overlay, numetabins, etabins, true, "x");
     vJetHist_rat_hi = GetDataOverMC (TString (Form ("zeeJetPtDataMCRatio_hi_iP%i", iP)), proj_hi, proj_mc_overlay, numetabins, etabins, true, "x");

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     vJetGraph_rat_sys = make_graph (vJetHist_rat);
     CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
     if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
     vJetGraph_rat_sys->SetFillColor (kBlack);
     vJetGraph_rat_sys->SetFillStyle (3001);

     vJetHist_rat->SetXTitle ("Jet #eta");
     vJetHist_rat->SetYTitle ("Data / MC");
     vJetHist_rat->SetAxisRange (0.85, 1.15, "Y");
     //vJetHist_rat->SetAxisRange (0.75, 1.35, "Y");
     vJetHist_rat->SetMarkerStyle (dataStyle);
     vJetHist_rat->SetMarkerColor (kBlack);
     vJetHist_rat->SetLineColor (kBlack);
     vJetHist_rat->GetYaxis ()->SetNdivisions (405);
     vJetHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTitleOffset (1);
     vJetHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
     vJetHist_rat->GetYaxis ()->CenterTitle (true);
     vJetHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTickLength (0.08);

     if (skipOldInsitu || iYear == 0) vJetHist_rat->DrawCopy ("e1 x0");
     else vJetHist_rat->DrawCopy ("same e1 x0");
     ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     for (TLine* line : zetalines) line->Draw ();

     //if (vJetHist_rat->Integral () != 0.) {
     // // add systematic errors
     // for (int ix = 1; ix < vJetHist_rat->GetNbinsX (); ix++) {
     //  const double sys_err = max (vJetGraph_rat_sys->GetErrorYhigh (ix-1), vJetGraph_rat_sys->GetErrorYlow (ix-1));
     //  vJetHist_rat->SetBinError (ix, sqrt (pow (vJetHist_rat->GetBinError (ix), 2) + pow (sys_err, 2)));
     // }

     // const TString func = "[0] + [1]*((x-[3])/[4]) + [2]*(2*((x-[3])/[4])^2-1)";
     // TF1* vJetRatioFit = new TF1 ("vJetPtDataMCRatioFit", func, etabins[0], etabins[numetabins]);
     // vJetRatioFit->SetParameter (0, 1);
     // vJetRatioFit->SetParameter (1, 0);
     // vJetRatioFit->SetParameter (2, 0);
     // vJetRatioFit->FixParameter (3, 0.5 * (etabins[numetabins] + etabins[0]));
     // vJetRatioFit->FixParameter (4, 0.5 * (etabins[numetabins] - etabins[0]));
     // vJetHist_rat->Fit (vJetRatioFit, "RN0Q");

     // vJetRatioCI = new TH1D ("vJetPtDataMCRatioCI", "", numetabins, etabins);
     // (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (vJetRatioCI);

     // vJetRatioFit->SetLineColor (2);
     // ( (TF1*)vJetRatioFit->Clone ())->Draw ("same");
     // vJetRatioCI->SetMarkerStyle (kDot);
     // vJetRatioCI->SetFillColorAlpha (2, 0.4);
     // vJetRatioCI->DrawCopy ("e3 same");

     // if (vJetRatioFit) { delete vJetRatioFit; vJetRatioFit = NULL; }
     // if (vJetRatioCI) { delete vJetRatioCI; vJetRatioCI = NULL; }
     //}

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
    } // end loop of insitu configurations

    if (iP < numpbins) plotName = Form ("z_ee_jet_iP%i.pdf", iP);
    else plotName = Form ("z_ee_jet_iP_combined.pdf");
    canvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));

   } // end loop over pT bins
  } // end loop over beam periods

  return;
}

} // end namespace
