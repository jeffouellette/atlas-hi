#include "ZGammaJetCrossCheckHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <GlobalParams.h>
#include <ArrayTemplates.h>

#include <TVirtualFitter.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TVectorT.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

namespace JetCalibration {

double min (TGraph* g) {
  double min = 1e300, x, y;
  for (int ix = 0; ix < g->GetN (); ix++) {
   g->GetPoint (ix, x, y);
   if (y < min) min = y;
  }
  return min;
}

double max (TGraph* g) {
  double max = -1e300, x, y;
  for (int ix = 0; ix < g->GetN (); ix++) {
   g->GetPoint (ix, x, y);
   if (y > max) max = y;
  }
  return max;
}


void AddWeighted (TGraphAsymmErrors** num, TGraphAsymmErrors** den, const TGraphAsymmErrors* n, const TGraphAsymmErrors* nsys_lo, const TGraphAsymmErrors* nsys_hi, const TGraphAsymmErrors* data) {
  double x, y, nx, ny, ei, ne, w;
  for (int ix = 0; ix < n->GetN (); ix++) {
   n->GetPoint (ix, nx, ny);

   if (n->GetErrorY (ix) <= 0 || data->GetErrorY (ix) == 0)
    continue;

   w = 1. / pow (n->GetErrorY (ix), 2);

   nsys_lo->GetPoint (ix, nx, ny);
   num[0]->GetPoint (ix, x, y);
   num[0]->SetPoint (ix, nx, y + w*ny);

   den[0]->GetPoint (ix, x, y);
   den[0]->SetPoint (ix, nx, y + w); // set new denominator

   ne = n->GetErrorY (ix);
   ei = num[1]->GetErrorY (ix);
   n->GetPoint (ix, nx, ny);
   num[1]->GetPoint (ix, x, y);
   num[1]->SetPoint (ix, nx, y + w*ny); // set new numerator
   num[1]->SetPointError (ix, n->GetErrorXlow (ix), n->GetErrorXhigh (ix), sqrt (ei*ei + w*w*ne*ne), sqrt (ei*ei + w*w*ne*ne)); // set new error in numerator
 
   den[1]->GetPoint (ix, x, y);
   den[1]->SetPoint (ix, nx, y + 1./(ne*ne)); // set new denominator
   den[1]->SetPointError (ix, n->GetErrorXlow (ix), n->GetErrorXhigh (ix), 0, 0); // set error in denominator

   nsys_hi->GetPoint (ix, nx, ny);
   num[2]->GetPoint (ix, x, y);
   num[2]->SetPoint (ix, nx, y + w*ny);

   den[2]->GetPoint (ix, x, y);
   den[2]->SetPoint (ix, nx, y + w); // set new denominator
  }
  return;
}


TGraphAsymmErrors* Divide (const TGraphAsymmErrors* num, const TGraphAsymmErrors* den) {
  TGraphAsymmErrors* rat = new TGraphAsymmErrors ();

  double nx, ny, dx, dy, ne, de, rx, ry, re;
  for (int ix = 0; ix < num->GetN (); ix++) {
   num->GetPoint (ix, nx, ny);
   den->GetPoint (ix, dx, dy);
   ne = num->GetErrorY (ix);
   de = den->GetErrorY (ix);

   rx = nx;
   if (dy > 0) {
    ry = ny / dy;
    re = sqrt (pow (ne / dy, 2) + pow (ry * de / dy, 2));
   }
   else {
    ry = 0;
    re = 0;
   }

   rat->SetPoint (ix, rx, ry);
   rat->SetPointError (ix, num->GetErrorXlow (ix), num->GetErrorXhigh (ix), re, re);
  }
  return rat;
}


void ZGammaJetCrossCheckHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("ZGammaJetCrossCheck/", "JetCalibration/");

  TH3D***** vJetHists = Get4DArray <TH3D*> (3, 3, 3, 4);
  TH2D**** vJetCounts = Get3DArray <TH2D*> (3, 3, 4);

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  for (short iData = 0; iData < 3; iData++) { // iData is 0 for data, 1 for MC
   if (iData == 2 && skipSignalMC)
    continue; // ignore signal MC if not desired
   const TString dataType = (iData == 0 ? "data" : (iData == 1 ? "mc_overlay" : "mc_signalonly"));

   for (short iPer = 0; iPer < 3; iPer++) {
    const TString period = (iPer == 0 ? "periodA" : (iPer == 1 ? "periodB" : "periodAB"));

    for (short iSpc = 0; iSpc < 4; iSpc++) {
     const TString spc = (iSpc == 0 ? "g" : (iSpc == 1 ? "zee" : (iSpc == 2 ? "zmumu" : "v")));

     for (short iErr = 0; iErr < 3; iErr++) {
      if (iErr != 1 && iData != 0)
       continue; // ignore systematics for MC
      const TString error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

      const TString name = Form ("%sJetPtRatio_2016_%s_%s_%s", spc.Data (), dataType.Data (), error.Data (), period.Data ());
      vJetHists[iPer][iData][iErr][iSpc] = (TH3D*)inFile->Get (name);
     }
 
     const TString name = Form ("%sJetCounts_2016_%s_%s", spc.Data (), dataType.Data (), period.Data ());
     vJetCounts[iPer][iData][iSpc] = (TH2D*)inFile->Get (name);
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
  TCanvas* topCanvas = new TCanvas ("topCanvas", "", 800, 600);
  TCanvas* bottomCanvas = new TCanvas ("bottomCanvas", "", 800, 600);
  const double dPadY = 1.0;
  const double uPadY = 1.0;
  topCanvas->SetBottomMargin (0.14);
  topCanvas->SetTopMargin (0.04);
  topCanvas->SetRightMargin (0.04);
  topCanvas->SetLeftMargin (0.14);
  bottomCanvas->SetBottomMargin (0.14);
  bottomCanvas->SetTopMargin (0.04);
  bottomCanvas->SetRightMargin (0.04);
  bottomCanvas->SetLeftMargin (0.14);
  //const double padRatio = 1.3; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  //const double dPadY = 1.0/ (padRatio+1.0);
  //const double uPadY = 1.0 - dPadY;
  //TPad* topPad = new TPad ("topPad", "", 0, dPadY, 1, 1);
  //TPad* bottomPad = new TPad ("bottomPad", "", 0, 0, 1, dPadY);
  //topPad->SetBottomMargin (0);
  //topPad->SetLeftMargin (-0.20);
  //bottomPad->SetTopMargin (0);
  //bottomPad->SetBottomMargin (0.30);
  //bottomPad->SetLeftMargin (-0.20);
  //topPad->Draw ();
  //bottomPad->Draw ();


  /**** Define local histograms, graphs, etc. ****/
  TH1D *vJetHist = NULL, *vJetHist_mc_overlay = NULL, *vJetHist_mc_signal = NULL, *vJetHist_lo = NULL, *vJetHist_hi = NULL, *vJetHist_rat = NULL, *vJetHist_rat_lo = NULL, *vJetHist_rat_hi = NULL, *vJetRatioCI = NULL;
  TH2D *proj = NULL, *proj_mc_overlay = NULL, *proj_mc_signal = NULL, *proj_lo = NULL, *proj_hi = NULL;
  TGraphAsymmErrors *vJetGraph = NULL, *vJetGraph_mc_overlay = NULL, *vJetGraph_mc_signal = NULL, *vJetGraph_sys = NULL, *vJetGraph_rat = NULL, *vJetGraph_rat_lo = NULL, *vJetGraph_rat_hi = NULL, *vJetGraph_rat_sys = NULL;
  TF1* vJetRatioFit = NULL;//, *polyN = NULL, *logShift = NULL;

  TGraphAsymmErrors** combinedNum = NULL; // 0 = sys_lo, 1 = stat, 2 = sys_hi
  TGraphAsymmErrors** combinedDen = NULL;
  TGraphAsymmErrors** combinedRatio = NULL; // 0 = sys_lo, 1 = stat, 2 = sys_hi

  char* plotName;

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
   const TString period = (iPer == 0 ? "Period A" : (iPer == 1 ? "Period B" : "Period A+B"));
   const TString perType = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));

   for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta

    const int eta_lo = (iEta != numetabins ? iEta+1 : eta_lo_comb);
    const int eta_hi = (iEta != numetabins ? iEta+1 : eta_hi_comb);

    /**** Plot combined VJet info ****/
    combinedNum = Get1DArray <TGraphAsymmErrors*> (3); // 0 = sys_lo, 1 = stat, 2 = sys_hi
    combinedDen = Get1DArray <TGraphAsymmErrors*> (3);
    combinedRatio = Get1DArray <TGraphAsymmErrors*> (2); // 0 = sys_lo, 1 = stat, 2 = sys_hi

    for (short iErr = 0; iErr < 3; iErr++) {
     combinedNum[iErr] = new TGraphAsymmErrors ();
     combinedDen[iErr] = new TGraphAsymmErrors ();
    }

    for (short iSpc = 0; iSpc < 3; iSpc++) {
     const TString species = (iSpc == 0 ? "#gamma + Jets" : (iSpc == 1 ? "Z (ee) + Jets" : (iSpc == 2 ? "Z (#mu#mu) + Jets" : "Combined V + Jets")));

     const Color_t dataColor = (iSpc == 0 ? 9 : (iSpc == 1 ? 8 : kMagenta+1));
     const Style_t dataStyle = 20;
     const Style_t mcOverlayStyle = 33;
     //const Style_t mcSignalStyle = 34;

     topCanvas->cd ();
     topCanvas->SetLogx ();

     proj = Project2D ("", vJetHists[1][iPer][0][1][iSpc], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist = GetProfileX ("vJetHist", proj, numpbins, pbins, true);

     double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;

     vJetGraph = make_graph (vJetHist);
     deltaize (vJetGraph, pow (1.012, iSpc-2), true);
     vJetGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{ref} #left[GeV#right]");
     vJetGraph->GetYaxis ()->SetTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
     vJetGraph->GetYaxis ()->SetRangeUser (middle - 0.35, middle + 0.35);
     vJetGraph->SetMarkerStyle (dataStyle);
     if (iSpc < 3)
      vJetGraph->SetMarkerSize (0.5);
     vJetGraph->SetMarkerColor (dataColor);
     vJetGraph->SetLineColor (dataColor);
     vJetGraph->SetLineWidth (2);
     vJetGraph->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetGraph->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetGraph->GetXaxis ()->SetTitleSize (0.04/uPadY);
     vJetGraph->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetGraph->GetXaxis ()->SetTitleOffset (1.5*uPadY);
     vJetGraph->GetYaxis ()->SetTitleOffset (1.5*uPadY);

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", vJetHists[1][iPer][0][0][iSpc], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numpbins, pbins, true);

     proj_hi = Project2D ("", vJetHists[1][iPer][0][2][iSpc], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numpbins, pbins, true);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     vJetGraph_sys->SetFillColor (dataColor);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", vJetHists[1][iPer][1][1][iSpc], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numpbins, pbins, true);

     vJetGraph_mc_overlay = make_graph (vJetHist_mc_overlay);
     vJetGraph_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetGraph_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetGraph_mc_overlay->SetLineColor (mcOverlayColor);

     //if (iPer != 1 && !skipSignalMC) {
     // proj_mc_signal = Project2D ("", vJetHists[iPer][2][1][iSpc], "x", "z", eta_lo, eta_hi);
     // vJetHist_mc_signal = GetProfileX ("vJetHist_mc_signal", proj_mc_signal, numpbins, pbins, true);

     // vJetGraph_mc_signal = make_graph (vJetHist_mc_signal);
     // vJetGraph_mc_signal->SetMarkerStyle (mcSignalStyle);
     // vJetGraph_mc_signal->SetMarkerColor (mcSignalColor);
     // vJetGraph_mc_signal->SetLineColor (mcSignalColor);
     //}

     if (iSpc == 0)
      ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("ap");
     else
      ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("p");
     if (iSpc == 2)
      ( (TGraphAsymmErrors*)vJetGraph_mc_overlay->Clone ())->Draw ("p");
     
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     int countsData = 0, countsMC = 0;//, countsMC_sig = 0;
     if (exclusive && iEta == numetabins) {
      countsData = vJetCounts[iPer][0][iSpc]->Integral () - vJetCounts[iPer][0][iSpc]->Integral (1, vJetCounts[iPer][0][iSpc]->GetNbinsX(), eta_lo, eta_hi);
      countsMC = vJetCounts[iPer][1][iSpc]->Integral () - vJetCounts[iPer][1][iSpc]->Integral (1, vJetCounts[iPer][1][iSpc]->GetNbinsX(), eta_lo, eta_hi);
      //if (!skipSignalMC)
      // countsMC_sig = vJetCounts[iPer][2][iSpc]->Integral () - vJetCounts[iPer][2][iSpc]->Integral (1, vJetCounts[iPer][2][iSpc]->GetNbinsX(), eta_lo, eta_hi);
     }
     else {
      countsData = vJetCounts[iPer][0][iSpc]->Integral (1, vJetCounts[iPer][0][iSpc]->GetNbinsX(), eta_lo, eta_hi);
      countsMC = vJetCounts[iPer][1][iSpc]->Integral (1, vJetCounts[iPer][1][iSpc]->GetNbinsX(), eta_lo, eta_hi);
      //if (!skipSignalMC)
      // countsMC_sig = vJetCounts[iPer][2][iSpc]->Integral (1, vJetCounts[iPer][2][iSpc]->GetNbinsX(), eta_lo, eta_hi);
     }

     myMarkerText (0.175, 0.12+0.07*(3-iSpc), dataColor, kFullCircle, Form ("%s (%i events)", species.Data (), countsData), 1.25, 0.04/uPadY);
     if (iSpc == 3)
      myMarkerText (0.175, 0.05, mcOverlayColor, kFullDiamond, Form ("Combined Pythia8 with Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);

     myText (0.66, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
     if (eta_lo != 1 || eta_hi != numetabins) {
      if (!exclusive || iEta != numetabins)
       myText (0.66, 0.805, kBlack, Form ("%g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04/uPadY);
      else
       myText (0.66, 0.805, kBlack, Form ("%g < #left|#eta_{det}^{Jet}#right|", etabins[eta_hi]), 0.04/uPadY);
      myText (0.66, 0.74, kBlack, period.Data (), 0.04/uPadY);
     }
     else
      myText (0.66, 0.74, kBlack, period.Data (), 0.04/uPadY);

     bottomCanvas->cd ();
     bottomCanvas->SetLogx ();

     for (int iMC = 1; iMC < 3; iMC++) { // loops over overlay then signal MC
      if (iMC == 2 && (skipSignalMC || iPer == 1)) continue; // no signal only for period B!

      const TString mcType = (iMC == 1 ? "overlay":"signal");
      TH2D* proj_mc = (iMC == 1 ? proj_mc_overlay:proj_mc_signal);

      vJetHist_rat = GetDataOverMC (TString ("vJetPtDataMCRatio"), proj, proj_mc, numpbins, pbins, false, "x");
      vJetHist_rat_lo = GetDataOverMC (TString ("vJetPtDataMCRatio_lo"), proj_lo, proj_mc, numpbins, pbins, false, "x");
      vJetHist_rat_hi = GetDataOverMC (TString ("vJetPtDataMCRatio_hi"), proj_hi, proj_mc, numpbins, pbins, false, "x");

      vJetGraph_rat = make_graph (vJetHist_rat);
      vJetGraph_rat_lo = make_graph (vJetHist_rat_lo);
      vJetGraph_rat_hi = make_graph (vJetHist_rat_hi);

      vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
      CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
      AddWeighted (combinedNum, combinedDen, vJetGraph_rat, vJetGraph_rat_lo, vJetGraph_rat_hi, vJetGraph);

      if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
      if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
      if (vJetGraph_rat_lo) { delete vJetGraph_rat_lo; vJetGraph_rat_lo = NULL; }
      if (vJetGraph_rat_hi) { delete vJetGraph_rat_hi; vJetGraph_rat_hi = NULL; }

      vJetGraph_rat_sys->SetFillColor (dataColor);
      vJetGraph_rat_sys->SetFillStyle (3001);

      deltaize (vJetGraph_rat, pow (1.012, iSpc-2), true);
      vJetGraph_rat->GetXaxis ()->SetTitle ("#it{p}_{T}^{V} #left[GeV#right]");
      vJetGraph_rat->GetYaxis ()->SetTitle ("<x_{J}^{ref}>_{Data} / <x_{J}^{ref}>_{MC}");
      vJetGraph_rat->GetYaxis ()->SetRangeUser (0.88, 1.12);
      //vJetGraph_rat->GetYaxis ()->SetRangeUser (0.75, 1.35);
      vJetGraph_rat->SetMarkerStyle (dataStyle);
      vJetGraph_rat->SetMarkerColor (dataColor);
      vJetGraph_rat->SetLineColor (dataColor);
      vJetGraph_rat->SetLineWidth (2);
      vJetGraph_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
      vJetGraph_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
      vJetGraph_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
      vJetGraph_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
      vJetGraph_rat->GetXaxis ()->SetTitleOffset (1.5*dPadY);
      vJetGraph_rat->GetYaxis ()->SetTitleOffset (1.5*dPadY);
      vJetGraph_rat->GetYaxis ()->CenterTitle (true);
      vJetGraph_rat->GetXaxis ()->SetTickLength (0.08);

      if (iMC == 1 && iSpc == 0)
       ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("ap");
      else
       ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("p");
      ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");

      for (TLine* line : glines) line->Draw ();

      myMarkerText (0.205, 0.12+0.06*(3-iSpc), dataColor, kFullCircle, Form ("%s (%i events)", species.Data (), countsData), 1.25, 0.04/uPadY);
 
      if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
      if (vJetGraph_rat) { delete vJetGraph_rat; vJetGraph_rat = NULL; }
      if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
     }

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_mc_signal) { delete proj_mc_signal; proj_mc_signal = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetHist_mc_signal) { delete vJetHist_mc_signal; vJetHist_mc_signal = NULL; }
     if (vJetGraph) { delete vJetGraph; vJetGraph = NULL; }
     if (vJetGraph_mc_overlay) { delete vJetGraph_mc_overlay; vJetGraph_mc_overlay = NULL; }
     if (vJetGraph_mc_signal) { delete vJetGraph_mc_signal; vJetGraph_mc_signal = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
    } // end loop over species

    combinedRatio[1] = Divide (combinedNum[1], combinedDen[1]);
    combinedRatio[0] = new TGraphAsymmErrors ();
    TGraphAsymmErrors* temp1 = Divide (combinedNum[2], combinedDen[2]);
    TGraphAsymmErrors* temp2 = Divide (combinedNum[0], combinedDen[0]);
    CalcSystematics (combinedRatio[0], combinedRatio[1], temp1, temp2);
    delete temp1;
    delete temp2;

    for (short iErr = 0; iErr < 2; iErr++) { 
     deltaize (combinedRatio[iErr], pow (1.012, 1), true);
     combinedRatio[iErr]->SetMarkerStyle (20);
     combinedRatio[iErr]->SetLineColor (kBlack);
     combinedRatio[iErr]->SetMarkerColor (kBlack);
     combinedRatio[iErr]->SetLineWidth (2);
    }
    combinedRatio[0]->SetFillColor (kBlack);
    combinedRatio[0]->SetFillStyle (3001);
    ( (TGraphAsymmErrors*)(combinedRatio[0])->Clone ())->Draw ("2");
    combinedRatio[0]->Draw ("2");
    ( (TGraphAsymmErrors*)(combinedRatio[1])->Clone ())->Draw ("p");
    combinedRatio[1]->Draw ("p");

    myMarkerText (0.205, 0.36, kBlack, kFullCircle, "Combined Z/#gamma + Jet", 1.25, 0.04/uPadY);
    myText (0.66, 0.33, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
    if (eta_lo != 1 || eta_hi != numetabins) {
     if (!exclusive || iEta != numetabins)
      myText (0.66, 0.255, kBlack, Form ("%g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04/uPadY);
     else
      myText (0.66, 0.255, kBlack, Form ("%g < #left|#eta_{det}^{Jet}#right|", etabins[eta_hi]), 0.04/uPadY);
     myText (0.66, 0.19, kBlack, period.Data (), 0.04/uPadY);
    }
    else
     myText (0.66, 0.19, kBlack, period.Data (), 0.04/uPadY);

    TH1D* sysCombinedRatioHist = new TH1D ("sysCombinedRatioHist", "", numpbins, pbins);
    for (int ix = 0; ix < combinedRatio[0]->GetN (); ix++) {
     double x, y;
     combinedRatio[0]->GetPoint (ix, x, y);
     sysCombinedRatioHist->SetBinContent (ix+1, y);
     sysCombinedRatioHist->SetBinError (ix+1, sqrt (pow (combinedRatio[0]->GetErrorY (ix), 2) + pow (combinedRatio[1]->GetErrorY (ix), 2))); 
    }

    //polyN = new TF1 ("polyN", "[0]+[1]*x+[2]*(2*x^2-1)+[3]*(4*x^3-3*x)+[4]*(8*x^4-8*x^2+1)", pbins[0], pbins[numpbins]);
    //polyN->SetParameter (0, 1);
    //polyN->SetParameter (1, 0);
    //polyN->SetParameter (2, 0);
    //polyN->SetParameter (3, 0);
    //polyN->SetParameter (4, 0);
    //logShift = new TF1 ("logShift", "(log(x)-[5])/[6]", pbins[0], pbins[numpbins]);
    //logShift->FixParameter (5, 0.5 * (log (pbins[numpbins]) + log (pbins[0])));
    //logShift->FixParameter (6, 0.5 * (log (pbins[numpbins]) - log (pbins[0])));

    const TString func = "[0] + [1]*((log(x)-[5])/[6]) + [2]*(2*((log(x)-[5])/[6])^2-1) + [3]*(4*((log(x)-[5])/[6])^3-3*(log(x)-[5])/[6]) + [4]*(8*((log(x)-[5])/[6])^4-8*((log(x)-[5])/[6])^2+1)";
    vJetRatioFit = new TF1 ("vJetPtDataMCRatioFit", func, pbins[0], pbins[numpbins]);
    vJetRatioFit->SetParameter (0, 1);
    vJetRatioFit->SetParameter (1, 0);
    vJetRatioFit->SetParameter (2, 0);
    vJetRatioFit->SetParameter (3, 0);
    vJetRatioFit->SetParameter (4, 0);
    //vJetRatioFit->FixParameter (5, 0.5 * (log (pbins[numpbins]) + log (pbins[0])));
    //vJetRatioFit->FixParameter (6, 0.5 * (log (pbins[numpbins]) - log (pbins[0])));
    sysCombinedRatioHist->Fit (vJetRatioFit, "RN0Q");

    vJetRatioCI = new TH1D ("vJetPtDataMCRatioCI", "", numpbins, pbins);
    (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (vJetRatioCI, 0.68);

    vJetRatioFit->SetLineColor (kGray+2);
    ( (TF1*)vJetRatioFit->Clone ())->Draw ("same");
    vJetRatioCI->SetMarkerStyle (kDot);
    vJetRatioCI->SetFillColorAlpha (kGray+2, 0.3);
    vJetRatioCI->DrawCopy ("e3 same");

    //if (polyN) { delete polyN; polyN = NULL; }
    //if (logShift) { delete logShift; logShift = NULL; }
    if (vJetRatioFit) { delete vJetRatioFit; vJetRatioFit = NULL; }
    if (vJetRatioCI) { delete vJetRatioCI; vJetRatioCI = NULL; }
    if (sysCombinedRatioHist) { delete sysCombinedRatioHist; sysCombinedRatioHist = NULL; }

    Delete1DArray (combinedNum, 3);
    Delete1DArray (combinedDen, 3);
    Delete1DArray (combinedRatio, 2);

    //} // end loop over insitu configurations
     
    if (iEta < numetabins)
     plotName = Form ("v_jet_xjref_iEta%i.pdf", iEta);
    else
     plotName = Form ("v_jet_xjref_iEta_combined.pdf");

    topCanvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));

    if (iEta < numetabins)
     plotName = Form ("v_jet_dataOverMC_iEta%i.pdf", iEta);
    else
     plotName = Form ("v_jet_dataOverMC_iEta_combined.pdf");

    bottomCanvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));

   } // end loop over etabins


   /**** Now loop over pT bins and plot response as function of eta^jet ****/
   for (short iP = 0; iP <= numpbins; iP++) {

    int p_lo = iP+1;
    int p_hi = iP+1;
    if (iP == numpbins) {
     p_lo = 5;//1;
     p_hi = numpbins;
    }

    /**** Plot combined VJet info ****/
    for (short iSpc = 0; iSpc < 3; iSpc++) {
     const TString species = (iSpc == 0 ? "#gamma + Jets" : (iSpc == 1 ? "Z (ee) + Jets" : (iSpc == 2 ? "Z (#mu#mu) + Jets" : "Combined V + Jets")));

     const Color_t dataColor = (iSpc == 0 ? 9 : (iSpc == 1 ? 8 : (iSpc == 2 ? kMagenta+1 : kBlack)));
     const Style_t dataStyle = 20;
     const Style_t mcOverlayStyle = 33;
     //const Style_t mcSignalStyle = 34;

     topCanvas->cd ();
     topCanvas->SetLogx (false);

     proj = Project2D ("", vJetHists[1][iPer][0][1][iSpc], "y", "z", p_lo, p_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numetabins, etabins, false);

     double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;

     vJetGraph = make_graph (vJetHist);
     deltaize (vJetGraph, 0.02*(iSpc-2), false);
     vJetGraph->GetXaxis ()->SetTitle ("Jet #eta");
     vJetGraph->GetYaxis ()->SetTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
     vJetGraph->GetYaxis ()->SetRangeUser (middle - 0.35, middle + 0.35);
     vJetGraph->SetMarkerStyle (dataStyle);
     vJetGraph->SetMarkerColor (dataColor);
     vJetGraph->SetLineColor (dataColor);
     vJetGraph->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetGraph->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetGraph->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetGraph->GetXaxis ()->SetTitleOffset (1.5*uPadY);
     vJetGraph->GetYaxis ()->SetTitleOffset (1.5*uPadY);

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", vJetHists[1][iPer][0][0][iSpc], "y", "z", p_lo, p_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj, numetabins, etabins, false);

     proj_hi = Project2D ("", vJetHists[1][iPer][0][2][iSpc], "y", "z", p_lo, p_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, false);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }
     vJetGraph_sys->SetFillColor (dataColor);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", vJetHists[1][iPer][1][1][iSpc], "y", "z", p_lo, p_hi);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numetabins, etabins, false);

     vJetGraph_mc_overlay = make_graph (vJetHist_mc_overlay);
     vJetGraph_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetGraph_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetGraph_mc_overlay->SetLineColor (mcOverlayColor);

     //if (iPer != 1 && !skipSignalMC) {
     // proj_mc_signal = Project2D ("", vJetHists[iPer][2][1][iSpc], "y", "z", p_lo, p_hi);
     // vJetHist_mc_signal = GetProfileX ("vJetHist_mc_signal", proj_mc_signal, numpbins, pbins, true);

     // vJetGraph_mc_signal = make_graph (vJetHist_mc_signal);
     // vJetGraph_mc_signal->SetMarkerStyle (mcSignalStyle);
     // vJetGraph_mc_signal->SetMarkerColor (mcSignalColor);
     // vJetGraph_mc_signal->SetLineColor (mcSignalColor);
     //}

     if (iSpc == 0)
      ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("ap");
     else
      ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("p");
     if (iSpc == 3)
      ( (TGraphAsymmErrors*)vJetGraph_mc_overlay->Clone ())->Draw ("p");
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     int countsData = 0, countsMC = 0;//, countsMC_sig = 0;
     countsData = vJetCounts[iPer][0][iSpc]->Integral (p_lo, p_hi, 1, vJetCounts[iPer][0][iSpc]->GetNbinsY());
     countsMC = vJetCounts[iPer][1][iSpc]->Integral (p_lo, p_hi, 1, vJetCounts[iPer][1][iSpc]->GetNbinsY());
     //if (!skipSignalMC)
     // countsMC_sig = vJetCounts[iPer][2][iSpc]->Integral (p_lo, p_hi, 1, vJetCounts[iPer][2][iSpc]->GetNbinsY());

     myMarkerText (0.175, 0.89-0.07*iSpc, dataColor, kFullCircle, Form ("%s (%i events)", species.Data (), countsData), 1.25, 0.04/uPadY);
     if (iSpc == 3)
      myMarkerText (0.175, 0.61, mcOverlayColor, kFullDiamond, Form ("Combined Pythia8 with Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);

     //myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countsData), 1.25, 0.04/uPadY);
     //myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
     //if (iPer != 1 && !skipSignalMC)
     // myMarkerText (0.175, 0.74, mcSignalColor, kFullCross, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", countsMC_sig), 1.25, 0.04/uPadY);
     myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
     if (p_lo != 1 || p_hi != numpbins)
      myText (0.155, 0.15, dataColor, Form ("%g < #it{p}_{T}^{V} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
     //else
     // myText (0.155, 0.15, dataColor, "V + Jet", 0.04/uPadY);
     myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);
     

     bottomCanvas->cd ();
     bottomCanvas->SetLogx (false);

     for (int iMC = 1; iMC < 3; iMC++) { // loops over overlay then signal MC
      if (iMC == 2 && (skipSignalMC || iPer == 1)) continue; // no signal only for period B!

      const TString mcType = (iMC == 1 ? "overlay":"signal");
      TH2D* proj_mc = (iMC == 1 ? proj_mc_overlay:proj_mc_signal);

      vJetHist_rat = GetDataOverMC (TString (Form ("vJetPtDataMCRatio_iP%i", iP)), proj, proj_mc, numetabins, etabins, false, "x");
      vJetHist_rat_lo = GetDataOverMC (TString (Form ("vJetPtDataMCRatio_lo_iP%i", iP)), proj_lo, proj_mc, numetabins, etabins, false, "x");
      vJetHist_rat_hi = GetDataOverMC (TString (Form ("vJetPtDataMCRatio_hi_iP%i", iP)), proj_hi, proj_mc, numetabins, etabins, false, "x");

      vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
      CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
      if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
      if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
      vJetGraph_rat_sys->SetFillColor (dataColor);
      vJetGraph_rat_sys->SetFillStyle (3001);

      vJetGraph_rat = make_graph (vJetHist_rat);
      deltaize (vJetGraph_rat, 0.02*(iSpc-2), false);
      vJetGraph_rat->GetXaxis ()->SetTitle ("Jet #eta");
      vJetGraph_rat->GetYaxis ()->SetTitle ("<x_{J}^{ref}>_{Data} / <x_{J}^{ref}>_{MC}");
      vJetGraph_rat->GetYaxis ()->SetRangeUser (0.91, 1.09);
      //vJetGraph_rat->GetYaxis ()->SetRangeUser (0.75, 1.35);
      vJetGraph_rat->SetMarkerStyle (dataStyle);
      if (iSpc < 3)
       vJetGraph_rat->SetMarkerSize (0.5);
      vJetGraph_rat->SetMarkerColor (dataColor);
      vJetGraph_rat->SetLineColor (dataColor);
      vJetGraph_rat->GetYaxis ()->SetNdivisions (405);
      vJetGraph_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
      vJetGraph_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
      vJetGraph_rat->GetXaxis ()->SetTitleOffset (1.5);
      vJetGraph_rat->GetYaxis ()->SetTitleOffset (1.5*dPadY);
      vJetGraph_rat->GetYaxis ()->CenterTitle (true);
      vJetGraph_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
      vJetGraph_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
      vJetGraph_rat->GetXaxis ()->SetTickLength (0.08);

      //TF1* midConst = new TF1 ("midConst", "[0]", -4.9, 4.9);
      //vJetHist_rat->Fit (midConst, "NQ0R");
      //constFit->SetPoint (0, 0, midConst->GetParameter (0));
      //constFit->SetPointError (0, 4.9, midConst->GetParError (0));

      if (iMC == 1 && iSpc == 0)
       ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("ap");
      else
       ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("p");
      ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");

      for (TLine* line : glines) line->Draw ();
      //for (TLine* line : dplines_bottom) line->Draw ();
      //constFit->Draw ("2");

      if (vJetHist_rat->Integral () != 0. && iSpc == 3) {
       // add systematic errors
       for (int ix = 1; ix < vJetHist_rat->GetNbinsX (); ix++) {
        const double sys_err = std::max (vJetGraph_rat_sys->GetErrorYhigh (ix-1), vJetGraph_rat_sys->GetErrorYlow (ix-1));
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

       vJetRatioFit->SetLineColor (kGray+2);
       ( (TF1*)vJetRatioFit->Clone ())->Draw ("same");
       vJetRatioCI->SetMarkerStyle (kDot);
       vJetRatioCI->SetFillColorAlpha (kGray+2, 0.3);
       vJetRatioCI->DrawCopy ("e3 same");

       if (vJetRatioFit) { delete vJetRatioFit; vJetRatioFit = NULL; }
       if (vJetRatioCI) { delete vJetRatioCI; vJetRatioCI = NULL; }
      }

      if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
      if (vJetGraph_rat) { delete vJetGraph_rat; vJetGraph_rat = NULL; }
      if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
     }

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_mc_signal) { delete proj_mc_signal; proj_mc_signal = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetHist_mc_signal) { delete vJetHist_mc_signal; vJetHist_mc_signal = NULL; }
     if (vJetGraph) { delete vJetGraph; vJetGraph = NULL; }
     if (vJetGraph_mc_overlay) { delete vJetGraph_mc_overlay; vJetGraph_mc_overlay = NULL; }
     if (vJetGraph_mc_signal) { delete vJetGraph_mc_signal; vJetGraph_mc_signal = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
    } // end loop over species

    if (iP < numpbins)
     plotName = Form ("v_jet_xjref_iP%i.pdf", iP);
    else
     plotName = Form ("v_jet_xjref_iP_combined.pdf");

    topCanvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));

    if (iP < numpbins)
     plotName = Form ("v_jet_dataOverMC_iP%i.pdf", iP);
    else
     plotName = Form ("v_jet_dataOverMC_iP_combined.pdf");

    bottomCanvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));

   } // end loop over pT bins
  } // end loop over beam periods


  return;
}

} // end namespace
