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
    continue;
    //ry = 0;
    //re = 0;
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

  //////////////////////////////////////////////////////////////////////////////
  // Create desired histograms
  //////////////////////////////////////////////////////////////////////////////
  TH3D***** vJetHists = Get4DArray <TH3D*> (3, 2, 3, 3); // iPer, iData, iErr, iSpc
  TH3D**** vJetCounts = Get3DArray <TH3D*> (3, 2, 3); //  iPer, iData, iSpc

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  for (short iPer = 0; iPer < 3; iPer++) {
    const char* per = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));

    for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
      const char* data = (iData == 0 ? "data" : "mc");

      for (short iSpc = 0; iSpc < 2; iSpc++) {
        const char* spc = (iSpc == 0 ? "photon" : "z");

        for (short iErr = 0; iErr < 3; iErr++) {
          if (iErr != 1 && iData != 0)
            continue; // ignore systematics for MC
          const char* error = (iErr == 0 ? "syslo" : (iErr == 1 ? "stat" : "syshi"));

          const char* name = Form ("%sJetPtRatio_%s_%s_%s%s", spc, per, data, error, (iSpc==0?"_signal":""));
          vJetHists[iPer][iData][iSpc][iErr] = (TH3D*)inFile->Get (name);
        }
        const char* name = Form ("%sJetCounts_%s_%s%s", spc, per, data, (iSpc==0?"_tight_unweighted":""));
        vJetCounts[iPer][iData][iSpc] = (TH3D*)inFile->Get (name);
      }

      for (short iErr = 0; iErr < 3; iErr++) {
        if (iErr != 1 && iData != 0)
          continue; // ignore systematics for MC
        const char* error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

        const char* name = Form ("vJetPtRatio_%s_%s_%s", data, error, per);
        vJetHists[iPer][iData][2][iErr] = new TH3D (name, "", numpbins, pbins, numetabins, etabins, 400, linspace (0, 2.0, 400));
        vJetHists[iPer][iData][2][iErr]->Sumw2 ();
      }
      const char* name = Form ("vJetCounts_%s_%s", data, per);
      vJetCounts[iPer][iData][2] = new TH3D (name, "", numpbins, pbins, numetabins, etabins, numphibins, phibins);
      vJetCounts[iPer][iData][2]->Sumw2 ();
      
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // Plotting elements 
  //////////////////////////////////////////////////////////////////////////////
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


  //////////////////////////////////////////////////////////////////////////////
  // Canvas definitions
  //////////////////////////////////////////////////////////////////////////////
  TCanvas* rawDistCanvas = new TCanvas ("rawDistCanvas", "", 1600, 1200);
  rawDistCanvas->SetLeftMargin (-0.12);
  rawDistCanvas->SetBottomMargin (-0.12);
  //rawDistCanvas->Divide (8, 6);

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
  TH1D *vJetHist = NULL, *vJetHist_mc = NULL, *vJetHist_lo = NULL, *vJetHist_hi = NULL, *vJetHist_rat = NULL, *vJetHist_rat_lo = NULL, *vJetHist_rat_hi = NULL, *vJetRatioCI = NULL;
  TH2D *proj = NULL, *proj_mc = NULL, *proj_lo = NULL, *proj_hi = NULL;
  TGraphAsymmErrors *vJetGraph = NULL, *vJetGraph_mc = NULL, *vJetGraph_sys = NULL, *vJetGraph_rat = NULL, *vJetGraph_rat_lo = NULL, *vJetGraph_rat_hi = NULL, *vJetGraph_rat_sys = NULL;
  TF1* vJetRatioFit = NULL;//, *polyN = NULL, *logShift = NULL;

  TGraphAsymmErrors** combinedNum = NULL; // 0 = sys_lo, 1 = stat, 2 = sys_hi
  TGraphAsymmErrors** combinedDen = NULL;
  TGraphAsymmErrors** combinedRatio = NULL; // 0 = sys, 1 = stat

  char* plotName;

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
    const char* period = (iPer == 0 ? "Period A" : (iPer == 1 ? "Period B" : "Period A+B"));

    for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta

      const int eta_lo = (iEta != numetabins ? iEta+1 : eta_lo_comb);
      const int eta_hi = (iEta != numetabins ? iEta+1 : eta_hi_comb);

      /**** Plot combined VJet info ****/
      combinedNum = Get1DArray <TGraphAsymmErrors*> (3); // 0 = sys_lo, 1 = stat, 2 = sys_hi
      combinedDen = Get1DArray <TGraphAsymmErrors*> (3);
      combinedRatio = Get1DArray <TGraphAsymmErrors*> (2); // 0 = sys, 1 = stat
      for (short iErr = 0; iErr < 3; iErr++) {
        combinedNum[iErr] = new TGraphAsymmErrors ();
        combinedDen[iErr] = new TGraphAsymmErrors ();
      }

      for (short iSpc = 0; iSpc < 2; iSpc++) {
        const TString species = (iSpc == 0 ? "#gamma + Jets" : (iSpc == 1 ? "Z + Jets" : "Combined V + Jets"));

        const Color_t dataColor = (iSpc == 0 ? 9 : (iSpc == 1 ? kMagenta+1 : 8));
        const Color_t mcColor = dataColor;
        const Style_t dataStyle = 20;
        const Style_t mcStyle = 24;

        topCanvas->cd ();
        topCanvas->SetLogx ();

        proj = Project2D ("", vJetHists[iPer][0][iSpc][1], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
        proj->RebinY (rebinFactor);
        vJetHist = GetProfileX ("vJetHist", proj, numpbins, pbins, true);

        double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
        if (10 * middle != floor (10*middle)) middle += 0.05;
        vJetGraph = make_graph (vJetHist);
        deltaize (vJetGraph, pow (1.05, (iSpc-0.5)), true);
        vJetGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{ref} #left[GeV#right]");
        vJetGraph->GetYaxis ()->SetTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
        vJetGraph->GetYaxis ()->SetRangeUser (middle - 0.35, middle + 0.35);
        vJetGraph->SetMarkerStyle (dataStyle);
        //if (iSpc < 3)
        //  vJetGraph->SetMarkerSize (0.5);
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
        proj_lo = Project2D ("", vJetHists[iPer][0][iSpc][0], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
        proj_lo->RebinY (rebinFactor);
        vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numpbins, pbins, true);

        proj_hi = Project2D ("", vJetHists[iPer][0][iSpc][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
        proj_hi->RebinY (rebinFactor);
        vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numpbins, pbins, true);

        vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
        CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
        if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
        if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

        vJetGraph_sys->SetFillColor (dataColor);
        vJetGraph_sys->SetFillStyle (3001);

        proj_mc = Project2D ("", vJetHists[iPer][1][iSpc][1], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
        proj_mc->RebinY (rebinFactor);
        vJetHist_mc = GetProfileX ("vJetHist_mc", proj_mc, numpbins, pbins, true);

        vJetGraph_mc = make_graph (vJetHist_mc);
        deltaize (vJetGraph_mc, pow (1.05, 0.5*iSpc-0.25), true);
        vJetGraph_mc->SetMarkerStyle (mcStyle);
        vJetGraph_mc->SetMarkerColor (dataColor);
        vJetGraph_mc->SetLineColor (dataColor);
        vJetGraph_mc->SetLineWidth (2);

        if (iSpc == 0)
          ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("ap");
        else
          ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("p");
        ( (TGraphAsymmErrors*)vJetGraph_mc->Clone ())->Draw ("p");
        
        ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

        int countsData = 0, countsMC = 0;
        if (exclusive && iEta == numetabins) {
          countsData = vJetCounts[iPer][0][iSpc]->Integral () - vJetCounts[iPer][0][iSpc]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
          countsMC = vJetCounts[iPer][1][iSpc]->Integral () - vJetCounts[iPer][1][iSpc]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
        }
        else {
          countsData = vJetCounts[iPer][0][iSpc]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
          countsMC = vJetCounts[iPer][1][iSpc]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
        }

        myMarkerText (0.205, 0.36-0.12*(1-iSpc), dataColor, dataStyle, Form ("%s, 2016 Data (%i events)", species.Data (), countsData), 1.25, 0.04/uPadY);
        myMarkerText (0.205, 0.30-0.12*(1-iSpc), mcColor, mcStyle, Form ("%s, Pythia8 MC (%i events)", species.Data (), countsMC), 1.25, 0.04/uPadY);
        //if (iSpc == 2)
        //  myMarkerText (0.205, 0.22, dataColor, mcStyle, Form ("Combined Pythia8 MC (%i events)", countsMC), 1.25, 0.04/uPadY);

        myText (0.66, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
        if (eta_lo != 1 || eta_hi != numetabins) {
          if (!exclusive || iEta != numetabins)
            myText (0.66, 0.815, kBlack, Form ("%g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04/uPadY);
          else
            myText (0.66, 0.815, kBlack, Form ("%g < #left|#eta_{det}^{Jet}#right|", etabins[eta_hi]), 0.04/uPadY);
          myText (0.66, 0.76, kBlack, period, 0.04/uPadY);
        }
        else
          myText (0.66, 0.76, kBlack, period, 0.04/uPadY);

        bottomCanvas->cd ();
        bottomCanvas->SetLogx ();

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

        deltaize (vJetGraph_rat, pow (1.05, iSpc-1), true);
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

        if (iSpc == 0)
          ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("ap");
        else
          ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("p");
        ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");

        for (TLine* line : glines) line->Draw ();

        myMarkerText (0.205, 0.22+0.06*(1-iSpc), dataColor, kFullCircle, Form ("%s (%i events)", species.Data (), countsData), 1.25, 0.04/uPadY);
 
        if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
        if (vJetGraph_rat) { delete vJetGraph_rat; vJetGraph_rat = NULL; }
        if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }

        if (proj) { delete proj; proj = NULL; }
        if (proj_mc) { delete proj_mc; proj_mc = NULL; }
        if (proj_lo) { delete proj_lo; proj_lo = NULL; }
        if (proj_hi) { delete proj_hi; proj_hi = NULL; }

        if (vJetHist) { delete vJetHist; vJetHist = NULL; }
        if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
        if (vJetGraph) { delete vJetGraph; vJetGraph = NULL; }
        if (vJetGraph_mc) { delete vJetGraph_mc; vJetGraph_mc = NULL; }
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
        deltaize (combinedRatio[iErr], pow (1.05, 1), true);
        combinedRatio[iErr]->SetMarkerStyle (20);
        combinedRatio[iErr]->SetLineColor (kBlack);
        combinedRatio[iErr]->SetMarkerColor (kBlack);
        combinedRatio[iErr]->SetLineWidth (2);
      }
      combinedRatio[0]->SetFillColor (kBlack);
      combinedRatio[0]->SetFillStyle (3001);
      combinedRatio[0]->Draw ("2");
      combinedRatio[1]->Draw ("p");

      myMarkerText (0.205, 0.34, kBlack, kFullCircle, "Combined Z/#gamma + Jet", 1.25, 0.04/uPadY);
      myText (0.66, 0.34, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      if (eta_lo != 1 || eta_hi != numetabins) {
        if (!exclusive || iEta != numetabins)
          myText (0.66, 0.275, kBlack, Form ("%g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04/uPadY);
        else
          myText (0.66, 0.275, kBlack, Form ("%g < #left|#eta_{det}^{Jet}#right|", etabins[eta_hi]), 0.04/uPadY);
        myText (0.66, 0.22, kBlack, period, 0.04/uPadY);
      }
      else
        myText (0.66, 0.22, kBlack, period, 0.04/uPadY);

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
      vJetRatioFit->FixParameter (5, 0.5 * (log (pbins[numpbins]) + log (pbins[0])));
      vJetRatioFit->FixParameter (6, 0.5 * (log (pbins[numpbins]) - log (pbins[0])));
      //sysCombinedRatioHist->Fit (vJetRatioFit, "RN0Q");
      combinedRatio[1]->Fit (vJetRatioFit, "RN0Q");

      const float chisq = vJetRatioFit->GetChisquare ();
      myText (0.175, 0.4, kBlack, Form ("#chi^{2} = %g", chisq), 0.04/dPadY);

      vJetRatioCI = new TH1D ("vJetPtDataMCRatioCI", "", numpbins, pbins);
      (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (vJetRatioCI, 0.68);

      vJetRatioFit->SetLineColor (kGray+2);
      ( (TF1*)vJetRatioFit->Clone ())->Draw ("same");
      vJetRatioCI->SetMarkerStyle (kDot);
      vJetRatioCI->SetFillColorAlpha (kGray+2, 0.3);
      vJetRatioCI->DrawCopy ("e3 same");

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

      //if (polyN) { delete polyN; polyN = NULL; }
      //if (logShift) { delete logShift; logShift = NULL; }
      if (vJetRatioFit) { delete vJetRatioFit; vJetRatioFit = NULL; }
      if (vJetRatioCI) { delete vJetRatioCI; vJetRatioCI = NULL; }
      if (sysCombinedRatioHist) { delete sysCombinedRatioHist; sysCombinedRatioHist = NULL; }

      Delete1DArray (combinedNum, 3);
      Delete1DArray (combinedDen, 3);
      Delete1DArray (combinedRatio, 2);

    } // end loop over etabins


    /**** Now loop over pT bins and plot response as function of eta^jet ****/
    for (short iP = 0; iP <= numpbins; iP++) {

      const int p_lo = (iP != numpbins ? iP+1 : p_lo_comb);
      const int p_hi = (iP != numpbins ? iP+1 : p_hi_comb);

      /**** Plot combined VJet info ****/
      combinedNum = Get1DArray <TGraphAsymmErrors*> (3); // 0 = sys_lo, 1 = stat, 2 = sys_hi
      combinedDen = Get1DArray <TGraphAsymmErrors*> (3);
      combinedRatio = Get1DArray <TGraphAsymmErrors*> (2); // 0 = sys_lo, 1 = stat, 2 = sys_hi
      for (short iErr = 0; iErr < 3; iErr++) {
        combinedNum[iErr] = new TGraphAsymmErrors ();
        combinedDen[iErr] = new TGraphAsymmErrors ();
      }

      for (short iSpc = 0; iSpc < 2; iSpc++) {
        const TString species = (iSpc == 0 ? "#gamma + Jets" : (iSpc == 1 ? "Z + Jets" : "Combined V + Jets"));

        const Color_t dataColor = (iSpc == 0 ? 9 : (iSpc == 1 ? kMagenta+1 : 8));
        const Color_t mcColor = dataColor;
        const Style_t dataStyle = 20;
        const Style_t mcStyle = 24;

        topCanvas->cd ();
        topCanvas->SetLogx (false);

        proj = Project2D ("", vJetHists[iPer][0][iSpc][1], "y", "z", p_lo, p_hi);
        proj->RebinY (rebinFactor);
        vJetHist = GetProfileX ("vJetHist", proj, numetabins, etabins, false);

        double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
        if (10 * middle != floor (10*middle)) middle += 0.05;

        vJetGraph = make_graph (vJetHist);
        deltaize (vJetGraph, 0.03*(iSpc-0.5), false);
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
        proj_lo = Project2D ("", vJetHists[iPer][0][iSpc][0], "y", "z", p_lo, p_hi);
        proj_lo->RebinY (rebinFactor);
        vJetHist_lo = GetProfileX ("vJetHist_lo", proj, numetabins, etabins, false);

        proj_hi = Project2D ("", vJetHists[iPer][0][iSpc][2], "y", "z", p_lo, p_hi);
        proj_hi->RebinY (rebinFactor);
        vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, false);

        vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
        CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
        if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
        if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }
        vJetGraph_sys->SetFillColor (dataColor);
        vJetGraph_sys->SetFillStyle (3001);

        proj_mc = Project2D ("", vJetHists[iPer][1][iSpc][1], "y", "z", p_lo, p_hi);
        proj_mc->RebinY (rebinFactor);
        vJetHist_mc = GetProfileX ("vJetHist_mc", proj_mc, numetabins, etabins, false);

        vJetGraph_mc = make_graph (vJetHist_mc);
        deltaize (vJetGraph, 0.03*(0.5*iSpc-0.25), false);
        vJetGraph_mc->SetMarkerStyle (mcStyle);
        vJetGraph_mc->SetMarkerColor (mcColor);
        vJetGraph_mc->SetLineColor (mcColor);
        vJetGraph_mc->SetLineWidth (2);

        if (iSpc == 0)
          ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("ap");
        else
          ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("p");
        ( (TGraphAsymmErrors*)vJetGraph_mc->Clone ())->Draw ("p");
        ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

        int countsData = 0, countsMC = 0;
        countsData = vJetCounts[iPer][0][iSpc]->Integral (p_lo, p_hi, 1, numetabins, 1, numphibins);
        countsMC = vJetCounts[iPer][1][iSpc]->Integral (p_lo, p_hi, 1, numetabins, 1, numphibins);

        myMarkerText (0.205, 0.36-0.12*(1-iSpc), dataColor, dataStyle, Form ("%s, 2016 Data (%i events)", species.Data (), countsData), 1.25, 0.04/uPadY);
        myMarkerText (0.205, 0.30-0.12*(1-iSpc), mcColor, mcStyle, Form ("%s, Pythia8 MC (%i events)", species.Data (), countsMC), 1.25, 0.04/uPadY);

        myText (0.66, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
        if (p_lo != 1 || p_hi != numpbins)
          myText (0.66, 0.815, kBlack, Form ("%g < #it{p}_{T}^{V} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
        myText (0.66, 0.76, kBlack, period, 0.04/uPadY);
        

        bottomCanvas->cd ();
        bottomCanvas->SetLogx (false);

        vJetHist_rat = GetDataOverMC (TString (Form ("vJetPtDataMCRatio_iP%i", iP)), proj, proj_mc, numetabins, etabins, false, "x");
        vJetHist_rat_lo = GetDataOverMC (TString (Form ("vJetPtDataMCRatio_lo_iP%i", iP)), proj_lo, proj_mc, numetabins, etabins, false, "x");
        vJetHist_rat_hi = GetDataOverMC (TString (Form ("vJetPtDataMCRatio_hi_iP%i", iP)), proj_hi, proj_mc, numetabins, etabins, false, "x");

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

        deltaize (vJetGraph_rat, 0.02*(iSpc-2), false);
        vJetGraph_rat->GetXaxis ()->SetTitle ("Jet #eta");
        vJetGraph_rat->GetYaxis ()->SetTitle ("<x_{J}^{ref}>_{Data} / <x_{J}^{ref}>_{MC}");
        vJetGraph_rat->GetYaxis ()->SetRangeUser (0.91, 1.09);
        //vJetGraph_rat->GetYaxis ()->SetRangeUser (0.75, 1.35);
        vJetGraph_rat->SetMarkerStyle (dataStyle);
        //if (iSpc < 3)
        //  vJetGraph_rat->SetMarkerSize (0.5);
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

        if (iSpc == 0)
          ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("ap");
        else
          ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("p");
        ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");

        for (TLine* line : glines) line->Draw ();

        myMarkerText (0.205, 0.22+0.06*(1-iSpc), dataColor, kFullCircle, Form ("%s (%i events)", species.Data (), countsData), 1.25, 0.04/uPadY);

        //if (vJetHist_rat->Integral () != 0. && iSpc == 2) {
        //  // add systematic errors
        //  for (int ix = 1; ix < vJetHist_rat->GetNbinsX (); ix++) {
        //    const double sys_err = std::max (vJetGraph_rat_sys->GetErrorYhigh (ix-1), vJetGraph_rat_sys->GetErrorYlow (ix-1));
        //    vJetHist_rat->SetBinError (ix, sqrt (pow (vJetHist_rat->GetBinError (ix), 2) + pow (sys_err, 2)));
        //  }

        //  const TString func = "[0] + [1]*((x-[3])/[4]) + [2]*(2*((x-[3])/[4])^2-1)";
        //  TF1* vJetRatioFit = new TF1 ("vJetPtDataMCRatioFit", func, etabins[0], etabins[numetabins]);
        //  vJetRatioFit->SetParameter (0, 1);
        //  vJetRatioFit->SetParameter (1, 0);
        //  vJetRatioFit->SetParameter (2, 0);
        //  vJetRatioFit->FixParameter (3, 0.5 * (etabins[numetabins] + etabins[0]));
        //  vJetRatioFit->FixParameter (4, 0.5 * (etabins[numetabins] - etabins[0]));
        //  vJetHist_rat->Fit (vJetRatioFit, "RN0Q");

        //  const float chisq = vJetRatioFit->GetChisquare ();
        //  myText (0.155, 0.4, kBlack, Form ("#chi^{2} = %g", chisq), 0.04/dPadY);

        //  vJetRatioCI = new TH1D ("vJetPtDataMCRatioCI", "", numetabins, etabins);
        //  (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (vJetRatioCI);

        //  vJetRatioFit->SetLineColor (kGray+2);
        //  ( (TF1*)vJetRatioFit->Clone ())->Draw ("same");
        //  vJetRatioCI->SetMarkerStyle (kDot);
        //  vJetRatioCI->SetFillColorAlpha (kGray+2, 0.3);
        //  vJetRatioCI->DrawCopy ("e3 same");

        //  if (vJetRatioFit) { delete vJetRatioFit; vJetRatioFit = NULL; }
        //  if (vJetRatioCI) { delete vJetRatioCI; vJetRatioCI = NULL; }
        //}

        if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
        if (vJetGraph_rat) { delete vJetGraph_rat; vJetGraph_rat = NULL; }
        if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }

        if (proj) { delete proj; proj = NULL; }
        if (proj_mc) { delete proj_mc; proj_mc = NULL; }
        if (proj_lo) { delete proj_lo; proj_lo = NULL; }
        if (proj_hi) { delete proj_hi; proj_hi = NULL; }

        if (vJetHist) { delete vJetHist; vJetHist = NULL; }
        if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
        if (vJetGraph) { delete vJetGraph; vJetGraph = NULL; }
        if (vJetGraph_mc) { delete vJetGraph_mc; vJetGraph_mc = NULL; }
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
        deltaize (combinedRatio[iErr], pow (1.05, 1), true);
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

      myMarkerText (0.205, 0.34, kBlack, kFullCircle, "Combined Z/#gamma + Jet", 1.25, 0.04/uPadY);
      myText (0.66, 0.34, kBlack, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      if (p_lo != 1 || p_hi != numpbins) {
        if (!exclusive || iP != numpbins)
          myText (0.66, 0.275, kBlack, Form ("%g < #it{p}_{T}^{V} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
        else
          myText (0.66, 0.275, kBlack, Form ("%g < #left|#it{p}_{T}^{V}#right|", pbins[p_hi]), 0.04/uPadY);
        myText (0.66, 0.22, kBlack, period, 0.04/uPadY);
      }
      else
        myText (0.66, 0.22, kBlack, period, 0.04/uPadY);

      TH1D* sysCombinedRatioHist = new TH1D ("sysCombinedRatioHist", "", numetabins, etabins);
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

      const TString func = "[0] + [1]*((x-[4])/[5]) + [2]*(2*((x-[4])/[5])^2-1) + [3]*(4*((x-[4])/[5])^3-3*(x-[4])/[5])";// + [4]*(8*((x-[5])/[6])^4-8*((x-[5])/[6])^2+1)";
      vJetRatioFit = new TF1 ("vJetPtDataMCRatioFit", func, etabins[0], etabins[numetabins]);
      vJetRatioFit->SetParameter (0, 1);
      vJetRatioFit->SetParameter (1, 0);
      vJetRatioFit->SetParameter (2, 0);
      vJetRatioFit->SetParameter (3, 0);
      //vJetRatioFit->SetParameter (4, 0);
      vJetRatioFit->FixParameter (4, 0.5 * (etabins[numetabins] + etabins[0]));
      vJetRatioFit->FixParameter (5, 0.5 * (etabins[numetabins] - etabins[0]));
      sysCombinedRatioHist->Fit (vJetRatioFit, "RN0Q");

      const float chisq = vJetRatioFit->GetChisquare ();
      myText (0.175, 0.4, kBlack, Form ("#chi^{2} = %g", chisq), 0.04/dPadY);

      vJetRatioCI = new TH1D ("vJetPtDataMCRatioCI", "", numetabins, etabins);
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


    if (!plot_xjref) continue;

    /**** Plots xjref distributions, binned by ptref and eta ****/
    rawDistCanvas->cd ();
    gPad->SetLogx (false);
    for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta
      const int eta_lo = (iEta != numetabins ? iEta+1 : eta_lo_comb);
      const int eta_hi = (iEta != numetabins ? iEta+1 : eta_hi_comb);

      for (short iP = 0; iP <= numpbins; iP++) {
        const int p_lo = (iP != numpbins ? iP+1 : p_lo_comb);
        const int p_hi =  (iP != numpbins ? iP+1 : p_hi_comb);

        const Style_t dataStyle = 20;
        const Style_t mcStyle = 24;

        for (short iSpc = 0; iSpc < 2; iSpc++) {
          const Color_t dataColor = (iSpc == 0 ? 9 : (iSpc == 1 ? kMagenta+1 : 8));
          const Color_t mcColor = dataColor;
          // get data hist
          proj = Project2D ("", vJetHists[iPer][0][iSpc][0], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
          vJetHist = proj->ProjectionY ("vJetProjection", p_lo, p_hi);
          // rebin & scale data hist
          vJetHist->Rebin (rebinFactor);
          if (vJetHist->Integral () != 0) vJetHist->Scale (1./vJetHist->Integral ());
          // format data hist
          vJetHist->SetXTitle ("#it{x}_{J}^{ref}");
          vJetHist->SetYTitle ("Counts / Total");
          vJetHist->SetMarkerStyle (dataStyle);
          vJetHist->SetMarkerColor (dataColor);
          vJetHist->SetLineColor (dataColor);

          // get MC hist
          proj_mc = Project2D ("", vJetHists[iPer][1][iSpc][1], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
          vJetHist_mc = proj_mc->ProjectionY ("vJetProjection_mc", p_lo, p_hi);
          // rebin & scale MC hist
          vJetHist_mc->Rebin (rebinFactor);
          if (vJetHist_mc->Integral () != 0) vJetHist_mc->Scale (1./vJetHist_mc->Integral ()); 
          // format MC hist
          vJetHist_mc->SetMarkerStyle (mcStyle);
          vJetHist_mc->SetMarkerColor (mcColor);
          vJetHist_mc->SetLineColor (mcColor);

          const float max = std::max (vJetHist->GetMaximum (), vJetHist_mc->GetMaximum ());

          vJetHist->GetYaxis ()->SetRangeUser (0., 1.3*max);
          vJetHist_mc->GetYaxis ()->SetRangeUser (0., 1.3*max);
          
          vJetHist->GetXaxis ()->SetLabelSize (0.04);
          vJetHist->GetYaxis ()->SetLabelSize (0.04);
          vJetHist->GetXaxis ()->SetTitleSize (0.04);
          vJetHist->GetYaxis ()->SetTitleSize (0.04);
          vJetHist->GetYaxis ()->SetTitleOffset (1.1);

          vJetHist_mc->SetMarkerStyle (mcStyle);
          vJetHist_mc->SetMarkerColor (mcColor);
          vJetHist_mc->SetLineColor (mcColor);

          if (iSpc == 0)
            vJetHist_mc->DrawCopy ("e1 x0");
          else
            vJetHist_mc->DrawCopy ("same e1 x0");
          vJetHist->DrawCopy ("same e1 x0");

          float mean, mean_err, mean_mc, mean_mc_err, stddev;
          TF1 *fit_data = NULL, *fit_mc = NULL; 
          if (useGaussian) {
            mean = vJetHist->GetMean ();
            stddev = vJetHist->GetStdDev ();
            fit_data = new TF1 ("fit_data", "gaus (0)", mean - 2.8*stddev, mean + 2.8*stddev);
            //fit_data = new TF1 ("fit_data", "[0]*x^([1])*e^(-x/[2])", 0, 2);// mean - 1*stddev, mean + 2*stddev);
            //fit_data = new TF1 ("fit_data", "[0]/(sqrt(2*pi)*[2])  *e^(-(x-[1])^2/(2*[2]^2))  *(1+erf([3]/sqrt(2)*(x-[1])/[2]))", mean-3*stddev, mean+3*stddev);
            //fit_data->SetParameter (0, 1);
            //fit_data->SetParameter (1, mean);
            //fit_data->SetParameter (2, stddev);
            //fit_data->SetParameter (3, 0);
            vJetHist->Fit (fit_data, "Q0R");

            mean = vJetHist_mc->GetMean ();
            stddev = vJetHist_mc->GetStdDev ();
            fit_mc = new TF1 ("fit_mc", "gaus (0)", mean - 2.8*stddev, mean + 2.8*stddev);
            //fit_mc = new TF1 ("fit_mc", "[0]*x^([1])*e^(-x/[2])", 0, 2);//, mean - 1*stddev, mean + 2*stddev);
            //fit_mc = new TF1 ("fit_mc", "[0]/(sqrt(2*pi)*[2]) * e^(-(x-[1])^2/(2*[2]^2)) * (1+erf([3]/sqrt(2)*(x-[1])/[2]))", mean-3*stddev, mean+3*stddev);
            //fit_mc->SetParameter (0, 1);
            //fit_mc->SetParameter (1, mean);
            //fit_mc->SetParameter (2, stddev);
            //fit_mc->SetParameter (3, 0);
            vJetHist_mc->Fit (fit_mc, "Q0R");

            mean = fit_data->GetParameter (1);
            mean_err = fit_data->GetParError (1);
            //float m = fit_data->GetParameter (1);
            //float s = fit_data->GetParameter (2);
            //float a = fit_data->GetParameter (3);
            //float me = fit_data->GetParError (1);
            //float se = fit_data->GetParError (2);
            //float ae = fit_data->GetParError (3);
            //mean = m + sqrt (2/pi)*s*a / sqrt (1+a*a);
            //mean_err = sqrt (me*me + (2/pi)*se*se*a*a / (1+a*a) + (2/pi)*s*s*ae*ae / pow (1+a*a, 3));

            mean_mc = fit_mc->GetParameter (1);
            mean_mc_err = fit_mc->GetParError (1);
            //m = fit_mc->GetParameter (1);
            //s = fit_mc->GetParameter (2);
            //a = fit_mc->GetParameter (3);
            //me = fit_mc->GetParError (1);
            //se = fit_mc->GetParError (2);
            //ae = fit_mc->GetParError (3);
            //mean_mc = m + sqrt (2/pi)*s*a / sqrt (1+a*a);
            //mean_mc_err = sqrt (me*me + (2/pi)*se*se*a*a / (1+a*a) + (2/pi)*s*s*ae*ae / pow (1+a*a, 3));
          }
          else {
            mean = vJetHist->GetMean ();
            mean_err = vJetHist->GetMeanError ();
            mean_mc = vJetHist_mc->GetMean ();
            mean_mc_err = vJetHist_mc->GetMeanError ();
          }

          fit_data->SetLineColor (dataColor);
          fit_mc->SetLineColor (mcColor);
          ( (TF1*)fit_mc->Clone ())->Draw ("same");
          ( (TF1*)fit_data->Clone ())->Draw ("same");

          if (fit_data) delete fit_data;
          if (fit_mc) delete fit_mc;

          int countsData = 0, countsMC = 0;
          if (exclusive && iEta == numetabins) {
            countsData = vJetCounts[iPer][0][iSpc]->Integral () - vJetCounts[iPer][0][iSpc]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
            countsMC = vJetCounts[iPer][1][iSpc]->Integral () - vJetCounts[iPer][1][iSpc]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
          }
          else {
            countsData = vJetCounts[iPer][0][iSpc]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
            countsMC = vJetCounts[iPer][1][iSpc]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
          }

          myMarkerText (0.175, 0.88-0.12*iSpc, dataColor, dataStyle, Form ("%s + Jets, 2016 Data (%i events)", (iSpc==0?"#gamma":"Z"), countsData), 1.25, 0.04);
          myMarkerText (0.175, 0.82-0.12*iSpc, mcColor, mcStyle, Form ("%s + Jets, Pythia8 MC (%i events)", (iSpc==0?"#gamma":"Z"), countsMC), 1.25, 0.04);

          myText (0.655, 0.88-0.12*iSpc, dataColor, Form ("#mu_{Data} = %s", FormatMeasurement (mean, mean_err, 1)), 0.04);
          myText (0.655, 0.82-0.12*iSpc, dataColor, Form ("#mu_{MC} = %s", FormatMeasurement (mean_mc, mean_mc_err, 1)), 0.04);

          if (vJetHist) { delete vJetHist; vJetHist = NULL; }
          if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }

          if (proj) { delete proj; proj = NULL; }
          if (proj_mc) { delete proj_mc; proj_mc = NULL; }
        }

        myText (0.155, 0.43, kBlack, Form ("#gamma + Jet, %s", period), 0.04);
        myText (0.155, 0.37, kBlack, Form ("%g < #it{p}_{T}^{ref} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04);
        myText (0.155, 0.31, kBlack, Form ("%g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04);
        myText (0.155, 0.25, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);

        //bottomPad->cd ();
        //bottomPad->SetLogx (false);
        //vJetHist->Divide (vJetHist_mc);

        //vJetHist->SetYTitle ("Data / MC");
        //vJetHist->SetAxisRange (0.45, 1.65, "Y");
        //vJetHist->GetYaxis ()->SetNdivisions (605);
        //vJetHist->GetXaxis ()->SetTitleSize (0.04/dPadY);
        //vJetHist->GetYaxis ()->SetTitleSize (0.04/dPadY);
        //vJetHist->GetXaxis ()->SetTitleOffset (1);
        //vJetHist->GetYaxis ()->SetTitleOffset (1.1*dPadY);
        //vJetHist->GetYaxis ()->CenterTitle (true);
        //vJetHist->GetXaxis ()->SetLabelSize (0.04/dPadY);
        //vJetHist->GetYaxis ()->SetLabelSize (0.04/dPadY);
        //vJetHist->GetXaxis ()->SetTickLength (0.08);

        //vJetHist->DrawCopy ("e1 x0"); 
        //for (TLine* line : xlines) line->Draw ();

        plotName = Form ("xjref_dists/gamma_jet_iEta%i_iP%i.pdf", iEta, iP);
        rawDistCanvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));
        
      } // end loop over pT bins

    } // end loop over eta bins
  } // end loop over beam periods


  return;
}

} // end namespace
