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

  TH3D***** gJetHists = Get4DArray <TH3D*> (3, 3, 3, 3); // iPer, iData, iErr, iPCut
  TH2D***** gJetCounts = Get4DArray <TH2D*> (3, 3, 3, 2); // iPer, iData, iPCut, iWgt

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

      const TString name = Form ("gJetPtRatio_2016_%s_%s_%s_%s", dataType.Data (), error.Data (), period.Data (), pCut.Data ());
      gJetHists[iPer][iData][iErr][iPCut] = (TH3D*)inFile->Get (name);
     }

     for (short iWgt = 0; iWgt < 2; iWgt++) {
      const TString weight = (iWgt == 0 ? "unweighted" : "weighted");

      const TString name = Form ("gJetCounts_2016_%s_%s_%s_%s", dataType.Data (), period.Data (), pCut.Data (), weight.Data ());
      gJetCounts[iPer][iData][iPCut][iWgt] = (TH2D*)inFile->Get (name);
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
  TCanvas* perCombCanvas = new TCanvas ("perCombCanvas", "", 800, 600);

  TCanvas* prefCanvas = new TCanvas ("prefCanvas", "", 800, 600);
  prefCanvas->SetLeftMargin (-0.20);
//  prefCanvas->SetRightMargin(

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
  TH1D *vJetHist = NULL, *vJetHist_mc_overlay = NULL, *vJetHist_mc_signal = NULL, *vJetHist_lo = NULL, *vJetHist_hi = NULL, *vJetHist_rat = NULL, *vJetHist_rat_lo = NULL, *vJetHist_rat_hi = NULL;
  TH2D *proj = NULL, *proj_mc_overlay = NULL, *proj_mc_signal = NULL, *proj_lo = NULL, *proj_hi = NULL;
  TGraphAsymmErrors *vJetGraph = NULL, *vJetGraph_mc_overlay = NULL, *vJetGraph_mc_signal = NULL, *vJetGraph_sys = NULL, *vJetGraph_rat = NULL, *vJetGraph_rat_sys = NULL;
  char* plotName;

  TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");

  for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta

   const int eta_lo = (iEta != numetabins ? iEta+1 : eta_lo_comb);
   const int eta_hi = (iEta != numetabins ? iEta+1 : eta_hi_comb);

   for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
    const TString period = (iPer == 0 ? "Period A" : (iPer == 1 ? "Period B" : "Period A+B"));
    const TString perType = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));

    /**** Plots GammaJet info as a function of p_T^ref****/
    const Color_t dataColor = kBlack;
    const Style_t dataStyle = 20;
    const Style_t mcOverlayStyle = 33;
    //const Style_t mcSignalStyle = 34;

    topPad->cd ();
    topPad->SetLogx ();

    proj = Project2D ("", gJetHists[iPer][0][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
    TString name = Form ("gJetHist_%s_data_stat_signal_iEta%i", perType.Data (), iEta);
    vJetHist = GetProfileX (name, proj, numpbins, pbins, true);
    vJetHist->Write ();

    double middle = 0.05 * floor (20*vJetHist->GetMean (2)); //* 0.5 * (vJetHist->GetMaximum () + std::max (vJetHist->GetMinimum (), 0)));
    if (10 * middle != floor (10*middle)) middle += 0.05;
    double stddev = 0.05 * floor (20*vJetHist->GetStdDev (2));
    if (10 * stddev != floor (10*stddev)) stddev += 0.05;

    vJetGraph = make_graph (vJetHist);
    vJetGraph->GetYaxis ()->SetTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
    vJetGraph->GetYaxis ()->SetRangeUser (0.6, 1.4);
    //vJetGraph->GetYaxis ()->SetRangeUser (middle - 2*stddev, middle + 2*stddev);
    vJetGraph->SetMarkerStyle (dataStyle);
    vJetGraph->SetMarkerColor (dataColor);
    vJetGraph->SetLineColor (dataColor);
    vJetGraph->GetXaxis ()->SetLabelSize (0.04/uPadY);
    vJetGraph->GetYaxis ()->SetLabelSize (0.04/uPadY);
    vJetGraph->GetYaxis ()->SetTitleSize (0.04/uPadY);
    vJetGraph->GetYaxis ()->SetTitleOffset (uPadY);

    // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
    proj_lo = Project2D ("", gJetHists[iPer][0][0][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
    name = Form ("gJetHist_%s_data_sys_lo_signal_iEta%i", perType.Data (), iEta);
    vJetHist_lo = GetProfileX (name, proj_lo, numpbins, pbins, true);
    vJetHist_lo->Write ();

    proj_hi = Project2D ("", gJetHists[iPer][0][2][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
    name = Form ("gJetHist_%s_data_sys_hi_signal_iEta%i", perType.Data (), iEta);
    vJetHist_hi = GetProfileX (name, proj_hi, numpbins, pbins, true);
    vJetHist_hi->Write ();

    vJetGraph_sys = make_graph (vJetHist); // for plotting systematics
    CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    vJetGraph_sys->SetFillColor (dataColor);
    vJetGraph_sys->SetFillStyle (3001);

    proj_mc_overlay = Project2D ("", gJetHists[iPer][1][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
    name = Form ("gJetHist_%s_mc_stat_signal_iEta%i", perType.Data (), iEta);
    vJetHist_mc_overlay = GetProfileX (name, proj_mc_overlay, numpbins, pbins, true);
    vJetHist_mc_overlay->Write ();

    vJetGraph_mc_overlay = make_graph (vJetHist_mc_overlay);
    vJetGraph_mc_overlay->SetMarkerStyle (mcOverlayStyle);
    vJetGraph_mc_overlay->SetMarkerColor (mcOverlayColor);
    vJetGraph_mc_overlay->SetLineColor (mcOverlayColor);

    //if (iPer != 1 && !skipSignalMC) {
    // proj_mc_signal = Project2D ("", gJetHists[iPer][2][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
    // vJetHist_mc_signal = GetProfileX ("vJetHist_mc_signal", proj_mc_signal, numpbins, pbins, true);

    // vJetGraph_mc_signal = make_graph (vJetHist_mc_signal);
    // vJetGraph_mc_signal->SetMarkerStyle (mcSignalStyle);
    // vJetGraph_mc_signal->SetMarkerColor (mcSignalColor);
    // vJetGraph_mc_signal->SetLineColor (mcSignalColor);
    //}
    if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
    if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

    topPad->cd ();
    ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("ap");
    ( (TGraphAsymmErrors*)vJetGraph_mc_overlay->Clone ())->Draw ("p");
    ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

    if (iPer < 2) {
     perCombCanvas->cd ();
     perCombCanvas->SetLogx ();
     vJetGraph->SetMarkerStyle (iPer == 0 ? 20 : 24);
     vJetGraph_mc_overlay->SetMarkerStyle (iPer == 0 ? 33 : 27);
     deltaize (vJetGraph, 0.97+(0.04*iPer), true);
     deltaize (vJetGraph_mc_overlay, 0.99+(0.04*iPer), true);
     vJetGraph->GetYaxis ()->SetRangeUser (0.91, 1.08);
     if (iPer == 0)
      ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("ap");
     else
      ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("p");
     ( (TGraphAsymmErrors*)vJetGraph_mc_overlay->Clone ())->Draw ("p");
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");
    }

    int countsData = 0, countsMC = 0;//, countsMC_sig = 0;
    if (exclusive && iEta == numetabins) {
     countsData = gJetCounts[iPer][0][0][0]->Integral () - gJetCounts[iPer][0][0][0]->Integral (1, gJetCounts[iPer][0][0][0]->GetNbinsX(), eta_lo, eta_hi);
     countsMC = gJetCounts[iPer][1][0][0]->Integral () - gJetCounts[iPer][1][0][0]->Integral (1, gJetCounts[iPer][1][0][0]->GetNbinsX(), eta_lo, eta_hi);
     //if (!skipSignalMC)
     // countsMC_sig = gJetCounts[iPer][2][0][0]->Integral () - gJetCounts[iPer][2][0][0]->Integral (1, gJetCounts[iPer][2][0][0]->GetNbinsX(), eta_lo, eta_hi);
    }
    else {
     countsData = gJetCounts[iPer][0][0][0]->Integral (1, gJetCounts[iPer][0][0][0]->GetNbinsX(), eta_lo, eta_hi);
     countsMC = gJetCounts[iPer][1][0][0]->Integral (1, gJetCounts[iPer][1][0][0]->GetNbinsX(), eta_lo, eta_hi);
     //if (!skipSignalMC)
     // countsMC_sig = gJetCounts[iPer][2][0][0]->Integral (1, gJetCounts[iPer][2][0][0]->GetNbinsX(), eta_lo, eta_hi);
    }

    topPad->cd ();
    myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countsData), 1.25, 0.04/uPadY);
    myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
    //if (iPer != 1 && !skipSignalMC)
    // myMarkerText (0.175, 0.74, mcSignalColor, kFullCross, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", countsMC_sig), 1.25, 0.04/uPadY);
    myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
    if (eta_lo != 1 || eta_hi != numetabins)
     if (!exclusive || iEta != numetabins)
      myText (0.155, 0.15, dataColor, Form ("#gamma + Jet, %g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04/uPadY);
     else
      myText (0.155, 0.15, dataColor, Form ("#gamma + Jet, %g < #left|#eta_{det}^{Jet}#right|", etabins[eta_hi]), 0.04/uPadY);
    else
     myText (0.155, 0.15, dataColor, "#gamma + Jet", 0.04/uPadY);
    myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);

    if (iPer < 2) {
     perCombCanvas->cd ();
     myMarkerText (0.58, 0.44-0.05*iPer, dataColor, (iPer == 0 ? 20 : 24), Form ("%s Data (%i events)", period.Data (), countsData), 1.25, 0.04);
     myMarkerText (0.58, 0.34-0.05*iPer, mcOverlayColor, (iPer == 0 ? 33 : 27), Form ("%s MC (%i events)", period.Data (), countsMC), 1.25, 0.04);
     if (iPer == 0) {
      myText (0.175, 0.27, dataColor, "#bf{#it{ATLAS}} Internal", 0.04);
      if (eta_lo != 1 || eta_hi != numetabins)
       if (!exclusive || iEta != numetabins)
        myText (0.175, 0.22, dataColor, Form ("#gamma + Jet, %g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04);
       else
        myText (0.175, 0.22, dataColor, Form ("#gamma + Jet, %g < #left|#eta_{det}^{Jet}#right|", etabins[eta_hi]), 0.04);
      else
       myText (0.175, 0.22, dataColor, "#gamma + Jet", 0.04);
     }
    }

    bottomPad->cd ();
    bottomPad->SetLogx ();

    for (int iMC = 1; iMC < 3; iMC++) { // loops over overlay then signal MC
     if (iMC == 2 && (skipSignalMC || iPer == 1)) continue; // no signal only for period B!

     const TString mcType = (iMC == 1 ? "overlay":"signal");
     TH2D* proj_mc = (iMC == 1 ? proj_mc_overlay : proj_mc_signal);

     vJetHist_rat = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_%s_iEta%i", mcType.Data (), iEta)), proj, proj_mc, numpbins, pbins, true, "x");
     vJetHist_rat_lo = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_lo_%s_iEta%i", mcType.Data (), iEta)), proj_lo, proj_mc, numpbins, pbins, true, "x");
     vJetHist_rat_hi = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_hi_%s_iEta%i", mcType.Data (), iEta)), proj_hi, proj_mc, numpbins, pbins, true, "x");

     vJetGraph_rat_sys = make_graph (vJetHist_rat);
     CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
     if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
     //vJetGraph_rat_sys->SetFillColor (iMC == 1 ? dataColor:mcSignalColor);
     vJetGraph_rat_sys->SetFillColor (dataColor);
     vJetGraph_rat_sys->SetFillStyle (3001);

     vJetGraph_rat = make_graph (vJetHist_rat);
     vJetGraph_rat->GetXaxis ()->SetTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
     vJetGraph_rat->GetYaxis ()->SetTitle ("Data / MC");
     vJetGraph_rat->GetYaxis ()->SetRangeUser (0.91, 1.09);
     vJetGraph_rat->SetMarkerStyle (dataStyle);
     vJetGraph_rat->SetMarkerColor (dataColor);
     vJetGraph_rat->SetLineColor (iMC == 1 ? dataColor:mcSignalColor);
     vJetGraph_rat->GetYaxis ()->SetNdivisions (405);
     vJetGraph_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     vJetGraph_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     vJetGraph_rat->GetXaxis ()->SetTitleOffset (1);
     vJetGraph_rat->GetYaxis ()->SetTitleOffset (dPadY);
     vJetGraph_rat->GetYaxis ()->CenterTitle (true);
     vJetGraph_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     vJetGraph_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     vJetGraph_rat->GetXaxis ()->SetTickLength (0.08);

     if (iMC == 1)
      ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("ap");
     else
      ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("p");
     ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     for (TLine* line : glines) line->Draw ();

     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat) { delete vJetGraph_rat; vJetGraph_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
    }
    if (vJetHist) { delete vJetHist; vJetHist = NULL; }
    if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
    if (vJetHist_mc_signal) { delete vJetHist_mc_signal; vJetHist_mc_signal = NULL; }
    if (vJetGraph) { delete vJetGraph; vJetGraph = NULL; }
    if (vJetGraph_mc_overlay) { delete vJetGraph_mc_overlay; vJetGraph_mc_overlay = NULL; }
    if (vJetGraph_mc_signal) { delete vJetGraph_mc_signal; vJetGraph_mc_signal = NULL; }
    if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }

    if (proj) { delete proj; proj = NULL; }
    if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
    if (proj_mc_signal) { delete proj_mc_signal; proj_mc_signal = NULL; }
    if (proj_lo) { delete proj_lo; proj_lo = NULL; }
    if (proj_hi) { delete proj_hi; proj_hi = NULL; }

    if (iEta < numetabins) plotName = Form ("gamma_jet_iEta%i.pdf", iEta);
    else plotName = Form ("gamma_jet_iEta_combined.pdf");

    canvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));

    if (!plot_xjref || iEta == numetabins) continue;

    /**** Plots xjref distributions, binned by ptref ****/
    for (short iP = 0; iP < numpbins; iP++) {
     const double pref_lo = pbins[iP];
     const double pref_hi =  pbins[iP+1];

     const Color_t dataColor = kBlack;
     const Style_t dataStyle = 20;
     const Style_t mcOverlayStyle = 33;
     //const Style_t mcSignalStyle = 34;

     prefCanvas->cd ();
     gPad->SetLogx (false);

     // get data hist
     proj = Project2D ("", gJetHists[iPer][0][1][0], "x", "z", iEta+1, iEta+1);
     vJetHist = proj->ProjectionY ("vJetProjection", iP+1, iP+1);
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
     proj_mc_overlay = Project2D ("", gJetHists[iPer][1][1][0], "x", "z", iEta+1, iEta+1);
     vJetHist_mc_overlay = proj_mc_overlay->ProjectionY ("vJetProjection_mc", iP+1, iP+1);
     // rebin & scale MC hist
     vJetHist_mc_overlay->Rebin (rebinFactor);
     if (vJetHist_mc_overlay->Integral () != 0) vJetHist_mc_overlay->Scale (1./vJetHist_mc_overlay->Integral ()); 
     // format MC hist
     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

     const float max = std::max (vJetHist->GetMaximum (), vJetHist_mc_overlay->GetMaximum ());

     vJetHist->GetYaxis ()->SetRangeUser (0., 1.3*max);
     vJetHist_mc_overlay->GetYaxis ()->SetRangeUser (0., 1.3*max);
     
     vJetHist->GetXaxis ()->SetLabelSize (0.04);
     vJetHist->GetYaxis ()->SetLabelSize (0.04);
     vJetHist->GetXaxis ()->SetTitleSize (0.04);
     vJetHist->GetYaxis ()->SetTitleSize (0.04);
     vJetHist->GetYaxis ()->SetTitleOffset (1.1);

     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

     vJetHist_mc_overlay->DrawCopy ("e1 x0");
     vJetHist->DrawCopy ("same e1 x0");

     float mean, mean_err, mean_mc, mean_mc_err;//, stddev;
     TF1 *fit_data = NULL, *fit_mc = NULL; 
     if (useGaussian) {
      //mean = vJetHist->GetMean ();
      //stddev = vJetHist->GetStdDev ();
      fit_data = new TF1 ("fit_data", "gaus (0)", 0.5, 1.5);//mean - 1*stddev, mean + 2*stddev);
      //fit_data = new TF1 ("fit_data", "[0]*x^([1])*e^(-x/[2])", 0, 2);// mean - 1*stddev, mean + 2*stddev);
      //fit_data = new TF1 ("fit_data", "[0]/(sqrt(2*pi)*[2])  *e^(-(x-[1])^2/(2*[2]^2))  *(1+erf([3]/sqrt(2)*(x-[1])/[2]))", mean-3*stddev, mean+3*stddev);
      //fit_data->SetParameter (0, 1);
      //fit_data->SetParameter (1, mean);
      //fit_data->SetParameter (2, stddev);
      //fit_data->SetParameter (3, 0);
      vJetHist->Fit (fit_data, "Q0R");

      //mean = vJetHist_mc_overlay->GetMean ();
      //stddev = vJetHist_mc_overlay->GetStdDev ();
      fit_mc = new TF1 ("fit_mc", "gaus (0)", 0.5, 1.5);//mean - 1*stddev, mean + 2*stddev);
      //fit_mc = new TF1 ("fit_mc", "[0]*x^([1])*e^(-x/[2])", 0, 2);//, mean - 1*stddev, mean + 2*stddev);
      //fit_mc = new TF1 ("fit_mc", "[0]/(sqrt(2*pi)*[2]) * e^(-(x-[1])^2/(2*[2]^2)) * (1+erf([3]/sqrt(2)*(x-[1])/[2]))", mean-3*stddev, mean+3*stddev);
      //fit_mc->SetParameter (0, 1);
      //fit_mc->SetParameter (1, mean);
      //fit_mc->SetParameter (2, stddev);
      //fit_mc->SetParameter (3, 0);
      vJetHist_mc_overlay->Fit (fit_mc, "Q0R");

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
      mean_mc = vJetHist_mc_overlay->GetMean ();
      mean_mc_err = vJetHist_mc_overlay->GetMeanError ();
     }

     fit_data->SetLineColor (dataColor);
     fit_mc->SetLineColor (mcOverlayColor);
     ( (TF1*)fit_mc->Clone ())->Draw ("same");
     ( (TF1*)fit_data->Clone ())->Draw ("same");

     if (fit_data) delete fit_data;
     if (fit_mc) delete fit_mc;

     const int countData = gJetCounts[iPer][0][0][0]->GetBinContent (iP+1, iEta+1);
     const int countMC = gJetCounts[iPer][1][0][0]->GetBinContent (iP+1, iEta+1);

     myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countData), 1.25, 0.04);
     myMarkerText (0.175, 0.82, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV %s (%i events)", (runValidation ? "":"with #it{p}-Pb Overlay"), countMC), 1.25, 0.04);

     myText (0.155, 0.72, dataColor, Form ("<#it{x}_{J}^{ref}>^{data} = %.2f #pm %.2f", mean, mean_err), 0.04);
     myText (0.155, 0.66, dataColor, Form ("<#it{x}_{J}^{ref}>^{MC} = %.2f #pm %.2f", mean_mc, mean_mc_err), 0.04);

     myText (0.155, 0.43, dataColor, Form ("#gamma + Jet, %s", period.Data ()), 0.04);
     myText (0.155, 0.37, dataColor, Form ("%g < #it{p}_{T}^{ref} < %g", pbins[iP], pbins[iP+1]), 0.04);
     myText (0.155, 0.31, dataColor, Form ("%g < #eta_{det}^{Jet} < %g", etabins[iEta], etabins[iEta+1]), 0.04);
     myText (0.155, 0.25, dataColor, "#bf{#it{ATLAS}} Internal", 0.04);
     

     //bottomPad->cd ();
     //bottomPad->SetLogx (false);
     //vJetHist->Divide (vJetHist_mc_overlay);

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

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     plotName = Form ("pref_slices/gamma_jet_iEta%i_iP%i.pdf", iEta, iP);
     prefCanvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));
      
    } // end loop over pT bins

   } // end loop over periods
   if (iEta < numetabins) plotName = Form ("gamma_jet_iEta%i.pdf", iEta);
   else plotName = Form ("gamma_jet_iEta_combined.pdf");
   perCombCanvas->SaveAs (Form ("%s/PeriodSuperimposed/%s.pdf", plotPath.Data (), plotName));
  } // end loop over etabins


  /**** Now loop over pT bins and plot response as function of eta^jet ****/
  for (short iP = 0; iP <= numpbins; iP++) {

   const int p_lo = (iP != numpbins ? iP+1 : p_lo_comb);
   const int p_hi = (iP != numpbins ? iP+1 : p_hi_comb);

   for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
    const TString period = (iPer == 0 ? "Period A" : (iPer == 1 ? "Period B" : "Period A+B"));
    const TString perType = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));

    /**** Plots GammaJet info as a function of eta^jet****/
    const Color_t dataColor = kBlack;
    const Style_t dataStyle = 20;
    const Style_t mcOverlayStyle = 33;
    //const Style_t mcSignalStyle = 34;

    topPad->cd ();
    topPad->SetLogx (false);

    proj = Project2D ("", gJetHists[iPer][0][1][2], "y", "z", p_lo, p_hi);
    vJetHist = GetProfileX ("vJetHist", proj, numetabins, etabins, true);

    double middle = 0.05 * floor (20 * 0.5 * (vJetHist->GetMaximum () + vJetHist->GetMinimum ()));
    if (10 * middle != floor (10*middle)) middle += 0.05;

    vJetGraph = make_graph (vJetHist);
    vJetGraph->GetYaxis ()->SetTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
    vJetGraph->GetYaxis ()->SetRangeUser (middle - 0.35, middle + 0.35);
    vJetGraph->SetMarkerStyle (dataStyle);
    vJetGraph->SetMarkerColor (dataColor);
    vJetGraph->SetLineColor (dataColor);
    vJetGraph->GetXaxis ()->SetLabelSize (0.04/uPadY);
    vJetGraph->GetYaxis ()->SetLabelSize (0.04/uPadY);
    vJetGraph->GetYaxis ()->SetTitleSize (0.04/uPadY);
    vJetGraph->GetYaxis ()->SetTitleOffset (uPadY);

    // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
    proj_lo = Project2D ("", gJetHists[iPer][0][0][2], "y", "z", p_lo, p_hi);
    vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numetabins, etabins, true);

    proj_hi = Project2D ("", gJetHists[iPer][0][2][2], "y", "z", p_lo, p_hi);
    vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, true);

    vJetGraph_sys = make_graph (vJetHist); // for plotting systematics
    CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    vJetGraph_sys->SetFillColor (dataColor);
    vJetGraph_sys->SetFillStyle (3001);

    proj_mc_overlay = Project2D ("", gJetHists[iPer][1][1][2], "y", "z", p_lo, p_hi);
    vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numetabins, etabins, true);

    vJetGraph_mc_overlay = make_graph (vJetHist_mc_overlay);
    vJetGraph_mc_overlay->SetMarkerStyle (mcOverlayStyle);
    vJetGraph_mc_overlay->SetMarkerColor (mcOverlayColor);
    vJetGraph_mc_overlay->SetLineColor (mcOverlayColor);

    //if (iPer != 1 && !skipSignalMC) {
    // proj_mc_signal = Project2D ("", gJetHists[iPer][2][1][2], "y", "z", p_lo, p_hi);
    // vJetHist_mc_signal = GetProfileX ("vJetHist_mc_signal", proj_mc_signal, numetabins, etabins, true);

    // vJetGraph_mc_signal = make_graph (vJetHist_mc_signal);
    // vJetGraph_mc_signal->SetMarkerStyle (mcSignalStyle);
    // vJetGraph_mc_signal->SetMarkerColor (mcSignalColor);
    // vJetGraph_mc_signal->SetLineColor (mcSignalColor);
    //}
    if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
    if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

    //if (iPer != 1 && !skipSignalMC) vJetHist_mc_signal->DrawCopy ("same e1 x0");
    ( (TGraphAsymmErrors*)vJetGraph->Clone ())->Draw ("ap");
    ( (TGraphAsymmErrors*)vJetGraph_mc_overlay->Clone ())->Draw ("p");
    ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

    int countsData = 0, countsMC = 0, countsMC_sig = 0;
    countsData = gJetCounts[iPer][0][0][0]->Integral (p_lo, p_hi, 1, gJetCounts[iPer][0][0][0]->GetNbinsY());
    countsMC = gJetCounts[iPer][1][0][0]->Integral (p_lo, p_hi, 1, gJetCounts[iPer][1][0][0]->GetNbinsY());
    if (!skipSignalMC)
     countsMC_sig = gJetCounts[iPer][2][0][0]->Integral (p_lo, p_hi, 1, gJetCounts[iPer][2][0][0]->GetNbinsY());

    myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countsData), 1.25, 0.04/uPadY);
    myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
    //if (iPer != 1 && !skipSignalMC)
    // myMarkerText (0.175, 0.74, mcSignalColor, kFullCross, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", countsMC_sig), 1.25, 0.04/uPadY);
    //if (!skipOldInsitu) {
    // myMarkerText (0.65, 0.66, dataColor, kFullCircle, "2016 Insitu Factors", 1.25, 0.04/uPadY);
    // myMarkerText (0.65, 0.59, dataColor, kOpenCircle, "2015 Insitu Factors", 1.25, 0.04/uPadY);
    //}
    myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
    if (p_lo != 1 || p_hi != numpbins)
     myText (0.155, 0.15, dataColor, Form ("#gamma + Jet, %g < #it{p}_{T}^{#gamma} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
    else
     myText (0.155, 0.15, dataColor, "#gamma + Jet", 0.04/uPadY);
    myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);
    

    bottomPad->cd ();
    bottomPad->SetLogx (false);

    for (int iMC = 1; iMC < 3; iMC++) { // loops over overlay then signal MC
     if (iMC == 2 && (skipSignalMC || iPer == 1)) continue; // no signal only for period B!

     const TString mcType = (iMC == 1 ? "overlay":"signal");
     TH2D* proj_mc = (iMC == 1 ? proj_mc_overlay:proj_mc_signal);

     vJetHist_rat = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_iP%i", iP)), proj, proj_mc, numetabins, etabins, true, "x");
     vJetHist_rat_lo = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_lo_iP%i", iP)), proj_lo, proj_mc, numetabins, etabins, true, "x");
     vJetHist_rat_hi = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_hi_iP%i", iP)), proj_hi, proj_mc, numetabins, etabins, true, "x");

     vJetGraph_rat_sys = make_graph (vJetHist_rat);
     CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
     if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
     //vJetGraph_rat_sys->SetFillColor (iMC == 1 ? dataColor:mcSignalColor);
     vJetGraph_rat_sys->SetFillColor (dataColor);
     vJetGraph_rat_sys->SetFillStyle (3001);

     vJetGraph_rat = make_graph (vJetHist_rat);
     vJetGraph_rat->GetXaxis ()->SetTitle ("Jet #eta");
     vJetGraph_rat->GetYaxis ()->SetTitle ("Data / MC");
     vJetGraph_rat->GetYaxis ()->SetRangeUser (0.91, 1.09);
     vJetGraph_rat->SetMarkerStyle (dataStyle);
     vJetGraph_rat->SetMarkerColor (dataColor);
     vJetGraph_rat->SetLineColor (dataColor);
     vJetGraph_rat->GetYaxis ()->SetNdivisions (405);
     vJetGraph_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     vJetGraph_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     vJetGraph_rat->GetXaxis ()->SetTitleOffset (1);
     vJetGraph_rat->GetYaxis ()->SetTitleOffset (dPadY);
     vJetGraph_rat->GetYaxis ()->CenterTitle (true);
     vJetGraph_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     vJetGraph_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     vJetGraph_rat->GetXaxis ()->SetTickLength (0.08);

     if (iMC == 1)
      ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("ap");
     else
      ( (TGraphAsymmErrors*)vJetGraph_rat->Clone ())->Draw ("p");
     ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     for (TLine* line : getalines) line->Draw ();

     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat) { delete vJetGraph_rat; vJetGraph_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
    } // end loop over MC types
    if (vJetHist) { delete vJetHist; vJetHist = NULL; }
    if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
    if (vJetHist_mc_signal) { delete vJetHist_mc_signal; vJetHist_mc_signal = NULL; }
    if (vJetGraph) { delete vJetGraph; vJetGraph = NULL; }
    if (vJetGraph_mc_overlay) { delete vJetGraph_mc_overlay; vJetGraph_mc_overlay = NULL; }
    if (vJetGraph_mc_signal) { delete vJetGraph_mc_signal; vJetGraph_mc_signal = NULL; }
    if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }

    if (proj) { delete proj; proj = NULL; }
    if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
    if (proj_mc_signal) { delete proj_mc_signal; proj_mc_signal = NULL; }
    if (proj_lo) { delete proj_lo; proj_lo = NULL; }
    if (proj_hi) { delete proj_hi; proj_hi = NULL; }

    if (iP < numpbins) plotName = Form ("gamma_jet_iP%i.pdf", iP);
    else plotName = Form ("gamma_jet_iP_combined.pdf");

    canvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));

   } // end loop over beam periods
  } // end loop over pT bins

  outFile->Close ();

  return;
}

} // end namespace
