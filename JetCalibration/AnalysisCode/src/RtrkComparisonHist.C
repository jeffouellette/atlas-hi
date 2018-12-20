#include "RtrkComparisonHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TTree.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TVectorT.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

namespace JetCalibration {

void RtrkComparisonHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("RtrkComparison/", "JetCalibration/");

  TH3D**** jetRtrkHists = Get3DArray <TH3D*> (3, 3, 3); // iPer, iData, iErr
  TH2D*** jetRtrkCounts = Get2DArray <TH2D*> (3, 3);

  for (short iPer = 0; iPer < 3; iPer++) {
    const char* period = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));
 
    for (short iData = 0; iData < 3; iData++) { // iData is 0 for data, 1 for MC overlay, 2 for MC signal
      const char* data = (iData == 0 ? "data" : (iData == 1 ? "mc" : "mc_signal"));

      for (short iErr = 0; iErr < 3; iErr++) {
        const char* error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

        jetRtrkHists[iPer][iData][iErr] = new TH3D (Form ("jetRtrkDist_%s_%s_%s", data, error, period), "", numpbins, pbins, numtrketabins, trketabins, numrtrkbins, rtrkbins);
        jetRtrkHists[iPer][iData][iErr]->Sumw2 ();
      }
  
      jetRtrkCounts[iPer][iData] = new TH2D (Form ("jetRtrkCounts_%s_%s", data, period), "", numpbins, pbins, numtrketabins, trketabins);
      jetRtrkCounts[iPer][iData]->Sumw2 ();

    }
  }
  

  //////////////////////////////////////////////////////////////////////////////
  // Load analyzed TTrees
  //////////////////////////////////////////////////////////////////////////////
  float jpt = 0, jeta = 0, jphi = 0, je = 0, jpterr = 0, jtrk500 = 0, jtrk1000 = 0, ppt = 0, peta = 0, pphi = 0, purityFactor = 0, fcalScaleFactor = 0;
  double evtWeight = 0;
  bool isMC = false, isSignalOnly = false, isPeriodA = false, isTightPhoton = false;

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");
  TTree* inTree = (TTree*)inFile->Get ("RtrkTree");

  inTree->SetBranchAddress ("evt_weight", &evtWeight);
  inTree->SetBranchAddress ("purity_factor", &purityFactor);
  inTree->SetBranchAddress ("jet_pt", &jpt);
  inTree->SetBranchAddress ("jet_eta", &jeta);
  inTree->SetBranchAddress ("jet_phi", &jphi);
  inTree->SetBranchAddress ("jet_e", &je);
  inTree->SetBranchAddress ("jet_pt_sys", &jpterr);
  inTree->SetBranchAddress ("jet_SumPtTrkPt500", &jtrk500);
  inTree->SetBranchAddress ("jet_SumPtTrkPt1000", &jtrk1000);
  inTree->SetBranchAddress ("photon_pt", &ppt);
  inTree->SetBranchAddress ("photon_eta", &peta);
  inTree->SetBranchAddress ("photon_phi", &pphi);
  inTree->SetBranchAddress ("photon_tight", &isTightPhoton);
  inTree->SetBranchAddress ("fcalScaleFactor", &fcalScaleFactor);
  inTree->SetBranchAddress ("isMC", &isMC);
  inTree->SetBranchAddress ("isSignalOnly", &isSignalOnly);
  inTree->SetBranchAddress ("isPeriodA", &isPeriodA);


  //////////////////////////////////////////////////////////////////////////////
  // Fill desired histograms
  //////////////////////////////////////////////////////////////////////////////
  const long nJets = inTree->GetEntries ();
  for (long iJet = 0; iJet < nJets; iJet++) {
    inTree->GetEntry (iJet);

    if (ppt < 40)
      continue;
    if (!isTightPhoton)
      continue;

    const short iPer = isPeriodA ? 0 : 1;
    const short iMC = isMC ? (isSignalOnly ? 2 : 1) : 0;

    const float jtrk = jtrk500;

    // for Rtrk with jets
    jetRtrkHists[iPer][iMC][1]->Fill (jpt, jeta, jtrk/jpt, evtWeight*fcalScaleFactor);
    jetRtrkHists[2][iMC][1]->Fill (jpt, jeta, jtrk/jpt, evtWeight*fcalScaleFactor);
    jetRtrkCounts[iPer][iMC]->Fill (jpt, jeta);
    jetRtrkCounts[2][iMC]->Fill (jpt, jeta);

    jetRtrkHists[iPer][iMC][0]->Fill (jpt-jpterr, jeta, jtrk/jpt, evtWeight*fcalScaleFactor);
    jetRtrkHists[2][iMC][0]->Fill (jpt-jpterr, jeta, jtrk/jpt, evtWeight*fcalScaleFactor);
    jetRtrkHists[iPer][iMC][2]->Fill (jpt+jpterr, jeta, jtrk/jpt, evtWeight*fcalScaleFactor);
    jetRtrkHists[2][iMC][2]->Fill (jpt+jpterr, jeta, jtrk/jpt, evtWeight*fcalScaleFactor);

  }

  //////////////////////////////////////////////////////////////////////////////
  // Save histograms for interactive access
  //////////////////////////////////////////////////////////////////////////////
  TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");
  for (short iPer = 0; iPer < 3; iPer++) {
    for (short iData = 0; iData < 3; iData++) {
      // for Rtrk with jets
      for (short iErr = 0; iErr < 3; iErr++) {
        if (iErr != 1 && iData != 0)
          continue;
        jetRtrkHists[iPer][iData][iErr]->Write ();
      }
      jetRtrkCounts[iPer][iData]->Write ();
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // Plotting elements 
  //////////////////////////////////////////////////////////////////////////////
  TLine* ptlines[3] = {};
  TLine* etalines[3] = {};
  TLine* rtrklines[3] = {};
  for (short i = 0; i < 3; i++) {
    const float dpt = 0.05;
    const float deta = 0.05;
    const float drtrk = 0.05;

    ptlines[i] = new TLine (pbins[0], 1.0-dpt+dpt*i, pbins[numpbins], 1.0-dpt+dpt*i);
    etalines[i] = new TLine (trketabins[0], 1.0-deta+deta*i, trketabins[numtrketabins], 1.0-deta+deta*i);
    rtrklines[i] = new TLine (rtrkbins[0], 1.0-drtrk+drtrk*i, rtrkbins[numrtrkbins], 1.0-drtrk+drtrk*i);

    if (1.0-dpt+dpt*i == 1) ptlines[i]->SetLineStyle (1);
    else ptlines[i]->SetLineStyle (3);
    if (1.0-deta+deta*i == 1) etalines[i]->SetLineStyle (1);
    else etalines[i]->SetLineStyle (3);
    if (1.0-drtrk+drtrk*i == 1) rtrklines[i]->SetLineStyle (1);
    else rtrklines[i]->SetLineStyle (3);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Canvas definitions
  //////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas = new TCanvas ("canvas", "", 800, 800);
  const double padRatio = 1.2; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = 1.0/ (padRatio+1.0);
  const double uPadY = 1.0 - dPadY;
  TPad* topPad = new TPad ("topPad", "", 0, dPadY, 1, 1);
  TPad* bottomPad = new TPad ("bottomPad", "", 0, 0, 1, dPadY);
  topPad->SetTopMargin (0.04);
  topPad->SetBottomMargin (0);
  topPad->SetLeftMargin (0.15);
  topPad->SetRightMargin (0.04);
  bottomPad->SetTopMargin (0);
  bottomPad->SetBottomMargin (0.20);
  bottomPad->SetLeftMargin (0.15);
  bottomPad->SetRightMargin (0.04);
  topPad->Draw ();
  bottomPad->Draw ();

  //////////////////////////////////////////////////////////////////////////////
  // Define local histograms, graphs, etc.
  //////////////////////////////////////////////////////////////////////////////
  TH1D *jetRtrkHist = NULL, *jetRtrkHist_mc = NULL, *jetRtrkHist_mc_sig = NULL, *jetRtrkHist_lo = NULL, *jetRtrkHist_hi = NULL, *jetRtrkHist_rat = NULL, *jetRtrkHist_rat_sig = NULL, *jetRtrkHist_rat_lo = NULL, *jetRtrkHist_rat_hi = NULL;
  TH2D *proj2d = NULL, *proj2d_mc = NULL, *proj2d_mc_sig = NULL, *proj2d_lo = NULL, *proj2d_hi = NULL;
  TGraphAsymmErrors *jetRtrkGraph_sys = NULL, *jetRtrkGraph_rat_sys = NULL, *jetRtrkGraph_rat_sys_sig = NULL;

  const double* rtrk_los = linspace (0.2, 0.2, std::max (numpbins, numetabins));
  const double* rtrk_his = linspace (0.9, 0.9, std::max (numpbins, numetabins));

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
    const char* period = (iPer == 0 ? "Period A" : (iPer == 1 ? "Period B" : "Period A+B"));
    const char* per = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));

    for (short iEta = 0; iEta <= numtrketabins; iEta++) { // loop over bins in eta
      const short eta_lo = (iEta != numtrketabins ? iEta+1 : 1);
      const short eta_hi = (iEta != numtrketabins ? iEta+1 : numtrketabins);

      const Color_t dataColor = kBlack;
      const Style_t markerStyle = 20;

      topPad->cd ();
      topPad->SetLogx ();

      proj2d = Project2D ("", jetRtrkHists[iPer][0][1], "x", "z", eta_lo, eta_hi, exclusive && iEta == numtrketabins);
      proj2d->RebinY (rebinFactor);
      jetRtrkHist = GetProfileX ("jetRtrk_Hist", proj2d, numpbins, pbins, true, rtrk_los, rtrk_his);

      jetRtrkHist->SetMarkerStyle (markerStyle);
      jetRtrkHist->SetMarkerColor (dataColor);
      jetRtrkHist->SetLineColor (dataColor);

      // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
      proj2d_lo = Project2D ("", jetRtrkHists[iPer][0][0], "x", "z", eta_lo, eta_hi, exclusive && iEta == numtrketabins);
      proj2d_lo->RebinY (rebinFactor);
      jetRtrkHist_lo = GetProfileX ("jetRtrk_Hist_lo", proj2d_lo, numpbins, pbins, true, rtrk_los, rtrk_his);

      proj2d_hi = Project2D ("", jetRtrkHists[iPer][0][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numtrketabins);
      proj2d_hi->RebinY (rebinFactor);
      jetRtrkHist_hi = GetProfileX ("jetRtrk_Hist_hi", proj2d_hi, numpbins, pbins, true, rtrk_los, rtrk_his);

      jetRtrkGraph_sys = new TGraphAsymmErrors (jetRtrkHist); // for plotting systematics
      CalcSystematics (jetRtrkGraph_sys, jetRtrkHist, jetRtrkHist_hi, jetRtrkHist_lo);
      if (jetRtrkHist_lo) delete jetRtrkHist_lo;
      if (jetRtrkHist_hi) delete jetRtrkHist_hi;

      jetRtrkGraph_sys->SetFillColor (kBlack);
      jetRtrkGraph_sys->SetFillStyle (3001);

      proj2d_mc = Project2D ("", jetRtrkHists[iPer][1][1], "x", "z", eta_lo, eta_hi, exclusive && iEta == numtrketabins);
      proj2d_mc->RebinY (rebinFactor);
      jetRtrkHist_mc = GetProfileX ("jetRtrk_Hist_mc", proj2d_mc, numpbins, pbins, true, rtrk_los, rtrk_his);

      if (iPer != 1 && !skipSignalMC) {
        proj2d_mc_sig = Project2D ("", jetRtrkHists[iPer][2][1], "x", "z", eta_lo, eta_hi, exclusive && iEta == numtrketabins);
        proj2d_mc_sig->RebinY (rebinFactor);
        jetRtrkHist_mc_sig = GetProfileX ("jetRtrk_Hist_mc_sig", proj2d_mc_sig, numpbins, pbins, true, rtrk_los, rtrk_his);

        jetRtrkHist_mc_sig->SetMarkerStyle (markerStyle);
        jetRtrkHist_mc_sig->SetMarkerColor (mcSignalColor);
        jetRtrkHist_mc_sig->SetLineColor (mcSignalColor);
      }

      jetRtrkHist_mc->GetYaxis ()->SetTitle ("<#Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{Jet}>");
      jetRtrkHist_mc->GetYaxis ()->SetRangeUser (0., 2.0);
      jetRtrkHist_mc->SetMarkerStyle (markerStyle);
      jetRtrkHist_mc->SetMarkerColor (mcOverlayColor);
      jetRtrkHist_mc->SetLineColor (mcOverlayColor);
      jetRtrkHist_mc->GetXaxis ()->SetLabelSize (0.032/uPadY);
      jetRtrkHist_mc->GetYaxis ()->SetLabelSize (0.032/uPadY);
      jetRtrkHist_mc->GetYaxis ()->SetTitleSize (0.032/uPadY);
      jetRtrkHist_mc->GetYaxis ()->SetTitleOffset (1.5*uPadY);

      jetRtrkHist_mc->Draw ("e1 x0");
      if (iPer != 1 && !skipSignalMC)
        jetRtrkHist_mc_sig->Draw ("same e1 x0");
      jetRtrkHist->Draw ("same e1 x0");
      jetRtrkGraph_sys->Draw ("2");

      int nJetData = 0, nJetMC = 0, nJetMCSig = 0;

      if (exclusive && iEta == numtrketabins) {
        nJetData = jetRtrkCounts[iPer][0]->Integral () - jetRtrkCounts[iPer][0]->Integral (1, numpbins, eta_lo, eta_hi);
        nJetMC = jetRtrkCounts[iPer][1]->Integral () - jetRtrkCounts[iPer][1]->Integral (1, numpbins, eta_lo, eta_hi);
        nJetMCSig = jetRtrkCounts[iPer][2]->Integral () - jetRtrkCounts[iPer][2]->Integral (1, numpbins, eta_lo, eta_hi);
      }
      else {
        nJetData = jetRtrkCounts[iPer][0]->Integral (1, numpbins, eta_lo, eta_hi);
        nJetMC = jetRtrkCounts[iPer][1]->Integral (1, numpbins, eta_lo, eta_hi);
        nJetMCSig = jetRtrkCounts[iPer][2]->Integral (1, numpbins, eta_lo, eta_hi);
      }
      int ntext = 0;
      myMarkerText (0.225, 0.88-(ntext++)*0.07, dataColor, markerStyle, Form ("2016 Data (%i events)", nJetData), 1.25, 0.032/uPadY);
      myMarkerText (0.225, 0.88-(ntext++)*0.07, mcOverlayColor, markerStyle, Form ("Pythia8 MC + Overlay (%i events)", nJetMC), 1.25, 0.032/uPadY);
      if (iPer != 1 && !skipSignalMC)
        myMarkerText (0.225, 0.88-(ntext++)*0.07, mcSignalColor, markerStyle, Form ("Pythia8 MC (%i events)", nJetMCSig), 1.25, 0.032/uPadY);
      if (eta_lo != 1 || eta_hi != numtrketabins) {
        if (exclusive && iEta == numtrketabins) {
          if (fabs (trketabins[eta_lo-1]) == fabs (trketabins[eta_hi]))
            myText (0.195, 0.88-(ntext++)*0.07, kBlack, Form ("#left|#eta_{det}^{jet}#right| > %g", trketabins[eta_hi]), 0.032/uPadY);
          else
            myText (0.195, 0.88-(ntext++)*0.07, kBlack, Form ("#eta_{det}^{jet} < %g #union #eta_{det}^{jet} > %g", trketabins[eta_lo-1], trketabins[eta_hi]), 0.032/uPadY);
        }
        else
          myText (0.195, 0.88-(ntext++)*0.07, kBlack, Form ("%g < #eta_{det}^{jet} < %g", trketabins[eta_lo-1], trketabins[eta_hi]), 0.032/uPadY);
      }
      myText (0.195, 0.88-(ntext++)*0.07, kBlack, period, 0.032/uPadY);

      bottomPad->cd ();
      bottomPad->SetLogx ();

      jetRtrkHist_rat = GetDataOverMC (Form ("jetRtrk_DataMCRatio_%s_iEta%i", per, iEta), proj2d, proj2d_mc, numpbins, pbins, true, "x", rtrk_los, rtrk_his);
      jetRtrkHist_rat_lo = GetDataOverMC (Form ("jetRtrk_DataMCRatio_lo_%s_iEta%i", per, iEta), proj2d_lo, proj2d_mc, numpbins, pbins, true, "x", rtrk_los, rtrk_his);
      jetRtrkHist_rat_hi = GetDataOverMC (Form ("jetRtrk_DataMCRatio_hi_%s_iEta%i", per, iEta), proj2d_hi, proj2d_mc, numpbins, pbins, true, "x", rtrk_los, rtrk_his);

      jetRtrkGraph_rat_sys = new TGraphAsymmErrors (jetRtrkHist_rat);
      jetRtrkGraph_rat_sys->SetName (Form ("jetRtrk_DataMCRatio_sys_%s_iEta%i", per, iEta));
      CalcSystematics (jetRtrkGraph_rat_sys, jetRtrkHist_rat, jetRtrkHist_rat_hi, jetRtrkHist_rat_lo);
      if (jetRtrkHist_rat_lo) delete jetRtrkHist_rat_lo;
      if (jetRtrkHist_rat_hi) delete jetRtrkHist_rat_hi;

      jetRtrkGraph_rat_sys->SetFillColor (kBlack);
      jetRtrkGraph_rat_sys->SetFillStyle (3001);

      jetRtrkHist_rat->SetMarkerStyle (markerStyle);
      jetRtrkHist_rat->SetLineColor (dataColor);
      jetRtrkHist_rat->SetMarkerColor (dataColor);

      if (iPer != 1 && !skipSignalMC) {
        jetRtrkHist_rat_sig = GetDataOverMC (Form ("jetRtrk_DataMCRatio_Signal_%s_iEta%i", per, iEta), proj2d, proj2d_mc_sig, numpbins, pbins, true, "x", rtrk_los, rtrk_his);
        jetRtrkHist_rat_lo = GetDataOverMC (Form ("jetRtrk_DataMCRatio_Signal_lo_%s_iEta%i", per, iEta), proj2d_lo, proj2d_mc_sig, numpbins, pbins, true, "x", rtrk_los, rtrk_his);
        jetRtrkHist_rat_hi = GetDataOverMC (Form ("jetRtrk_DataMCRatio_Signal_hi_%s_iEta%i", per, iEta), proj2d_hi, proj2d_mc_sig, numpbins, pbins, true, "x", rtrk_los, rtrk_his);

        jetRtrkGraph_rat_sys_sig = new TGraphAsymmErrors (jetRtrkHist_rat_sig);
        jetRtrkGraph_rat_sys_sig->SetName (Form ("jetRtrk_DataMCRatio_Signal_sys_%s_iEta%i", per, iEta));
        CalcSystematics (jetRtrkGraph_rat_sys_sig, jetRtrkHist_rat_sig, jetRtrkHist_rat_hi, jetRtrkHist_rat_lo);
        if (jetRtrkHist_rat_lo) delete jetRtrkHist_rat_lo;
        if (jetRtrkHist_rat_hi) delete jetRtrkHist_rat_hi;

        jetRtrkGraph_rat_sys_sig->SetFillColor (mcSignalColor);
        jetRtrkGraph_rat_sys_sig->SetFillStyle (3001);

        jetRtrkHist_rat_sig->SetMarkerStyle (markerStyle);
        jetRtrkHist_rat_sig->SetLineColor (mcSignalColor);
        jetRtrkHist_rat_sig->SetMarkerColor (mcSignalColor);
      }

      jetRtrkHist_rat->GetXaxis ()->SetTitle ("#it{p}_{T}^{Jet} #left[GeV#right]");
      jetRtrkHist_rat->GetYaxis ()->SetTitle ("Data / MC");
      jetRtrkHist_rat->GetYaxis ()->SetRangeUser (0.94, 1.06);
      //jetRtrkHist_rat->GetYaxis ()->SetNdivisions (405);
      jetRtrkHist_rat->GetXaxis ()->SetTitleSize (0.032/dPadY);
      jetRtrkHist_rat->GetYaxis ()->SetTitleSize (0.032/dPadY);
      jetRtrkHist_rat->GetXaxis ()->SetTitleOffset (1);
      jetRtrkHist_rat->GetYaxis ()->SetTitleOffset (1.5*dPadY);
      jetRtrkHist_rat->GetYaxis ()->CenterTitle (true);
      jetRtrkHist_rat->GetXaxis ()->SetLabelSize (0.032/dPadY);
      jetRtrkHist_rat->GetYaxis ()->SetLabelSize (0.032/dPadY);
      jetRtrkHist_rat->GetXaxis ()->SetTickLength (0.08);

      jetRtrkHist_rat->Draw ("e1 x0");
      jetRtrkGraph_rat_sys->Draw ("2");
      if (iPer != 1 && !skipSignalMC) {
        jetRtrkHist_rat_sig->Draw ("same e1 x0");
        jetRtrkGraph_rat_sys_sig->Draw ("2");
      }
      for (TLine* line : ptlines) line->Draw ("same");

      TF1* constFit = new TF1 ("constFit", "[0]", pbins[0], pbins[numpbins]);
      jetRtrkHist_rat->Fit (constFit, "R0QN");
      myText (0.195, 0.33, kBlack, Form ("#mu = %.3f #pm %.3f", constFit->GetParameter (0), constFit->GetParError (0)), 0.032/dPadY);
      if (constFit) { delete constFit; constFit = NULL; }

      char* plotName;
      if (iEta < numtrketabins) plotName = Form ("jet_rtrk_iEta%i.pdf", iEta);
      else plotName = Form ("jet_rtrk_iEta_combined.pdf");

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
      }

      if (jetRtrkHist_rat) jetRtrkHist_rat->Write ();
      if (jetRtrkGraph_rat_sys) jetRtrkGraph_rat_sys->Write ();
      if (jetRtrkHist_rat_sig) jetRtrkHist_rat_sig->Write ();
      if (jetRtrkGraph_rat_sys_sig) jetRtrkGraph_rat_sys_sig->Write ();

      if (proj2d) { delete proj2d; proj2d = NULL; }
      if (proj2d_lo) { delete proj2d_lo; proj2d_lo = NULL; }
      if (proj2d_hi) { delete proj2d_hi; proj2d_hi = NULL; }
      if (proj2d_mc) { delete proj2d_mc; proj2d_mc = NULL; }
      if (proj2d_mc_sig) { delete proj2d_mc_sig; proj2d_mc_sig = NULL; }

      if (jetRtrkHist) { delete jetRtrkHist; jetRtrkHist = NULL; }
      if (jetRtrkHist_mc) { delete jetRtrkHist_mc; jetRtrkHist_mc = NULL; }
      if (jetRtrkHist_mc_sig) { delete jetRtrkHist_mc_sig; jetRtrkHist_mc_sig = NULL; }
      if (jetRtrkGraph_sys) { delete jetRtrkGraph_sys; jetRtrkGraph_sys = NULL; }
      if (jetRtrkHist_rat) { delete jetRtrkHist_rat; jetRtrkHist_rat = NULL; }
      if (jetRtrkHist_rat_sig) { delete jetRtrkHist_rat_sig; jetRtrkHist_rat_sig = NULL; }
      if (jetRtrkGraph_rat_sys) { delete jetRtrkGraph_rat_sys; jetRtrkGraph_rat_sys = NULL; }
      if (jetRtrkGraph_rat_sys_sig) { delete jetRtrkGraph_rat_sys_sig; jetRtrkGraph_rat_sys_sig = NULL; }
    } // end loop over eta bins

    for (short iP = 0; iP <= numpbins; iP++) { // loop over bins in p
      const short p_lo = (iP != numpbins ? iP+1 : p_lo_comb);
      const short p_hi = (iP != numpbins ? iP+1 : p_hi_comb);

      const Color_t dataColor = kBlack;
      const Style_t markerStyle = 20;

      topPad->cd ();
      topPad->SetLogx (false);

      proj2d = Project2D ("", jetRtrkHists[iPer][0][1], "y", "z", p_lo, p_hi);
      proj2d->RebinY (rebinFactor);
      jetRtrkHist = GetProfileX ("jetRtrk_Hist", proj2d, numtrketabins, trketabins, true, rtrk_los, rtrk_his);

      jetRtrkHist->SetMarkerStyle (markerStyle);
      jetRtrkHist->SetMarkerColor (dataColor);
      jetRtrkHist->SetLineColor (dataColor);

      // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
      proj2d_lo = Project2D ("", jetRtrkHists[iPer][0][0], "y", "z", p_lo, p_hi);
      proj2d_lo->RebinY (rebinFactor);
      jetRtrkHist_lo = GetProfileX ("jetRtrk_Hist_lo", proj2d_lo, numtrketabins, pbins, true, rtrk_los, rtrk_his);

      proj2d_hi = Project2D ("", jetRtrkHists[iPer][0][2], "y", "z", p_lo, p_hi);
      proj2d_hi->RebinY (rebinFactor);
      jetRtrkHist_hi = GetProfileX ("jetRtrk_Hist_hi", proj2d_hi, numtrketabins, pbins, true, rtrk_los, rtrk_his);

      jetRtrkGraph_sys = new TGraphAsymmErrors (jetRtrkHist); // for plotting systematics
      CalcSystematics (jetRtrkGraph_sys, jetRtrkHist, jetRtrkHist_hi, jetRtrkHist_lo);
      if (jetRtrkHist_lo) delete jetRtrkHist_lo;
      if (jetRtrkHist_hi) delete jetRtrkHist_hi;

      jetRtrkGraph_sys->SetFillColor (kBlack);
      jetRtrkGraph_sys->SetFillStyle (3001);

      proj2d_mc = Project2D ("", jetRtrkHists[iPer][1][1], "y", "z", p_lo, p_hi);
      proj2d_mc->RebinY (rebinFactor);
      jetRtrkHist_mc = GetProfileX ("jetRtrk_Hist_mc", proj2d_mc, numtrketabins, trketabins, true, rtrk_los, rtrk_his);

      jetRtrkHist_mc->SetMarkerStyle (markerStyle);
      jetRtrkHist_mc->SetMarkerColor (mcOverlayColor);
      jetRtrkHist_mc->SetLineColor (mcOverlayColor);

      proj2d_mc_sig = Project2D ("", jetRtrkHists[iPer][2][1], "y", "z", p_lo, p_hi);
      proj2d_mc_sig->RebinY (rebinFactor);
      jetRtrkHist_mc_sig = GetProfileX ("jetRtrk_Hist_mc_sig", proj2d_mc_sig, numtrketabins, trketabins, true, rtrk_los, rtrk_his);

      jetRtrkHist_mc_sig->SetMarkerStyle (markerStyle);
      jetRtrkHist_mc_sig->SetMarkerColor (mcSignalColor);
      jetRtrkHist_mc_sig->SetLineColor (mcSignalColor);

      jetRtrkHist_mc->GetYaxis ()->SetTitle ("<#Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{Jet}>");
      jetRtrkHist_mc->GetYaxis ()->SetRangeUser (0., 2.0);
      jetRtrkHist_mc->GetXaxis ()->SetLabelSize (0.032/uPadY);
      jetRtrkHist_mc->GetYaxis ()->SetLabelSize (0.032/uPadY);
      jetRtrkHist_mc->GetYaxis ()->SetTitleSize (0.032/uPadY);
      jetRtrkHist_mc->GetYaxis ()->SetTitleOffset (1.5*uPadY);

      jetRtrkHist_mc->Draw ("e1 x0");
      if (iPer != 1 && !skipSignalMC)
        jetRtrkHist_mc_sig->Draw ("same e1 x0");
      jetRtrkHist->Draw ("same e1 x0");
      jetRtrkGraph_sys->Draw ("2");

      const int nJetData = jetRtrkCounts[iPer][0]->Integral (p_lo, p_hi, 1, numtrketabins);
      const int nJetMC = jetRtrkCounts[iPer][1]->Integral (p_lo, p_hi, 1, numtrketabins);
      const int nJetMCSig = jetRtrkCounts[iPer][2]->Integral (p_lo, p_hi, 1, numtrketabins);
      int ntext = 0;
      myMarkerText (0.225, 0.88-(ntext++)*0.07, dataColor, markerStyle, Form ("2016 Data (%i events)", nJetData), 1.25, 0.032/uPadY);
      myMarkerText (0.225, 0.88-(ntext++)*0.07, mcOverlayColor, markerStyle, Form ("Pythia8 MC + Overlay (%i events)", nJetMC), 1.25, 0.032/uPadY);
      if (iPer != 1 && !skipSignalMC)
        myMarkerText (0.225, 0.88-(ntext++)*0.07, mcSignalColor, markerStyle, Form ("Pythia8 MC (%i events)", nJetMCSig), 1.25, 0.032/uPadY);
      if (p_lo != 1 || p_hi != numpbins)
        myText (0.195, 0.88-(ntext++)*0.07, kBlack, Form ("%g < #it{p}_{T}^{Jet} < %g", pbins[p_lo-1], pbins[p_hi]), 0.032/uPadY);
      myText (0.195, 0.88-(ntext++)*0.07, kBlack, period, 0.032/uPadY);

      bottomPad->cd ();
      bottomPad->SetLogx (false);

      jetRtrkHist_rat = GetDataOverMC (Form ("jetRtrk_DataMCRatio_%s_iP%i", per, iP), proj2d, proj2d_mc, numtrketabins, trketabins, true, "x", rtrk_los, rtrk_his);
      jetRtrkHist_rat_lo = GetDataOverMC (Form ("jetRtrk_DataMCRatio_lo_%s_iP%i", per, iP), proj2d_lo, proj2d_mc, numtrketabins, trketabins, true, "x", rtrk_los, rtrk_his);
      jetRtrkHist_rat_hi = GetDataOverMC (Form ("jetRtrk_DataMCRatio_hi_%s_iP%i", per, iP), proj2d_hi, proj2d_mc, numtrketabins, trketabins, true, "x", rtrk_los, rtrk_his);

      jetRtrkGraph_rat_sys = new TGraphAsymmErrors (jetRtrkHist_rat);
      jetRtrkGraph_rat_sys->SetName (Form ("jetRtrk_DataMCRatio_sys_%s_iP%i", per, iP));
      CalcSystematics (jetRtrkGraph_rat_sys, jetRtrkHist_rat, jetRtrkHist_rat_hi, jetRtrkHist_rat_lo);
      if (jetRtrkHist_rat_lo) delete jetRtrkHist_rat_lo;
      if (jetRtrkHist_rat_hi) delete jetRtrkHist_rat_hi;

      jetRtrkGraph_rat_sys->SetFillColor (dataColor);
      jetRtrkGraph_rat_sys->SetFillStyle (3001);

      jetRtrkHist_rat->SetMarkerStyle (markerStyle);
      jetRtrkHist_rat->SetLineColor (dataColor);
      jetRtrkHist_rat->SetMarkerColor (dataColor);

      if (iPer != 1 && !skipSignalMC) {
        jetRtrkHist_rat_sig = GetDataOverMC (Form ("jetRtrk_DataMCRatio_Signal_%s_iP%i", per, iP), proj2d, proj2d_mc_sig, numtrketabins, trketabins, true, "x", rtrk_los, rtrk_his);
        jetRtrkHist_rat_lo = GetDataOverMC (Form ("jetRtrk_DataMCRatio_Signal_lo_%s_iP%i", per, iP), proj2d_lo, proj2d_mc_sig, numtrketabins, trketabins, true, "x", rtrk_los, rtrk_his);
        jetRtrkHist_rat_hi = GetDataOverMC (Form ("jetRtrk_DataMCRatio_Signal_hi_%s_iP%i", per, iP), proj2d_hi, proj2d_mc_sig, numtrketabins, trketabins, true, "x", rtrk_los, rtrk_his);

        jetRtrkGraph_rat_sys_sig = new TGraphAsymmErrors (jetRtrkHist_rat_sig);
        jetRtrkGraph_rat_sys_sig->SetName (Form ("jetRtrk_DataMCRatio_Signal_sys_%s_iP%i", per, iP));
        CalcSystematics (jetRtrkGraph_rat_sys_sig, jetRtrkHist_rat_sig, jetRtrkHist_rat_hi, jetRtrkHist_rat_lo);
        if (jetRtrkHist_rat_lo) delete jetRtrkHist_rat_lo;
        if (jetRtrkHist_rat_hi) delete jetRtrkHist_rat_hi;

        jetRtrkGraph_rat_sys_sig->SetFillColor (mcSignalColor);
        jetRtrkGraph_rat_sys_sig->SetFillStyle (3001);

        jetRtrkHist_rat_sig->SetMarkerStyle (markerStyle);
        jetRtrkHist_rat_sig->SetLineColor (mcSignalColor);
        jetRtrkHist_rat_sig->SetMarkerColor (mcSignalColor);
      }

      jetRtrkHist_rat->GetXaxis ()->SetTitle ("#eta_{det}^{Jet}");
      jetRtrkHist_rat->GetYaxis ()->SetTitle ("Data / MC");
      jetRtrkHist_rat->GetYaxis ()->SetRangeUser (0.94, 1.06);
      //jetRtrkHist_rat->GetYaxis ()->SetNdivisions (405);
      jetRtrkHist_rat->GetXaxis ()->SetTitleSize (0.032/dPadY);
      jetRtrkHist_rat->GetYaxis ()->SetTitleSize (0.032/dPadY);
      jetRtrkHist_rat->GetXaxis ()->SetTitleOffset (1);
      jetRtrkHist_rat->GetYaxis ()->SetTitleOffset (1.5*dPadY);
      jetRtrkHist_rat->GetYaxis ()->CenterTitle (true);
      jetRtrkHist_rat->GetXaxis ()->SetLabelSize (0.032/dPadY);
      jetRtrkHist_rat->GetYaxis ()->SetLabelSize (0.032/dPadY);
      jetRtrkHist_rat->GetXaxis ()->SetTickLength (0.08);

      jetRtrkHist_rat->Draw ("e1 x0");
      jetRtrkGraph_rat_sys->Draw ("2");
      if (iPer != 1 && !skipSignalMC) {
        jetRtrkHist_rat_sig->Draw ("same e1 x0");
        jetRtrkGraph_rat_sys_sig->Draw ("2");
      }
      for (TLine* line : etalines) line->Draw ("same");

      TF1* constFit = new TF1 ("constFit", "[0]", etabins[0], etabins[numetabins]);
      jetRtrkHist_rat->Fit (constFit, "R0QN");
      myText (0.195, 0.33, kBlack, Form ("#mu = %.3f #pm %.3f", constFit->GetParameter (0), constFit->GetParError (0)), 0.032/dPadY);
      if (constFit) { delete constFit; constFit = NULL; }

      char* plotName;
      if (iP < numpbins) plotName = Form ("jet_rtrk_iP%i.pdf", iP);
      else plotName = Form ("jet_rtrk_iP_combined.pdf");

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
      }

      if (jetRtrkHist_rat) jetRtrkHist_rat->Write ();
      if (jetRtrkGraph_rat_sys) jetRtrkGraph_rat_sys->Write ();
      if (jetRtrkHist_rat_sig) jetRtrkHist_rat_sig->Write ();
      if (jetRtrkGraph_rat_sys_sig) jetRtrkGraph_rat_sys_sig->Write ();

      if (proj2d) { delete proj2d; proj2d = NULL; }
      if (proj2d_lo) { delete proj2d_lo; proj2d_lo = NULL; }
      if (proj2d_hi) { delete proj2d_hi; proj2d_hi = NULL; }
      if (proj2d_mc) { delete proj2d_mc; proj2d_mc = NULL; }
      if (proj2d_mc_sig) { delete proj2d_mc_sig; proj2d_mc_sig = NULL; }

      if (jetRtrkHist) { delete jetRtrkHist; jetRtrkHist = NULL; }
      if (jetRtrkHist_mc) { delete jetRtrkHist_mc; jetRtrkHist_mc = NULL; }
      if (jetRtrkHist_mc_sig) { delete jetRtrkHist_mc_sig; jetRtrkHist_mc_sig = NULL; }
      if (jetRtrkGraph_sys) { delete jetRtrkGraph_sys; jetRtrkGraph_sys = NULL; }
      if (jetRtrkHist_rat) { delete jetRtrkHist_rat; jetRtrkHist_rat = NULL; }
      if (jetRtrkHist_rat_sig) { delete jetRtrkHist_rat_sig; jetRtrkHist_rat_sig = NULL; }
      if (jetRtrkGraph_rat_sys) { delete jetRtrkGraph_rat_sys; jetRtrkGraph_rat_sys = NULL; }
      if (jetRtrkGraph_rat_sys_sig) { delete jetRtrkGraph_rat_sys_sig; jetRtrkGraph_rat_sys_sig = NULL; }
    } // end loop over pT bins


    /**** Plots rtrk distributions, binned by eta^jet ****/
    for (short iEta = 0; iEta <= numtrketabins; iEta++) {
      const short eta_lo = (iEta != numtrketabins ? iEta+1 : 1);
      const short eta_hi = (iEta != numtrketabins ? iEta+1 : 4);

      for (short iP = 0; iP <= numpbins; iP++) {
        const int p_lo = (iP != numpbins ? iP+1 : p_lo_comb);
        const int p_hi =  (iP != numpbins ? iP+1 : p_hi_comb);

        const Color_t dataColor = kBlack;
        const Style_t markerStyle = 20;

        topPad->cd ();
        topPad->SetLogx (0);

        proj2d = Project2D ("", jetRtrkHists[iPer][0][1], "y", "z", p_lo, p_hi);
        proj2d->RebinY (rebinFactor);
        jetRtrkHist = proj2d->ProjectionY ("jetRtrkHist", eta_lo, eta_hi);

        if (jetRtrkHist->Integral () != 0) jetRtrkHist->Scale (1./jetRtrkHist->Integral ());

        proj2d_mc = Project2D ("", jetRtrkHists[iPer][1][1], "y", "z", p_lo, p_hi);
        proj2d_mc->RebinY (rebinFactor);
        jetRtrkHist_mc = proj2d_mc->ProjectionY ("vJetProjection_mc", eta_lo, eta_hi);

        if (jetRtrkHist_mc->Integral () != 0) jetRtrkHist_mc->Scale (1./jetRtrkHist_mc->Integral ()); 

        jetRtrkHist_mc->SetMarkerStyle (markerStyle);
        jetRtrkHist_mc->SetMarkerColor (mcOverlayColor);
        jetRtrkHist_mc->SetLineColor (mcOverlayColor);

        jetRtrkHist_mc->GetYaxis ()->SetTitle ("Counts / Total");
        jetRtrkHist_mc->GetYaxis ()->SetRangeUser (0., 1.6*std::max (jetRtrkHist->GetMaximum (), jetRtrkHist_mc->GetMaximum ()));
        jetRtrkHist_mc->GetXaxis ()->SetLabelSize (0.032/uPadY);
        jetRtrkHist_mc->GetYaxis ()->SetLabelSize (0.032/uPadY);
        jetRtrkHist_mc->GetYaxis ()->SetTitleSize (0.032/uPadY);
        jetRtrkHist_mc->GetXaxis ()->SetTitleOffset (1.);
        jetRtrkHist_mc->GetYaxis ()->SetTitleOffset (1.5*uPadY);

        float mean, mean_err, mean_mc, mean_mc_err;

        TF1 *gaus_data = 0, *gaus_mc = 0;
  
        mean = jetRtrkHist->GetMean ();
        float stddev = jetRtrkHist->GetStdDev ();
        gaus_data = new TF1 ("gaus_data", "gaus (0)", rtrk_los[iP], rtrk_his[iP]);
        jetRtrkHist->Fit (gaus_data, "Q0R");

        mean = jetRtrkHist_mc->GetMean ();
        stddev = jetRtrkHist_mc->GetStdDev ();
        gaus_mc = new TF1 ("gaus_mc", "gaus (0)", rtrk_los[iP], rtrk_his[iP]);
        jetRtrkHist_mc->Fit (gaus_mc, "Q0R");

        mean = gaus_data->GetParameter (1);
        mean_err = gaus_data->GetParError (1);
        mean_mc = gaus_mc->GetParameter (1);
        mean_mc_err = gaus_mc->GetParError (1);

        jetRtrkHist->Scale (1./gaus_data->Integral (gaus_data->GetXmin (), gaus_data->GetXmax ()));
        jetRtrkHist_mc->Scale (1./gaus_mc->Integral (gaus_mc->GetXmin (), gaus_mc->GetXmax ()));
        gaus_data->SetParameter (0, gaus_data->GetParameter (0) / gaus_data->Integral (gaus_data->GetXmin (), gaus_data->GetXmax ()));
        gaus_mc->SetParameter (0, gaus_mc->GetParameter (0) / gaus_mc->Integral (gaus_mc->GetXmin (), gaus_mc->GetXmax ()));

        jetRtrkHist_mc->GetYaxis ()->SetRangeUser (0., 1.6* std::max (jetRtrkHist->GetMaximum (), jetRtrkHist_mc->GetMaximum ()));

        jetRtrkHist_mc->DrawCopy ("e1 x0");
        jetRtrkHist->DrawCopy ("same e1 x0");

        gaus_data->SetLineColor (dataColor);
        gaus_mc->SetLineColor (mcOverlayColor);
        gaus_data->Draw ("same");
        gaus_mc->Draw ("same");

        const int countsData = jetRtrkCounts[iPer][0]->Integral (p_lo, p_hi, eta_lo, eta_hi);
        const int countsMC = jetRtrkCounts[iPer][1]->Integral (p_lo, p_hi, eta_lo, eta_hi);

        myMarkerText (0.225, 0.88, dataColor, markerStyle, Form ("2016 Data (%i events)", countsData), 1.25, 0.032/uPadY);
        myMarkerText (0.225, 0.80, mcOverlayColor, markerStyle, Form ("Pythia8 MC (%i events)", countsMC), 1.25, 0.032/uPadY);

        myText (0.65, 0.88, dataColor, Form ("#mu_{data} = %s", FormatMeasurement (mean, mean_err)), 0.032/uPadY);
        myText (0.65, 0.80, dataColor, Form ("#mu_{MC} = %s", FormatMeasurement (mean_mc, mean_mc_err)), 0.032/uPadY);

        myText (0.68, 0.34, dataColor, "#bf{#it{ATLAS}} Internal", 0.032/uPadY);
        myText (0.68, 0.25, dataColor, period, 0.032/uPadY);
        myText (0.68, 0.16, dataColor, Form ("%g < #it{p}_{T}^{Jet} < %g", pbins[p_lo-1], pbins[p_hi]), 0.032/uPadY);
        myText (0.68, 0.08, dataColor, Form ("%g < #eta_{det}^{Jet} < %g", trketabins[eta_lo-1], trketabins[eta_hi]), 0.032/uPadY);

        bottomPad->cd ();
        bottomPad->SetLogx (false);
        jetRtrkHist->Divide (jetRtrkHist_mc);

        jetRtrkHist->GetXaxis ()->SetTitle ("r_{trk} = #Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{Jet}");
        jetRtrkHist->GetYaxis ()->SetTitle ("Data / MC");
        jetRtrkHist->GetYaxis ()->SetRangeUser (0.45, 1.65);
        jetRtrkHist->GetYaxis ()->SetNdivisions (605);
        jetRtrkHist->GetXaxis ()->SetTitleSize (0.032/dPadY);
        jetRtrkHist->GetYaxis ()->SetTitleSize (0.032/dPadY);
        jetRtrkHist->GetXaxis ()->SetTitleOffset (1);
        jetRtrkHist->GetYaxis ()->SetTitleOffset (1.5*dPadY);
        jetRtrkHist->GetYaxis ()->CenterTitle (true);
        jetRtrkHist->GetXaxis ()->SetLabelSize (0.032/dPadY);
        jetRtrkHist->GetYaxis ()->SetLabelSize (0.032/dPadY);
        jetRtrkHist->GetXaxis ()->SetTickLength (0.08);

        jetRtrkHist->DrawCopy ("e1 x0"); 
        for (TLine* line : rtrklines) line->Draw ();

        char* plotName = Form ("rtrk_dists/rtrk_iEta%i_iP%i.pdf", iEta, iP);
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

        if (gaus_data) delete gaus_data;
        if (gaus_mc) delete gaus_mc;

        if (jetRtrkHist) { delete jetRtrkHist; jetRtrkHist = NULL; }
        if (jetRtrkHist_mc) { delete jetRtrkHist_mc; jetRtrkHist_mc = NULL; }

        if (proj2d) { delete proj2d; proj2d = NULL; }
        if (proj2d_mc) { delete proj2d_mc; proj2d_mc = NULL; }

      } // end loop over pT bins
    } // end loop over eta bins
  } // end loop over periods

  outFile->Close ();
  if (outFile) { delete outFile; outFile = NULL; }

  return;
}

} // end namespace
