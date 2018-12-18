#include "PhotonRtrkComparisonHist.h"
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

void PhotonRtrkComparisonHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("PhotonRtrkComparison/", "JetCalibration/");

  TH3D***** photonRtrkHists = Get4DArray <TH3D*> (3, 3, 3, 3); // iPer, iData, iErr, iPCut
  TH3D***** photonRtrkCounts = Get4DArray <TH3D*> (3, 3, 3, 2); // iPer, iData, iPCut, iWgt

  for (short iPer = 0; iPer < 3; iPer++) {
    const char* period = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));
 
    for (short iData = 0; iData < 3; iData++) { // iData is 0 for data, 1 for MC
      const char* data = (iData == 0 ? "data" : (iData == 1 ? "mc" : "mc_signal"));

      for (short iPCut = 0; iPCut < 3; iPCut++) {
        const char* pCut = (iPCut == 0 ? "tight" : (iPCut == 1 ? "lnt" : "signal"));
  
        for (short iErr = 0; iErr < 3; iErr++) {
          const char* error = (iErr == 0 ? "sys_lo" : (iErr == 1 ? "stat" : "sys_hi"));

          photonRtrkHists[iPer][iData][iErr][iPCut] = new TH3D (Form ("photonRtrkDist_%s_%s_%s_%s", data, error, period, pCut), "", numpbins, pbins, numtrketabins, trketabins, numrtrkbins, rtrkbins);
          photonRtrkHists[iPer][iData][iErr][iPCut]->Sumw2 ();
        }
  
        for (short iWgt = 0; iWgt < 2; iWgt++) {
          const char* weight = (iWgt == 0 ? "unweighted" : "weighted");

          photonRtrkCounts[iPer][iData][iPCut][iWgt] = new TH3D (Form ("photonRtrkCounts_%s_%s_%s_%s", data, period, pCut, weight), "", numpbins, pbins, numtrketabins, trketabins, numphibins, phibins);
          photonRtrkCounts[iPer][iData][iPCut][iWgt]->Sumw2 ();
        }
      }
    }
  }
  

  //////////////////////////////////////////////////////////////////////////////
  // Load analyzed TTrees
  //////////////////////////////////////////////////////////////////////////////
  float /*jpt = 0, jeta = 0, jphi = 0, je = 0,*/ jpterr = 0, jtrk500 = 0, jtrk1000 = 0, ppt = 0, peta = 0, pphi = 0, purityFactor = 0, fcalScaleFactor = 0;
  double evtWeight = 0;
  bool isMC = false, isSignalOnly = false, isPeriodA = false, isTightPhoton = false;

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");
  TTree* inTree = (TTree*)inFile->Get ("RtrkTree");

  inTree->SetBranchAddress ("evt_weight", &evtWeight);
  inTree->SetBranchAddress ("purity_factor", &purityFactor);
  //inTree->SetBranchAddress ("jet_pt", &jpt);
  //inTree->SetBranchAddress ("jet_eta", &jeta);
  //inTree->SetBranchAddress ("jet_phi", &jphi);
  //inTree->SetBranchAddress ("jet_e", &je);
  //inTree->SetBranchAddress ("jet_pt_sys", &jpterr);
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

    const short iPer = isPeriodA ? 0 : 1;
    const short iMC = isMC ? (isSignalOnly ? 2 : 1) : 0;

    const float jtrk = jtrk1000;

    const short iPCut = isTightPhoton ? 0 : 1;
    photonRtrkHists[iPer][iMC][1][iPCut]->Fill (ppt, peta, jtrk/ppt, evtWeight*purityFactor);
    photonRtrkHists[2][iMC][1][iPCut]->Fill (ppt, peta, jtrk/ppt, evtWeight*purityFactor);

    photonRtrkCounts[iPer][iMC][iPCut][0]->Fill (ppt, peta, pphi);
    photonRtrkCounts[2][iMC][iPCut][0]->Fill (ppt, peta, pphi);
    photonRtrkCounts[iPer][iMC][iPCut][1]->Fill (ppt, peta, pphi, evtWeight);
    photonRtrkCounts[2][iMC][iPCut][1]->Fill (ppt, peta, pphi, evtWeight);

    photonRtrkHists[iPer][iMC][0][iPCut]->Fill (ppt-jpterr, peta, jtrk/ppt, evtWeight*purityFactor);
    photonRtrkHists[2][iMC][0][iPCut]->Fill (ppt-jpterr, peta, jtrk/ppt, evtWeight*purityFactor);
    photonRtrkHists[iPer][iMC][2][iPCut]->Fill (ppt+jpterr, peta, jtrk/ppt, evtWeight*purityFactor);
    photonRtrkHists[2][iMC][2][iPCut]->Fill (ppt+jpterr, peta, jtrk/ppt, evtWeight*purityFactor);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Subtract hadronic background from tight distributions
  //////////////////////////////////////////////////////////////////////////////
  for (short iPer = 0; iPer < 3; iPer++) {
    for (short iErr = 0; iErr < 3; iErr++) { // subtract for data
      if (subtractBackground) {
        SubtractLooseNonTight (photonRtrkHists[iPer][0][iErr][2], // signal
                               photonRtrkHists[iPer][0][iErr][0], // tight
                               photonRtrkHists[iPer][0][iErr][1], // loose non tight
                               photonRtrkCounts[iPer][0][0][1], // tight weighted counts
                               photonRtrkCounts[iPer][0][1][1]); // loose non tight weighted counts
      }
      else {
        delete photonRtrkHists[iPer][0][iErr][2];
        photonRtrkHists[iPer][0][iErr][2] = photonRtrkHists[iPer][0][iErr][0];
      }
    }
    // subtract for MC
    if (subtractBackground) {
      SubtractLooseNonTight (photonRtrkHists[iPer][1][1][2],
                             photonRtrkHists[iPer][1][1][0],
                             photonRtrkHists[iPer][1][1][1],
                             photonRtrkCounts[iPer][1][0][1],
                             photonRtrkCounts[iPer][1][1][1]);
    }
    else {
      delete photonRtrkHists[iPer][1][1][2];
      photonRtrkHists[iPer][1][1][2] = photonRtrkHists[iPer][1][1][0];
    }
    // subtract for MC signal
    if (iPer != 1 && !skipSignalMC) {
      if (subtractBackground) {
        SubtractLooseNonTight (photonRtrkHists[iPer][2][1][2],
                               photonRtrkHists[iPer][2][1][0],
                               photonRtrkHists[iPer][2][1][1],
                               photonRtrkCounts[iPer][2][0][1],
                               photonRtrkCounts[iPer][2][1][1]);
      }
      else {
        delete photonRtrkHists[iPer][2][1][2];
        photonRtrkHists[iPer][2][1][2] = photonRtrkHists[iPer][2][1][0];
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Save histograms for interactive access
  //////////////////////////////////////////////////////////////////////////////
  TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");
  for (short iPer = 0; iPer < 3; iPer++) {
    for (short iData = 0; iData < 3; iData++) {
      for (short iPCut = 0; iPCut < 3; iPCut++) {
        for (short iErr = 0; iErr < 3; iErr++) {
          if (iErr != 1 && iData != 0)
            continue;
          photonRtrkHists[iPer][iData][iErr][iPCut]->Write ();
        }
      
        for (short iWgt = 0; iWgt < 2; iWgt++) {  
          photonRtrkCounts[iPer][iData][iPCut][iWgt]->Write ();
        }
      }
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
  TCanvas* canvas = new TCanvas ("canvas", "", 800, 600);
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
  TH1D *photonRtrkHist = NULL, *photonRtrkHist_mc = NULL, *photonRtrkHist_mc_sig = NULL, *photonRtrkHist_lo = NULL, *photonRtrkHist_hi = NULL, *photonRtrkHist_rat = NULL, *photonRtrkHist_rat_sig = NULL, *photonRtrkHist_rat_lo = NULL, *photonRtrkHist_rat_hi = NULL;
  TH2D *proj2d = NULL, *proj2d_mc = NULL, *proj2d_mc_sig = NULL, *proj2d_lo = NULL, *proj2d_hi = NULL;
  TGraphAsymmErrors *photonRtrkGraph_sys = NULL, *photonRtrkGraph_rat_sys = NULL, *photonRtrkGraph_rat_sys_sig = NULL;

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

      proj2d = Project2D ("", photonRtrkHists[iPer][0][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numtrketabins);
      proj2d->RebinY (rebinFactor);
      photonRtrkHist = GetProfileX ("photonRtrk_Hist", proj2d, numpbins, pbins, false);

      photonRtrkHist->SetMarkerStyle (markerStyle);
      photonRtrkHist->SetMarkerColor (dataColor);
      photonRtrkHist->SetLineColor (dataColor);

      // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
      proj2d_lo = Project2D ("", photonRtrkHists[iPer][0][0][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numtrketabins);
      proj2d_lo->RebinY (rebinFactor);
      photonRtrkHist_lo = GetProfileX ("photonRtrk_Hist_lo", proj2d_lo, numpbins, pbins, false);

      proj2d_hi = Project2D ("", photonRtrkHists[iPer][0][2][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numtrketabins);
      proj2d_hi->RebinY (rebinFactor);
      photonRtrkHist_hi = GetProfileX ("photonRtrk_Hist_hi", proj2d_hi, numpbins, pbins, false);

      photonRtrkGraph_sys = new TGraphAsymmErrors (photonRtrkHist); // for plotting systematics
      CalcSystematics (photonRtrkGraph_sys, photonRtrkHist, photonRtrkHist_hi, photonRtrkHist_lo);
      if (photonRtrkHist_lo) delete photonRtrkHist_lo;
      if (photonRtrkHist_hi) delete photonRtrkHist_hi;

      photonRtrkGraph_sys->SetFillColor (kBlack);
      photonRtrkGraph_sys->SetFillStyle (3001);

      proj2d_mc = Project2D ("", photonRtrkHists[iPer][1][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numtrketabins);
      proj2d_mc->RebinY (rebinFactor);
      photonRtrkHist_mc = GetProfileX ("photonRtrk_Hist_mc", proj2d_mc, numpbins, pbins, false);

      photonRtrkHist_mc->SetMarkerStyle (markerStyle);
      photonRtrkHist_mc->SetMarkerColor (mcOverlayColor);
      photonRtrkHist_mc->SetLineColor (mcOverlayColor);

      if (iPer != 1 && !skipSignalMC) {
        proj2d_mc_sig = Project2D ("", photonRtrkHists[iPer][2][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numtrketabins);
        proj2d_mc_sig->RebinY (rebinFactor);
        photonRtrkHist_mc_sig = GetProfileX ("photonRtrk_Hist_mc_sig", proj2d_mc_sig, numpbins, pbins, false);

        photonRtrkHist_mc_sig->SetMarkerStyle (markerStyle);
        photonRtrkHist_mc_sig->SetMarkerColor (mcSignalColor);
        photonRtrkHist_mc_sig->SetLineColor (mcSignalColor);
      }

      photonRtrkHist_mc->GetYaxis ()->SetTitle ("< #Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{#gamma}>");
      photonRtrkHist_mc->GetYaxis ()->SetRangeUser (0., 2.0);
      photonRtrkHist_mc->GetXaxis ()->SetLabelSize (0.04/uPadY);
      photonRtrkHist_mc->GetYaxis ()->SetLabelSize (0.04/uPadY);
      photonRtrkHist_mc->GetYaxis ()->SetTitleSize (0.04/uPadY);
      photonRtrkHist_mc->GetYaxis ()->SetTitleOffset (1.5*uPadY);

      photonRtrkHist_mc->Draw ("e1 x0");
      if (iPer != 1 && !skipSignalMC)
        photonRtrkHist_mc_sig->Draw ("same e1 x0");
      photonRtrkHist->Draw ("same e1 x0");
      photonRtrkGraph_sys->Draw ("2");

      int countsData = 0, countsMC = 0, countsMCSig = 0;
      if (exclusive && iEta == numtrketabins) {
        countsData = photonRtrkCounts[iPer][0][0][0]->Integral () - photonRtrkCounts[iPer][0][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
        countsMC = photonRtrkCounts[iPer][1][0][0]->Integral () - photonRtrkCounts[iPer][1][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
        if (iPer != 1 && !skipSignalMC)
          countsMCSig = photonRtrkCounts[iPer][2][0][0]->Integral () - photonRtrkCounts[iPer][2][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
      }
      else {
        countsData = photonRtrkCounts[iPer][0][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
        countsMC = photonRtrkCounts[iPer][1][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
        if (iPer != 1 && !skipSignalMC)
          countsMCSig = photonRtrkCounts[iPer][2][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
      }
      myMarkerText (0.225, 0.88, dataColor, markerStyle, Form ("2016 Data (%i events)", countsData), 1.25, 0.04/uPadY);
      myMarkerText (0.225, 0.81, mcOverlayColor, markerStyle, Form ("Pythia8 MC + Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
      if (iPer != 1 && !skipSignalMC)
        myMarkerText (0.225, 0.74, mcSignalColor, markerStyle, Form ("Pythia8 MC (%i events)", countsMCSig), 1.25, 0.04/uPadY);
      if (eta_lo != 1 || eta_hi != numtrketabins) {
        if (exclusive && iEta == numtrketabins) {
          if (fabs (trketabins[eta_lo-1]) == fabs (trketabins[eta_hi]))
            myText (0.195, 0.58, kBlack, Form ("#left|#eta_{det}^{#gamma}#right| > %g", trketabins[eta_hi]), 0.04/uPadY);
          else
            myText (0.195, 0.58, kBlack, Form ("#eta_{det}^{#gamma} < %g #union #eta_{det}^{#gamma} > %g", trketabins[eta_lo-1], trketabins[eta_hi]), 0.04/uPadY);
        }
        else
          myText (0.195, 0.58, kBlack, Form ("%g < #eta_{det}^{#gamma} < %g", trketabins[eta_lo-1], trketabins[eta_hi]), 0.04/uPadY);
      }
      myText (0.195, 0.66, kBlack, period, 0.04/uPadY);

      bottomPad->cd ();
      bottomPad->SetLogx ();

      photonRtrkHist_rat = GetDataOverMC (Form ("photonRtrk_DataMCRatio_%s_iEta%i", per, iEta), proj2d, proj2d_mc, numpbins, pbins, false, "x");
      photonRtrkHist_rat_lo = GetDataOverMC (Form ("photonRtrk_DataMCRatio_lo_%s_iEta%i", per, iEta), proj2d_lo, proj2d_mc, numpbins, pbins, false, "x");
      photonRtrkHist_rat_hi = GetDataOverMC (Form ("photonRtrk_DataMCRatio_hi_%s_iEta%i", per, iEta), proj2d_hi, proj2d_mc, numpbins, pbins, false, "x");

      photonRtrkGraph_rat_sys = new TGraphAsymmErrors (photonRtrkHist_rat);
      photonRtrkGraph_rat_sys->SetName (Form ("photonRtrk_DataMCRatio_sys_%s_iEta%i", per, iEta));
      CalcSystematics (photonRtrkGraph_rat_sys, photonRtrkHist_rat, photonRtrkHist_rat_hi, photonRtrkHist_rat_lo);
      if (photonRtrkHist_rat_lo) delete photonRtrkHist_rat_lo;
      if (photonRtrkHist_rat_hi) delete photonRtrkHist_rat_hi;

      photonRtrkGraph_rat_sys->SetFillColor (dataColor);
      photonRtrkGraph_rat_sys->SetFillStyle (3001);

      photonRtrkHist_rat->SetMarkerStyle (markerStyle);
      photonRtrkHist_rat->SetLineColor (dataColor);
      photonRtrkHist_rat->SetMarkerColor (dataColor);

      if (iPer != 1 && !skipSignalMC) {
        photonRtrkHist_rat_sig = GetDataOverMC (Form ("photonRtrk_DataMCRatio_Signal_%s_iEta%i", per, iEta), proj2d, proj2d_mc_sig, numpbins, pbins, false, "x");
        photonRtrkHist_rat_lo = GetDataOverMC (Form ("photonRtrk_DataMCRatio_Signal_lo_%s_iEta%i", per, iEta), proj2d_lo, proj2d_mc_sig, numpbins, pbins, false, "x");
        photonRtrkHist_rat_hi = GetDataOverMC (Form ("photonRtrk_DataMCRatio_Signal_hi_%s_iEta%i", per, iEta), proj2d_hi, proj2d_mc_sig, numpbins, pbins, false, "x");

        photonRtrkGraph_rat_sys_sig = new TGraphAsymmErrors (photonRtrkHist_rat_sig);
        photonRtrkGraph_rat_sys_sig->SetName (Form ("photonRtrk_DataMCRatio_Signal_sys_%s_iEta%i", per, iEta));
        CalcSystematics (photonRtrkGraph_rat_sys_sig, photonRtrkHist_rat_sig, photonRtrkHist_rat_hi, photonRtrkHist_rat_lo);
        if (photonRtrkHist_rat_lo) delete photonRtrkHist_rat_lo;
        if (photonRtrkHist_rat_hi) delete photonRtrkHist_rat_hi;

        photonRtrkGraph_rat_sys_sig->SetFillColor (mcSignalColor);
        photonRtrkGraph_rat_sys_sig->SetFillStyle (3001);

        photonRtrkHist_rat_sig->SetMarkerStyle (markerStyle);
        photonRtrkHist_rat_sig->SetLineColor (mcSignalColor);
        photonRtrkHist_rat_sig->SetMarkerColor (mcSignalColor);
      }

      photonRtrkHist_rat->GetXaxis ()->SetTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
      photonRtrkHist_rat->GetYaxis ()->SetTitle ("Data / MC");
      photonRtrkHist_rat->GetYaxis ()->SetRangeUser (0.94, 1.06);
      photonRtrkHist_rat->GetYaxis ()->SetNdivisions (405);
      photonRtrkHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
      photonRtrkHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
      photonRtrkHist_rat->GetXaxis ()->SetTitleOffset (1);
      photonRtrkHist_rat->GetYaxis ()->SetTitleOffset (1.5*dPadY);
      photonRtrkHist_rat->GetYaxis ()->CenterTitle (true);
      photonRtrkHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
      photonRtrkHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
      photonRtrkHist_rat->GetXaxis ()->SetTickLength (0.08);

      photonRtrkHist_rat->Draw ("e1 x0");
      photonRtrkGraph_rat_sys->Draw ("2");
      if (iPer != 1 && !skipSignalMC) {
        photonRtrkHist_rat_sig->Draw ("same e1 x0");
        photonRtrkGraph_rat_sys_sig->Draw ("2");
      }
      for (TLine* line : ptlines) line->Draw ("same");

      TF1* constFit = new TF1 ("constFit", "[0]", pbins[0], pbins[numpbins]);
      photonRtrkHist_rat->Fit (constFit, "R0QN");
      myText (0.195, 0.33, kBlack, Form ("#mu = %.3f #pm %.3f", constFit->GetParameter (0), constFit->GetParError (0)), 0.04/dPadY);
      if (constFit) { delete constFit; constFit = NULL; }

      char* plotName;
      if (iEta < numtrketabins) plotName = Form ("photon_rtrk_iEta%i.pdf", iEta);
      else plotName = Form ("photon_rtrk_iEta_combined.pdf");

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

      if (photonRtrkHist_rat) photonRtrkHist_rat->Write ();
      if (photonRtrkGraph_rat_sys) photonRtrkGraph_rat_sys->Write ();
      if (photonRtrkHist_rat_sig) photonRtrkHist_rat_sig->Write ();
      if (photonRtrkGraph_rat_sys_sig) photonRtrkGraph_rat_sys_sig->Write ();

      if (proj2d) { delete proj2d; proj2d = NULL; }
      if (proj2d_lo) { delete proj2d_lo; proj2d_lo = NULL; }
      if (proj2d_hi) { delete proj2d_hi; proj2d_hi = NULL; }
      if (proj2d_mc) { delete proj2d_mc; proj2d_mc = NULL; }
      if (proj2d_mc_sig) { delete proj2d_mc_sig; proj2d_mc_sig = NULL; }

      if (photonRtrkHist) { delete photonRtrkHist; photonRtrkHist = NULL; }
      if (photonRtrkHist_mc) { delete photonRtrkHist_mc; photonRtrkHist_mc = NULL; }
      if (photonRtrkHist_mc_sig) { delete photonRtrkHist_mc_sig; photonRtrkHist_mc_sig = NULL; }
      if (photonRtrkGraph_sys) { delete photonRtrkGraph_sys; photonRtrkGraph_sys = NULL; }
      if (photonRtrkHist_rat) { delete photonRtrkHist_rat; photonRtrkHist_rat = NULL; }
      if (photonRtrkHist_rat_sig) { delete photonRtrkHist_rat_sig; photonRtrkHist_rat_sig = NULL; }
      if (photonRtrkGraph_rat_sys) { delete photonRtrkGraph_rat_sys; photonRtrkGraph_rat_sys = NULL; }
      if (photonRtrkGraph_rat_sys_sig) { delete photonRtrkGraph_rat_sys_sig; photonRtrkGraph_rat_sys_sig = NULL; }
    } // end loop over eta bins

    for (short iP = 0; iP <= numpbins; iP++) { // loop over bins in p
      const short p_lo = (iP != numpbins ? iP+1 : p_lo_comb);
      const short p_hi = (iP != numpbins ? iP+1 : p_hi_comb);

      const Color_t dataColor = kBlack;
      const Style_t markerStyle = 20;

      topPad->cd ();
      topPad->SetLogx (false);

      proj2d = Project2D ("", photonRtrkHists[iPer][0][1][2], "y", "z", p_lo, p_hi);
      proj2d->RebinY (rebinFactor);
      photonRtrkHist = GetProfileX ("photonRtrk_Hist", proj2d, numtrketabins, trketabins, false);

      photonRtrkHist->SetMarkerStyle (markerStyle);
      photonRtrkHist->SetMarkerColor (dataColor);
      photonRtrkHist->SetLineColor (dataColor);

      // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
      proj2d_lo = Project2D ("", photonRtrkHists[iPer][0][0][2], "y", "z", p_lo, p_hi);
      proj2d_lo->RebinY (rebinFactor);
      photonRtrkHist_lo = GetProfileX ("photonRtrk_Hist_lo", proj2d_lo, numtrketabins, pbins, false);

      proj2d_hi = Project2D ("", photonRtrkHists[iPer][0][2][2], "y", "z", p_lo, p_hi);
      proj2d_hi->RebinY (rebinFactor);
      photonRtrkHist_hi = GetProfileX ("photonRtrk_Hist_hi", proj2d_hi, numtrketabins, pbins, false);

      photonRtrkGraph_sys = new TGraphAsymmErrors (photonRtrkHist); // for plotting systematics
      CalcSystematics (photonRtrkGraph_sys, photonRtrkHist, photonRtrkHist_hi, photonRtrkHist_lo);
      if (photonRtrkHist_lo) delete photonRtrkHist_lo;
      if (photonRtrkHist_hi) delete photonRtrkHist_hi;

      photonRtrkGraph_sys->SetFillColor (kBlack);
      photonRtrkGraph_sys->SetFillStyle (3001);

      proj2d_mc = Project2D ("", photonRtrkHists[iPer][1][1][2], "y", "z", p_lo, p_hi);
      proj2d_mc->RebinY (rebinFactor);
      photonRtrkHist_mc = GetProfileX ("photonRtrk_Hist_mc", proj2d_mc, numtrketabins, trketabins, false);

      photonRtrkHist_mc->SetMarkerStyle (markerStyle);
      photonRtrkHist_mc->SetMarkerColor (mcOverlayColor);
      photonRtrkHist_mc->SetLineColor (mcOverlayColor);

      if (iPer != 1 && !skipSignalMC) {
        proj2d_mc_sig = Project2D ("", photonRtrkHists[iPer][2][1][2], "y", "z", p_lo, p_hi);
        proj2d_mc_sig->RebinY (rebinFactor);
        photonRtrkHist_mc_sig = GetProfileX ("photonRtrk_Hist_mc_sig", proj2d_mc_sig, numtrketabins, trketabins, false);

        photonRtrkHist_mc_sig->SetMarkerStyle (markerStyle);
        photonRtrkHist_mc_sig->SetMarkerColor (mcSignalColor);
        photonRtrkHist_mc_sig->SetLineColor (mcSignalColor);
      }

      photonRtrkHist_mc->GetYaxis ()->SetTitle ("< #Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{#gamma}>");
      photonRtrkHist_mc->GetYaxis ()->SetRangeUser (0., 2.0);
      photonRtrkHist_mc->GetXaxis ()->SetLabelSize (0.04/uPadY);
      photonRtrkHist_mc->GetYaxis ()->SetLabelSize (0.04/uPadY);
      photonRtrkHist_mc->GetYaxis ()->SetTitleSize (0.04/uPadY);
      photonRtrkHist_mc->GetYaxis ()->SetTitleOffset (1.5*uPadY);

      photonRtrkHist_mc->Draw ("e1 x0");
      if (iPer != 1 && !skipSignalMC)
        photonRtrkHist_mc_sig->Draw ("same e1 x0");
      photonRtrkHist->Draw ("same e1 x0");
      photonRtrkGraph_sys->Draw ("2");

      const int countsData = photonRtrkCounts[iPer][0][0][0]->Integral (p_lo, p_hi, 1, numtrketabins, 1, numphibins);
      const int countsMC = photonRtrkCounts[iPer][1][0][0]->Integral (p_lo, p_hi, 1, numtrketabins, 1, numphibins);
      const int countsMCSig = photonRtrkCounts[iPer][2][0][0]->Integral (p_lo, p_hi, 1, numtrketabins, 1, numphibins);
      myMarkerText (0.225, 0.88, dataColor, markerStyle, Form ("2016 Data (%i events)", countsData), 1.25, 0.04/uPadY);
      myMarkerText (0.225, 0.81, mcOverlayColor, markerStyle, Form ("Pythia8 MC + Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
      if (iPer != 1 && !skipSignalMC)
        myMarkerText (0.225, 0.74, mcSignalColor, markerStyle, Form ("Pythia8 MC (%i events)", countsMCSig), 1.25, 0.04/uPadY);
      if (p_lo != 1 || p_hi != numpbins)
       myText (0.195, 0.58, kBlack, Form ("%g < #it{p}_{T}^{#gamma} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
      myText (0.195, 0.66, kBlack, period, 0.04/uPadY);

      bottomPad->cd ();
      bottomPad->SetLogx (false);

      photonRtrkHist_rat = GetDataOverMC (Form ("photonRtrk_DataMCRatio_%s_iP%i", per, iP), proj2d, proj2d_mc, numtrketabins, trketabins, false, "x");
      photonRtrkHist_rat_lo = GetDataOverMC (Form ("photonRtrk_DataMCRatio_lo_%s_iP%i", per, iP), proj2d_lo, proj2d_mc, numtrketabins, trketabins, false, "x");
      photonRtrkHist_rat_hi = GetDataOverMC (Form ("photonRtrk_DataMCRatio_hi_%s_iP%i", per, iP), proj2d_hi, proj2d_mc, numtrketabins, trketabins, false, "x");

      photonRtrkGraph_rat_sys = new TGraphAsymmErrors (photonRtrkHist_rat);
      photonRtrkGraph_rat_sys->SetName (Form ("photonRtrk_DataMCRatio_sys_%s_iP%i", per, iP));
      CalcSystematics (photonRtrkGraph_rat_sys, photonRtrkHist_rat, photonRtrkHist_rat_hi, photonRtrkHist_rat_lo);
      if (photonRtrkHist_rat_lo) delete photonRtrkHist_rat_lo;
      if (photonRtrkHist_rat_hi) delete photonRtrkHist_rat_hi;

      photonRtrkGraph_rat_sys->SetFillColor (dataColor);
      photonRtrkGraph_rat_sys->SetFillStyle (3001);

      photonRtrkHist_rat->SetMarkerStyle (markerStyle);
      photonRtrkHist_rat->SetLineColor (dataColor);
      photonRtrkHist_rat->SetMarkerColor (dataColor);

      if (iPer != 1 && !skipSignalMC) {
        photonRtrkHist_rat_sig = GetDataOverMC (Form ("photonRtrk_DataMCRatio_Signal_%s_iP%i", per, iP), proj2d, proj2d_mc_sig, numtrketabins, trketabins, false, "x");
        photonRtrkHist_rat_lo = GetDataOverMC (Form ("photonRtrk_DataMCRatio_Signal_lo_%s_iP%i", per, iP), proj2d_lo, proj2d_mc_sig, numtrketabins, trketabins, false, "x");
        photonRtrkHist_rat_hi = GetDataOverMC (Form ("photonRtrk_DataMCRatio_Signal_hi_%s_iP%i", per, iP), proj2d_hi, proj2d_mc_sig, numtrketabins, trketabins, false, "x");

        photonRtrkGraph_rat_sys_sig = new TGraphAsymmErrors (photonRtrkHist_rat_sig);
        photonRtrkGraph_rat_sys_sig->SetName (Form ("photonRtrk_DataMCRatio_Signal_sys_%s_iP%i", per, iP));
        CalcSystematics (photonRtrkGraph_rat_sys_sig, photonRtrkHist_rat_sig, photonRtrkHist_rat_hi, photonRtrkHist_rat_lo);
        if (photonRtrkHist_rat_lo) delete photonRtrkHist_rat_lo;
        if (photonRtrkHist_rat_hi) delete photonRtrkHist_rat_hi;

        photonRtrkGraph_rat_sys_sig->SetFillColor (mcSignalColor);
        photonRtrkGraph_rat_sys_sig->SetFillStyle (3001);

        photonRtrkHist_rat_sig->SetMarkerStyle (markerStyle);
        photonRtrkHist_rat_sig->SetLineColor (mcSignalColor);
        photonRtrkHist_rat_sig->SetMarkerColor (mcSignalColor);
      }

      photonRtrkHist_rat->GetXaxis ()->SetTitle ("#eta_{det}^{#gamma}");
      photonRtrkHist_rat->GetYaxis ()->SetTitle ("Data / MC");
      photonRtrkHist_rat->GetYaxis ()->SetRangeUser (0.94, 1.06);
      photonRtrkHist_rat->GetYaxis ()->SetNdivisions (405);
      photonRtrkHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
      photonRtrkHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
      photonRtrkHist_rat->GetXaxis ()->SetTitleOffset (1);
      photonRtrkHist_rat->GetYaxis ()->SetTitleOffset (1.2*dPadY);
      photonRtrkHist_rat->GetYaxis ()->CenterTitle (true);
      photonRtrkHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
      photonRtrkHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
      photonRtrkHist_rat->GetXaxis ()->SetTickLength (0.08);

      photonRtrkHist_rat->Draw ("e1 x0");
      photonRtrkGraph_rat_sys->Draw ("2");
      if (iPer != 1 && !skipSignalMC) {
        photonRtrkHist_rat_sig->Draw ("same e1 x0");
        photonRtrkGraph_rat_sys_sig->Draw ("2");
      }
      for (TLine* line : etalines) line->Draw ("same");

      TF1* constFit = new TF1 ("constFit", "[0]", etabins[0], etabins[numetabins]);
      photonRtrkHist_rat->Fit (constFit, "R0QN");
      myText (0.195, 0.33, kBlack, Form ("#mu = %.3f #pm %.3f", constFit->GetParameter (0), constFit->GetParError (0)), 0.04/dPadY);
      if (constFit) { delete constFit; constFit = NULL; }

      char* plotName;
      if (iP < numpbins) plotName = Form ("photon_rtrk_iP%i.pdf", iP);
      else plotName = Form ("photon_rtrk_iP_combined.pdf");

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

      if (photonRtrkHist_rat) photonRtrkHist_rat->Write ();
      if (photonRtrkGraph_rat_sys) photonRtrkGraph_rat_sys->Write ();
      if (photonRtrkHist_rat_sig) photonRtrkHist_rat_sig->Write ();
      if (photonRtrkGraph_rat_sys_sig) photonRtrkGraph_rat_sys_sig->Write ();

      if (proj2d) { delete proj2d; proj2d = NULL; }
      if (proj2d_lo) { delete proj2d_lo; proj2d_lo = NULL; }
      if (proj2d_hi) { delete proj2d_hi; proj2d_hi = NULL; }
      if (proj2d_mc) { delete proj2d_mc; proj2d_mc = NULL; }
      if (proj2d_mc_sig) { delete proj2d_mc_sig; proj2d_mc_sig = NULL; }

      if (photonRtrkHist) { delete photonRtrkHist; photonRtrkHist = NULL; }
      if (photonRtrkHist_mc) { delete photonRtrkHist_mc; photonRtrkHist_mc = NULL; }
      if (photonRtrkHist_mc_sig) { delete photonRtrkHist_mc_sig; photonRtrkHist_mc_sig = NULL; }
      if (photonRtrkGraph_sys) { delete photonRtrkGraph_sys; photonRtrkGraph_sys = NULL; }
      if (photonRtrkHist_rat) { delete photonRtrkHist_rat; photonRtrkHist_rat = NULL; }
      if (photonRtrkHist_rat_sig) { delete photonRtrkHist_rat_sig; photonRtrkHist_rat_sig = NULL; }
      if (photonRtrkGraph_rat_sys) { delete photonRtrkGraph_rat_sys; photonRtrkGraph_rat_sys = NULL; }
      if (photonRtrkGraph_rat_sys_sig) { delete photonRtrkGraph_rat_sys_sig; photonRtrkGraph_rat_sys_sig = NULL; }
    } // end loop over pT bins


    /**** Plots rtrk distributions, binned by eta^photon ****/
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

        proj2d = Project2D ("", photonRtrkHists[iPer][0][1][2], "y", "z", p_lo, p_hi);
        proj2d->RebinY (rebinFactor);
        photonRtrkHist = proj2d->ProjectionY ("photonRtrkHist", eta_lo, eta_hi);

        if (photonRtrkHist->Integral () != 0) photonRtrkHist->Scale (1./photonRtrkHist->Integral ());

        proj2d_mc = Project2D ("", photonRtrkHists[iPer][1][1][2], "y", "z", p_lo, p_hi);
        proj2d_mc->RebinY (rebinFactor);
        photonRtrkHist_mc = proj2d_mc->ProjectionY ("vJetProjection_mc", eta_lo, eta_hi);

        if (photonRtrkHist_mc->Integral () != 0) photonRtrkHist_mc->Scale (1./photonRtrkHist_mc->Integral ()); 

        photonRtrkHist_mc->GetXaxis ()->SetTitle ("r_{trk} = #Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{#gamma}");
        photonRtrkHist_mc->GetYaxis ()->SetTitle ("Counts / Total");
        photonRtrkHist_mc->GetYaxis ()->SetRangeUser (0., 1.6*std::max (photonRtrkHist->GetMaximum (), photonRtrkHist_mc->GetMaximum ()));
        photonRtrkHist_mc->SetMarkerStyle (markerStyle);
        photonRtrkHist_mc->SetMarkerColor (mcOverlayColor);
        photonRtrkHist_mc->SetLineColor (mcOverlayColor);

        photonRtrkHist_mc->GetXaxis ()->SetLabelSize (0.04/uPadY);
        photonRtrkHist_mc->GetYaxis ()->SetLabelSize (0.04/uPadY);
        photonRtrkHist_mc->GetYaxis ()->SetTitleSize (0.04/uPadY);
        photonRtrkHist_mc->GetYaxis ()->SetTitleOffset (1.1*uPadY);

        float mean, mean_err, mean_mc, mean_mc_err;

        TF1 *gaus_data = 0, *gaus_mc = 0;
  
        mean = photonRtrkHist->GetMean ();
        float stddev = photonRtrkHist->GetStdDev ();
        gaus_data = new TF1 ("gaus_data", "gaus (0)", 0.1, 2);
        photonRtrkHist->Fit (gaus_data, "Q0R");

        mean = photonRtrkHist_mc->GetMean ();
        stddev = photonRtrkHist_mc->GetStdDev ();
        gaus_mc = new TF1 ("gaus_mc", "gaus (0)", 0.1, 2);
        photonRtrkHist_mc->Fit (gaus_mc, "Q0R");

        mean = gaus_data->GetParameter (1);
        mean_err = gaus_data->GetParError (1);
        mean_mc = gaus_mc->GetParameter (1);
        mean_mc_err = gaus_mc->GetParError (1);

        photonRtrkHist->Scale (1./gaus_data->Integral (gaus_data->GetXmin (), gaus_data->GetXmax ()));
        photonRtrkHist_mc->Scale (1./gaus_mc->Integral (gaus_mc->GetXmin (), gaus_mc->GetXmax ()));
        gaus_data->SetParameter (0, gaus_data->GetParameter (0) / gaus_data->Integral (gaus_data->GetXmin (), gaus_data->GetXmax ()));
        gaus_mc->SetParameter (0, gaus_mc->GetParameter (0) / gaus_mc->Integral (gaus_mc->GetXmin (), gaus_mc->GetXmax ()));

        photonRtrkHist_mc->GetYaxis ()->SetRangeUser (0., 1.6* std::max (photonRtrkHist->GetMaximum (), photonRtrkHist_mc->GetMaximum ()));

        photonRtrkHist_mc->DrawCopy ("e1 x0");
        photonRtrkHist->DrawCopy ("same e1 x0");

        gaus_data->SetLineColor (dataColor);
        gaus_mc->SetLineColor (mcOverlayColor);
        gaus_data->Draw ("same");
        gaus_mc->Draw ("same");

        const int countsData = photonRtrkCounts[iPer][0][0][0]->Integral (p_lo, p_hi, eta_lo, eta_hi, 1, numphibins);
        const int countsMC = photonRtrkCounts[iPer][1][0][0]->Integral (p_lo, p_hi, eta_lo, eta_hi, 1, numphibins);

        myMarkerText (0.225, 0.88, dataColor, markerStyle, Form ("2016 Data (%i events)", countsData), 1.25, 0.04/uPadY);
        myMarkerText (0.225, 0.80, mcOverlayColor, markerStyle, Form ("Pythia8 MC (%i events)", countsMC), 1.25, 0.04/uPadY);

        myText (0.65, 0.88, dataColor, Form ("#mu_{data} = %s", FormatMeasurement (mean, mean_err)), 0.04/uPadY);
        myText (0.65, 0.80, dataColor, Form ("#mu_{MC} = %s", FormatMeasurement (mean_mc, mean_mc_err)), 0.04/uPadY);

        myText (0.68, 0.34, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
        myText (0.68, 0.25, dataColor, period, 0.04/uPadY);
        myText (0.68, 0.16, dataColor, Form ("%g < #it{p}_{T}^{#gamma} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
        myText (0.68, 0.08, dataColor, Form ("%g < #eta_{det}^{#gamma} < %g", trketabins[eta_lo-1], trketabins[eta_hi]), 0.04/uPadY);

        bottomPad->cd ();
        bottomPad->SetLogx (false);
        photonRtrkHist->Divide (photonRtrkHist_mc);

        photonRtrkHist->SetYTitle ("Data / MC");
        photonRtrkHist->SetAxisRange (0.45, 1.65, "Y");
        photonRtrkHist->GetYaxis ()->SetNdivisions (605);
        photonRtrkHist->GetXaxis ()->SetTitleSize (0.04/dPadY);
        photonRtrkHist->GetYaxis ()->SetTitleSize (0.04/dPadY);
        photonRtrkHist->GetXaxis ()->SetTitleOffset (1);
        photonRtrkHist->GetYaxis ()->SetTitleOffset (1.1*dPadY);
        photonRtrkHist->GetYaxis ()->CenterTitle (true);
        photonRtrkHist->GetXaxis ()->SetLabelSize (0.04/dPadY);
        photonRtrkHist->GetYaxis ()->SetLabelSize (0.04/dPadY);
        photonRtrkHist->GetXaxis ()->SetTickLength (0.08);

        photonRtrkHist->DrawCopy ("e1 x0"); 
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

        if (photonRtrkHist) { delete photonRtrkHist; photonRtrkHist = NULL; }
        if (photonRtrkHist_mc) { delete photonRtrkHist_mc; photonRtrkHist_mc = NULL; }

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
