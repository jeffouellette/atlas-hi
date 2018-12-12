#include "PhotonJetsHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TTree.h>
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


void PhotonJetsHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("PhotonJets/", "JetCalibration/");

  //////////////////////////////////////////////////////////////////////////////
  // Create desired histograms
  //////////////////////////////////////////////////////////////////////////////
  TH3D***** photonJetHists_pt_eta = Get4DArray <TH3D*> (3, 2, 3, 3); // iPer, iData, iErr, iPCut
  TH3D***** photonJetCounts = Get4DArray <TH3D*> (3, 2, 3, 2); // iPer, iData, iPCut, iWgt

  for (short iPer = 0; iPer < 3; iPer++) {
    const char* per = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));

    for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
      const char* data = (iData == 0 ? "data" : "mc");

      for (short iPCut = 0; iPCut < 3; iPCut++) {
        const char* pCut = (iPCut == 0 ? "tight" : (iPCut == 1 ? "lnt" : "signal")); // lnt = loose non-tight

        for (short iErr = 0; iErr < 3; iErr++) {
          if (iErr != 1 && iData != 0)
           continue; // ignore systematics for MC
          const char* error = (iErr == 0 ? "syslo" : (iErr == 1 ? "stat" : "syshi"));

          const char* name = Form ("photonJetPtRatio_%s_%s_%s_%s", per, data, error, pCut);
          photonJetHists_pt_eta[iPer][iData][iErr][iPCut] = new TH3D (name, "", numpbins, pbins, numetabins, etabins, 1000, linspace (0, 2.0, 1000));
          photonJetHists_pt_eta[iPer][iData][iErr][iPCut]->Sumw2 ();
        }

        for (short iWgt = 0; iWgt < 2; iWgt++) {
          const char* weight = (iWgt == 0 ? "unweighted" : "weighted");

          const char* name = Form ("photonJetCounts_%s_%s_%s_%s", per, data, pCut, weight);
          photonJetCounts[iPer][iData][iPCut][iWgt] = new TH3D (name, "", numpbins, pbins, numetabins, etabins, numphibins, phibins);
          photonJetCounts[iPer][iData][iPCut][iWgt]->Sumw2 ();
        }
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // Load analyzed TTrees
  //////////////////////////////////////////////////////////////////////////////
  float purityFactor = 0, ppt = 0, peta = 0, pphi = 0, jpt = 0, jeta = 0, jphi = 0, je = 0, jpterr = 0, dPhi = 0;
  double evtWeight = 0;
  bool ptight = false, isMC = false, isPeriodA = false;

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");
  TTree* inTree = (TTree*)inFile->Get ("PhotonJetTree");

  inTree->SetBranchAddress ("evt_weight", &evtWeight);
  inTree->SetBranchAddress ("purity_factor", &purityFactor);
  inTree->SetBranchAddress ("photon_pt", &ppt);
  //inTree->SetBranchAddress ("photon_eta", &peta);
  //inTree->SetBranchAddress ("photon_phi", &pphi);
  inTree->SetBranchAddress ("jet_pt", &jpt);
  inTree->SetBranchAddress ("jet_eta", &jeta);
  inTree->SetBranchAddress ("jet_phi", &jphi);
  //inTree->SetBranchAddress ("jet_e", &je);
  inTree->SetBranchAddress ("delta_phi", &dPhi);
  inTree->SetBranchAddress ("jet_pt_sys", &jpterr);
  inTree->SetBranchAddress ("tight_photon", &ptight);
  inTree->SetBranchAddress ("isMC", &isMC);
  inTree->SetBranchAddress ("isPeriodA", &isPeriodA);


  //////////////////////////////////////////////////////////////////////////////
  // Fill desired histograms
  //////////////////////////////////////////////////////////////////////////////
  const int nPhotonJets = inTree->GetEntries ();
  for (int iPhotonJet = 0; iPhotonJet < nPhotonJets; iPhotonJet++) {
    inTree->GetEntry (iPhotonJet);

    if (ppt > 0) {
      double xjref = jpt / (ppt*cos(pi-dPhi));
      double xjreferr = jpterr / (ppt*cos(pi-dPhi));

      photonJetHists_pt_eta[isPeriodA?0:1][isMC?1:0][1][ptight?0:1]->Fill (ppt, jeta, xjref, evtWeight*purityFactor);
      photonJetHists_pt_eta[2][isMC?1:0][1][ptight?0:1]->Fill (ppt, jeta, xjref, evtWeight*purityFactor);

      if (!isMC) {
        photonJetHists_pt_eta[isPeriodA?0:1][0][0][ptight?0:1]->Fill (ppt, jeta, xjref-xjreferr, evtWeight*purityFactor);
        photonJetHists_pt_eta[isPeriodA?0:1][0][2][ptight?0:1]->Fill (ppt, jeta, xjref+xjreferr, evtWeight*purityFactor);
        photonJetHists_pt_eta[2][0][0][ptight?0:1]->Fill (ppt, jeta, xjref-xjreferr, evtWeight*purityFactor);
        photonJetHists_pt_eta[2][0][2][ptight?0:1]->Fill (ppt, jeta, xjref+xjreferr, evtWeight*purityFactor);
      }
    }

    photonJetCounts[isPeriodA?0:1][isMC?1:0][ptight?0:1][0]->Fill (ppt, jeta, jphi);
    photonJetCounts[2][isMC?1:0][ptight?0:1][0]->Fill (ppt, jeta, jphi);
    photonJetCounts[isPeriodA?0:1][isMC?1:0][ptight?0:1][1]->Fill (ppt, jeta, jphi, evtWeight);
    photonJetCounts[2][isMC?1:0][ptight?0:1][1]->Fill (ppt, jeta, jphi, evtWeight);
  }
  inTree = NULL;
  inFile->Close ();
  if (inFile) { delete inFile; inFile = NULL; }


  //////////////////////////////////////////////////////////////////////////////
  // Subtract hadronic background from tight distributions
  //////////////////////////////////////////////////////////////////////////////
  for (short iPer = 0; iPer < 3; iPer++) {
    for (short iErr = 0; iErr < 3; iErr++) { // subtract for data
      if (subtractBackground) {
        SubtractLooseNonTight (photonJetHists_pt_eta[iPer][0][iErr][2], // signal
                               photonJetHists_pt_eta[iPer][0][iErr][0], // tight
                               photonJetHists_pt_eta[iPer][0][iErr][1], // loose non tight
                               photonJetCounts[iPer][0][0][1], // tight weighted counts
                               photonJetCounts[iPer][0][1][1]); // loose non tight weighted counts
      }
      else {
        delete photonJetHists_pt_eta[iPer][0][iErr][2];
        photonJetHists_pt_eta[iPer][0][iErr][2] = photonJetHists_pt_eta[iPer][0][iErr][0];
      }
    }
    // subtract for MC
    if (subtractBackground) {
      SubtractLooseNonTight (photonJetHists_pt_eta[iPer][1][1][2],
                             photonJetHists_pt_eta[iPer][1][1][0],
                             photonJetHists_pt_eta[iPer][1][1][1],
                             photonJetCounts[iPer][1][0][1],
                             photonJetCounts[iPer][1][1][1]);
    }
    else {
      delete photonJetHists_pt_eta[iPer][1][1][2];
      photonJetHists_pt_eta[iPer][1][1][2] = photonJetHists_pt_eta[iPer][1][1][0];
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // Save histograms for interactive access
  //////////////////////////////////////////////////////////////////////////////
  TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");
  for (short iPer = 0; iPer < 3; iPer++) {
    for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
      for (short iPCut = 0; iPCut < 3; iPCut++) {
        for (short iErr = 0; iErr < 3; iErr++) {
          if (iErr != 1 && iData != 0)
            continue;
          photonJetHists_pt_eta[iPer][iData][iErr][iPCut]->Write ();
        }
      
      
        for (short iWgt = 0; iWgt < 2; iWgt++) {
          photonJetCounts[iPer][iData][iPCut][iWgt]->Write ();
        }
      }
    }
  }
  outFile->Close ();
  if (outFile) { delete outFile; outFile = NULL; }


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

  TCanvas* canvas = new TCanvas ("canvas", "", 800, 1000);
  const double padRatio = 1.2; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = 1.0/ (padRatio+1.0);
  const double uPadY = 1.0 - dPadY;
  TPad* topPad = new TPad ("topPad", "", 0, dPadY, 1, 1);
  TPad* bottomPad = new TPad ("bottomPad", "", 0, 0, 1, dPadY);
  topPad->SetBottomMargin (0);
  topPad->SetLeftMargin (-0.12);
  bottomPad->SetTopMargin (0);
  bottomPad->SetBottomMargin (0.30);
  bottomPad->SetLeftMargin (-0.12);
  topPad->Draw ();
  bottomPad->Draw ();


  /**** Define local histograms, graphs, etc. ****/
  TH1D *vJetHist = NULL, *vJetHist_mc = NULL, *vJetHist_lo = NULL, *vJetHist_hi = NULL, *vJetHist_rat = NULL, *vJetHist_rat_lo = NULL, *vJetHist_rat_hi = NULL;
  TH2D *proj = NULL, *proj_mc = NULL, *proj_lo = NULL, *proj_hi = NULL;
  TGraphAsymmErrors *vJetGraph = NULL, *vJetGraph_mc = NULL, *vJetGraph_sys = NULL, *vJetGraph_rat = NULL, *vJetGraph_rat_sys = NULL;
  char* plotName;

  //TFile* outFile = new TFile (Form ("%s/histograms.root", rootPath.Data ()), "recreate");

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
    const char* period = (iPer == 0 ? "Period A" : (iPer == 1 ? "Period B" : "Period A+B"));
    const char* per = (iPer == 0 ? "pA" : (iPer == 1 ? "pB" : "pAB"));

    for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta
      const int eta_lo = (iEta != numetabins ? iEta+1 : eta_lo_comb);
      const int eta_hi = (iEta != numetabins ? iEta+1 : eta_hi_comb);

      const Color_t dataColor = kBlack;
      const Color_t mcColor = kRed;
      const Style_t dataStyle = 20;
      const Style_t mcStyle = 20;

      topPad->cd ();
      topPad->SetLogx ();

      proj = Project2D ("", photonJetHists_pt_eta[iPer][0][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
      proj->RebinY (rebinFactor);
      TString name = Form ("photonJetHist_%s_data_stat_signal_iEta%i", per, iEta);
      vJetHist = GetProfileX (name, proj, numpbins, pbins, true);

      vJetGraph = make_graph (vJetHist);
      vJetGraph->SetMarkerStyle (dataStyle);
      vJetGraph->SetMarkerColor (dataColor);
      vJetGraph->SetLineColor (dataColor);
      vJetGraph->SetLineWidth (2);

      // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
      proj_lo = Project2D ("", photonJetHists_pt_eta[iPer][0][0][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
      proj_lo->RebinY (rebinFactor);
      name = Form ("photonJetHist_%s_data_sys_lo_signal_iEta%i", per, iEta);
      vJetHist_lo = GetProfileX (name, proj_lo, numpbins, pbins, true);

      proj_hi = Project2D ("", photonJetHists_pt_eta[iPer][0][2][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
      proj_hi->RebinY (rebinFactor);
      name = Form ("photonJetHist_%s_data_sys_hi_signal_iEta%i", per, iEta);
      vJetHist_hi = GetProfileX (name, proj_hi, numpbins, pbins, true);

      vJetGraph_sys = make_graph (vJetHist); // for plotting systematics
      CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
      if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
      if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

      vJetGraph_sys->SetFillColor (dataColor);
      vJetGraph_sys->SetFillStyle (3001);

      proj_mc = Project2D ("", photonJetHists_pt_eta[iPer][1][1][2], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
      proj_mc->RebinY (rebinFactor);
      name = Form ("photonJetHist_%s_mc_stat_signal_iEta%i", per, iEta);
      vJetHist_mc = GetProfileX (name, proj_mc, numpbins, pbins, true);

      double middle = 0.05 * floor (20. * vJetHist_mc->Integral () / vJetHist_mc->GetNbinsX ()); // gets mean along y
      if (10 * middle != floor (10*middle)) middle += 0.05;

      vJetGraph_mc = make_graph (vJetHist_mc);
      vJetGraph_mc->GetYaxis ()->SetTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
      vJetGraph_mc->GetYaxis ()->SetRangeUser (middle - 0.35, middle + 0.35);
      vJetGraph_mc->SetMarkerStyle (mcStyle);
      vJetGraph_mc->SetMarkerColor (mcColor);
      vJetGraph_mc->SetLineColor (mcColor);
      vJetGraph_mc->SetLineWidth (2);
      vJetGraph_mc->GetXaxis ()->SetLabelSize (0.032/uPadY);
      vJetGraph_mc->GetYaxis ()->SetLabelSize (0.032/uPadY);
      vJetGraph_mc->GetYaxis ()->SetTitleSize (0.032/uPadY);
      vJetGraph_mc->GetYaxis ()->SetTitleOffset (1.2*uPadY);

      vJetGraph_mc->Draw ("ap");
      vJetGraph->Draw ("p");
      vJetGraph_sys->Draw ("2");

      int countsData = 0, countsMC = 0;
      if (exclusive && iEta == numetabins) {
        countsData = photonJetCounts[iPer][0][0][0]->Integral () - photonJetCounts[iPer][0][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
        countsMC = photonJetCounts[iPer][1][0][0]->Integral () - photonJetCounts[iPer][1][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
      }
      else {
        countsData = photonJetCounts[iPer][0][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
        countsMC = photonJetCounts[iPer][1][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
      }

      topPad->cd ();
      myMarkerText (0.175, 0.88, dataColor, dataStyle, Form ("2016 Data (%i events)", countsData), 1.25, 0.032/uPadY);
      myMarkerText (0.175, 0.82, mcColor, mcStyle, Form ("MC (%i events)", countsMC), 1.25, 0.032/uPadY);
      myText (0.155, 0.20, dataColor, "#bf{#it{ATLAS}} Internal", 0.032/uPadY);
      if (eta_lo != 1 || eta_hi != numetabins)
        if (!exclusive || iEta != numetabins)
          myText (0.155, 0.135, dataColor, Form ("#gamma + Jet, %g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.032/uPadY);
        else
          myText (0.155, 0.135, dataColor, Form ("#gamma + Jet, %g < #left|#eta_{det}^{Jet}#right|", etabins[eta_hi]), 0.032/uPadY);
      else
        myText (0.155, 0.14, dataColor, "#gamma + Jet", 0.032/uPadY);
      myText (0.155, 0.08, dataColor, period, 0.032/uPadY);


      bottomPad->cd ();
      bottomPad->SetLogx ();

      vJetHist_rat = GetDataOverMC (Form ("photonJetPtDataMCRatio_mc_iEta%i", iEta), proj, proj_mc, numpbins, pbins, true, "x");
      vJetHist_rat_lo = GetDataOverMC (Form ("photonJetPtDataMCRatio_lo_mc_iEta%i", iEta), proj_lo, proj_mc, numpbins, pbins, true, "x");
      vJetHist_rat_hi = GetDataOverMC (Form ("photonJetPtDataMCRatio_hi_mc_iEta%i", iEta), proj_hi, proj_mc, numpbins, pbins, true, "x");

      vJetGraph_rat_sys = make_graph (vJetHist_rat);
      CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
      if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
      if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
      vJetGraph_rat_sys->SetFillColor (dataColor);
      vJetGraph_rat_sys->SetFillStyle (3001);

      vJetGraph_rat = make_graph (vJetHist_rat);
      vJetGraph_rat->GetXaxis ()->SetTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
      vJetGraph_rat->GetYaxis ()->SetTitle ("Data / MC");
      vJetGraph_rat->GetYaxis ()->SetRangeUser (0.91, 1.09);
      vJetGraph_rat->SetMarkerStyle (dataStyle);
      vJetGraph_rat->SetMarkerColor (dataColor);
      vJetGraph_rat->SetLineColor (dataColor);
      vJetGraph_rat->SetLineWidth (2);
      vJetGraph_rat->GetYaxis ()->SetNdivisions (405);
      vJetGraph_rat->GetXaxis ()->SetTitleSize (0.032/dPadY);
      vJetGraph_rat->GetYaxis ()->SetTitleSize (0.032/dPadY);
      vJetGraph_rat->GetXaxis ()->SetTitleOffset (1);
      vJetGraph_rat->GetYaxis ()->SetTitleOffset (1.2*dPadY);
      vJetGraph_rat->GetYaxis ()->CenterTitle (true);
      vJetGraph_rat->GetXaxis ()->SetLabelSize (0.032/dPadY);
      vJetGraph_rat->GetYaxis ()->SetLabelSize (0.032/dPadY);
      vJetGraph_rat->GetXaxis ()->SetTickLength (0.08);

      vJetGraph_rat->Draw ("ap");
      vJetGraph_rat_sys->Draw ("2");
      for (TLine* line : glines) line->Draw ();

      if (iEta < numetabins) plotName = Form ("gamma_jet_iEta%i.pdf", iEta);
      else plotName = Form ("gamma_jet_iEta_combined.pdf");

      canvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));

      if (vJetHist) { delete vJetHist; vJetHist = NULL; }
      if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
      if (vJetGraph) { delete vJetGraph; vJetGraph = NULL; }
      if (vJetGraph_mc) { delete vJetGraph_mc; vJetGraph_mc = NULL; }
      if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }

      if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
      if (vJetGraph_rat) { delete vJetGraph_rat; vJetGraph_rat = NULL; }
      if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }

      if (proj) { delete proj; proj = NULL; }
      if (proj_mc) { delete proj_mc; proj_mc = NULL; }
      if (proj_lo) { delete proj_lo; proj_lo = NULL; }
      if (proj_hi) { delete proj_hi; proj_hi = NULL; }

    } // end loop over etabins


    /**** Now loop over pT bins and plot response as function of eta^jet ****/
    for (short iP = 0; iP <= numpbins; iP++) {

      const int p_lo = (iP != numpbins ? iP+1 : p_lo_comb);
      const int p_hi = (iP != numpbins ? iP+1 : p_hi_comb);

      /**** Plots PhotonJet info as a function of eta^jet****/
      const Color_t dataColor = kBlack;
      const Style_t dataStyle = 20;
      const Style_t mcStyle = 20;

      topPad->cd ();
      topPad->SetLogx (false);

      proj = Project2D ("", photonJetHists_pt_eta[iPer][0][1][2], "y", "z", p_lo, p_hi);
      proj->RebinY (rebinFactor);
      vJetHist = GetProfileX ("vJetHist", proj, numetabins, etabins, true);

      vJetGraph = make_graph (vJetHist);
      vJetGraph->SetMarkerStyle (dataStyle);
      vJetGraph->SetMarkerColor (dataColor);
      vJetGraph->SetLineColor (dataColor);
      vJetGraph->SetLineWidth (2);

      // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
      proj_lo = Project2D ("", photonJetHists_pt_eta[iPer][0][0][2], "y", "z", p_lo, p_hi);
      proj_lo->RebinY (rebinFactor);
      vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numetabins, etabins, true);

      proj_hi = Project2D ("", photonJetHists_pt_eta[iPer][0][2][2], "y", "z", p_lo, p_hi);
      proj_hi->RebinY (rebinFactor);
      vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, true);

      vJetGraph_sys = make_graph (vJetHist); // for plotting systematics
      CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
      if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
      if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

      vJetGraph_sys->SetFillColor (dataColor);
      vJetGraph_sys->SetFillStyle (3001);

      proj_mc = Project2D ("", photonJetHists_pt_eta[iPer][1][1][2], "y", "z", p_lo, p_hi);
      proj_mc->RebinY (rebinFactor);
      vJetHist_mc = GetProfileX ("vJetHist_mc", proj_mc, numetabins, etabins, true);

      double middle = 0.05 * floor (20. * vJetHist_mc->Integral () / vJetHist_mc->GetNbinsX ()); // gets mean along y
      if (10 * middle != floor (10*middle)) middle += 0.05;

      vJetGraph_mc = make_graph (vJetHist_mc);
      vJetGraph_mc->GetYaxis ()->SetTitle ("<#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
      vJetGraph_mc->GetYaxis ()->SetRangeUser (middle - 0.35, middle + 0.35);
      vJetGraph_mc->SetMarkerStyle (mcStyle);
      vJetGraph_mc->SetMarkerColor (mcColor);
      vJetGraph_mc->SetLineColor (mcColor);
      vJetGraph_mc->SetLineWidth (2);
      vJetGraph_mc->GetXaxis ()->SetLabelSize (0.032/uPadY);
      vJetGraph_mc->GetYaxis ()->SetLabelSize (0.032/uPadY);
      vJetGraph_mc->GetYaxis ()->SetTitleSize (0.032/uPadY);
      vJetGraph_mc->GetYaxis ()->SetTitleOffset (uPadY);

      vJetGraph_mc->Draw ("ap");
      vJetGraph->Draw ("p");
      vJetGraph_sys->Draw ("2");

      int countsData = 0, countsMC = 0;
      countsData = photonJetCounts[iPer][0][0][0]->Integral (p_lo, p_hi, 1, numetabins, 1, numphibins);
      countsMC = photonJetCounts[iPer][1][0][0]->Integral (p_lo, p_hi, 1, numetabins, 1, numphibins);

      myMarkerText (0.175, 0.88, dataColor, dataStyle, Form ("2016 Data (%i events)", countsData), 1.25, 0.032/uPadY);
      myMarkerText (0.175, 0.82, mcColor, mcStyle, Form ("Pythia8 MC (%i events)", countsMC), 1.25, 0.032/uPadY);
      myText (0.155, 0.2, dataColor, "#bf{#it{ATLAS}} Internal", 0.032/uPadY);
      if (p_lo != 1 || p_hi != numpbins)
        myText (0.155, 0.135, dataColor, Form ("#gamma + Jet, %g < #it{p}_{T}^{#gamma} < %g", pbins[p_lo-1], pbins[p_hi]), 0.032/uPadY);
      else
        myText (0.155, 0.14, dataColor, "#gamma + Jet", 0.032/uPadY);
      myText (0.155, 0.08, dataColor, period, 0.032/uPadY);
      

      bottomPad->cd ();
      bottomPad->SetLogx (false);

      vJetHist_rat = GetDataOverMC (Form ("photonJetPtDataMCRatio_iP%i", iP), proj, proj_mc, numetabins, etabins, true, "x");
      vJetHist_rat_lo = GetDataOverMC (Form ("photonJetPtDataMCRatio_lo_iP%i", iP), proj_lo, proj_mc, numetabins, etabins, true, "x");
      vJetHist_rat_hi = GetDataOverMC (Form ("photonJetPtDataMCRatio_hi_iP%i", iP), proj_hi, proj_mc, numetabins, etabins, true, "x");

      vJetGraph_rat_sys = make_graph (vJetHist_rat);
      CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
      if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
      if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
      vJetGraph_rat_sys->SetFillColor (dataColor);
      vJetGraph_rat_sys->SetFillStyle (3001);

      vJetGraph_rat = make_graph (vJetHist_rat);
      vJetGraph_rat->GetXaxis ()->SetTitle ("Jet #eta");
      vJetGraph_rat->GetYaxis ()->SetTitle ("Data / MC");
      vJetGraph_rat->GetYaxis ()->SetRangeUser (0.91, 1.09);
      vJetGraph_rat->SetMarkerStyle (dataStyle);
      vJetGraph_rat->SetMarkerColor (dataColor);
      vJetGraph_rat->SetLineColor (dataColor);
      vJetGraph_rat->SetLineWidth (2);
      vJetGraph_rat->GetYaxis ()->SetNdivisions (405);
      vJetGraph_rat->GetXaxis ()->SetTitleSize (0.032/dPadY);
      vJetGraph_rat->GetYaxis ()->SetTitleSize (0.032/dPadY);
      vJetGraph_rat->GetXaxis ()->SetTitleOffset (1);
      vJetGraph_rat->GetYaxis ()->SetTitleOffset (dPadY);
      vJetGraph_rat->GetYaxis ()->CenterTitle (true);
      vJetGraph_rat->GetXaxis ()->SetLabelSize (0.032/dPadY);
      vJetGraph_rat->GetYaxis ()->SetLabelSize (0.032/dPadY);
      vJetGraph_rat->GetXaxis ()->SetTickLength (0.08);

      vJetGraph_rat->Draw ("ap");
      vJetGraph_rat_sys->Draw ("2");
      for (TLine* line : getalines) line->Draw ();

      if (iP < numpbins) plotName = Form ("gamma_jet_iP%i.pdf", iP);
      else plotName = Form ("gamma_jet_iP_combined.pdf");

      canvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));

      if (vJetHist) { delete vJetHist; vJetHist = NULL; }
      if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
      if (vJetGraph) { delete vJetGraph; vJetGraph = NULL; }
      if (vJetGraph_mc) { delete vJetGraph_mc; vJetGraph_mc = NULL; }
      if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }

      if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
      if (vJetGraph_rat) { delete vJetGraph_rat; vJetGraph_rat = NULL; }
      if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }

      if (proj) { delete proj; proj = NULL; }
      if (proj_mc) { delete proj_mc; proj_mc = NULL; }
      if (proj_lo) { delete proj_lo; proj_lo = NULL; }
      if (proj_hi) { delete proj_hi; proj_hi = NULL; }

    } // end loop over pT bins


    if (!plot_xjref) continue;

    for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta
      const int eta_lo = (iEta != numetabins ? iEta+1 : eta_lo_comb);
      const int eta_hi = (iEta != numetabins ? iEta+1 : eta_hi_comb);

      /**** Plots xjref distributions, binned by ptref ****/
      for (short iP = 0; iP <= numpbins; iP++) {
        const int p_lo = (iP != numpbins ? iP+1 : p_lo_comb);
        const int p_hi =  (iP != numpbins ? iP+1 : p_hi_comb);

        const Color_t dataColor = kBlack;
        const Style_t dataStyle = 20;
        const Style_t mcStyle = 20;

        rawDistCanvas->cd ();
        gPad->SetLogx (false);

        // get data hist
        proj = Project2D ("", photonJetHists_pt_eta[iPer][0][1][0], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
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
        proj_mc = Project2D ("", photonJetHists_pt_eta[iPer][1][1][0], "x", "z", eta_lo, eta_hi, exclusive && iEta == numetabins);
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

        vJetHist_mc->DrawCopy ("e1 x0");
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
          countsData = photonJetCounts[iPer][0][0][0]->Integral () - photonJetCounts[iPer][0][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
          countsMC = photonJetCounts[iPer][1][0][0]->Integral () - photonJetCounts[iPer][1][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
        }
        else {
          countsData = photonJetCounts[iPer][0][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
          countsMC = photonJetCounts[iPer][1][0][0]->Integral (1, numpbins, eta_lo, eta_hi, 1, numphibins);
        }

        myMarkerText (0.175, 0.88, dataColor, dataStyle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countsData), 1.25, 0.04);
        myMarkerText (0.175, 0.82, mcColor, mcStyle, Form ("Pythia8 #it{pp} 8.16 TeV %s (%i events)", (runValidation ? "":"with #it{p}-Pb Overlay"), countsMC), 1.25, 0.04);

        myText (0.155, 0.72, dataColor, Form ("#mu_{Data} = %s", FormatMeasurement (mean, mean_err, 1)), 0.04);
        myText (0.155, 0.66, dataColor, Form ("#mu_{MC} = %s", FormatMeasurement (mean_mc, mean_mc_err, 1)), 0.04);

        myText (0.155, 0.43, dataColor, Form ("#gamma + Jet, %s", period), 0.04);
        myText (0.155, 0.37, dataColor, Form ("%g < #it{p}_{T}^{ref} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04);
        myText (0.155, 0.31, dataColor, Form ("%g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04);
        myText (0.155, 0.25, dataColor, "#bf{#it{ATLAS}} Internal", 0.04);

        plotName = Form ("xjref_dists/gamma_jet_iEta%i_iP%i.pdf", iEta, iP);
        rawDistCanvas->SaveAs (Form ("%s/Period%s/%s", plotPath.Data (), iPer == 0 ? "A" : (iPer == 1 ? "B" : "AB"), plotName));
        
        if (vJetHist) { delete vJetHist; vJetHist = NULL; }
        if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
        if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
        if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

        if (proj) { delete proj; proj = NULL; }
        if (proj_mc) { delete proj_mc; proj_mc = NULL; }
        if (proj_lo) { delete proj_lo; proj_lo = NULL; }
        if (proj_hi) { delete proj_hi; proj_hi = NULL; }
        
      } // end loop over pT bins

    } // end loop over eta bins

  } // end loop over periods

  return;
}

} // end namespace
