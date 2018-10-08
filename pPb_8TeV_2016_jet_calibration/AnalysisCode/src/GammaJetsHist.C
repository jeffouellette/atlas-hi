#include "GammaJetsHist.h"
#include "Params.h"
#include "Utils.h"

#include <ArrayTemplates.h>

#include <TF1.h>
#include <TH1D.h>
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

namespace pPb8TeV2016JetCalibration {

void GammaJetsHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("GammaJets/", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);
  vector<TString> gammaJetOverlaySampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetOverlaySampleIds.push_back (TString ("Pbp_Overlay_GammaJet_Slice") + to_string (i+1));
   gammaJetOverlaySampleIds.push_back (TString ("pPb_Overlay_GammaJet_Slice") + to_string (i+1));
  }
  vector<TString> gammaJetSignalSampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetSignalSampleIds.push_back (TString ("pPb_Signal_GammaJet_Slice") + to_string (i+1));
  }

  TH3D***** gJetHists = Get4DArray <TH3D*> (2, 3, 3, 3); // iYear, iPer, iData, iErr
  //TH2D* gJetHistsSys[3][numetabins+1][3][3];
  TH2D**** gJetCounts = Get3DArray <TH2D*> (2, 3, 3);

  for (short iYear = 0; iYear < 2; iYear++) {
   const TString year = (iYear == 0 ? "2015" : "2016");

   for (short iData = 0; iData < 3; iData++) { // iData is 0 for data, 1 for MC
    TString dataType = "data";
    if (iData == 1) dataType = "mc_overlay";
    else if (iData == 2) dataType = "mc_signalonly";

    for (short iPer = 0; iPer < 3; iPer++) {
     TString period = "periodA";
     if (iPer == 1) period = "periodB";
     else if (iPer == 2) period = "periodAB";

     for (short iErr = 0; iErr < 3; iErr++) {
      TString error = "sys_lo";
      if (iErr == 1) error = "stat";
      else if (iErr == 2) error = "sys_hi";

      gJetHists[iYear][iPer][iData][iErr] = new TH3D (Form ("gJetPtRatio_%s_%s_%s_%s", year.Data (), dataType.Data (), error.Data (), period.Data ()), "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
      gJetHists[iYear][iPer][iData][iErr]->Sumw2 ();

      //gJetHistsSys[iPer][iEta][iData][iErr] = new TH2D (Form ("gJetPtRatioSys_iEta%i_%s_%s_%s", iEta, dataType.Data (), error.Data (), period.Data ()), ";#it{p}_{T}^{jet} #left[GeV#right];#Delta#it{x}_{J}^{ref}#it{p}_{T}^{ref}/#it{p}_{T}^{jet}", numpbins, pbins, numSigmaBins, -maxSigma, maxSigma);
      //gJetHistsSys[iPer][iEta][iData][iErr]->Sumw2 ();
     }

     gJetCounts[iYear][iPer][iData] = new TH2D (Form ("gJetCounts_%s_%s_%s", year.Data (), dataType.Data (), period.Data()), "", numpbins, pbins, numetabins, etabins);
     gJetCounts[iYear][iPer][iData]->Sumw2 ();
    }
   }
  }

  //int***** nGammaJet = Get5DArray <int> (2, 3, 3, numetabins+1, numpbins+1);

  for (short iYear = 0; iYear < 2; iYear++) {
   if (skipOldInsitu && iYear == 0) continue;

   const TString path = (iYear == 0 ? rootPath + "/2015factors/" : rootPath);

   TSystemDirectory dir (path.Data (), path.Data ());
   TList* sysfiles = dir.GetListOfFiles ();
   if (!sysfiles) {
    cout << "Cannot get list of files! Exiting." << endl;
    return;
   }
   TSystemFile *sysfile;
   TString fname;
   TString histName;
   TIter next (sysfiles);
   int numFiles = 0;
   while ( (sysfile= (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
     if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data () << endl;

     // do this if file is data
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains (to_string (runNumber))) { // if data, do this
       numFiles++;
       cout << "Reading in " << path+fname << endl;
       TFile* thisFile = new TFile (path + fname, "READ");
       const short iPer = (runNumber < 313500 ? 0 : 1);

       for (short iErr = 0; iErr < 3; iErr++) {
        TString error = "sys_lo";
        if (iErr == 1) error = "stat";
        else if (iErr == 2) error = "sys_hi";

        TH3D* temp = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%i_data_%s", runNumber, error.Data ()));
        gJetHists[iYear][iPer][0][iErr]->Add (temp);
        gJetHists[iYear][2][0][iErr]->Add (temp);

        //if (iErr == 1) 
        // gJetHistsSys[iPer][iEta][0][iErr]->Add ( (TH2D*)thisFile->Get (Form ("gJetPtRatioSys_dataSet%i_iEta%i_data_%s", runNumber, iEta, error.Data ())));
        ////gJetHistsSys[2][iEta][0][iErr]->Add ( (TH2D*)thisFile->Get (Form ("gJetPtRatioSys_dataSet%i_iEta%i_data_%s", runNumber, act_iEta, error.Data ())));
       }

       TH2D* temp = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%i_data", runNumber));
       gJetCounts[iYear][iPer][0]->Add (temp);
       gJetCounts[iYear][2][0]->Add (temp);

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if gamma jet OVERLAY MC sample
     for (TString gammaJetOverlaySampleId : gammaJetOverlaySampleIds) { // check for gamma jet MC
      if (fname.Contains (gammaJetOverlaySampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << path+fname << endl;
       TFile* thisFile = new TFile (path + fname, "READ");
       const short iPer = (gammaJetOverlaySampleId.Contains ("pPb") ? 0 : 1);

       TH3D* temp3 = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%s_mc_overlay_stat", gammaJetOverlaySampleId.Data ()));
       gJetHists[iYear][iPer][1][1]->Add (temp3);
       gJetHists[iYear][2][1][1]->Add (temp3);

       TH2D* temp2 = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%s_mc_overlay", gammaJetOverlaySampleId.Data ()));
       gJetCounts[iYear][iPer][1]->Add (temp2);
       gJetCounts[iYear][2][1]->Add (temp2);

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if gamma jet SIGNAL MC sample
     for (TString gammaJetSignalSampleId : gammaJetSignalSampleIds) { // check for gamma jet MC
      if (fname.Contains (gammaJetSignalSampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << path+fname << endl;
       TFile* thisFile = new TFile (path + fname, "READ");
       const short iPer = (gammaJetSignalSampleId.Contains ("pPb") ? 0 : 1);

       TH3D* temp = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%s_mc_signal_stat", gammaJetSignalSampleId.Data ()));
       gJetHists[iYear][iPer][2][1]->Add (temp);
       gJetHists[iYear][2][2][1]->Add (temp);

       TH2D* temp2 = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%s_mc_signal", gammaJetSignalSampleId.Data ()));
       gJetCounts[iYear][iPer][2]->Add (temp2);
       gJetCounts[iYear][2][2]->Add (temp2);

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
    }
   }
   cout << numFiles << " files read in." << endl;
  }
  /**** End loop over input files ****/


  TLine* glines[5] = {};
  TLine* getalines[5] = {};
  TLine* xlines[5] = {};
  //TLine* dplines[5] = {};
  //TLine* dplines_bottom[5] = {};
  //float dpbounds[5] = {35, 50, 70, 140, 280};
  for (short i = 0; i < 5; i++) {
   const float dg = 0.05;
   const float dx = 0.2;

   glines[i] = new TLine (pbins[0], 1.0-1*dg+dg*i, pbins[numpbins], 1.0-1*dg+dg*i);
   getalines[i] = new TLine (etabins[0], 1.0-1*dg+dg*i, etabins[numetabins], 1.0-1*dg+dg*i);
   xlines[i] = new TLine (xjrefbins[0], 1.0-2*dx+dx*i, xjrefbins[numxjrefbins], 1.0-2*dx+dx*i);
   //dplines[i] = new TLine (dpbounds[i], 0.55, dpbounds[i], 1.85);
   //dplines_bottom[i] = new TLine (dpbounds[i], 0.91, dpbounds[i], 1.09);

   if (1.0-1*dg+dg*i == 1) glines[i]->SetLineStyle (1);
   else glines[i]->SetLineStyle (3);
   if (1.0-1*dg+dg*i == 1) getalines[i]->SetLineStyle (1);
   else getalines[i]->SetLineStyle (3);
   if (1.0-2*dx+dx*i == 1) xlines[i]->SetLineStyle (1);
   else xlines[i]->SetLineStyle (3);
   //dplines[i]->SetLineStyle (3);
   //dplines_bottom[i]->SetLineStyle (3);
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
  TH1D *vJetHist = NULL, *vJetHist_mc_overlay = NULL, *vJetHist_mc_signal = NULL, *vJetHist_lo = NULL, *vJetHist_hi = NULL, *vJetHist_rat = NULL, *vJetHist_rat_lo = NULL, *vJetHist_rat_hi = NULL;
  TH2D *proj = NULL, *proj_mc_overlay = NULL, *proj_mc_signal = NULL, *proj_lo = NULL, *proj_hi = NULL;
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

    int eta_lo = iEta+1;
    int eta_hi = iEta+1;
    if (iEta == numetabins) {
     eta_lo = 6;
     eta_hi = 9;
    }

    /**** Plots GammaJet info as a function of p_T^ref****/
    for (short iYear = 0; iYear < 2; iYear++) {
     if (skipOldInsitu && iYear == 0) continue;

     const Style_t dataStyle = (iYear == 0 ? 24 : 20);
     const Style_t mcOverlayStyle = 33;
     const Style_t mcSignalStyle = 34;

     topPad->cd ();
     topPad->SetLogx ();

     proj = Project2D ("", gJetHists[iYear][iPer][0][1], "x", "z", eta_lo, eta_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numpbins, pbins, true);

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
     proj_lo = Project2D ("", gJetHists[iYear][iPer][0][0], "x", "z", eta_lo, eta_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numpbins, pbins, true);

     proj_hi = Project2D ("", gJetHists[iYear][iPer][0][2], "x", "z", eta_lo, eta_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numpbins, pbins, true);

     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     vJetGraph_sys->SetFillColor (dataColor);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", gJetHists[iYear][iPer][1][1], "x", "z", eta_lo, eta_hi);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numpbins, pbins, true);

     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

     if (iPer != 1 && !skipSignalMC) {
      proj_mc_signal = Project2D ("", gJetHists[iYear][iPer][2][1], "x", "z", eta_lo, eta_hi);
      vJetHist_mc_signal = GetProfileX ("vJetHist_mc_signal", proj_mc_signal, numpbins, pbins, true);

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
      const int countsData = gJetCounts[iYear][iPer][0]->Integral (1, gJetCounts[iYear][iPer][0]->GetNbinsX(), eta_lo, eta_hi);
      const int countsMC = gJetCounts[iYear][iPer][1]->Integral (1, gJetCounts[iYear][iPer][1]->GetNbinsX(), eta_lo, eta_hi);
      const int countsMC_sig = gJetCounts[iYear][iPer][2]->Integral (1, gJetCounts[iYear][iPer][2]->GetNbinsX(), eta_lo, eta_hi);

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
       myText (0.155, 0.15, dataColor, Form ("#gamma + Jet, %g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04/uPadY);
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

      vJetHist_rat = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_%s_iEta%i", mcType.Data (), iEta)), proj, proj_mc, numpbins, pbins, true, "x");
      vJetHist_rat_lo = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_lo_%s_iEta%i", mcType.Data (), iEta)), proj_lo, proj_mc, numpbins, pbins, true, "x");
      vJetHist_rat_hi = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_hi_%s_iEta%i", mcType.Data (), iEta)), proj_hi, proj_mc, numpbins, pbins, true, "x");

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

      proj = Project2D ("", gJetHists[iYear][iPer][0][1], "x", "z", eta_lo, eta_hi);
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

      proj_mc_overlay = Project2D ("", gJetHists[iYear][iPer][1][1], "x", "z", eta_lo, eta_hi);
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
       const int countData = gJetCounts[iYear][iPer][0]->GetBinContent (iP, iEta);
       const int countMC = gJetCounts[iYear][iPer][1]->GetBinContent (iP, iEta);

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
     p_lo = 8;
     p_hi = numpbins;
    }


    /**** Plots GammaJet info as a function of eta^jet****/
    for (short iYear = 0; iYear < 2; iYear++) {
     if (skipOldInsitu && iYear == 0) continue;

     const Style_t dataStyle = (iYear == 0 ? 24 : 20);
     const Style_t mcOverlayStyle = 33;
     const Style_t mcSignalStyle = 34;

     topPad->cd ();
     topPad->SetLogx (false);

     proj = Project2D ("", gJetHists[iYear][iPer][0][1], "y", "z", p_lo, p_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numetabins, etabins, true);

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
     proj_lo = Project2D ("", gJetHists[iYear][iPer][0][0], "y", "z", p_lo, p_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numetabins, etabins, true);

     proj_hi = Project2D ("", gJetHists[iYear][iPer][0][2], "y", "z", p_lo, p_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, true);

     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     vJetGraph_sys->SetFillColor (dataColor);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", gJetHists[iYear][iPer][1][1], "y", "z", p_lo, p_hi);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numetabins, etabins, true);

     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

     if (iPer != 1 && !skipSignalMC) {
      proj_mc_signal = Project2D ("", gJetHists[iYear][iPer][2][1], "y", "z", p_lo, p_hi);
      vJetHist_mc_signal = GetProfileX ("vJetHist_mc_signal", proj_mc_signal, numetabins, etabins, true);

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
      const int countsData = gJetCounts[iYear][iPer][0]->Integral (p_lo, p_hi, 1, gJetCounts[iYear][iPer][0]->GetNbinsY());
      const int countsMC = gJetCounts[iYear][iPer][1]->Integral (p_lo, p_hi, 1, gJetCounts[iYear][iPer][1]->GetNbinsY());
      const int countsMC_sig = gJetCounts[iYear][iPer][2]->Integral (p_lo, p_hi, 1, gJetCounts[iYear][iPer][2]->GetNbinsY());

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

      vJetHist_rat = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_iP%i", iP)), proj, proj_mc, numetabins, etabins, true, "x");
      vJetHist_rat_lo = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_lo_iP%i", iP)), proj_lo, proj_mc, numetabins, etabins, true, "x");
      vJetHist_rat_hi = GetDataOverMC (TString (Form ("gJetPtDataMCRatio_hi_iP%i", iP)), proj_hi, proj_mc, numetabins, etabins, true, "x");

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
