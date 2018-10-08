#include "ZGammaJetCrossCheckHist.h"
#include "Params.h"
#include "Utils.h"

#include <GlobalParams.h>
#include <ArrayTemplates.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TFile.h>
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

/**
 * Equivalent to adding h2 to h1, except with the weight defined by a TF1 at each x value.
 * Could be modified to take in a different weight at each x, y value or x, y, z value.
 * To ensure different histograms are reweighted correctly, they are normalized by default.
 */
void AddWx (TF1* fx, TH3D* h1, TH3D* h2, const bool normalize = true) {
  if (h2->Integral() == 0) {
   return; // true if there is nothing to add
  }
  for (int ix = 1; ix <= h1->GetNbinsX(); ix++) {
   const double wx = fx->Eval (h1->GetXaxis ()->GetBinCenter (ix)) / h2->Integral(); // this cannot = 0
   for (int iy = 1; iy <= h1->GetNbinsY(); iy++) {
    for (int iz = 1; iz <= h1->GetNbinsZ(); iz++) {
     h1->SetBinContent (ix, iy, iz, h1->GetBinContent (ix, iy, iz) + wx * h2->GetBinContent (ix, iy, iz));
     h1->SetBinError (ix, iy, iz, sqrt (pow (h1->GetBinError (ix, iy, iz), 2) + pow (wx * h2->GetBinError (ix, iy, iz), 2)));
    }
   }
  }
  return;
}


void ZGammaJetCrossCheckHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("ZGammaJetCrossCheck/", "pPb_8TeV_2016_jet_calibration/");

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
  vector<TString> zeeJetSampleIds (0);
  zeeJetSampleIds.push_back ("Pbp_Overlay_ZeeJet");
  zeeJetSampleIds.push_back ("pPb_Overlay_ZeeJet");

  vector<TString> zmumuJetSampleIds (0);
  zmumuJetSampleIds.push_back ("Pbp_Overlay_ZmumuJet");
  zmumuJetSampleIds.push_back ("pPb_Overlay_ZmumuJet");

  TH3D***** vJetHists = Get4DArray <TH3D*> (2, 3, 3, 3);
  TH2D**** vJetCounts = Get3DArray <TH2D*> (2, 3, 3);

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

      vJetHists[iYear][iPer][iData][iErr] = new TH3D (Form ("vJetPtRatio_%s_%s_%s_%s", year.Data (), dataType.Data (), error.Data (), period.Data ()), "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
      vJetHists[iYear][iPer][iData][iErr]->Sumw2();
     }

     vJetCounts[iYear][iPer][iData] = new TH2D (Form ("vJetCounts_%s_%s_%s", year.Data (), dataType.Data (), period.Data ()), "", numpbins, pbins, numetabins, etabins);
     vJetCounts[iYear][iPer][iData]->Sumw2();
    }
   }
  }

  TF1* gWeight = new TF1 ("gWeight", "1 / (1 + exp (-[0] * (x - [1])))", 0, pbins[numpbins]);
  gWeight->SetParNames ("g");
  gWeight->SetParameter (0, 5); // guess- changes significantly over ~10 GeV?
  gWeight->SetParameter (1, 60); // guess- centered at 80 GeV?

  TF1* zWeight = new TF1 ("zWeight", "1-0.5*g", 0, pbins[numpbins]);

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
     if (debugStatements) cout << "Status: In ZGammaJetCrossCheckHist.C: Found " << fname.Data () << endl;

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

        TH3D* temp3 = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%i_data_%s", runNumber, error.Data ()));
        AddWx (gWeight, vJetHists[iYear][iPer][0][iErr], temp3);
        AddWx (gWeight, vJetHists[iYear][2][0][iErr], temp3);

        temp3 = (TH3D*)thisFile->Get (Form ("zeeJetPtRatio_dataSet%i_data_%s", runNumber, error.Data ()));
        AddWx (zWeight, vJetHists[iYear][iPer][0][iErr], temp3);
        AddWx (zWeight, vJetHists[iYear][2][0][iErr], temp3);

        temp3 = (TH3D*)thisFile->Get (Form ("zmumuJetPtRatio_dataSet%i_data_%s", runNumber, error.Data ()));
        AddWx (zWeight, vJetHists[iYear][iPer][0][iErr], temp3);
        AddWx (zWeight, vJetHists[iYear][2][0][iErr], temp3);
       }

       TH2D* temp2 = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%i_data", runNumber));
       vJetCounts[iYear][iPer][0]->Add (temp2);
       vJetCounts[iYear][2][0]->Add (temp2);

       temp2 = (TH2D*)thisFile->Get (Form ("zeeJetCounts_dataSet%i_data", runNumber));
       vJetCounts[iYear][iPer][0]->Add (temp2);
       vJetCounts[iYear][2][0]->Add (temp2);

       temp2 = (TH2D*)thisFile->Get (Form ("zmumuJetCounts_dataSet%i_data", runNumber));
       vJetCounts[iYear][iPer][0]->Add (temp2);
       vJetCounts[iYear][2][0]->Add (temp2);

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
       AddWx (gWeight, vJetHists[iYear][iPer][1][1], temp3);
       AddWx (gWeight, vJetHists[iYear][2][1][1], temp3);

       TH2D* temp2 = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%s_mc_overlay", gammaJetOverlaySampleId.Data ()));
       vJetCounts[iYear][iPer][1]->Add (temp2);
       vJetCounts[iYear][2][1]->Add (temp2);

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

       TH3D* temp3 = (TH3D*)thisFile->Get (Form ("gJetPtRatio_dataSet%s_mc_signal_stat", gammaJetSignalSampleId.Data ()));
       AddWx (gWeight, vJetHists[iYear][iPer][2][1], temp3);
       AddWx (gWeight, vJetHists[iYear][2][2][1], temp3);

       TH2D* temp2 = (TH2D*)thisFile->Get (Form ("gJetCounts_dataSet%s_mc_signal", gammaJetSignalSampleId.Data ()));
       vJetCounts[iYear][iPer][2]->Add (temp2);
       vJetCounts[iYear][2][2]->Add (temp2);

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if Z->ee MC sample
     for (TString zeeJetSampleId : zeeJetSampleIds) { // check for Z->ee MC
      if (fname.Contains (zeeJetSampleId)) { // if Z->ee MC do this
       numFiles++;
       cout << "Reading in " << path+fname << endl;
       TFile* thisFile = new TFile (path + fname, "READ");
       const short iPer = (zeeJetSampleId.Contains ("pPb") ? 0 : 1);

       TH3D* temp3 = (TH3D*)thisFile->Get (Form ("zeeJetPtRatio_dataSet%s_mc_overlay_stat", zeeJetSampleId.Data ()));
       AddWx (zWeight, vJetHists[iYear][iPer][1][1], temp3);
       AddWx (zWeight, vJetHists[iYear][2][1][1], temp3);

       TH2D* temp2 = (TH2D*)thisFile->Get (Form ("zeeJetCounts_dataSet%s_mc_overlay", zeeJetSampleId.Data ()));
       vJetCounts[iYear][iPer][1]->Add (temp2);
       vJetCounts[iYear][2][1]->Add (temp2);

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if Z->mumu sample
     for (TString zmumuJetSampleId : zmumuJetSampleIds) { // check for Z->mumu MC
      if (fname.Contains (zmumuJetSampleId)) { // if Z->mumu sample do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (zmumuJetSampleId.Contains ("pPb") ? 0 : 1);

       TH3D* temp3 = (TH3D*)thisFile->Get (Form ("zmumuJetPtRatio_dataSet%s_mc_overlay_stat", zmumuJetSampleId.Data ()));
       AddWx (zWeight, vJetHists[iYear][iPer][1][1], temp3);
       AddWx (zWeight, vJetHists[iYear][2][1][1], temp3);

       TH2D* temp2 = (TH2D*)thisFile->Get (Form ("zmumuJetCounts_dataSet%s_mc_overlay", zmumuJetSampleId.Data ()));
       vJetCounts[iYear][iPer][1]->Add (temp2);
       vJetCounts[iYear][2][1]->Add (temp2);

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


  TLine* zlines[5] = {};
  TLine* glines[5] = {};
  TLine* zetalines[5] = {};
  TLine* getalines[5] = {};
  TLine* xlines[5] = {};
  //TLine* dplines[5] = {};
  //TLine* dplines_bottom[5] = {};
  //float dpbounds[5] = {35, 50, 70, 140, 280};
  for (short i = 0; i < 5; i++) {
   const float dz = 0.05;
   const float dg = 0.05;
   const float dx = 0.2;

   zlines[i] = new TLine (pbins[0], 1.0-2*dz+dz*i, pbins[numpbins], 1.0-2*dz+dz*i);
   glines[i] = new TLine (pbins[0], 1.0-1*dg+dg*i, pbins[numpbins], 1.0-1*dg+dg*i);
   zetalines[i] = new TLine (etabins[0], 1.0-2*dz+dz*i, etabins[numetabins], 1.0-2*dz+dz*i);
   getalines[i] = new TLine (etabins[0], 1.0-1*dg+dg*i, etabins[numetabins], 1.0-1*dg+dg*i);
   xlines[i] = new TLine (xjrefbins[0], 1.0-2*dx+dx*i, xjrefbins[numxjrefbins], 1.0-2*dx+dx*i);
   //dplines[i] = new TLine (dpbounds[i], 0.55, dpbounds[i], 1.85);
   //dplines_bottom[i] = new TLine (dpbounds[i], 0.91, dpbounds[i], 1.09);

   if (1.0-2*dz+dz*i == 1) zlines[i]->SetLineStyle (1);
   else zlines[i]->SetLineStyle (3);
   if (1.0-1*dg+dg*i == 1) glines[i]->SetLineStyle (1);
   else glines[i]->SetLineStyle (3);
   if (1.0-2*dz+dz*i == 1) zetalines[i]->SetLineStyle (1);
   else zetalines[i]->SetLineStyle (3);
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
     eta_lo = 6;//1;
     eta_hi = 9;//numetabins;
    }

    /**** Plot combined VJet info ****/
    for (short iYear = 0; iYear < 2; iYear++) {
     if (skipOldInsitu && iYear == 0) continue;

     const Style_t dataStyle = (iYear == 0 ? 24 : 20);
     const Style_t mcOverlayStyle = 33;
     const Style_t mcSignalStyle = 34;

     topPad->cd ();
     topPad->SetLogx ();

     proj = Project2D ("", vJetHists[1][iPer][0][1], "x", "z", eta_lo, eta_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numpbins, pbins, false);

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

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", vJetHists[1][iPer][0][0], "x", "z", eta_lo, eta_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numpbins, pbins, false);

     proj_hi = Project2D ("", vJetHists[1][iPer][0][2], "x", "z", eta_lo, eta_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numpbins, pbins, false);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     vJetGraph_sys->SetFillColor (dataColor);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", vJetHists[1][iPer][1][1], "x", "z", eta_lo, eta_hi);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numpbins, pbins, false);

     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

     if (iPer != 1 && !skipSignalMC) {
      proj_mc_signal = Project2D ("", vJetHists[iYear][iPer][2][1], "x", "z", eta_lo, eta_hi);
      vJetHist_mc_signal = GetProfileX ("vJetHist_mc_signal", proj_mc_signal, numpbins, pbins, true);

      vJetHist_mc_signal->SetMarkerStyle (mcSignalStyle);
      vJetHist_mc_signal->SetMarkerColor (mcSignalColor);
      vJetHist_mc_signal->SetLineColor (mcSignalColor);
     }

     if (skipOldInsitu || iYear == 0) {
      vJetHist->DrawCopy ("e1 x0");
      vJetHist_mc_overlay->DrawCopy ("same e1 x0"); // insitu factors are not applied to MC
     }
     else {
      vJetHist->DrawCopy ("same e1 x0");
     }
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     const int countsData = vJetCounts[1][iPer][0]->Integral (1, vJetCounts[1][iPer][0]->GetNbinsX(), eta_lo, eta_hi);
     const int countsMC = vJetCounts[1][iPer][1]->Integral (1, vJetCounts[1][iPer][1]->GetNbinsX(), eta_lo, eta_hi);
     const int countsMC_sig = vJetCounts[1][iPer][2]->Integral (1, vJetCounts[1][iPer][2]->GetNbinsX(), eta_lo, eta_hi);

     myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countsData), 1.25, 0.04/uPadY);
     myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
     if (iPer != 1 && !skipSignalMC)
      myMarkerText (0.175, 0.74, mcSignalColor, kFullCross, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", countsMC_sig), 1.25, 0.04/uPadY);
     myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
     if (eta_lo != 1 || eta_hi != numetabins)
      myText (0.155, 0.15, dataColor, Form ("V + Jet, %g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04/uPadY);
     else
      myText (0.155, 0.15, dataColor, "V + Jet", 0.04/uPadY);
     myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);

     bottomPad->cd ();
     bottomPad->SetLogx ();

     for (int iMC = 1; iMC < 3; iMC++) { // loops over overlay then signal MC
      if (iMC == 2 && (skipSignalMC || iPer == 1)) continue; // no signal only for period B!

      const TString mcType = (iMC == 1 ? "overlay":"signal");
      TH2D* proj_mc = (iMC == 1 ? proj_mc_overlay:proj_mc_signal);

      vJetHist_rat = GetDataOverMC (TString (Form ("vJetPtDataMCRatio_iEta%i", iEta)), proj, proj_mc, numpbins, pbins, false, "x");
      vJetHist_rat_lo = GetDataOverMC (TString (Form ("vJetPtDataMCRatio_lo_iEta%i", iEta)), proj_lo, proj_mc, numpbins, pbins, false, "x");
      vJetHist_rat_hi = GetDataOverMC (TString (Form ("vJetPtDataMCRatio_hi_iEta%i", iEta)), proj_hi, proj_mc, numpbins, pbins, false, "x");

      vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
      CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
      if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
      if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }
      vJetGraph_rat_sys->SetFillColor (dataColor);
      vJetGraph_rat_sys->SetFillStyle (3001);

      vJetHist_rat->SetXTitle ("#it{p}_{T}^{Z} #left[GeV#right]");
      vJetHist_rat->SetYTitle ("Data / MC");
      vJetHist_rat->SetAxisRange (0.85, 1.15, "Y");
      //vJetHist_rat->SetAxisRange (0.75, 1.35, "Y");
      vJetHist_rat->SetMarkerStyle (dataStyle);
      vJetHist_rat->SetMarkerColor (dataColor);
      vJetHist_rat->SetLineColor (dataColor);
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

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_mc_signal) { delete proj_mc_signal; proj_mc_signal = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetHist_mc_signal) { delete vJetHist_mc_signal; vJetHist_mc_signal = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
    } // end loop over insitu configurations
     
    if (iEta < numetabins) plotName = Form ("v_jet_iEta%i.pdf", iEta);
    else plotName = Form ("v_jet_iEta_combined.pdf");

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
    

   } // end loop over etabins


   /**** Now loop over pT bins and plot response as function of eta^jet ****/
   for (short iP = 0; iP <= numpbins; iP++) {

    int p_lo = iP+1;
    int p_hi = iP+1;
    if (iP == numpbins) {
     p_lo = 1;
     p_hi = numpbins;
    }

    /**** Plot combined VJet info ****/
    for (short iYear = 0; iYear < 2; iYear++) {
     if (skipOldInsitu && iYear == 0) continue;

     const Style_t dataStyle = (iYear == 0 ? 24 : 20);
     const Style_t mcOverlayStyle = 33;
     const Style_t mcSignalStyle = 34;

     topPad->cd ();
     topPad->SetLogx (false);

     proj = Project2D ("", vJetHists[1][iPer][0][1], "y", "z", p_lo, p_hi);
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

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", vJetHists[1][iPer][0][0], "y", "z", p_lo, p_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj, numetabins, etabins, false);

     proj_hi = Project2D ("", vJetHists[1][iPer][0][2], "y", "z", p_lo, p_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, false);

     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }
     vJetGraph_sys->SetFillColor (dataColor);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", vJetHists[1][iPer][1][1], "y", "z", p_lo, p_hi);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numetabins, etabins, false);

     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

     if (iPer != 1 && !skipSignalMC) {
      proj_mc_signal = Project2D ("", vJetHists[iYear][iPer][2][1], "y", "z", p_lo, p_hi);
      vJetHist_mc_signal = GetProfileX ("vJetHist_mc_signal", proj_mc_signal, numpbins, pbins, true);

      vJetHist_mc_signal->SetMarkerStyle (mcSignalStyle);
      vJetHist_mc_signal->SetMarkerColor (mcSignalColor);
      vJetHist_mc_signal->SetLineColor (mcSignalColor);
     }

     if (skipOldInsitu || iYear == 0) {
      vJetHist->DrawCopy ("e1 x0");
      vJetHist_mc_overlay->DrawCopy ("same e1 x0"); // insitu factors are not applied to MC
     }
     else {
      vJetHist->DrawCopy ("same e1 x0");
     }
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     const int countsData = vJetCounts[1][iPer][0]->Integral (p_lo, p_hi, 1, vJetCounts[1][iPer][0]->GetNbinsY());
     const int countsMC = vJetCounts[1][iPer][1]->Integral (p_lo, p_hi, 1, vJetCounts[1][iPer][1]->GetNbinsY());
     const int countsMC_sig = vJetCounts[1][iPer][2]->Integral (p_lo, p_hi, 1, vJetCounts[1][iPer][2]->GetNbinsY());

     myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countsData), 1.25, 0.04/uPadY);
     myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", countsMC), 1.25, 0.04/uPadY);
     if (iPer != 1 && !skipSignalMC)
      myMarkerText (0.175, 0.74, mcSignalColor, kFullCross, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", countsMC_sig), 1.25, 0.04/uPadY);
     myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
     if (p_lo != 1 || p_hi != numpbins)
      myText (0.155, 0.15, dataColor, Form ("V + Jet, %g < #it{p}_{T}^{V} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
     else
      myText (0.155, 0.15, dataColor, "V + Jet", 0.04/uPadY);
     myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);

     bottomPad->cd ();
     bottomPad->SetLogx (false);

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

      vJetHist_rat->SetXTitle ("Jet #eta");
      vJetHist_rat->SetYTitle ("Data / MC");
      vJetHist_rat->SetAxisRange (0.85, 1.15, "Y");
      //vJetHist_rat->SetAxisRange (0.75, 1.35, "Y");
      vJetHist_rat->SetMarkerStyle (dataStyle);
      vJetHist_rat->SetMarkerColor (dataColor);
      vJetHist_rat->SetLineColor (dataColor);
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

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_mc_signal) { delete proj_mc_signal; proj_mc_signal = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetHist_mc_signal) { delete vJetHist_mc_signal; vJetHist_mc_signal = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
    } // end loop over insitu configurations

    if (iP < numpbins) plotName = Form ("v_jet_iP%i.pdf", iP);
    else plotName = Form ("v_jet_iP_combined.pdf");

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


  return;
}

} // end namespace
