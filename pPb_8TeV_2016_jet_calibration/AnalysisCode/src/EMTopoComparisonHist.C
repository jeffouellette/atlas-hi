#include "EMTopoComparisonHist.h"
#include "Params.h"
#include "Utils.h"

#include <ArrayTemplates.h>
#include <GlobalParams.h>

#include <TF1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TVectorT.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasUtils.h>
#include <AtlasStyle.h>

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {


void EMTopoComparisonHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("EMTopoComparison/", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);
  vector<TString> gammaJetSampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetSampleIds.push_back (TString ("Pbp_Overlay_GammaJet_Slice") + to_string (i+1));
   gammaJetSampleIds.push_back (TString ("pPb_Overlay_GammaJet_Slice") + to_string (i+1));
  }
  vector<TString> zeeJetSampleIds (0);
  zeeJetSampleIds.push_back ("Pbp_Overlay_ZeeJet");
  zeeJetSampleIds.push_back ("pPb_Overlay_ZeeJet");

  vector<TString> zmumuJetSampleIds (0);
  zmumuJetSampleIds.push_back ("Pbp_Overlay_ZmumuJet");
  zmumuJetSampleIds.push_back ("pPb_Overlay_ZmumuJet");

  vector<TString> dijetSampleIds (0);
  dijetSampleIds.push_back ("pPb_Signal_Dijet_Slice2");

  TH3D***** zeeJetHists = Get4DArray <TH3D*> (2, 3, 2, 3);
  TH3D***** zmumuJetHists = Get4DArray <TH3D*> (2, 3, 2, 3);
  TH3D***** gJetHists = Get4DArray <TH3D*> (2, 3, 2, 3);

  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
   for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
    const TString dataType = (iData == 0 ? "data":"mc");

    for (short iErr = 0; iErr < 3; iErr++) {
     TString error = "sys_lo";
     if (iErr == 1) error = "stat";
     else if (iErr == 2) error = "sys_hi";

     for (short iPer = 0; iPer < 3; iPer++) {
      TString period = "periodA";
      if (iPer == 1) period = "periodB";
      else if (iPer == 2) period = "periodAB";

      zeeJetHists[iAlgo][iPer][iData][iErr] = new TH3D (Form ("zeeJetPtRatio_%s_%s_%s_%s", algo.Data (), dataType.Data (), error.Data (), period.Data ()), "", numpzbins, pzbins, numetabins, etabins, numxjrefbins, xjrefbins);
      zeeJetHists[iAlgo][iPer][iData][iErr]->Sumw2 ();

      zmumuJetHists[iAlgo][iPer][iData][iErr] = new TH3D (Form ("zmumuJetPtRatio_%s_%s_%s_%s", algo.Data (), dataType.Data (), error.Data (), period.Data ()), "", numpzbins, pzbins, numetabins, etabins, numxjrefbins, xjrefbins);
      zmumuJetHists[iAlgo][iPer][iData][iErr]->Sumw2 ();

      gJetHists[iAlgo][iPer][iData][iErr] = new TH3D (Form ("gJetPtRatio_%s_%s_%s_%s", algo.Data (), dataType.Data (), error.Data (), period.Data ()), "", numpbins, pbins, numetabins, etabins, numxjrefbins, xjrefbins);
      gJetHists[iAlgo][iPer][iData][iErr]->Sumw2 ();
     }
    }
   }
  }

  int***** nZeeJet = Get5DArray <int> (2, 3, 2, numetabins+1, numpzbins+1);
  int***** nZmumuJet = Get5DArray <int> (2, 3, 2, numetabins+1, numpzbins+1);
  int***** nGammaJet = Get5DArray <int> (2, 3, 2, numetabins+1, numpbins+1);

  int nTotalJets[2] = {};
  int nCleanJets[2] = {};

  {
   TSystemDirectory dir (rootPath.Data (), rootPath.Data ());
   TList* sysfiles = dir.GetListOfFiles ();
   if (!sysfiles) {
    cout << "Cannot get list of files! Exiting." << endl;
    return;
   }
   TSystemFile *sysfile;
   TString fname;
   TString histName;
   TIter next (sysfiles);
   TVectorD *nZeeJetVec, *nZmumuJetVec, *nGammaJetVec, *jetCleaningVec;
   int numFiles = 0;
   while ( (sysfile= (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
     if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data () << endl;

     // do this if file is data
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains (to_string (runNumber))) { // if data, do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (runNumber < 313500 ? 0 : 1);
       //infoVec = (TVectorD*) thisFile->Get (Form ("infoVec_%i", runNumber));
       nZeeJetVec = (TVectorD*) thisFile->Get (Form ("nZeeJetVec_%i", runNumber));
       nZmumuJetVec = (TVectorD*) thisFile->Get (Form ("nZmumuJetVec_%i", runNumber));
       nGammaJetVec = (TVectorD*) thisFile->Get (Form ("nGammaJetVec_%i", runNumber));

       jetCleaningVec = (TVectorD*) thisFile->Get (Form ("jetCleaningVec_%i", runNumber));
       if (jetCleaningVec) {
        nTotalJets[0] += (int) (*jetCleaningVec)[0];
        nCleanJets[0] += (int) (*jetCleaningVec)[1];
       }

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iErr = 0; iErr < 3; iErr++) {
         TString error = "sys_lo";
         if (iErr == 1) error = "stat";
         else if (iErr == 2) error = "sys_hi";

         TH3D* temp = (TH3D*) thisFile->Get (Form ("zeeJetPtRatio_dataSet%i_%s_data_%s", runNumber, algo.Data (), error.Data ()));
         zeeJetHists[iAlgo][iPer][0][iErr]->Add (temp);
         zeeJetHists[iAlgo][2][0][iErr]->Add (temp);

         temp = (TH3D*) thisFile->Get (Form ("zmumuJetPtRatio_dataSet%i_%s_data_%s", runNumber, algo.Data (), error.Data ()));
         zmumuJetHists[iAlgo][iPer][0][iErr]->Add (temp);
         zmumuJetHists[iAlgo][2][0][iErr]->Add (temp);

         temp = (TH3D*) thisFile->Get (Form ("gJetPtRatio_dataSet%i_%s_data_%s", runNumber, algo.Data (), error.Data ()));
         gJetHists[iAlgo][iPer][0][iErr]->Add (temp);
         gJetHists[iAlgo][2][0][iErr]->Add (temp);

        }

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         //const bool flipEta = runNumber < 313500 && iEta < numetabins;
         const bool flipEta = false;
         const short iEta_flip = (flipEta ? (numetabins - iEta - 1) : iEta);

         for (short iP = 0; iP <= numpzbins; iP++) {
          nZeeJet[iAlgo][iPer][0][iEta][iP] += (*nZeeJetVec)[iAlgo + 2* (iEta + iP* (numetabins+1))];
          nZeeJet[iAlgo][2][0][iEta][iP] += (*nZeeJetVec)[iAlgo + 2* (iEta_flip + iP* (numetabins+1))];

          nZmumuJet[iAlgo][iPer][0][iEta][iP] += (*nZmumuJetVec)[iAlgo + 2* (iEta + iP* (numetabins+1))];
          nZmumuJet[iAlgo][2][0][iEta][iP] += (*nZmumuJetVec)[iAlgo + 2* (iEta_flip + iP* (numetabins+1))];
         }

         for (short iP = 0; iP <= numpbins; iP++) {
          nGammaJet[iAlgo][iPer][0][iEta][iP] += (*nGammaJetVec)[iAlgo + 2* (iEta + iP* (numetabins+1))];
          nGammaJet[iAlgo][2][0][iEta][iP] += (*nGammaJetVec)[iAlgo + 2* (iEta_flip + iP* (numetabins+1))];
         }
        }
       }

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if gamma jet MC sample
     for (TString gammaJetSampleId : gammaJetSampleIds) { // check for gamma jet MC
      if (fname.Contains (gammaJetSampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (gammaJetSampleId.Contains ("pPb") ? 0 : 1);
       nGammaJetVec = (TVectorD*) thisFile->Get (Form ("nGammaJetVec_%s", gammaJetSampleId.Data ()));

       jetCleaningVec = (TVectorD*) thisFile->Get (Form ("jetCleaningVec_%s", gammaJetSampleId.Data ()));
       if (jetCleaningVec) {
        nTotalJets[1] += (int) (*jetCleaningVec)[0];
        nCleanJets[1] += (int) (*jetCleaningVec)[1];
       }

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH3D* temp = (TH3D*) thisFile->Get (Form ("gJetPtRatio_dataSet%s_%s_mc_stat", gammaJetSampleId.Data (), algo.Data ()));
        gJetHists[iAlgo][iPer][1][1]->Add (temp);
        gJetHists[iAlgo][2][1][1]->Add (temp);

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         //const bool flipEta = gammaJetSampleId.Contains ("pPb") && iEta < numetabins;
         const bool flipEta = false;
         const short iEta_flip = (flipEta ? (numetabins - iEta - 1) : iEta); // period A condition

         for (short iP = 0; iP <= numpbins; iP++) {
          nGammaJet[iAlgo][iPer][1][iEta][iP] += (*nGammaJetVec)[iAlgo + 2* (iEta + iP* (numetabins+1))];
          nGammaJet[iAlgo][2][1][iEta][iP] += (*nGammaJetVec)[iAlgo + 2* (iEta_flip + iP* (numetabins+1))];
         }
        }
       }

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if Z->ee MC sample
     for (TString zeeJetSampleId : zeeJetSampleIds) { // check for Z->ee MC
      if (fname.Contains (zeeJetSampleId)) { // if Z->ee MC do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (zeeJetSampleId.Contains ("pPb") ? 0 : 1);
       nZeeJetVec = (TVectorD*) thisFile->Get (Form ("nZeeJetVec_%s", zeeJetSampleId.Data ()));

       jetCleaningVec = (TVectorD*) thisFile->Get (Form ("jetCleaningVec_%s", zeeJetSampleId.Data ()));
       if (jetCleaningVec) {
        nTotalJets[1] += (int) (*jetCleaningVec)[0];
        nCleanJets[1] += (int) (*jetCleaningVec)[1];
       }

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        TH3D* temp = (TH3D*) thisFile->Get (Form ("zeeJetPtRatio_dataSet%s_%s_mc_stat", zeeJetSampleId.Data (), algo.Data ()));
        zeeJetHists[iAlgo][iPer][1][1]->Add (temp);
        zeeJetHists[iAlgo][2][1][1]->Add (temp);

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         //const bool flipEta = zeeJetSampleId.Contains ("pPb") && iEta < numetabins;
         const bool flipEta = false;
         const short iEta_flip = (flipEta ? numetabins - iEta - 1 : iEta); // period A condition

         for (short iP = 0; iP <= numpzbins; iP++) {
          nZeeJet[iAlgo][iPer][1][iEta][iP] += (*nZeeJetVec)[iAlgo + 2* (iEta + iP* (numetabins+1))];
          nZeeJet[iAlgo][2][1][iEta][iP] += (*nZeeJetVec)[iAlgo + 2* (iEta_flip + iP* (numetabins+1))];
         }
        }
       }

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
       nZmumuJetVec = (TVectorD*) thisFile->Get (Form ("nZmumuJetVec_%s", zmumuJetSampleId.Data ()));

       jetCleaningVec = (TVectorD*) thisFile->Get (Form ("jetCleaningVec_%s", zmumuJetSampleId.Data ()));
       if (jetCleaningVec) {
        nTotalJets[1] += (int) (*jetCleaningVec)[0];
        nCleanJets[1] += (int) (*jetCleaningVec)[1];
       }

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        TH3D* temp = (TH3D*) thisFile->Get (Form ("zmumuJetPtRatio_dataSet%s_%s_mc_stat", zmumuJetSampleId.Data (), algo.Data ()));
        zmumuJetHists[iAlgo][iPer][1][1]->Add (temp);
        zmumuJetHists[iAlgo][2][1][1]->Add (temp);

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         //const bool flipEta = zmumuJetSampleId.Contains ("pPb") && iEta < numetabins;
         const bool flipEta = false;
         const short iEta_flip = (flipEta ? numetabins - iEta - 1 : iEta); // period A condition

         for (short iP = 0; iP <= numpzbins; iP++) {
          nZmumuJet[iAlgo][iPer][1][iEta][iP] += (*nZmumuJetVec)[iAlgo + 2* (iEta + iP* (numetabins+1))];
          nZmumuJet[iAlgo][2][1][iEta][iP] += (*nZmumuJetVec)[iAlgo + 2* (iEta_flip + iP* (numetabins+1))];
         }
        }
       }

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
  TLine* xlines[5] = {};
  TLine* dplines[5] = {};
  TLine* dplines_bottom[5] = {};
  float dpbounds[5] = {35, 50, 70, 140, 280};
  for (short i = 0; i < 5; i++) {
   const float dz = 0.1;
   const float dg = 0.05;
   const float dx = 0.2;

   zlines[i] = new TLine (pzbins[0], 1.0-2*dz+dz*i, pzbins[numpzbins], 1.0-2*dz+dz*i);
   glines[i] = new TLine (pbins[0], 1.0-1*dg+dg*i, pbins[numpbins], 1.0-1*dg+dg*i);
   xlines[i] = new TLine (xjrefbins[0], 1.0-2*dx+dx*i, xjrefbins[numxjrefbins], 1.0-2*dx+dx*i);
   dplines[i] = new TLine (dpbounds[i], 0.75, dpbounds[i], 2.15);
   dplines_bottom[i] = new TLine (dpbounds[i], 0.91, dpbounds[i], 1.09);

   if (1.0-2*dz+dz*i == 1) zlines[i]->SetLineStyle (1);
   else zlines[i]->SetLineStyle (3);
   if (1.0-1*dg+dg*i == 1) glines[i]->SetLineStyle (1);
   else glines[i]->SetLineStyle (3);
   if (1.0-2*dx+dx*i == 1) xlines[i]->SetLineStyle (1);
   else xlines[i]->SetLineStyle (3);
   dplines[i]->SetLineStyle (3);
   dplines_bottom[i]->SetLineStyle (3);
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
  TH1D *vJetHist, *vJetHist_mc, *vJetHist_lo, *vJetHist_hi, *vJetHist_rat, *vJetHist_rat_lo, *vJetHist_rat_hi;
  TH2D *proj, *proj_mc, *proj_lo, *proj_hi;
  TGraphAsymmErrors *vJetGraph_sys, *vJetGraph_rat_sys;
  char* plotName;

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
   TString period = "Period A";
   if (iPer == 1) period = "Period B";
   else if (iPer == 2) period = "Period A+B";

   TString perType = "pA";
   if (iPer == 1) perType = "pB";
   else if (iPer == 2) perType = "pAB";

   for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta

    const int eta_lo = (iEta == numetabins ? 1 : iEta+1);
    const int eta_hi = (iEta == numetabins ? numetabins : iEta+1);

    /**** Plot ZmumuJet info ****/
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) { // loop over jet algorithms
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     const Style_t markerStyle = (iAlgo == 0 ? 20 : 24);

     topPad->cd ();
     topPad->SetLogx ();

     proj = Project2D ("", zmumuJetHists[iAlgo][iPer][0][1], "x", "z", eta_lo, eta_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numpzbins, pzbins, false);

     vJetHist->SetYTitle ("<#it{x}_{J}^{ref}>");
     double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;
     vJetHist->SetAxisRange (middle - 0.35, middle + 0.35, "Y");
     vJetHist->SetMarkerStyle (markerStyle);
     vJetHist->SetMarkerColor (dataColor);
     vJetHist->SetLineColor (dataColor);
     vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", zmumuJetHists[iAlgo][iPer][0][0], "x", "z", eta_lo, eta_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numpzbins, pzbins, false);

     proj_hi = Project2D ("", zmumuJetHists[iAlgo][iPer][0][2], "x", "z", eta_lo, eta_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numpzbins, pzbins, false);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     vJetGraph_sys->SetFillColor (kBlack);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc = Project2D ("", zmumuJetHists[iAlgo][iPer][1][1], "x", "z", eta_lo, eta_hi);
     //vJetHist_mc = GetProfileX (Form ("zmumuJetHist_mc_%s_%s_iEta%i", algo.Data (), perType.Data (), iEta), zmumuJetHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, false);
     vJetHist_mc = GetProfileX ("vJetHist_mc", proj_mc, numpzbins, pzbins, false);

     vJetHist_mc->SetMarkerStyle (markerStyle);
     vJetHist_mc->SetMarkerColor (mcOverlayColor);
     vJetHist_mc->SetLineColor (mcOverlayColor);

     if (iAlgo == 0) vJetHist->DrawCopy ("e1 x0");
     else vJetHist->DrawCopy ("same e1 x0");
     vJetHist_mc->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     if (iAlgo == 0) {
      myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nZmumuJet[iAlgo][iPer][0][iEta][numpzbins]), 1.25, 0.04/uPadY);
      myMarkerText (0.175, 0.81, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", nZmumuJet[iAlgo][iPer][1][iEta][numpzbins]), 1.25, 0.04/uPadY);
      myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      if (iEta < numetabins) myText (0.155, 0.15, dataColor, Form ("Z (#mu#mu) + Jet, %g < #eta_{det}^{Jet} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
      else myText (0.155, 0.15, dataColor, "Z (#mu#mu) + Jet", 0.04/uPadY);
      myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);
      myMarkerText (0.65, 0.66, dataColor, kFullCircle, "HI Jets", 1.25, 0.04/uPadY);
      myMarkerText (0.65, 0.59, dataColor, kOpenCircle, "EMTopo Jets", 1.25, 0.04/uPadY);
     }

     bottomPad->cd ();
     bottomPad->SetLogx ();

     vJetHist_rat = GetDataOverMC ("vJetHist_rat", proj, proj_mc, numpzbins, pzbins, false, "x");
     vJetHist_rat_lo = GetDataOverMC ("vJetHist_rat_lo", proj_lo, proj_mc, numpzbins, pzbins, false, "x");
     vJetHist_rat_hi = GetDataOverMC ("vJetHist_rat_hi", proj_hi, proj_mc, numpzbins, pzbins, false, "x");

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc) { delete proj_mc; proj_mc = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
     CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
     if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }

     vJetGraph_rat_sys->SetFillColor (kBlack);
     vJetGraph_rat_sys->SetFillStyle (3001);

     vJetHist_rat->SetXTitle ("#it{p}_{T}^{#mu#mu} #left[GeV#right]");
     vJetHist_rat->SetYTitle ("Data / MC");
     vJetHist_rat->SetAxisRange (0.85, 1.15, "Y");
     vJetHist_rat->SetMarkerStyle (markerStyle);
     //vJetHist_rat->SetAxisRange (0.75, 1.35, "Y");
     vJetHist_rat->GetYaxis ()->SetNdivisions (405);
     vJetHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTitleOffset (1);
     vJetHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
     vJetHist_rat->GetYaxis ()->CenterTitle (true);
     vJetHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTickLength (0.08);

     if (iAlgo == 0) vJetHist_rat->DrawCopy ("e1 x0");
     else vJetHist_rat->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");

     if (iAlgo == 0) for (TLine* line : zlines) line->Draw ("same");

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
    }

    if (iEta < numetabins) plotName = Form ("z_mumu_jet_iEta%i.pdf", iEta);
    else plotName = Form ("z_mumu_jet_iEta_combined.pdf");

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


    /**** Plots ZeeJet info ****/
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     const Style_t markerStyle = (iAlgo == 0 ? 20 : 24);
     topPad->cd ();
     topPad->SetLogx ();

     proj = Project2D ("", zeeJetHists[iAlgo][iPer][0][1], "x", "z", eta_lo, eta_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numpzbins, pzbins, false);

     vJetHist->SetYTitle ("<#it{x}_{J}^{ref}>");
     double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;
     vJetHist->SetAxisRange (middle - 0.35, middle + 0.35, "Y");
     vJetHist->SetMarkerStyle (markerStyle);
     vJetHist->SetMarkerColor (dataColor);
     vJetHist->SetLineColor (dataColor);
     vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", zeeJetHists[iAlgo][iPer][0][0], "x", "z", eta_lo, eta_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numpzbins, pzbins, false);

     proj_hi = Project2D ("", zeeJetHists[iAlgo][iPer][0][2], "x", "z", eta_lo, eta_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numpzbins, pzbins, false);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     vJetGraph_sys->SetFillColor (kBlack);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc = Project2D ("", zeeJetHists[iAlgo][iPer][1][1], "x", "z", eta_lo, eta_hi);
     vJetHist_mc = GetProfileX ("vJetHist_mc", proj_mc, numpzbins, pzbins, false);

     vJetHist_mc->SetMarkerStyle (markerStyle);
     vJetHist_mc->SetMarkerColor (mcOverlayColor);
     vJetHist_mc->SetLineColor (mcOverlayColor);

     if (iAlgo == 0) vJetHist->DrawCopy ("e1 x0");
     else vJetHist->DrawCopy ("same e1 x0");
     vJetHist_mc->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     if (iAlgo == 0) {
      myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nZeeJet[iAlgo][iPer][0][iEta][numpzbins]), 1.25, 0.04/uPadY);
      myMarkerText (0.175, 0.81, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", nZeeJet[iAlgo][iPer][1][iEta][numpzbins]), 1.25, 0.04/uPadY);
      myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      if (iEta < numetabins) myText (0.155, 0.15, dataColor, Form ("Z (ee) + Jet, %g < #eta_{det}^{Jet} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
      else myText (0.155, 0.15, dataColor, "Z (ee) + Jet", 0.04/uPadY);
      myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);
      myMarkerText (0.65, 0.66, dataColor, kFullCircle, "HI Jets", 1.25, 0.04/uPadY);
      myMarkerText (0.65, 0.59, dataColor, kOpenCircle, "EMTopo Jets", 1.25, 0.04/uPadY);
     }

     bottomPad->cd ();
     bottomPad->SetLogx ();

     vJetHist_rat = GetDataOverMC ("vJetHist_rat", proj, proj_mc, numpzbins, pzbins, false, "x");
     vJetHist_rat_lo = GetDataOverMC ("vJetHist_rat_lo", proj_lo, proj_mc, numpzbins, pzbins, false, "x");
     vJetHist_rat_hi = GetDataOverMC ("vJetHist_rat_hi", proj_hi, proj_mc, numpzbins, pzbins, false, "x");

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc) { delete proj_mc; proj_mc = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
     CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
     if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }

     vJetGraph_rat_sys->SetFillColor (kBlack);
     vJetGraph_rat_sys->SetFillStyle (3001);

     vJetHist_rat->SetXTitle ("#it{p}_{T}^{#it{ee}} #left[GeV#right]");
     vJetHist_rat->SetYTitle ("Data / MC");
     vJetHist_rat->SetAxisRange (0.85, 1.15, "Y");
     vJetHist_rat->SetMarkerStyle (markerStyle);
     //vJetHist_rat->SetAxisRange (0.75, 1.35, "Y");
     vJetHist_rat->GetYaxis ()->SetNdivisions (405);
     vJetHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTitleOffset (1);
     vJetHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
     vJetHist_rat->GetYaxis ()->CenterTitle (true);
     vJetHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTickLength (0.08);

     if (iAlgo == 0) vJetHist_rat->DrawCopy ("e1 x0");
     else vJetHist_rat->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     if (iAlgo == 0) for (TLine* line : zlines) line->Draw ("same");

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
    }

    if (iEta < numetabins) plotName = Form ("z_ee_jet_iEta%i.pdf", iEta);
    else plotName = Form ("z_ee_jet_iEta_combined.pdf");
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


    /**** Plots GammaJet info as a function of p_T^ref****/
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     const Style_t markerStyle = (iAlgo == 0 ? 20 : 24);
     topPad->cd ();
     topPad->SetLogx ();

     proj = Project2D ("", gJetHists[iAlgo][iPer][0][1], "x", "z", eta_lo, eta_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numpbins, pbins, true);

     vJetHist->SetYTitle ("<#it{x}_{J}^{ref}>");
     double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;
     vJetHist->SetAxisRange (middle - 0.35, middle + 0.35, "Y");
     vJetHist->SetMarkerStyle (markerStyle);
     vJetHist->SetMarkerColor (dataColor);
     vJetHist->SetLineColor (dataColor);
     vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", gJetHists[iAlgo][iPer][0][0], "x", "z", eta_lo, eta_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numpbins, pbins, true);

     proj_hi = Project2D ("", gJetHists[iAlgo][iPer][0][2], "x", "z", eta_lo, eta_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numpbins, pbins, true);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     vJetGraph_sys->SetFillColor (kBlack);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc = Project2D ("", gJetHists[iAlgo][iPer][1][1], "x", "z", eta_lo, eta_hi);
     vJetHist_mc = GetProfileX ("vJetHist_mc", proj_mc, numpbins, pbins, true);

     vJetHist_mc->SetMarkerStyle (markerStyle);
     vJetHist_mc->SetMarkerColor (mcOverlayColor);
     vJetHist_mc->SetLineColor (mcOverlayColor);

     if (iAlgo == 0) vJetHist->DrawCopy ("e1 x0");
     else vJetHist->DrawCopy ("same e1 x0");
     vJetHist_mc->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");
     if (iAlgo == 0)
      for (TLine* line : dplines) line->Draw ("same");

     if (iAlgo == 0) {
      myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nGammaJet[iAlgo][iPer][0][iEta][numpbins]), 1.25, 0.04/uPadY);
      myMarkerText (0.175, 0.81, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", nGammaJet[iAlgo][iPer][1][iEta][numpbins]), 1.25, 0.04/uPadY);
      myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      if (iEta < numetabins) myText (0.155, 0.15, dataColor, Form ("#gamma + Jet, %g < #eta_{det}^{Jet} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
      else myText (0.155, 0.15, dataColor, "#gamma + Jet", 0.04/uPadY);
      myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);
      myMarkerText (0.65, 0.66, dataColor, kFullCircle, "HI Jets", 1.25, 0.04/uPadY);
      myMarkerText (0.65, 0.59, dataColor, kOpenCircle, "EMTopo Jets", 1.25, 0.04/uPadY);
     }

     bottomPad->cd ();
     bottomPad->SetLogx ();

     vJetHist_rat = GetDataOverMC ("vJetHist_rat", proj, proj_mc, numpbins, pbins, true, "x");
     vJetHist_rat_lo = GetDataOverMC ("vJetHist_rat_lo", proj_lo, proj_mc, numpbins, pbins, true, "x");
     vJetHist_rat_hi = GetDataOverMC ("vJetHist_rat_hi", proj_hi, proj_mc, numpbins, pbins, true, "x");

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc) { delete proj_mc; proj_mc = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
     CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
     if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }

     vJetGraph_rat_sys->SetFillColor (kBlack);
     vJetGraph_rat_sys->SetFillStyle (3001);

     vJetHist_rat->SetXTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
     vJetHist_rat->SetYTitle ("Data / MC");
     vJetHist_rat->SetAxisRange (0.91, 1.09, "Y");
     vJetHist_rat->SetMarkerStyle (markerStyle);
     vJetHist_rat->GetYaxis ()->SetNdivisions (405);
     vJetHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTitleOffset (1);
     vJetHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
     vJetHist_rat->GetYaxis ()->CenterTitle (true);
     vJetHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTickLength (0.08);

     if (iAlgo == 0) vJetHist_rat->DrawCopy ("e1 x0"); 
     else vJetHist_rat->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     if (iAlgo == 0) {
      for (TLine* line : glines) line->Draw ("same");
      for (TLine* line : dplines_bottom) line->Draw ("same");
     }

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
    }

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
    }
   }


   for (short iP = 0; iP <= numpzbins; iP++) { // loop over bins in pt

    const int p_lo = (iP == numpzbins ? 1 : iP+1);
    const int p_hi = (iP == numpzbins ? numpzbins : iP+1);

    /**** Plot ZmumuJet info ****/
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) { // loop over jet algorithms
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     const Style_t markerStyle = (iAlgo == 0 ? 20 : 24);

     topPad->cd ();
     topPad->SetLogx (false);

     proj = Project2D ("", zmumuJetHists[iAlgo][iPer][0][1], "y", "z", p_lo, p_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numetabins, etabins, false);

     vJetHist->SetYTitle ("<#it{x}_{J}^{ref}>");
     double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;
     vJetHist->SetAxisRange (middle - 0.35, middle + 0.35, "Y");
     vJetHist->SetMarkerStyle (markerStyle);
     vJetHist->SetMarkerColor (dataColor);
     vJetHist->SetLineColor (dataColor);
     vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", zmumuJetHists[iAlgo][iPer][0][0], "y", "z", p_lo, p_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numetabins, etabins, false);

     proj_hi = Project2D ("", zmumuJetHists[iAlgo][iPer][0][2], "y", "z", p_lo, p_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, false);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     vJetGraph_sys->SetFillColor (kBlack);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc = Project2D ("", zmumuJetHists[iAlgo][iPer][1][1], "y", "z", p_lo, p_hi);
     //vJetHist_mc = GetProfileX (Form ("zmumuJetHist_mc_%s_%s_iP%i", algo.Data (), perType.Data (), iP), zmumuJetHists[iAlgo][iPer][iP][1][1], numetabins, etabins, false);
     vJetHist_mc = GetProfileX ("vJetHist_mc", proj_mc, numetabins, etabins, false);

     vJetHist_mc->SetMarkerStyle (markerStyle);
     vJetHist_mc->SetMarkerColor (mcOverlayColor);
     vJetHist_mc->SetLineColor (mcOverlayColor);

     if (iAlgo == 0) vJetHist->DrawCopy ("e1 x0");
     else vJetHist->DrawCopy ("same e1 x0");
     vJetHist_mc->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     if (iAlgo == 0) {
      myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nZmumuJet[iAlgo][iPer][0][numetabins][iP]), 1.25, 0.04/uPadY);
      myMarkerText (0.175, 0.81, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", nZmumuJet[iAlgo][iPer][1][numetabins][iP]), 1.25, 0.04/uPadY);
      myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      if (iP < numpzbins) myText (0.155, 0.15, dataColor, Form ("Z (#mu#mu) + Jet, %g < #it{p}_{T}^{Jet} < %g", pzbins[iP], pzbins[iP+1]), 0.04/uPadY);
      else myText (0.155, 0.15, dataColor, "Z (#mu#mu) + Jet", 0.04/uPadY);
      myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);
      myMarkerText (0.65, 0.66, dataColor, kFullCircle, "HI Jets", 1.25, 0.04/uPadY);
      myMarkerText (0.65, 0.59, dataColor, kOpenCircle, "EMTopo Jets", 1.25, 0.04/uPadY);
     }

     bottomPad->cd ();
     bottomPad->SetLogx (false);

     vJetHist_rat = GetDataOverMC ("vJetHist_rat", proj, proj_mc, numetabins, etabins, false, "x");
     vJetHist_rat_lo = GetDataOverMC ("vJetHist_rat_lo", proj_lo, proj_mc, numetabins, etabins, false, "x");
     vJetHist_rat_hi = GetDataOverMC ("vJetHist_rat_hi", proj_hi, proj_mc, numetabins, etabins, false, "x");

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc) { delete proj_mc; proj_mc = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
     CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
     if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }

     vJetGraph_rat_sys->SetFillColor (kBlack);
     vJetGraph_rat_sys->SetFillStyle (3001);

     vJetHist_rat->SetXTitle ("#it{p}_{T}^{#mu#mu} #left[GeV#right]");
     vJetHist_rat->SetYTitle ("Data / MC");
     vJetHist_rat->SetAxisRange (0.85, 1.15, "Y");
     vJetHist_rat->SetMarkerStyle (markerStyle);
     //vJetHist_rat->SetAxisRange (0.75, 1.35, "Y");
     vJetHist_rat->GetYaxis ()->SetNdivisions (405);
     vJetHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTitleOffset (1);
     vJetHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
     vJetHist_rat->GetYaxis ()->CenterTitle (true);
     vJetHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTickLength (0.08);

     if (iAlgo == 0) vJetHist_rat->DrawCopy ("e1 x0");
     else vJetHist_rat->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     if (iAlgo == 0) for (TLine* line : zlines) line->Draw ("same");

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
    }

    if (iP < numpzbins) plotName = Form ("z_mumu_jet_iP%i.pdf", iP);
    else plotName = Form ("z_mumu_jet_iP_combined.pdf");

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


    /**** Plots ZeeJet info ****/
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     const Style_t markerStyle = (iAlgo == 0 ? 20 : 24);
     topPad->cd ();
     topPad->SetLogx (false);

     proj = Project2D ("", zeeJetHists[iAlgo][iPer][0][1], "y", "z", p_lo, p_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numetabins, etabins, false);

     vJetHist->SetYTitle ("<#it{x}_{J}^{ref}>");
     double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;
     vJetHist->SetAxisRange (middle - 0.35, middle + 0.35, "Y");
     vJetHist->SetMarkerStyle (markerStyle);
     vJetHist->SetMarkerColor (dataColor);
     vJetHist->SetLineColor (dataColor);
     vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", zeeJetHists[iAlgo][iPer][0][0], "y", "z", p_lo, p_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numetabins, etabins, false);

     proj_hi = Project2D ("", zeeJetHists[iAlgo][iPer][0][2], "y", "z", p_lo, p_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, false);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     vJetGraph_sys->SetFillColor (kBlack);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc = Project2D ("", zeeJetHists[iAlgo][iPer][1][1], "y", "z", p_lo, p_hi);
     vJetHist_mc = GetProfileX ("vJetHist_mc", proj_mc, numetabins, etabins, false);

     vJetHist_mc->SetMarkerStyle (markerStyle);
     vJetHist_mc->SetMarkerColor (mcOverlayColor);
     vJetHist_mc->SetLineColor (mcOverlayColor);

     if (iAlgo == 0) vJetHist->DrawCopy ("e1 x0");
     else vJetHist->DrawCopy ("same e1 x0");
     vJetHist_mc->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     if (iAlgo == 0) {
      myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nZeeJet[iAlgo][iPer][0][numetabins][iP]), 1.25, 0.04/uPadY);
      myMarkerText (0.175, 0.81, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", nZeeJet[iAlgo][iPer][1][numetabins][iP]), 1.25, 0.04/uPadY);
      myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      if (iP < numpzbins) myText (0.155, 0.15, dataColor, Form ("Z (ee) + Jet, %g < #it{p}_{T}^{Jet} < %g", pzbins[iP], pzbins[iP+1]), 0.04/uPadY);
      else myText (0.155, 0.15, dataColor, "Z (ee) + Jet", 0.04/uPadY);
      myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);
      myMarkerText (0.65, 0.66, dataColor, kFullCircle, "HI Jets", 1.25, 0.04/uPadY);
      myMarkerText (0.65, 0.59, dataColor, kOpenCircle, "EMTopo Jets", 1.25, 0.04/uPadY);
     }

     bottomPad->cd ();
     bottomPad->SetLogx (false);

     vJetHist_rat = GetDataOverMC ("vJetHist_rat", proj, proj_mc, numetabins, etabins, false, "x");
     vJetHist_rat_lo = GetDataOverMC ("vJetHist_rat_lo", proj_lo, proj_mc, numetabins, etabins, false, "x");
     vJetHist_rat_hi = GetDataOverMC ("vJetHist_rat_hi", proj_hi, proj_mc, numetabins, etabins, false, "x");

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc) { delete proj_mc; proj_mc = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
     CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
     if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }

     vJetGraph_rat_sys->SetFillColor (kBlack);
     vJetGraph_rat_sys->SetFillStyle (3001);

     vJetHist_rat->SetXTitle ("#it{p}_{T}^{#it{ee}} #left[GeV#right]");
     vJetHist_rat->SetYTitle ("Data / MC");
     vJetHist_rat->SetAxisRange (0.85, 1.15, "Y");
     vJetHist_rat->SetMarkerStyle (markerStyle);
     //vJetHist_rat->SetAxisRange (0.75, 1.35, "Y");
     vJetHist_rat->GetYaxis ()->SetNdivisions (405);
     vJetHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTitleOffset (1);
     vJetHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
     vJetHist_rat->GetYaxis ()->CenterTitle (true);
     vJetHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTickLength (0.08);

     if (iAlgo == 0) vJetHist_rat->DrawCopy ("e1 x0");
     else vJetHist_rat->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     if (iAlgo == 0) for (TLine* line : zlines) line->Draw ("same");

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
    }

    if (iP < numpzbins) plotName = Form ("z_ee_jet_iP%i.pdf", iP);
    else plotName = Form ("z_ee_jet_iP_combined.pdf");
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

   }

   for (short iP = 0; iP <= numpbins; iP++) {

    const int p_lo = (iP == numpbins ? 1 : iP+1);
    const int p_hi = (iP == numpbins ? numpbins : iP+1);

    /**** Plots GammaJet info as a function of p_T^ref****/
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     const Style_t markerStyle = (iAlgo == 0 ? 20 : 24);
     topPad->cd ();
     topPad->SetLogx (false);

     proj = Project2D ("", gJetHists[iAlgo][iPer][0][1], "y", "z", p_lo, p_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numetabins, etabins, true);

     vJetHist->SetYTitle ("<#it{x}_{J}^{ref}>");
     double middle = 0.05 * floor (20 * vJetHist->Integral () / vJetHist->GetNbinsX ()); // gets mean along y
     if (10 * middle != floor (10*middle)) middle += 0.05;
     vJetHist->SetAxisRange (middle - 0.35, middle + 0.35, "Y");
     vJetHist->SetMarkerStyle (markerStyle);
     vJetHist->SetMarkerColor (dataColor);
     vJetHist->SetLineColor (dataColor);
     vJetHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     vJetHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     proj_lo = Project2D ("", gJetHists[iAlgo][iPer][0][0], "y", "z", p_lo, p_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numetabins, etabins, true);

     proj_hi = Project2D ("", gJetHists[iAlgo][iPer][0][2], "y", "z", p_lo, p_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, true);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     vJetGraph_sys->SetFillColor (kBlack);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc = Project2D ("", gJetHists[iAlgo][iPer][1][1], "y", "z", p_lo, p_hi);
     vJetHist_mc = GetProfileX ("vJetHist_mc", proj_mc, numetabins, etabins, true);
     vJetHist_mc->SetMarkerStyle (markerStyle);
     vJetHist_mc->SetMarkerColor (mcOverlayColor);
     vJetHist_mc->SetLineColor (mcOverlayColor);

     if (iAlgo == 0) vJetHist->DrawCopy ("e1 x0");
     else vJetHist->DrawCopy ("same e1 x0");
     vJetHist_mc->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");
     for (TLine* line : dplines) line->Draw ("same");

     if (iAlgo == 0) {
      myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nGammaJet[iAlgo][iPer][0][numetabins][iP]), 1.25, 0.04/uPadY);
      myMarkerText (0.175, 0.81, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", nGammaJet[iAlgo][iPer][1][numetabins][iP]), 1.25, 0.04/uPadY);
      myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
      if (iP < numpbins) myText (0.155, 0.15, dataColor, Form ("#gamma + Jet, %g < #it{p}_{T}^{Jet} < %g", pbins[iP], pbins[iP+1]), 0.04/uPadY);
      else myText (0.155, 0.15, dataColor, "#gamma + Jet", 0.04/uPadY);
      myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);
      myMarkerText (0.65, 0.66, dataColor, kFullCircle, "HI Jets", 1.25, 0.04/uPadY);
      myMarkerText (0.65, 0.59, dataColor, kOpenCircle, "EMTopo Jets", 1.25, 0.04/uPadY);
     }

     bottomPad->cd ();
     bottomPad->SetLogx (false);

     vJetHist_rat = GetDataOverMC ("vJetHist_rat", proj, proj_mc, numetabins, etabins, true, "x");
     vJetHist_rat_lo = GetDataOverMC ("vJetHist_rat_lo", proj_lo, proj_mc, numetabins, etabins, true, "x");
     vJetHist_rat_hi = GetDataOverMC ("vJetHist_rat_hi", proj_hi, proj_mc, numetabins, etabins, true, "x");

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc) { delete proj_mc; proj_mc = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

     vJetGraph_rat_sys = new TGraphAsymmErrors (vJetHist_rat);
     CalcSystematics (vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) { delete vJetHist_rat_lo; vJetHist_rat_lo = NULL; }
     if (vJetHist_rat_hi) { delete vJetHist_rat_hi; vJetHist_rat_hi = NULL; }

     vJetGraph_rat_sys->SetFillColor (kBlack);
     vJetGraph_rat_sys->SetFillStyle (3001);

     vJetHist_rat->SetXTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
     vJetHist_rat->SetYTitle ("Data / MC");
     vJetHist_rat->SetAxisRange (0.91, 1.09, "Y");
     vJetHist_rat->SetMarkerStyle (markerStyle);
     vJetHist_rat->GetYaxis ()->SetNdivisions (405);
     vJetHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTitleOffset (1);
     vJetHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
     vJetHist_rat->GetYaxis ()->CenterTitle (true);
     vJetHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     vJetHist_rat->GetXaxis ()->SetTickLength (0.08);

     if (iAlgo == 0) vJetHist_rat->DrawCopy ("e1 x0"); 
     else vJetHist_rat->DrawCopy ("same e1 x0");
     ((TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     if (iAlgo == 0) {
      for (TLine* line : glines) line->Draw ("same");
      for (TLine* line : dplines_bottom) line->Draw ("same");
     }

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc) { delete vJetHist_mc; vJetHist_mc = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
    }

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
    }
   }
  }

  Delete5DArray (nZeeJet, 2, 3, 2, numetabins+1, numpzbins+1);
  Delete5DArray (nZmumuJet, 2, 3, 2, numetabins+1, numpzbins+1);
  Delete5DArray (nGammaJet, 2, 3, 2, numetabins+1, numpbins+1);

  Delete4DArray (zeeJetHists, 2, 3, 2, 3);
  Delete4DArray (zmumuJetHists, 2, 3, 2, 3);
  Delete4DArray (gJetHists, 2, 3, 2, 3);

  double num = (double)nCleanJets[0];
  double den = (double)nTotalJets[0];
  double frac = num / den;
  double fracErr = sqrt ( ( (num+1.)* (num+2.)) / ( (den+2.)* (den+3.)) - ( (num+1.)* (num+1.))/ ( (den+2.)* (den+2.)) );

  cout << "Clean / total jets in data = " << nCleanJets[0]
       << " / " << nTotalJets[0]
       << " = " << frac
       << " +/- " << fracErr
       << endl;

  num = (double)nCleanJets[1];
  den = (double)nTotalJets[1];
  frac = num / den;
  fracErr = sqrt ( ( (num+1.)* (num+2.)) / ( (den+2.)* (den+3.)) - ( (num+1.)* (num+1.))/ ( (den+2.)* (den+2.)) );

  cout << "Clean / total jets in MC   = " << nCleanJets[1]
       << " / " << nTotalJets[1]
       << " = " << frac
       << " +/- " << fracErr
       << endl;

  return;
}

} // end namespace
