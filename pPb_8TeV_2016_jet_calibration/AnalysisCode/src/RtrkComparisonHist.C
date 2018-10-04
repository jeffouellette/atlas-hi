#include "RtrkComparisonHist.h"
#include "Params.h"
#include "Utils.h"

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
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

void RtrkComparisonHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("RtrkComparison/", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);
  vector<TString> gammaJetSampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetSampleIds.push_back (TString ("Pbp") + (runValidation ? "_Signal":"_Overlay") + "_GammaJet_Slice" + to_string (i+1));
   gammaJetSampleIds.push_back (TString ("pPb") + (runValidation ? "_Signal":"_Overlay") + "_GammaJet_Slice" + to_string (i+1));
  }
  vector<TString> zeeJetSampleIds (0);
  zeeJetSampleIds.push_back ("Pbp_Overlay_ZeeJet");
  zeeJetSampleIds.push_back ("pPb_Overlay_ZeeJet");

  vector<TString> zmumuJetSampleIds (0);
  zmumuJetSampleIds.push_back ("Pbp_Overlay_ZmumuJet");
  zmumuJetSampleIds.push_back ("pPb_Overlay_ZmumuJet");

  vector<TString> dijetSampleIds (0);
  dijetSampleIds.push_back ("pPb_Signal_Dijet_Slice2");

  TH2D* jetRtrkHists[2][3][numetabins+1][2][3];

  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
   for (short iEta = 0; iEta <= numetabins; iEta++) {
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

       jetRtrkHists[iAlgo][iPer][iEta][iData][iErr] = new TH2D (Form ("jetRtrkDist_iEta%i_%s_%s_%s_%s", iEta, algo.Data (), dataType.Data (), error.Data (), period.Data ()), "", numpzbins, pzbins, numrtrkbins, rtrkbins);
       jetRtrkHists[iAlgo][iPer][iEta][iData][iErr]->Sumw2 ();
      }
     }
    }
   }
  }

  int* nJet[2][3][2] = {{{}, {}, {}}, {{}, {}, {}}};
  for (int i = 0; i < 2; i++) {
   for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 2; k++) {
     nJet[i][j][k] = new int[numetabins+1];
     for (int iEta = 0; iEta <= numetabins; iEta++) {
      nJet[i][j][k][iEta] = 0;
     }
    }
   }
  }
  

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
   TVectorD *nJetVec;
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
       //infoVec = (TVectorD*)thisFile->Get (Form ("infoVec_%i", runNumber));
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%i", runNumber));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iEta = 0; iEta < numetabins; iEta++) {
         nJet[iAlgo][iPer][0][iEta] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][2][0][iEta] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][iPer][0][numetabins] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][2][0][numetabins] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];

         for (short iErr = 0; iErr < 3; iErr++) {
          TString error = "sys_lo";
          if (iErr == 1) error = "stat";
          else if (iErr == 2) error = "sys_hi";

          TH2D* temp = (TH2D*)thisFile->Get (Form ("jetRtrkDist_dataSet%i_%s_iEta%i_data_%s", runNumber, algo.Data (), iEta, error.Data ()));
          jetRtrkHists[iAlgo][iPer][iEta][0][iErr]->Add (temp);
          jetRtrkHists[iAlgo][2][iEta][0][iErr]->Add (temp);
          jetRtrkHists[iAlgo][iPer][numetabins][0][iErr]->Add (temp);
          jetRtrkHists[iAlgo][2][numetabins][0][iErr]->Add (temp);
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
       //infoVec = (TVectorD*)thisFile->Get (Form ("infoVec_%s", gammaJetSampleId.Data ()));
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", gammaJetSampleId.Data ()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iEta = 0; iEta < numetabins; iEta++) {
         nJet[iAlgo][iPer][1][iEta] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][2][1][iEta] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][iPer][1][numetabins] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][2][1][numetabins] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];

         // Only add the statistical error plots for MC (don't need to consider systematics)
         TH2D* temp = (TH2D*)thisFile->Get (Form ("jetRtrkDist_dataSet%s_%s_iEta%i_mc_stat", gammaJetSampleId.Data (), algo.Data (), iEta));
         jetRtrkHists[iAlgo][iPer][iEta][1][1]->Add (temp);
         jetRtrkHists[iAlgo][2][iEta][1][1]->Add (temp);
         jetRtrkHists[iAlgo][iPer][numetabins][1][1]->Add (temp);
         jetRtrkHists[iAlgo][2][numetabins][1][1]->Add (temp);
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
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", zeeJetSampleId.Data ()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iEta = 0; iEta < numetabins; iEta++) {
         nJet[iAlgo][iPer][1][iEta] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][2][1][iEta] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][iPer][1][numetabins] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][2][1][numetabins] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];

         // Only add the statistical error plots for MC (don't need to consider systematics)
         TH2D* temp = (TH2D*)thisFile->Get (Form ("jetRtrkDist_dataSet%s_%s_iEta%i_mc_stat", zeeJetSampleId.Data (), algo.Data (), iEta));
         jetRtrkHists[iAlgo][iPer][iEta][1][1]->Add (temp);
         jetRtrkHists[iAlgo][2][iEta][1][1]->Add (temp);
         jetRtrkHists[iAlgo][iPer][numetabins][1][1]->Add (temp);
         jetRtrkHists[iAlgo][2][numetabins][1][1]->Add (temp);
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
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", zmumuJetSampleId.Data ()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iEta = 0; iEta < numetabins; iEta++) {
         nJet[iAlgo][iPer][1][iEta] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][2][1][iEta] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][iPer][1][numetabins] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][2][1][numetabins] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];

         // Only add the statistical error plots for MC (don't need to consider systematics)
         TH2D* temp = (TH2D*)thisFile->Get (Form ("jetRtrkDist_dataSet%s_%s_iEta%i_mc_stat", zmumuJetSampleId.Data (), algo.Data (), iEta));
         jetRtrkHists[iAlgo][iPer][iEta][1][1]->Add (temp);
         jetRtrkHists[iAlgo][2][iEta][1][1]->Add (temp);
         jetRtrkHists[iAlgo][iPer][numetabins][1][1]->Add (temp);
         jetRtrkHists[iAlgo][2][numetabins][1][1]->Add (temp);
        }
       }

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if a dijet MC sample
     for (TString dijetSampleId : dijetSampleIds) { // check for Z->ee MC
      if (fname.Contains (dijetSampleId)) { // if Z->ee MC do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (dijetSampleId.Contains ("pPb") ? 0 : 1);
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", dijetSampleId.Data ()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iEta = 0; iEta < numetabins; iEta++) {
         nJet[iAlgo][iPer][1][iEta] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][2][1][iEta] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][iPer][1][numetabins] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];
         nJet[iAlgo][2][1][numetabins] += (*nJetVec)[iEta + iAlgo* (numetabins+1)];

         // Only add the statistical error plots for MC (don't need to consider systematics)
         TH2D* temp = (TH2D*)thisFile->Get (Form ("jetRtrkDist_dataSet%s_%s_iEta%i_mc_stat", dijetSampleId.Data (), algo.Data (), iEta));
         jetRtrkHists[iAlgo][iPer][iEta][1][1]->Add (temp);
         jetRtrkHists[iAlgo][2][iEta][1][1]->Add (temp);
         jetRtrkHists[iAlgo][iPer][numetabins][1][1]->Add (temp);
         jetRtrkHists[iAlgo][2][numetabins][1][1]->Add (temp);
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
  TH1D *jetRtrkHist, *jetRtrkHist_mc, *jetRtrkHist_lo, *jetRtrkHist_hi, *jetRtrkHist_rat, *jetRtrkHist_rat_lo, *jetRtrkHist_rat_hi;
  TGraphAsymmErrors *jetRtrkGraph_sys, *jetRtrkGraph_rat_sys;

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
   TString period = "Period A";
   if (iPer == 1) period = "Period B";
   else if (iPer == 2) period = "Period A+B";

   TString perType = "pA";
   if (iPer == 1) perType = "pB";
   else if (iPer == 2) perType = "pAB";

   for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta

    for (short iAlgo = 0; iAlgo < 2; iAlgo++) { // loop over jet algorithms
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     const Style_t markerStyle = (iAlgo == 0 ? 20 : 24);
     topPad->cd ();
     topPad->SetLogx ();
     jetRtrkHist = GetProfileX (Form ("jetRtrk_Hist_%s_%s_iEta%i", algo.Data (), perType.Data (), iEta), jetRtrkHists[iAlgo][iPer][iEta][0][1], numpzbins, pzbins, true);
     jetRtrkGraph_sys = new TGraphAsymmErrors (jetRtrkHist); // for plotting systematics
     jetRtrkHist->SetYTitle ("< #Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{calo}>");
     jetRtrkHist->SetAxisRange (0., 2.0, "Y");
     jetRtrkHist->SetMarkerStyle (markerStyle);
     jetRtrkHist->SetMarkerColor (dataColor);
     jetRtrkHist->SetLineColor (dataColor);
     jetRtrkHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     jetRtrkHist_lo = GetProfileX ("jetRtrk_Hist_lo", jetRtrkHists[iAlgo][iPer][iEta][0][0], numpzbins, pzbins, true);
     jetRtrkHist_hi = GetProfileX ("jetRtrk_Hist_hi", jetRtrkHists[iAlgo][iPer][iEta][0][2], numpzbins, pzbins, true);
     CalcSystematics (jetRtrkGraph_sys, jetRtrkHist, jetRtrkHist_hi, jetRtrkHist_lo);
     if (jetRtrkHist_lo) delete jetRtrkHist_lo;
     if (jetRtrkHist_hi) delete jetRtrkHist_hi;
     jetRtrkGraph_sys->SetFillColor (kBlack);
     jetRtrkGraph_sys->SetFillStyle (3001);

     jetRtrkHist_mc = GetProfileX (Form ("jetRtrk_Hist_mc_%s_%s_iEta%i", algo.Data (), perType.Data (), iEta), jetRtrkHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, true);
     jetRtrkHist_mc->SetMarkerStyle (markerStyle);
     jetRtrkHist_mc->SetMarkerColor (mcOverlayColor);
     jetRtrkHist_mc->SetLineColor (mcOverlayColor);

     if (iAlgo == 0) jetRtrkHist->DrawCopy ("e1 x0");
     else jetRtrkHist->DrawCopy ("same e1 x0");
     jetRtrkHist_mc->DrawCopy ("same e1 x0");
     jetRtrkGraph_sys->Draw ("2");

     if (iAlgo == 0) {
      myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nJet[iAlgo][iPer][0][iEta]), 1.25, 0.04/uPadY);
      myMarkerText (0.175, 0.81, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", nJet[iAlgo][iPer][1][iEta]), 1.25, 0.04/uPadY);
      if (iEta < numetabins) {
       if (iPer == 2) myText (0.155, 0.65, kBlack, Form ("%g < #eta_{Lab}^{Calo} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
       else myText (0.155, 0.65, kBlack, Form ("%g < #eta_{Lab}^{Calo} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
      }
      myText (0.155, 0.73, kBlack, period.Data (), 0.04/uPadY);
     }

     bottomPad->cd ();
     bottomPad->SetLogx ();

     jetRtrkHist_rat = GetDataOverMC (TString (Form ("jetRtrk_DataMCRatio_%s_%s_iEta%i", algo.Data (), perType.Data (), iEta)), jetRtrkHists[iAlgo][iPer][iEta][0][1], jetRtrkHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, false, "x");
     jetRtrkGraph_rat_sys = new TGraphAsymmErrors (jetRtrkHist_rat);
     jetRtrkHist_rat_lo = GetDataOverMC (TString (Form ("jetRtrk_DataMCRatio_lo_%s_iEta%i", algo.Data (), iEta)), jetRtrkHists[iAlgo][iPer][iEta][0][0], jetRtrkHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, false, "x");
     jetRtrkHist_rat_hi = GetDataOverMC (TString (Form ("jetRtrk_DataMCRatio_hi_%s_iEta%i", algo.Data (), iEta)), jetRtrkHists[iAlgo][iPer][iEta][0][2], jetRtrkHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, false, "x");
     CalcSystematics (jetRtrkGraph_rat_sys, jetRtrkHist_rat, jetRtrkHist_rat_hi, jetRtrkHist_rat_lo);
     if (jetRtrkHist_rat_lo) delete jetRtrkHist_rat_lo;
     if (jetRtrkHist_rat_hi) delete jetRtrkHist_rat_hi;
     jetRtrkGraph_rat_sys->SetFillColor (kBlack);
     jetRtrkGraph_rat_sys->SetFillStyle (3001);

     jetRtrkHist_rat->SetXTitle ("#it{p}_{T}^{calo} #left[GeV#right]");
     jetRtrkHist_rat->SetYTitle ("Data / MC");
     jetRtrkHist_rat->SetAxisRange (0.65, 1.55, "Y");
     jetRtrkHist_rat->SetMarkerStyle (markerStyle);
     jetRtrkHist_rat->GetYaxis ()->SetNdivisions (405);
     jetRtrkHist_rat->GetXaxis ()->SetTitleSize (0.04/dPadY);
     jetRtrkHist_rat->GetYaxis ()->SetTitleSize (0.04/dPadY);
     jetRtrkHist_rat->GetXaxis ()->SetTitleOffset (1);
     jetRtrkHist_rat->GetYaxis ()->SetTitleOffset (dPadY);
     jetRtrkHist_rat->GetYaxis ()->CenterTitle (true);
     jetRtrkHist_rat->GetXaxis ()->SetLabelSize (0.04/dPadY);
     jetRtrkHist_rat->GetYaxis ()->SetLabelSize (0.04/dPadY);
     jetRtrkHist_rat->GetXaxis ()->SetTickLength (0.08);

     if (iAlgo == 0) jetRtrkHist_rat->Draw ("e1 x0");
     else jetRtrkHist_rat->Draw ("same e1 x0");
     jetRtrkGraph_rat_sys->Draw ("2");
     for (TLine* line : zlines) line->Draw ("same");
    }

    char* plotName;
    if (iEta < numetabins) plotName = Form ("jet_rtrk_eta%i.pdf", iEta);
    else plotName = Form ("jet_rtrk_combined.pdf");

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
    if (jetRtrkHist) delete jetRtrkHist;
    if (jetRtrkHist_mc) delete jetRtrkHist_mc;
    if (jetRtrkGraph_sys) delete jetRtrkGraph_sys;
    if (jetRtrkHist_rat) delete jetRtrkHist_rat;
    if (jetRtrkGraph_rat_sys) delete jetRtrkGraph_rat_sys;
   }
  }

  return;
}

} // end namespace
