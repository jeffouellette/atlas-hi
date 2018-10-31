#include "RtrkComparisonHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

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

  TH3D***** jetRtrkHists = Get4DArray <TH3D*> (2, 3, 2, 3); // iAlgo, iPer, iData, iErr
  TH2D**** jetRtrkCounts = Get3DArray <TH2D*> (2, 3, 2);

  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
   for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
    const TString dataType = (iData == 0 ? "data":"mc");

    for (short iPer = 0; iPer < 3; iPer++) {
     TString period = "periodA";
     if (iPer == 1) period = "periodB";
     else if (iPer == 2) period = "periodAB";

     for (short iErr = 0; iErr < 3; iErr++) {
      TString error = "sys_lo";
      if (iErr == 1) error = "stat";
      else if (iErr == 2) error = "sys_hi";

      jetRtrkHists[iAlgo][iPer][iData][iErr] = new TH3D (Form ("jetRtrkDist_%s_%s_%s_%s", algo.Data (), dataType.Data (), error.Data (), period.Data ()), "", numpbins, pbins, numetabins, etabins, numrtrkbins, rtrkbins);
      jetRtrkHists[iAlgo][iPer][iData][iErr]->Sumw2 ();
     }

     jetRtrkCounts[iAlgo][iPer][iData] = new TH2D (Form ("jetRtrkCounts_%s_%s_%s", algo.Data (), dataType.Data (), period.Data ()), "", numpbins, pbins, numetabins, etabins);
     jetRtrkCounts[iAlgo][iPer][iData]->Sumw2 ();
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
   TString fname, histName;
   TIter next (sysfiles);
   int numFiles = 0;
   while ( (sysfile= (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
     if (debugStatements) cout << "Status: In RtrkComparisonHist.C: Found " << fname.Data () << endl;

     // do this if file is data
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains (to_string (runNumber))) { // if data, do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       const short iPer = (runNumber < 313500 ? 0 : 1);

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iErr = 0; iErr < 3; iErr++) {
         TString error = "sys_lo";
         if (iErr == 1) error = "stat";
         else if (iErr == 2) error = "sys_hi";

         TH3D* temp3 = (TH3D*)thisFile->Get (Form ("jetRtrkDist_dataSet%i_%s_data_%s", runNumber, algo.Data (), error.Data ()));
         jetRtrkHists[iAlgo][iPer][0][iErr]->Add (temp3);
         jetRtrkHists[iAlgo][2][0][iErr]->Add (temp3);
        }
        TH2D* temp2 = (TH2D*)thisFile->Get (Form ("jetRtrkCounts_dataSet%i_%s_data", runNumber, algo.Data ()));
        jetRtrkCounts[iAlgo][iPer][0]->Add (temp2);
        jetRtrkCounts[iAlgo][2][0]->Add (temp2);
        
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

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH3D* temp3 = (TH3D*)thisFile->Get (Form ("jetRtrkDist_dataSet%s_%s_mc_stat", gammaJetSampleId.Data (), algo.Data ()));
        jetRtrkHists[iAlgo][iPer][1][1]->Add (temp3);
        jetRtrkHists[iAlgo][2][1][1]->Add (temp3);

        TH2D* temp2 = (TH2D*)thisFile->Get (Form ("jetRtrkCounts_dataSet%s_%s_mc", gammaJetSampleId.Data (), algo.Data ()));
        jetRtrkCounts[iAlgo][iPer][1]->Add (temp2);
        jetRtrkCounts[iAlgo][2][1]->Add (temp2);
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

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH3D* temp3 = (TH3D*)thisFile->Get (Form ("jetRtrkDist_dataSet%s_%s_mc_stat", zeeJetSampleId.Data (), algo.Data ()));
        jetRtrkHists[iAlgo][iPer][1][1]->Add (temp3);
        jetRtrkHists[iAlgo][2][1][1]->Add (temp3);

        TH2D* temp2 = (TH2D*)thisFile->Get (Form ("jetRtrkCounts_dataSet%s_%s_mc", zeeJetSampleId.Data (), algo.Data ()));
        jetRtrkCounts[iAlgo][iPer][1]->Add (temp2);
        jetRtrkCounts[iAlgo][2][1]->Add (temp2);
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

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH3D* temp3 = (TH3D*)thisFile->Get (Form ("jetRtrkDist_dataSet%s_%s_mc_stat", zmumuJetSampleId.Data (), algo.Data ()));
        jetRtrkHists[iAlgo][iPer][1][1]->Add (temp3);
        jetRtrkHists[iAlgo][2][1][1]->Add (temp3);

        TH2D* temp2 = (TH2D*)thisFile->Get (Form ("jetRtrkCounts_dataSet%s_%s_mc", zmumuJetSampleId.Data (), algo.Data ()));
        jetRtrkCounts[iAlgo][iPer][1]->Add (temp2);
        jetRtrkCounts[iAlgo][2][1]->Add (temp2);
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

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH3D* temp3 = (TH3D*)thisFile->Get (Form ("jetRtrkDist_dataSet%s_%s_mc_stat", dijetSampleId.Data (), algo.Data ()));
        jetRtrkHists[iAlgo][iPer][1][1]->Add (temp3);
        jetRtrkHists[iAlgo][2][1][1]->Add (temp3);

        TH2D* temp2 = (TH2D*)thisFile->Get (Form ("jetRtrkCounts_dataSet%s_%s_mc", dijetSampleId.Data (), algo.Data ()));
        jetRtrkCounts[iAlgo][iPer][1]->Add (temp2);
        jetRtrkCounts[iAlgo][2][1]->Add (temp2);
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
  TLine* getalines[5] = {};
  TLine* xlines[5] = {};
  for (short i = 0; i < 5; i++) {
   const float dz = 0.1;
   const float dg = 0.05;
   const float dx = 0.2;

   zlines[i] = new TLine (pbins[0], 1.0-2*dz+dz*i, pbins[numpbins], 1.0-2*dz+dz*i);
   glines[i] = new TLine (pbins[0], 1.0-1*dg+dg*i, pbins[numpbins], 1.0-1*dg+dg*i);
   getalines[i] = new TLine (etabins[0], 1.0-1*dg+dg*i, etabins[numetabins], 1.0-1*dg+dg*i);
   xlines[i] = new TLine (xjrefbins[0], 1.0-2*dx+dx*i, xjrefbins[numxjrefbins], 1.0-2*dx+dx*i);

   if (1.0-2*dz+dz*i == 1) zlines[i]->SetLineStyle (1);
   else zlines[i]->SetLineStyle (3);
   if (1.0-1*dg+dg*i == 1) glines[i]->SetLineStyle (1);
   else glines[i]->SetLineStyle (3);
   if (1.0-1*dg+dg*i == 1) getalines[i]->SetLineStyle (1);
   else getalines[i]->SetLineStyle (3);
   if (1.0-2*dx+dx*i == 1) xlines[i]->SetLineStyle (1);
   else xlines[i]->SetLineStyle (3);
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
  TH2D *proj2d, *proj2d_mc, *proj2d_lo, *proj2d_hi;
  TGraphAsymmErrors *jetRtrkGraph_sys, *jetRtrkGraph_rat_sys;

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
   TString period = "Period A";
   if (iPer == 1) period = "Period B";
   else if (iPer == 2) period = "Period A+B";

   TString perType = "pA";
   if (iPer == 1) perType = "pB";
   else if (iPer == 2) perType = "pAB";

   for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta
    const short eta_lo = (iEta != numetabins ? iEta+1 : 5);
    const short eta_hi = (iEta != numetabins ? iEta+1 : 10);

    for (short iAlgo = 0; iAlgo < 2; iAlgo++) { // loop over jet algorithms
     const Style_t markerStyle = kFullDotLarge;
     topPad->cd ();
     topPad->SetLogx ();

     proj2d = Project2D ("", jetRtrkHists[iAlgo][iPer][0][1], "x", "z", eta_lo, eta_hi);
     jetRtrkHist = GetProfileX ("jetRtrk_Hist", proj2d, numpbins, pbins, true);
     jetRtrkHist->SetYTitle ("< #Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{jet}>");
     jetRtrkHist->SetAxisRange (0., 2.0, "Y");
     jetRtrkHist->SetMarkerStyle (markerStyle);
     jetRtrkHist->SetMarkerColor (dataColor);
     jetRtrkHist->SetLineColor (dataColor);
     jetRtrkHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     proj2d_lo = Project2D ("", jetRtrkHists[iAlgo][iPer][0][0], "x", "z", eta_lo, eta_hi);
     jetRtrkHist_lo = GetProfileX ("jetRtrk_Hist_lo", proj2d_lo, numpbins, pbins, true);

     proj2d_hi = Project2D ("", jetRtrkHists[iAlgo][iPer][0][2], "x", "z", eta_lo, eta_hi);
     jetRtrkHist_hi = GetProfileX ("jetRtrk_Hist_hi", proj2d_hi, numpbins, pbins, true);

     jetRtrkGraph_sys = new TGraphAsymmErrors (jetRtrkHist); // for plotting systematics
     CalcSystematics (jetRtrkGraph_sys, jetRtrkHist, jetRtrkHist_hi, jetRtrkHist_lo);
     if (jetRtrkHist_lo) delete jetRtrkHist_lo;
     if (jetRtrkHist_hi) delete jetRtrkHist_hi;

     jetRtrkGraph_sys->SetFillColor (kBlack);
     jetRtrkGraph_sys->SetFillStyle (3001);

     proj2d_mc = Project2D ("", jetRtrkHists[iAlgo][iPer][1][1], "x", "z", eta_lo, eta_hi);
     jetRtrkHist_mc = GetProfileX ("jetRtrk_Hist_mc", proj2d_mc, numpbins, pbins, true);
     jetRtrkHist_mc->SetMarkerStyle (markerStyle);
     jetRtrkHist_mc->SetMarkerColor (mcOverlayColor);
     jetRtrkHist_mc->SetLineColor (mcOverlayColor);

     jetRtrkHist->Draw ("e1 x0");
     jetRtrkHist_mc->Draw ("same e1 x0");
     jetRtrkGraph_sys->Draw ("2");

     const int nJetData = jetRtrkCounts[iAlgo][iPer][0]->Integral (1, numpbins, iEta, iEta);
     myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nJetData), 1.25, 0.04/uPadY);
     const int nJetMC = jetRtrkCounts[iAlgo][iPer][1]->Integral (1, numpbins, iEta, iEta);
     myMarkerText (0.175, 0.81, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", nJetMC), 1.25, 0.04/uPadY);
     if (eta_lo != 1 || eta_hi != numetabins)
      myText (0.155, 0.65, kBlack, Form ("%g < #eta_{det}^{Jet} < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04/uPadY);
     myText (0.155, 0.73, kBlack, period.Data (), 0.04/uPadY);

     bottomPad->cd ();
     bottomPad->SetLogx ();

     jetRtrkHist_rat = GetDataOverMC (TString ("jetRtrk_DataMCRatio"), proj2d, proj2d_mc, numpbins, pbins, false, "x");
     jetRtrkHist_rat_lo = GetDataOverMC (TString ("jetRtrk_DataMCRatio_lo"), proj2d_lo, proj2d_mc, numpbins, pbins, false, "x");
     jetRtrkHist_rat_hi = GetDataOverMC (TString ("jetRtrk_DataMCRatio_hi"), proj2d_hi, proj2d_mc, numpbins, pbins, false, "x");

     jetRtrkGraph_rat_sys = new TGraphAsymmErrors (jetRtrkHist_rat);
     CalcSystematics (jetRtrkGraph_rat_sys, jetRtrkHist_rat, jetRtrkHist_rat_hi, jetRtrkHist_rat_lo);
     if (jetRtrkHist_rat_lo) delete jetRtrkHist_rat_lo;
     if (jetRtrkHist_rat_hi) delete jetRtrkHist_rat_hi;

     jetRtrkGraph_rat_sys->SetFillColor (kBlack);
     jetRtrkGraph_rat_sys->SetFillStyle (3001);

     jetRtrkHist_rat->SetXTitle ("#it{p}_{T}^{Jet} #left[GeV#right]");
     jetRtrkHist_rat->SetYTitle ("Data / MC");
     jetRtrkHist_rat->SetAxisRange (0.91, 1.09, "Y");
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

     jetRtrkHist_rat->Draw ("e1 x0");
     jetRtrkGraph_rat_sys->Draw ("2");
     for (TLine* line : zlines) line->Draw ("same");

     char* plotName;
     if (iEta < numetabins) plotName = Form ("jet_rtrk_iEta%i.pdf", iEta);
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
     if (proj2d) { delete proj2d; proj2d = NULL; }
     if (proj2d_lo) { delete proj2d_lo; proj2d_lo = NULL; }
     if (proj2d_hi) { delete proj2d_hi; proj2d_hi = NULL; }
     if (proj2d_mc) { delete proj2d_mc; proj2d_mc = NULL; }

     if (jetRtrkHist) { delete jetRtrkHist; jetRtrkHist = NULL; }
     if (jetRtrkHist_mc) { delete jetRtrkHist_mc; jetRtrkHist_mc = NULL; }
     if (jetRtrkGraph_sys) { delete jetRtrkGraph_sys; jetRtrkGraph_sys = NULL; }
     if (jetRtrkHist_rat) { delete jetRtrkHist_rat; jetRtrkHist_rat = NULL; }
     if (jetRtrkGraph_rat_sys) { delete jetRtrkGraph_rat_sys; jetRtrkGraph_rat_sys = NULL; }
    } // end loop over jet algorithms
   } // end loop over eta bins

   for (short iP = 0; iP <= numpbins; iP++) { // loop over bins in p
    const short p_lo = (iP != numpbins ? iP+1 : 7);
    const short p_hi = (iP != numpbins ? iP+1 : 10);

    for (short iAlgo = 0; iAlgo < 2; iAlgo++) { // loop over jet algorithms
     const Style_t markerStyle = kFullDotLarge;
     topPad->cd ();
     topPad->SetLogx (false);

     proj2d = Project2D ("", jetRtrkHists[iAlgo][iPer][0][1], "y", "z", p_lo, p_hi);
     jetRtrkHist = GetProfileX ("jetRtrk_Hist", proj2d, numetabins, etabins, true);
     jetRtrkHist->SetYTitle ("< #Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{jet}>");
     jetRtrkHist->SetAxisRange (0., 2.0, "Y");
     jetRtrkHist->SetMarkerStyle (markerStyle);
     jetRtrkHist->SetMarkerColor (dataColor);
     jetRtrkHist->SetLineColor (dataColor);
     jetRtrkHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetTitleOffset (uPadY);

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     proj2d_lo = Project2D ("", jetRtrkHists[iAlgo][iPer][0][0], "y", "z", p_lo, p_hi);
     jetRtrkHist_lo = GetProfileX ("jetRtrk_Hist_lo", proj2d_lo, numetabins, pbins, true);

     proj2d_hi = Project2D ("", jetRtrkHists[iAlgo][iPer][0][2], "y", "z", p_lo, p_hi);
     jetRtrkHist_hi = GetProfileX ("jetRtrk_Hist_hi", proj2d_hi, numetabins, pbins, true);

     jetRtrkGraph_sys = new TGraphAsymmErrors (jetRtrkHist); // for plotting systematics
     CalcSystematics (jetRtrkGraph_sys, jetRtrkHist, jetRtrkHist_hi, jetRtrkHist_lo);
     if (jetRtrkHist_lo) delete jetRtrkHist_lo;
     if (jetRtrkHist_hi) delete jetRtrkHist_hi;

     jetRtrkGraph_sys->SetFillColor (kBlack);
     jetRtrkGraph_sys->SetFillStyle (3001);

     proj2d_mc = Project2D ("", jetRtrkHists[iAlgo][iPer][1][1], "y", "z", p_lo, p_hi);
     jetRtrkHist_mc = GetProfileX ("jetRtrk_Hist_mc", proj2d_mc, numetabins, etabins, true);
     jetRtrkHist_mc->SetMarkerStyle (markerStyle);
     jetRtrkHist_mc->SetMarkerColor (mcOverlayColor);
     jetRtrkHist_mc->SetLineColor (mcOverlayColor);

     jetRtrkHist->Draw ("e1 x0");
     jetRtrkHist_mc->Draw ("same e1 x0");
     jetRtrkGraph_sys->Draw ("2");

     const int nJetData = jetRtrkCounts[iAlgo][iPer][0]->Integral (p_lo, p_hi, 1, numetabins);
     myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nJetData), 1.25, 0.04/uPadY);
     const int nJetMC = jetRtrkCounts[iAlgo][iPer][1]->Integral (p_lo, p_hi, 1, numetabins);
     myMarkerText (0.175, 0.81, mcOverlayColor, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", nJetMC), 1.25, 0.04/uPadY);
     if (p_lo != 1 || p_hi != numpbins)
      myText (0.155, 0.65, kBlack, Form ("%g < #it{p}_{T}^{Jet} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
     myText (0.155, 0.73, kBlack, period.Data (), 0.04/uPadY);

     bottomPad->cd ();
     bottomPad->SetLogx (false);

     jetRtrkHist_rat = GetDataOverMC (TString ("jetRtrk_DataMCRatio"), proj2d, proj2d_mc, numetabins, etabins, false, "x");
     jetRtrkHist_rat_lo = GetDataOverMC (TString ("jetRtrk_DataMCRatio_lo"), proj2d_lo, proj2d_mc, numetabins, etabins, false, "x");
     jetRtrkHist_rat_hi = GetDataOverMC (TString ("jetRtrk_DataMCRatio_hi"), proj2d_hi, proj2d_mc, numetabins, etabins, false, "x");

     jetRtrkGraph_rat_sys = new TGraphAsymmErrors (jetRtrkHist_rat);
     CalcSystematics (jetRtrkGraph_rat_sys, jetRtrkHist_rat, jetRtrkHist_rat_hi, jetRtrkHist_rat_lo);
     if (jetRtrkHist_rat_lo) delete jetRtrkHist_rat_lo;
     if (jetRtrkHist_rat_hi) delete jetRtrkHist_rat_hi;

     jetRtrkGraph_rat_sys->SetFillColor (kBlack);
     jetRtrkGraph_rat_sys->SetFillStyle (3001);

     jetRtrkHist_rat->SetXTitle ("#eta_{det}^{Jet}");
     jetRtrkHist_rat->SetYTitle ("Data / MC");
     jetRtrkHist_rat->SetAxisRange (0.91, 1.09, "Y");
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

     jetRtrkHist_rat->Draw ("e1 x0");
     jetRtrkGraph_rat_sys->Draw ("2");
     for (TLine* line : getalines) line->Draw ("same");

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
     if (proj2d) { delete proj2d; proj2d = NULL; }
     if (proj2d_lo) { delete proj2d_lo; proj2d_lo = NULL; }
     if (proj2d_hi) { delete proj2d_hi; proj2d_hi = NULL; }
     if (proj2d_mc) { delete proj2d_mc; proj2d_mc = NULL; }

     if (jetRtrkHist) { delete jetRtrkHist; jetRtrkHist = NULL; }
     if (jetRtrkHist_mc) { delete jetRtrkHist_mc; jetRtrkHist_mc = NULL; }
     if (jetRtrkGraph_sys) { delete jetRtrkGraph_sys; jetRtrkGraph_sys = NULL; }
     if (jetRtrkHist_rat) { delete jetRtrkHist_rat; jetRtrkHist_rat = NULL; }
     if (jetRtrkGraph_rat_sys) { delete jetRtrkGraph_rat_sys; jetRtrkGraph_rat_sys = NULL; }
    } // end loop over jet algos
   } // end loop over pT bins


   /**** Plots rtrk distributions, binned by eta^jet ****/
   for (short iEta = 0; iEta < numetabins; iEta++) {
    const int p_lo = 7;
    const int p_hi =  10;

    const Style_t dataStyle = kFullCircle;
    const Style_t mcStyle = 33;
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) { // loop over jet algorithms
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

     topPad->cd ();
     topPad->SetLogx (0);

     proj2d = Project2D ("", jetRtrkHists[iAlgo][iPer][0][1], "y", "z", p_lo, p_hi);
     jetRtrkHist = proj2d->ProjectionY ("jetRtrkHist", iEta, iEta);

     jetRtrkHist->Rebin (rebinFactor);
     if (jetRtrkHist->Integral () != 0) jetRtrkHist->Scale (1./jetRtrkHist->Integral ());
     jetRtrkHist->SetXTitle ("r_{trk} = #Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{Jet}");
     jetRtrkHist->SetYTitle ("Counts / Total");
     jetRtrkHist->SetMarkerStyle (dataStyle);
     jetRtrkHist->SetMarkerColor (dataColor);
     jetRtrkHist->SetLineColor (dataColor);
     //jetRtrkHist->GetXaxis ()->SetRangeUser (0., 2.);
     jetRtrkHist->GetYaxis ()->SetRangeUser (0., 0.6);//jetRtrkHist->GetYaxis ()->GetXmax ());

     jetRtrkHist->GetXaxis ()->SetLabelSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetLabelSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetTitleSize (0.04/uPadY);
     jetRtrkHist->GetYaxis ()->SetTitleOffset (1.1*uPadY);

     proj2d_mc = Project2D ("", jetRtrkHists[iAlgo][iPer][1][1], "y", "z", p_lo, p_hi);
     jetRtrkHist_mc = proj2d_mc->ProjectionY ("vJetProjection_mc", iEta, iEta);

     jetRtrkHist_mc->Rebin (rebinFactor);
     if (jetRtrkHist_mc->Integral () != 0) jetRtrkHist_mc->Scale (1./jetRtrkHist_mc->Integral ()); 

     jetRtrkHist_mc->SetMarkerStyle (mcStyle);
     jetRtrkHist_mc->SetMarkerColor (mcOverlayColor);
     jetRtrkHist_mc->SetLineColor (mcOverlayColor);

     jetRtrkHist->DrawCopy ("e1 x0");
     jetRtrkHist_mc->DrawCopy ("same e1 x0"); // insitu factors are not applied to MC

     float mean, mean_err, mean_mc, mean_mc_err;
  
     if (useGaussian) {
      TF1* gaus_data = new TF1 ("gaus_data", "gaus (0)", 0, 4.0);
      jetRtrkHist->Fit (gaus_data, "Q0R");
      TF1* gaus_mc = new TF1 ("gaus_mc", "gaus (0)", 0, 4.0);
      jetRtrkHist_mc->Fit (gaus_mc, "Q0R");
      mean = gaus_data->GetParameter (1);
      mean_err = gaus_data->GetParError (1);
      mean_mc = gaus_mc->GetParameter (1);
      mean_mc_err = gaus_mc->GetParError (1);
      if (gaus_data) delete gaus_data;
      if (gaus_mc) delete gaus_mc;
     }
     else {
      mean = jetRtrkHist->GetMean ();
      mean_err = jetRtrkHist->GetMeanError ();
      mean_mc = jetRtrkHist_mc->GetMean ();
      mean_mc_err = jetRtrkHist_mc->GetMeanError ();
     }

     const int countData = jetRtrkCounts[iAlgo][iPer][0]->Integral (p_lo, p_hi, iEta, iEta);
     const int countMC = jetRtrkCounts[iAlgo][iPer][1]->Integral (p_lo, p_hi, iEta, iEta);

     myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", countData), 1.25, 0.04/uPadY);
     myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV (%i events)", countMC), 1.25, 0.04/uPadY);

     myText (0.155, 0.73, dataColor, Form ("<#it{r}_{trk}>^{data} = %.2f #pm %.2f", mean, mean_err), 0.04/uPadY);
     myText (0.155, 0.64, dataColor, Form ("<#it{r}_{trk}>^{MC} = %.2f #pm %.2f", mean_mc, mean_mc_err), 0.04/uPadY);

     myText (0.68, 0.34, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
     myText (0.68, 0.25, dataColor, period.Data (), 0.04/uPadY);
     myText (0.68, 0.16, dataColor, Form ("%g < #it{p}_{T}^{Jet} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04/uPadY);
     myText (0.68, 0.08, dataColor, Form ("%g < #eta_{det}^{Jet} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);

     bottomPad->cd ();
     bottomPad->SetLogx (false);
     jetRtrkHist->Divide (jetRtrkHist_mc);

     jetRtrkHist->SetYTitle ("Data / MC");
     jetRtrkHist->SetAxisRange (0.45, 1.65, "Y");
     jetRtrkHist->GetYaxis ()->SetNdivisions (605);
     jetRtrkHist->GetXaxis ()->SetTitleSize (0.04/dPadY);
     jetRtrkHist->GetYaxis ()->SetTitleSize (0.04/dPadY);
     jetRtrkHist->GetXaxis ()->SetTitleOffset (1);
     jetRtrkHist->GetYaxis ()->SetTitleOffset (1.1*dPadY);
     jetRtrkHist->GetYaxis ()->CenterTitle (true);
     jetRtrkHist->GetXaxis ()->SetLabelSize (0.04/dPadY);
     jetRtrkHist->GetYaxis ()->SetLabelSize (0.04/dPadY);
     jetRtrkHist->GetXaxis ()->SetTickLength (0.08);

     jetRtrkHist->DrawCopy ("e1 x0"); 
     for (TLine* line : xlines) line->Draw ();

     if (jetRtrkHist) { delete jetRtrkHist; jetRtrkHist = NULL; }
     if (jetRtrkHist_mc) { delete jetRtrkHist_mc; jetRtrkHist_mc = NULL; }

     if (proj2d) { delete proj2d; proj2d = NULL; }
     if (proj2d_mc) { delete proj2d_mc; proj2d_mc = NULL; }
      // end loop over insitu configurations

     char* plotName = Form ("rtrk_dists/%s_jets_rtrk_iEta%i.pdf", algo.Data (), iEta);
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
    } // end loop over jet algos 
   } // end loop over eta bins
  } // end loop over periods

  return;
}

} // end namespace
