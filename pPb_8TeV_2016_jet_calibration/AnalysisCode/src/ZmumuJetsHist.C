#include "ZmumuJetsHist.h"
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


void ZmumuJetsHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("ZmumuHists/", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);

  vector<TString> zmumuJetSampleIds (0);
  zmumuJetSampleIds.push_back ("Pbp_Overlay_ZmumuJet");
  zmumuJetSampleIds.push_back ("pPb_Overlay_ZmumuJet");

  TH3D* zmumuJetHists[2][3][3][3];
  //TH2D* zmumuJetHistsSys[3][numetabins][3][3];

  for (short iYear = 0; iYear < 2; iYear++) {
   const TString year = (iYear == 0 ? "2015" : "2016");

   for (short iData = 0; iData < 3; iData++) { // iData is 0 for data, 1 for MC
    TString dataType = "data";
    if (iData == 1) dataType = "mc_overlay";
    else if (iData == 2) dataType = "mc_signalonly";

    for (short iErr = 0; iErr < 3; iErr++) {
     TString error = "sys_lo";
     if (iErr == 1) error = "stat";
     else if (iErr == 2) error = "sys_hi";

     for (short iPer = 0; iPer < 3; iPer++) {
      TString period = "periodA";
      if (iPer == 1) period = "periodB";
      else if (iPer == 2) period = "periodAB";

      zmumuJetHists[iYear][iPer][iData][iErr] = new TH3D (Form ("zmumuJetPtRatio_%s_%s_%s_%s", year.Data (), dataType.Data (), error.Data (), period.Data ()), "", numpzbins, pzbins, numetabins, etabins, numxjrefbins, xjrefbins);
      zmumuJetHists[iYear][iPer][iData][iErr]->Sumw2 ();
     }
    }
   }
  }

  int***** nZmumuJet = Get5DArray <int> (2, 3, 2, numetabins+1, numpzbins+1);

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
   TVectorD *nZmumuJetVec;
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
       //infoVec = (TVectorD*)thisFile->Get (Form ("infoVec_%i", runNumber));
       nZmumuJetVec = (TVectorD*)thisFile->Get (Form ("nZmumuJetVec_%i", runNumber));

       for (short iErr = 0; iErr < 3; iErr++) {
        TString error = "sys_lo";
        if (iErr == 1) error = "stat";
        else if (iErr == 2) error = "sys_hi";

        TH3D* temp = (TH3D*)thisFile->Get (Form ("zmumuJetPtRatio_dataSet%i_data_%s", runNumber, error.Data ()));
        zmumuJetHists[iYear][iPer][0][iErr]->Add (temp);
        zmumuJetHists[iYear][2][0][iErr]->Add (temp);
       }

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        //const bool flipEta = runNumber < 313500 && iEta < numetabins;
        const bool flipEta = false;
        const short act_iEta = (flipEta ? (numetabins - iEta - 1) : iEta);

        for (short iP = 0; iP <= numpzbins; iP++) {
         nZmumuJet[iYear][iPer][0][iEta][iP] += (*nZmumuJetVec)[iEta + iP* (numetabins+1)];
         nZmumuJet[iYear][2][0][iEta][iP] += (*nZmumuJetVec)[act_iEta + iP* (numetabins+1)];
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
       nZmumuJetVec = (TVectorD*)thisFile->Get (Form ("nZmumuJetVec_%s", zmumuJetSampleId.Data ()));

       TH3D* temp = (TH3D*)thisFile->Get (Form ("zmumuJetPtRatio_dataSet%s_mc_overlay_stat", zmumuJetSampleId.Data ()));
       zmumuJetHists[iYear][iPer][1][1]->Add (temp);
       zmumuJetHists[iYear][2][1][1]->Add (temp);

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        //const bool flipEta = zmumuJetSampleId.Contains ("pPb") && iEta < numetabins;
        const bool flipEta = false;
        const short act_iEta = (flipEta ? numetabins - iEta - 1 : iEta); // period A condition

        for (short iP = 0; iP <= numpzbins; iP++) {
         nZmumuJet[iYear][iPer][1][iEta][iP] += (*nZmumuJetVec)[iEta + iP* (numetabins+1)];
         nZmumuJet[iYear][2][1][iEta][iP] += (*nZmumuJetVec)[act_iEta + iP* (numetabins+1)];

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
  TLine* zetalines[5] = {};
  for (short i = 0; i < 5; i++) {
   const float dz = 0.05;

   zlines[i] = new TLine (pzbins[0], 1.0-2*dz+dz*i, pzbins[numpzbins], 1.0-2*dz+dz*i);
   zetalines[i] = new TLine (etabins[0], 1.0-2*dz+dz*i, etabins[numetabins], 1.0-2*dz+dz*i);

   if (1.0-2*dz+dz*i == 1) zlines[i]->SetLineStyle (1);
   else zlines[i]->SetLineStyle (3);
   if (1.0-2*dz+dz*i == 1) zetalines[i]->SetLineStyle (1);
   else zetalines[i]->SetLineStyle (3);
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
  TH1D *vJetHist = NULL, *vJetHist_mc_overlay = NULL, *vJetHist_lo = NULL, *vJetHist_hi = NULL, *vJetHist_rat = NULL, *vJetHist_rat_lo = NULL, *vJetHist_rat_hi = NULL;
  TH2D *proj = NULL, *proj_mc_overlay = NULL, *proj_lo = NULL, *proj_hi = NULL;
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
     eta_lo = 1;
     eta_hi = numetabins;
    }

    /**** Plot ZmumuJet info ****/
    {
     const Style_t dataStyle = 20;
     const Style_t mcOverlayStyle = 33;

     topPad->cd ();
     topPad->SetLogx ();

     proj = Project2D ("", zmumuJetHists[1][iPer][0][1], "x", "z", eta_lo, eta_hi);
     vJetHist = GetProfileX ("vJetHist", proj, numpzbins, pzbins, false);

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
     proj_lo = Project2D ("", zmumuJetHists[1][iPer][0][0], "x", "z", eta_lo, eta_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj_lo, numpzbins, pzbins, false);

     proj_hi = Project2D ("", zmumuJetHists[1][iPer][0][2], "x", "z", eta_lo, eta_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numpzbins, pzbins, false);

     vJetGraph_sys = new TGraphAsymmErrors (vJetHist); // for plotting systematics
     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }

     vJetGraph_sys->SetFillColor (dataColor);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", zmumuJetHists[1][iPer][1][1], "x", "z", eta_lo, eta_hi);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numpzbins, pzbins, false);

     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

     vJetHist->DrawCopy ("e1 x0");
     vJetHist_mc_overlay->DrawCopy ("same e1 x0"); // insitu factors are not applied to MC
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", nZmumuJet[1][iPer][0][iEta][numpzbins]), 1.25, 0.04/uPadY);
     myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", nZmumuJet[1][iPer][1][iEta][numpzbins]), 1.25, 0.04/uPadY);
     myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
     if (iEta < numetabins) myText (0.155, 0.15, dataColor, Form ("Z (#mu#mu) + Jet, %g < #eta_{det}^{Jet} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
     else myText (0.155, 0.15, dataColor, "Z (#mu#mu) + Jet", 0.04/uPadY);
     myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);

     bottomPad->cd ();
     bottomPad->SetLogx ();

     vJetHist_rat = GetDataOverMC (TString (Form ("zmumuJetPtDataMCRatio_iEta%i", iEta)), proj, proj_mc_overlay, numpzbins, pzbins, false, "x");
     vJetHist_rat_lo = GetDataOverMC (TString (Form ("zmumuJetPtDataMCRatio_lo_iEta%i", iEta)), proj_lo, proj_mc_overlay, numpzbins, pzbins, false, "x");
     vJetHist_rat_hi = GetDataOverMC (TString (Form ("zmumuJetPtDataMCRatio_hi_iEta%i", iEta)), proj_hi, proj_mc_overlay, numpzbins, pzbins, false, "x");

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

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

     vJetHist_rat->DrawCopy ("e1 x0");
     ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     for (TLine* line : zlines) line->Draw ();

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }
     

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
    }
   }

   /**** Now loop over pT bins and plot response as function of eta^jet ****/
   for (short iP = 0; iP <= numpzbins; iP++) {

    int p_lo = iP+1;
    int p_hi = iP+1;
    if (iP == numpzbins) {
     p_lo = 1;
     p_hi = numpzbins;
    }

    /**** Plot ZmumuJet info ****/
    {
     const Style_t dataStyle = 20;
     const Style_t mcOverlayStyle = 33;
     //const Style_t mcSignalStyle = 34;

     topPad->cd ();
     topPad->SetLogx (false);

     proj = Project2D ("", zmumuJetHists[1][iPer][0][1], "y", "z", p_lo, p_hi);
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
     proj_lo = Project2D ("", zmumuJetHists[1][iPer][0][0], "y", "z", p_lo, p_hi);
     vJetHist_lo = GetProfileX ("vJetHist_lo", proj, numetabins, etabins, false);

     proj_hi = Project2D ("", zmumuJetHists[1][iPer][0][2], "y", "z", p_lo, p_hi);
     vJetHist_hi = GetProfileX ("vJetHist_hi", proj_hi, numetabins, etabins, false);

     CalcSystematics (vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) { delete vJetHist_lo; vJetHist_lo = NULL; }
     if (vJetHist_hi) { delete vJetHist_hi; vJetHist_hi = NULL; }
     vJetGraph_sys->SetFillColor (dataColor);
     vJetGraph_sys->SetFillStyle (3001);

     proj_mc_overlay = Project2D ("", zmumuJetHists[1][iPer][1][1], "y", "z", p_lo, p_hi);
     vJetHist_mc_overlay = GetProfileX ("vJetHist_mc_overlay", proj_mc_overlay, numetabins, etabins, false);

     vJetHist_mc_overlay->SetMarkerStyle (mcOverlayStyle);
     vJetHist_mc_overlay->SetMarkerColor (mcOverlayColor);
     vJetHist_mc_overlay->SetLineColor (mcOverlayColor);

     vJetHist->DrawCopy ("e1 x0");
     vJetHist_mc_overlay->DrawCopy ("same e1 x0"); // insitu factors are not applied to MC
     ( (TGraphAsymmErrors*)vJetGraph_sys->Clone ())->Draw ("2");

     myMarkerText (0.175, 0.88, dataColor, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV with Insitu Corrections (%i events)", nZmumuJet[1][iPer][0][numetabins][iP]), 1.25, 0.04/uPadY);
     myMarkerText (0.175, 0.81, mcOverlayColor, kFullDiamond, Form ("Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay (%i events)", nZmumuJet[1][iPer][1][numetabins][iP]), 1.25, 0.04/uPadY);
     myText (0.155, 0.22, dataColor, "#bf{#it{ATLAS}} Internal", 0.04/uPadY);
     if (iP < numpzbins) myText (0.155, 0.15, dataColor, Form ("Z (#mu#mu) + Jet, %g < #it{p}_{T}^{#mu#mu} < %g", pbins[iP], pbins[iP+1]), 0.04/uPadY);
     else myText (0.155, 0.15, dataColor, "Z (#mu#mu) + Jet", 0.04/uPadY);
     myText (0.155, 0.08, dataColor, period.Data (), 0.04/uPadY);

     bottomPad->cd ();
     bottomPad->SetLogx (false);

     vJetHist_rat = GetDataOverMC (TString (Form ("zmumuJetPtDataMCRatio_iP%i", iP)), proj, proj_mc_overlay, numetabins, etabins, false, "x");
     vJetHist_rat_lo = GetDataOverMC (TString (Form ("zmumuJetPtDataMCRatio_lo_iP%i", iP)), proj_lo, proj_mc_overlay, numetabins, etabins, false, "x");
     vJetHist_rat_hi = GetDataOverMC (TString (Form ("zmumuJetPtDataMCRatio_hi_iP%i", iP)), proj_hi, proj_mc_overlay, numetabins, etabins, false, "x");

     if (proj) { delete proj; proj = NULL; }
     if (proj_mc_overlay) { delete proj_mc_overlay; proj_mc_overlay = NULL; }
     if (proj_lo) { delete proj_lo; proj_lo = NULL; }
     if (proj_hi) { delete proj_hi; proj_hi = NULL; }

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

     vJetHist_rat->DrawCopy ("e1 x0");
     ( (TGraphAsymmErrors*)vJetGraph_rat_sys->Clone ())->Draw ("2");
     for (TLine* line : zetalines) line->Draw ();

     if (vJetHist) { delete vJetHist; vJetHist = NULL; }
     if (vJetHist_mc_overlay) { delete vJetHist_mc_overlay; vJetHist_mc_overlay = NULL; }
     if (vJetGraph_sys) { delete vJetGraph_sys; vJetGraph_sys = NULL; }
     if (vJetHist_rat) { delete vJetHist_rat; vJetHist_rat = NULL; }
     if (vJetGraph_rat_sys) { delete vJetGraph_rat_sys; vJetGraph_rat_sys = NULL; }

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
    }

   }
  }

  return;
}

} // end namespace
