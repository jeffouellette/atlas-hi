#include "IntercalibrationHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utils.h>
#include <ArrayTemplates.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

using namespace atlashi;

namespace JetCalibration {

TString GetMCType (const short iMC) {
  if (iMC == 0) return "gammajets_overlay";
  else if (iMC == 1) return "gammajets_signal";
  else if (iMC == 2) return "dijets_signal";
  else if (iMC == 3) return "dijets_overlay_valid";
  else return "";
}


void IntercalibrationHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("Intercalibration/", "JetCalibration/");

  // Setup list of data and lists of MC samples
  vector<TString> gammaJetOverlaySampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetOverlaySampleIds.push_back (TString ("Pbp_Overlay_GammaJet_Slice") + to_string (i+1));
   gammaJetOverlaySampleIds.push_back (TString ("pPb_Overlay_GammaJet_Slice") + to_string (i+1));
  }

  vector<TString> gammaJetSignalSampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetSignalSampleIds.push_back (TString ("pPb_Signal_GammaJet_Slice") + to_string (i+1));
  }

  vector<TString> dijetSignalSampleIds (0);
  dijetSignalSampleIds.push_back ("pPb_Signal_Dijet_Slice2");

  vector<TString> dijetOverlaySampleIds (0);
  dijetOverlaySampleIds.push_back ("pPb_Overlay_Dijet_Slice2");
  dijetOverlaySampleIds.push_back ("Pbp_Overlay_Dijet_Slice2");


  /**** Initialize histograms ****/
  TH3D** jetCalibRespHists = Get1DArray <TH3D*> (4); // MC type
  TH3D** jetRecoRespHists = Get1DArray <TH3D*> (4);
  TH2D** jetCalibRespCounts = Get1DArray <TH2D*> (4);
  TH2D** jetRecoRespCounts = Get1DArray <TH2D*> (4);

  for (short iMC = 0; iMC < 4; iMC++) {
   const TString mcType = GetMCType (iMC);

   jetCalibRespHists[iMC] = new TH3D (Form ("jetCalibRespHists_%s", mcType.Data ()), "", numpbins, pbins, numetabins, etabins, numclosurebins, closurebins);
   jetCalibRespHists[iMC]->Sumw2 ();
   jetRecoRespHists[iMC] = new TH3D (Form ("jetRecoRespHists_%s", mcType.Data ()), "", numpbins, pbins, numetabins, etabins, numclosurebins, closurebins);
   jetRecoRespHists[iMC]->Sumw2 ();

   jetCalibRespCounts[iMC] = new TH2D (Form ("jetCalibRespCounts_%s", mcType.Data ()), "", numpbins, pbins, numetabins, etabins);
   jetCalibRespCounts[iMC]->Sumw2();
   jetRecoRespCounts[iMC] = new TH2D (Form ("jetRecoRespCounts_%s", mcType.Data ()), "", numpbins, pbins, numetabins, etabins);
   jetRecoRespCounts[iMC]->Sumw2();
  }

  TH1D* truthPhotonsPerEvent = new TH1D ("truthPhotonsPerEvent", "", 9, -0.5, 8.5);
  TH1D* truthRatePerEvent = new TH1D ("truthRatePerEvent", "", 50, 0, 1);
  truthPhotonsPerEvent->Sumw2 ();
  truthRatePerEvent->Sumw2();

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
   int numFiles = 0;
   while ( (sysfile= (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
     if (debugStatements) cout << "Status: In IntercalibrationHist.C: Found " << fname.Data () << endl;

     // do this if gamma jet MC sample (OVERLAY)
     for (TString gammaJetOverlaySampleId : gammaJetOverlaySampleIds) { // check for gamma jet MC
      if (fname.Contains (gammaJetOverlaySampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");

       jetCalibRespHists[0]->Add ( (TH3D*)thisFile->Get (Form ("jetCalibResp_dataSet%s", gammaJetOverlaySampleId.Data ())));
       jetRecoRespHists[0]->Add ( (TH3D*)thisFile->Get (Form ("jetRecoResp_dataSet%s", gammaJetOverlaySampleId.Data ())));
       jetCalibRespCounts[0]->Add ( (TH2D*)thisFile->Get (Form ("jetCalibCounts_dataSet%s", gammaJetOverlaySampleId.Data ())));
       jetRecoRespCounts[0]->Add ( (TH2D*)thisFile->Get (Form ("jetRecoCounts_dataSet%s", gammaJetOverlaySampleId.Data ())));

       truthPhotonsPerEvent->Add ( (TH1D*)thisFile->Get (Form ("truthPhotonsPerEvent_%s", gammaJetOverlaySampleId.Data ())));
       truthRatePerEvent->Add ( (TH1D*)thisFile->Get (Form ("truthRatePerEvent_%s", gammaJetOverlaySampleId.Data ())));

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if gamma jet MC sample (SIGNAL)
     for (TString gammaJetSignalSampleId : gammaJetSignalSampleIds) { // check for gamma jet MC
      if (fname.Contains (gammaJetSignalSampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");

       jetCalibRespHists[1]->Add ( (TH3D*)thisFile->Get (Form ("jetCalibResp_dataSet%s", gammaJetSignalSampleId.Data ())));
       jetRecoRespHists[1]->Add ( (TH3D*)thisFile->Get (Form ("jetRecoResp_dataSet%s", gammaJetSignalSampleId.Data ())));
       jetCalibRespCounts[1]->Add ( (TH2D*)thisFile->Get (Form ("jetCalibCounts_dataSet%s", gammaJetSignalSampleId.Data ())));
       jetRecoRespCounts[1]->Add ( (TH2D*)thisFile->Get (Form ("jetRecoCounts_dataSet%s", gammaJetSignalSampleId.Data ())));

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if dijet MC sample (SIGNAL)
     for (TString dijetSignalSampleId : dijetSignalSampleIds) { // check for gamma jet MC
      if (fname.Contains (dijetSignalSampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");

       jetCalibRespHists[2]->Add ( (TH3D*)thisFile->Get (Form ("jetCalibResp_dataSet%s", dijetSignalSampleId.Data ())));
       jetRecoRespHists[2]->Add ( (TH3D*)thisFile->Get (Form ("jetRecoResp_dataSet%s", dijetSignalSampleId.Data ())));
       jetCalibRespCounts[2]->Add ( (TH2D*)thisFile->Get (Form ("jetCalibCounts_dataSet%s", dijetSignalSampleId.Data ())));
       jetRecoRespCounts[2]->Add ( (TH2D*)thisFile->Get (Form ("jetRecoCounts_dataSet%s", dijetSignalSampleId.Data ())));

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     // do this if dijet MC sample (OVERLAY)
     for (TString dijetOverlaySampleId : dijetOverlaySampleIds) { // check for gamma jet MC
      if (fname.Contains (dijetOverlaySampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");

       jetCalibRespHists[3]->Add ( (TH3D*)thisFile->Get (Form ("jetCalibResp_dataSet%s", dijetOverlaySampleId.Data ())));
       jetRecoRespHists[3]->Add ( (TH3D*)thisFile->Get (Form ("jetRecoResp_dataSet%s", dijetOverlaySampleId.Data ())));
       jetCalibRespCounts[3]->Add ( (TH2D*)thisFile->Get (Form ("jetCalibCounts_dataSet%s", dijetOverlaySampleId.Data ())));
       jetRecoRespCounts[3]->Add ( (TH2D*)thisFile->Get (Form ("jetRecoCounts_dataSet%s", dijetOverlaySampleId.Data ())));

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



  TH1D* thisHist = NULL;
  TH2D* proj2d = NULL;
  TGraphAsymmErrors* thisGraph = NULL;
  Color_t colors[15] = {kGray+2, kAzure, kViolet, kMagenta, kPink+10, kPink, kOrange+10, kOrange, kSpring+10, kSpring, kTeal+10, kTeal, kAzure+10, kGray, kBlack};


  /**** Plots the jet energy response ****/

  TH1D*** jetCalibEnergyScale = Get2DArray <TH1D*> (4, numetabins+1);
  TH1D*** jetCalibEnergyRes = Get2DArray <TH1D*> (4, numetabins+1);
  TH1D*** jetRecoEnergyScale = Get2DArray <TH1D*> (4, numetabins+1);
  TH1D*** jetRecoEnergyRes = Get2DArray <TH1D*> (4, numetabins+1);

  TCanvas** jesCombCanvas = Get1DArray <TCanvas*> (numetabins);
  for (short iEta = 0; iEta < numetabins; iEta++) {
   jesCombCanvas[iEta] = new TCanvas (Form ("jetEnergyScaleCombinedCanvas_iEta%i", iEta), "", 800, 600);
  }

  for (short iMC = 0; iMC < 4; iMC++) {
   const TString mcType = GetMCType (iMC);

   for (int iEta = 0; iEta <= numetabins; iEta++) {
    const int eta_lo = (iEta != numetabins ? iEta+1 : 1);
    const int eta_hi = (iEta != numetabins ? iEta+1 : numetabins);

    TCanvas* jesCanvas = new TCanvas (Form ("jetEnergyScalePad_%i", iEta), "", 2400, 1800);

    jetCalibEnergyScale[iMC][iEta] = new TH1D (Form ("%s_jetCalibEnergyScale_eta%i", mcType.Data (), iEta), "", numpbins, pbins);
    jetCalibEnergyRes[iMC][iEta] = new TH1D (Form ("%s_jetCalibEnergyRes_eta%i", mcType.Data (), iEta), "", numpbins, pbins);
    jetRecoEnergyScale[iMC][iEta] = new TH1D (Form ("%s_jetRecoEnergyScale_eta%i", mcType.Data (), iEta), "", numpbins, pbins);
    jetRecoEnergyRes[iMC][iEta] = new TH1D (Form ("%s_jetRecoEnergyRes_eta%i", mcType.Data (), iEta), "", numpbins, pbins);

    jesCanvas->cd ();
    jesCanvas->Divide (4, 3);

    TF1** recoFits = Get1DArray<TF1*> (12);
    TF1** calibFits = Get1DArray<TF1*> (12);
    for (int iP = 0; iP <= numpbins; iP++) {
     const int p_lo = (iP != numpbins ? iP+1 : 1);
     const int p_hi = (iP != numpbins ? iP+1 : numpbins);

     if (iP <= 11)
      jesCanvas->cd (iP+1);

     proj2d = Project2D ("", jetRecoRespHists[iMC], "x", "z", eta_lo, eta_hi);
     thisHist = proj2d->ProjectionY ("_py", p_lo, p_hi);
     if (thisHist->Integral () != 0)
      thisHist->Scale (1./thisHist->Integral ());

     int nJets = jetRecoRespCounts[iMC]->Integral (p_lo, p_hi, eta_lo, eta_hi);

     TF1* recoFit = new TF1 (Form ("recoFit_iP%i", iP), "gaus (0)", 0, 2.0);
     if (nJets > 0) {
      recoFit->SetParameter (1, thisHist->GetMean ());
      recoFit->SetParameter (2, thisHist->GetStdDev ());
      thisHist->Fit (recoFit, "RN0Q");
      double m = recoFit->GetParameter (1);
      double s = recoFit->GetParameter (2);
      if (recoFit) delete recoFit;

      recoFit = new TF1 (Form ("recoFit2_iP%i", iP), "gaus (0)", m - 2.0*s, m + 2.0*s);
      recoFit->SetParameter (1, thisHist->GetMean ());
      recoFit->SetParameter (2, thisHist->GetStdDev ());
      thisHist->Fit (recoFit, "RN0Q");
     }

     if (iP <= 11) {
      thisHist->SetMarkerStyle (kFullDotMedium);
      thisHist->SetLineColor (kBlack);
      thisHist->SetMarkerColor (kBlack);
      thisHist->GetYaxis ()->SetTitle ("Counts / Total");
      thisHist->DrawCopy ("e1 x0");

      recoFit->SetLineColor (kBlack);
      recoFit->Draw ("same");
      recoFits[iP] = recoFit;
     }

     jetRecoEnergyScale[iMC][iEta]->SetBinContent (iP+1, recoFit->GetParameter (1));
     jetRecoEnergyScale[iMC][iEta]->SetBinError (iP+1, recoFit->GetParError (1));
     jetRecoEnergyRes[iMC][iEta]->SetBinContent (iP+1, recoFit->GetParameter (2));
     jetRecoEnergyRes[iMC][iEta]->SetBinError (iP+1, recoFit->GetParError (2));

     if (iP == 8 && iEta != numetabins && (iMC == 0 || iMC == 2)) {
      jesCombCanvas[iEta]->cd();
      thisHist->SetMarkerColor (iMC+1);
      thisHist->SetLineColor (iMC+1);
      recoFit->SetLineColor (iMC+1);
      if (iMC == 0) thisHist->DrawCopy ("e1 x0");
      else thisHist->DrawCopy ("same e1 x0");
      ( (TF1*)recoFit->Clone())->Draw ("same");
      jesCanvas->cd (iP+1);
     }
     if (proj2d) { delete proj2d; proj2d = NULL; }

     proj2d = Project2D ("", jetCalibRespHists[iMC], "x", "z", eta_lo, eta_hi);
     thisHist = proj2d->ProjectionY ("_py", p_lo, p_hi);
     if (thisHist->Integral () != 0) {
      thisHist->Scale (1./thisHist->Integral ());
     }

     nJets = jetCalibRespCounts[iMC]->Integral (p_lo, p_hi, eta_lo, eta_hi);

     TF1* calibFit = new TF1 (Form ("calibFit_iP%i", iP), "gaus (0)", 0, 2.0);
     if (nJets > 0) {
      calibFit->SetParameter (1, thisHist->GetMean ());
      calibFit->SetParameter (2, thisHist->GetStdDev ());
      thisHist->Fit (calibFit, "RN0Q");
      double m = calibFit->GetParameter (1);
      double s = calibFit->GetParameter (2);
      if (calibFit) delete calibFit;

      calibFit = new TF1 (Form ("calibFit2_iP%i", iP), "gaus (0)", m - 2.0*s, m + 2.0*s);
      calibFit->SetParameter (1, thisHist->GetMean ());
      calibFit->SetParameter (2, thisHist->GetStdDev ());
      thisHist->Fit (calibFit, "RN0Q");
     }

     if (iP <= 11) {
      thisHist->SetMarkerStyle (kFullDotMedium);
      thisHist->SetLineColor (kBlue);
      thisHist->SetMarkerColor (kBlue);
      thisHist->DrawCopy ("same e1 x0");

      calibFit->SetLineColor (kBlue);
      calibFit->Draw ("same");
      calibFits[iP] = calibFit;
     }

     jetCalibEnergyScale[iMC][iEta]->SetBinContent (iP+1, calibFit->GetParameter (1));
     jetCalibEnergyScale[iMC][iEta]->SetBinError (iP+1, calibFit->GetParError (1));
     jetCalibEnergyRes[iMC][iEta]->SetBinContent (iP+1, calibFit->GetParameter (2));
     jetCalibEnergyRes[iMC][iEta]->SetBinError (iP+1, calibFit->GetParError (2));

     if (iP == 8 && iEta != numetabins && (iMC == 0 || iMC == 2)) {
      jesCombCanvas[iEta]->cd();
      thisHist->SetMarkerColor (iMC+2);
      thisHist->SetLineColor (iMC+2);
      calibFit->SetLineColor (iMC+2);
      thisHist->DrawCopy ("same e1 x0");
      ( (TF1*)calibFit->Clone())->Draw ("same");
      jesCanvas->cd (iP+1);
     }
     if (proj2d) { delete proj2d; proj2d = NULL; }
     if (thisHist) { delete thisHist; thisHist = NULL; }


     if (iP <= 11) {
      myText (0.18, 0.90, kBlack, Form ("#mu = %s", FormatMeasurement (calibFit->GetParameter (1), calibFit->GetParError (1), 1)), 0.04 * 2);
      myText (0.18, 0.80, kBlack, Form ("#sigma = %s", FormatMeasurement (calibFit->GetParameter (2), calibFit->GetParError (2), 1)), 0.04 * 2);
      if (calcPtClosure)
       myText (0.18, 0.70, kBlack, Form ("%g < #it{p}_{T}^{true} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04 * 2);
      else
       myText (0.18, 0.70, kBlack, Form ("%g < #it{E}_{true} < %g", pbins[p_lo-1], pbins[p_hi]), 0.04 * 2);
     }
     else {
      if (recoFit) delete recoFit;
      if (calibFit) delete calibFit;
     }
    } // end loop over pT bins
    

    jesCanvas->cd (1);
    myText (0.6, 0.9, kBlack, Form ("%g < #eta < %g", etabins[eta_lo-1], etabins[eta_hi]), 0.04*2);

    jesCanvas->SaveAs (Form ("%s/jetEnergyResponse/%s_iEta%i.pdf", plotPath.Data (), mcType.Data (), iEta));

    Delete1DArray (recoFits, 12);
    Delete1DArray (calibFits, 12);
    
    if (jesCanvas) delete jesCanvas;
   } // end loop over eta bins
  }

  for (short iEta = 0; iEta < numetabins; iEta++) {
   jesCombCanvas[iEta]->cd ();
   if (calcPtClosure)
    myText (0.18, 0.90, kBlack, Form ("%g < #it{p}_{T}^{true} < %g", pbins[8], pbins[9]), 0.04);
   else
    myText (0.18, 0.90, kBlack, Form ("%g < #it{E}_{true} < %g", pbins[8], pbins[9]), 0.04);

   myText (0.18, 0.84, kBlack, Form ("%g < #eta_{Det}^{Jet} < %g", etabins[iEta], etabins[iEta+1]), 0.04);
   myMarkerText (0.6, 0.9, 1, kFullCircle, "#gamma-tagged EM scale", 1.25, 0.04);
   myMarkerText (0.6, 0.84, 2, kFullCircle, "#gamma-tagged EtaJES scale", 1.25, 0.04);
   myMarkerText (0.6, 0.78, 3, kFullCircle, "Inclusive EM scale", 1.25, 0.04);
   myMarkerText (0.6, 0.72, 4, kFullCircle, "Inclusive EtaJES scale", 1.25, 0.04);
   jesCombCanvas[iEta]->SaveAs (Form ("%s/jetEnergyResponse/combined_iEta%i.pdf", plotPath.Data (), iEta));
  }
  Delete1DArray (jesCombCanvas, numetabins);


  for (short iMC = 0; iMC < 4; iMC++) {
   const TString mcType = GetMCType (iMC);

   TCanvas* jesSummaryCanvas = new TCanvas ("jesSummaryCanvas", "", 800, 600);

   jesSummaryCanvas->cd ();
   gPad->SetLogx ();
   gStyle->SetErrorX (0.5);

   const double delta = 0.005;

   for (short iJet = 0; iJet < 2; iJet++) {

    TGraphAsymmErrors** hArr = Get1DArray <TGraphAsymmErrors*> (numetabins+1);
    for (int iEta = 0; iEta <= numetabins; iEta++) {
     thisHist = (iJet == 0 ? jetCalibEnergyScale[iMC][iEta] : jetRecoEnergyScale[iMC][iEta]);

     thisGraph = make_graph (thisHist);
     double t_delta = 1.;
     if (iEta < numetabins) t_delta += delta * (iEta - (numetabins/2) + (int) (iEta>= (numetabins/2)));
     deltaize (thisGraph, t_delta, true);

     if (calcPtClosure) {
      thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{true} #left[GeV#right]");
      thisGraph->GetYaxis ()->SetTitle (Form ("<#it{p}_{T}^{%s} / #it{p}_{T}^{true}>", iJet == 0 ? "calib":"reco"));
     }
     else {
      thisGraph->GetXaxis ()->SetTitle ("#it{E}_{true} #left[GeV#right]");
      thisGraph->GetYaxis ()->SetTitle (Form ("<#it{E}_{%s} / #it{E}_{true}>", iJet == 0 ? "calib":"reco"));
     }

     if (iMC == 2 || iMC == 3) {
      if (calcPtClosure)
       thisGraph->GetXaxis ()->SetRangeUser (50, 150);
      else
       thisGraph->GetXaxis ()->SetRangeUser (50, pbins[numpbins]);
     }
     else {
      thisGraph->GetXaxis ()->SetRangeUser (pbins[0], pbins[numpbins]);
     }

     if (iJet == 0)
      thisGraph->GetYaxis ()->SetRangeUser (0.85, 1.20);
     else 
      thisGraph->GetYaxis ()->SetRangeUser (0.35, 1.10);
     thisGraph->SetLineColor (colors[iEta]);
     thisGraph->SetMarkerColor (colors[iEta]);
     thisGraph->SetMarkerStyle (kFullDotLarge);
     if (iEta < numetabins)
      thisGraph->SetMarkerSize (0.5);
     else
      thisGraph->SetLineWidth (2);

     if (iEta == 0)
      thisGraph->Draw ("AP");
     else
      thisGraph->Draw ("P");
     hArr[iEta] = thisGraph;

     if (iEta == numetabins)
      myMarkerText (0.22, 0.49, colors[iEta], kFullCircle, "#eta-integrated", 1.25, 0.03);
     else
      myMarkerText (0.22+ (iEta>= (numetabins/2))*0.18, 0.45-0.04*iEta+0.04* (numetabins/2)* (iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
    }

    if (iMC == 0 || iMC == 3)
     myText (0.20, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
    else
     myText (0.20, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV", 0.04);
    const int counts = (int)((iJet == 0 ? jetCalibRespCounts[iMC]:jetRecoRespCounts[iMC])->Integral ());
    myText (0.20, 0.78, kBlack, Form ("%i jets", counts), 0.04);
    myText (0.20, 0.90, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);

    jesSummaryCanvas->SaveAs (Form ("%s/%s_%s_JetEnergyScale.pdf", plotPath.Data (), mcType.Data (), iJet==0 ? "calib":"reco"));
    Delete1DArray (hArr, numetabins+1);

    hArr = Get1DArray <TGraphAsymmErrors*> (numetabins+1);
    for (int iEta = 0; iEta <= numetabins; iEta++) {
     thisHist = (iJet == 0 ? jetCalibEnergyRes[iMC][iEta] : jetRecoEnergyRes[iMC][iEta]);

     thisGraph = make_graph (thisHist);
     double t_delta = 1;
     if (iEta < numetabins) t_delta += delta * (iEta - (numetabins/2) + (int) (iEta>= (numetabins/2)));
     deltaize (thisGraph, t_delta, true);

     if (calcPtClosure) {
      thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{true} #left[GeV#right]");
      thisGraph->GetYaxis ()->SetTitle (Form ("#sigma#left[<#it{p}_{T}^{%s} / #it{p}_{T}^{true}>#right]", iJet == 0 ? "calib":"reco"));
     }
     else {
      thisGraph->GetXaxis ()->SetTitle ("#it{E}_{true} #left[GeV#right]");
      thisGraph->GetYaxis ()->SetTitle (Form ("#sigma#left[<#it{E}_{%s} / #it{E}_{true}>#right]", iJet == 0 ? "calib":"reco"));
     }

     if (iMC == 2 || iMC == 3) {
      if (calcPtClosure)
       thisGraph->GetXaxis ()->SetRangeUser (50, 150);
      else
       thisGraph->GetXaxis ()->SetRangeUser (50, pbins[numpbins]);
     }
     else {
      thisGraph->GetXaxis ()->SetRangeUser (pbins[0], pbins[numpbins]);
     }

     thisGraph->GetYaxis ()->SetRangeUser (0, 0.3);

     thisGraph->SetLineColor (colors[iEta]);
     thisGraph->SetMarkerColor (colors[iEta]);
     thisGraph->SetMarkerStyle (kFullDotLarge);
     if (iEta < numetabins)
      thisGraph->SetMarkerSize (0.5);
     else
      thisGraph->SetLineWidth (2);

     if (iEta == 0)
      thisGraph->Draw ("AP");
     else
      thisGraph->Draw ("P");
     hArr[iEta] = thisGraph;

     if (iEta == numetabins)
      myMarkerText (0.57, 0.89, colors[iEta], kFullCircle, "#eta-integrated", 1.25, 0.03);
     else
      myMarkerText (0.57+ (iEta>= (numetabins/2))*0.18, 0.85-0.04*iEta+0.04* (numetabins/2)* (iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
    }

    if (iMC == 0 || iMC == 3)
     myText (0.20, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
    else
     myText (0.20, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV", 0.04);
    myText (0.20, 0.76, kBlack, Form ("%i jets", counts), 0.04);
    myText (0.20, 0.90, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);

    jesSummaryCanvas->SaveAs (Form ("%s/%s_%s_JetEnergyResolution.pdf", plotPath.Data (), mcType.Data (), iJet==0 ? "calib":"reco"));
    Delete1DArray (hArr, numetabins+1);
   }

   if (jesSummaryCanvas) delete jesSummaryCanvas;
  }

  Delete2DArray (jetCalibEnergyScale, 4, numetabins+1);
  Delete2DArray (jetRecoEnergyScale, 4, numetabins+1);
  Delete2DArray (jetCalibEnergyRes, 4, numetabins+1);
  Delete2DArray (jetRecoEnergyRes, 4, numetabins+1);

  Delete1DArray (jetCalibRespHists, 4);
  Delete1DArray (jetRecoRespHists, 4);
  Delete1DArray (jetCalibRespCounts, 4);
  Delete1DArray (jetRecoRespCounts, 4);

  TCanvas* truthMatchedPhotonCanvas = new TCanvas ("truthMatchedPhotonCanvas", "", 800, 600);
  truthMatchedPhotonCanvas->cd();
  truthPhotonsPerEvent->Draw ("hist");
  truthMatchedPhotonCanvas->SaveAs (Form ("%s/truthPhotonsPerEvent.pdf", plotPath.Data ()));

  truthRatePerEvent->Draw ("hist");
  truthMatchedPhotonCanvas->SaveAs (Form ("%s/truthRatePerEvent.pdf", plotPath.Data ()));

  return;
}

}
