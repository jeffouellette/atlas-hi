#include "MCValidatorHist.h"
#include "Params.h"
#include "CalibUtils.h"

#include <Utilities.h>
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


void MCValidatorHist () {

  const int numetabins_mc = 26;
  const double etabins_mc[27] = {-4.5, -4.0, -3.6, -3.2, -2.8, -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.8, 3.2, 3.6, 4.0, 4.5};

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("MCValidator/", "JetCalibration/");

  vector<TString> dijetOverlaySampleIds (0);
  dijetOverlaySampleIds.push_back ("pPb_Overlay_Dijet_Slice2");
  dijetOverlaySampleIds.push_back ("Pbp_Overlay_Dijet_Slice2");


  /**** Initialize histograms ****/
  TH3D* jetCalibRespHist = new TH3D ("jetCalibRespHist", "", numpbins, pbins, numetabins, etabins, numclosurebins, closurebins);
  jetCalibRespHist->Sumw2 ();
  TH3D* jetRecoRespHist = new TH3D ("jetRecoRespHist", "", numpbins, pbins, numetabins, etabins, numclosurebins, closurebins);
  jetRecoRespHist->Sumw2 ();

  TH2D* jetCalibRespCounts = new TH2D ("jetCalibRespCounts", "", numpbins, pbins, numetabins, etabins);
  jetCalibRespCounts->Sumw2();
  TH2D* jetRecoRespCounts = new TH2D ("jetRecoRespCounts", "", numpbins, pbins, numetabins, etabins);
  jetRecoRespCounts->Sumw2();

  TH3D** jetSamplingHist = Get1DArray <TH3D*> (28);
  for (short iJS = 0; iJS < 28; iJS++) {
    jetSamplingHist[iJS] = new TH3D (Form ("jetSamplingHist_iJS%i", iJS), "", numetabins_mc, etabins_mc, numphibins, phibins, numpbins, pbins);
    jetSamplingHist[iJS]->Sumw2 ();
  }

  TH2D* jetEtaPhiCorr = new TH2D ("jetEtaPhiCorr", "", numetabins_mc, etabins_mc, 48, -pi, pi);
  jetEtaPhiCorr->Sumw2 ();

  TH2D* jetPtSpectrum = new TH2D ("jetPtSpectrum", "", numpbins, pbins, numetabins, etabins);
  jetPtSpectrum->Sumw2 ();

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
        if (debugStatements) cout << "Status: In MCValidatorHist.C: Found " << fname.Data () << endl;

        // do this if dijet MC sample (OVERLAY)
        for (TString dijetOverlaySampleId : dijetOverlaySampleIds) { // check for gamma jet MC
          if (fname.Contains (dijetOverlaySampleId)) { // if gamma jet MC sample
            numFiles++;
            cout << "Reading in " << rootPath+fname << endl;
            TFile* thisFile = new TFile (rootPath + fname, "READ");

            jetCalibRespHist->Add ( (TH3D*)thisFile->Get (Form ("jetCalibResp_dataSet%s", dijetOverlaySampleId.Data ())));
            jetRecoRespHist->Add ( (TH3D*)thisFile->Get (Form ("jetRecoResp_dataSet%s", dijetOverlaySampleId.Data ())));
            jetCalibRespCounts->Add ( (TH2D*)thisFile->Get (Form ("jetCalibCounts_dataSet%s", dijetOverlaySampleId.Data ())));
            jetRecoRespCounts->Add ( (TH2D*)thisFile->Get (Form ("jetRecoCounts_dataSet%s", dijetOverlaySampleId.Data ())));

            for (short iJS = 0; iJS < 28; iJS++) {
              jetSamplingHist[iJS]->Add ( (TH3D*)thisFile->Get (Form ("jetSamplingHist_iJS%i_dataSet%s", iJS, dijetOverlaySampleId.Data ())));
            }

            jetEtaPhiCorr->Add ( (TH2D*)thisFile->Get (Form ("jetEtaPhiCorr_dataSet%s", dijetOverlaySampleId.Data ())));

            jetPtSpectrum->Add ( (TH2D*)thisFile->Get (Form ("jetPtSpectrum_dataSet%s", dijetOverlaySampleId.Data ())));

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

  TH1D** jetCalibEnergyScale = Get1DArray <TH1D*> (numetabins+1);
  TH1D** jetCalibEnergyRes = Get1DArray <TH1D*> (numetabins+1);
  TH1D** jetRecoEnergyScale = Get1DArray <TH1D*> (numetabins+1);
  TH1D** jetRecoEnergyRes = Get1DArray <TH1D*> (numetabins+1);

  for (int iEta = 0; iEta <= numetabins; iEta++) {
    const int eta_lo = (iEta != numetabins ? iEta+1 : 1);
    const int eta_hi = (iEta != numetabins ? iEta+1 : numetabins);

    TCanvas* jesCanvas = new TCanvas (Form ("jetEnergyScalePad_%i", iEta), "", 2400, 1800);

    jetCalibEnergyScale[iEta] = new TH1D (Form ("jetCalibEnergyScale_eta%i", iEta), "", numpbins, pbins);
    jetCalibEnergyRes[iEta] = new TH1D (Form ("jetCalibEnergyRes_eta%i", iEta), "", numpbins, pbins);
    jetRecoEnergyScale[iEta] = new TH1D (Form ("jetRecoEnergyScale_eta%i", iEta), "", numpbins, pbins);
    jetRecoEnergyRes[iEta] = new TH1D (Form ("jetRecoEnergyRes_eta%i", iEta), "", numpbins, pbins);

    jesCanvas->cd ();
    jesCanvas->Divide (4, 3);

    TF1** recoFits = Get1DArray<TF1*> (12);
    TF1** calibFits = Get1DArray<TF1*> (12);
    for (int iP = 0; iP <= numpbins; iP++) {
      const int p_lo = (iP != numpbins ? iP+1 : 1);
      const int p_hi = (iP != numpbins ? iP+1 : numpbins);

      if (iP <= 11)
        jesCanvas->cd (iP+1);

      proj2d = Project2D ("", jetRecoRespHist, "x", "z", eta_lo, eta_hi);
      thisHist = proj2d->ProjectionY ("_py", p_lo, p_hi);
      if (thisHist->Integral () != 0)
        thisHist->Scale (1./thisHist->Integral ());

      int nJets = jetRecoRespCounts->Integral (p_lo, p_hi, eta_lo, eta_hi);

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

      jetRecoEnergyScale[iEta]->SetBinContent (iP+1, recoFit->GetParameter (1));
      jetRecoEnergyScale[iEta]->SetBinError (iP+1, recoFit->GetParError (1));
      jetRecoEnergyRes[iEta]->SetBinContent (iP+1, recoFit->GetParameter (2));
      jetRecoEnergyRes[iEta]->SetBinError (iP+1, recoFit->GetParError (2));

      if (proj2d) { delete proj2d; proj2d = NULL; }

      proj2d = Project2D ("", jetCalibRespHist, "x", "z", eta_lo, eta_hi);
      thisHist = proj2d->ProjectionY ("_py", p_lo, p_hi);
      if (thisHist->Integral () != 0) {
        thisHist->Scale (1./thisHist->Integral ());
      }

      nJets = jetCalibRespCounts->Integral (p_lo, p_hi, eta_lo, eta_hi);

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

      jetCalibEnergyScale[iEta]->SetBinContent (iP+1, calibFit->GetParameter (1));
      jetCalibEnergyScale[iEta]->SetBinError (iP+1, calibFit->GetParError (1));
      jetCalibEnergyRes[iEta]->SetBinContent (iP+1, calibFit->GetParameter (2));
      jetCalibEnergyRes[iEta]->SetBinError (iP+1, calibFit->GetParError (2));

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

    jesCanvas->SaveAs (Form ("%s/jetEnergyResponse/iEta%i.pdf", plotPath.Data (), iEta));

    Delete1DArray (recoFits, 12);
    Delete1DArray (calibFits, 12);
    
    if (jesCanvas) delete jesCanvas;
  } // end loop over eta bins
  



  TCanvas* jesSummaryCanvas = new TCanvas ("jesSummaryCanvas", "", 800, 600);

  jesSummaryCanvas->cd ();
  gPad->SetLogx ();
  gStyle->SetErrorX (0.5);

  const double delta = 0.005;

  for (short iJet = 0; iJet < 2; iJet++) {

    TGraphAsymmErrors** hArr = Get1DArray <TGraphAsymmErrors*> (numetabins+1);
    for (int iEta = 0; iEta <= numetabins; iEta++) {
      thisHist = (iJet == 0 ? jetCalibEnergyScale[iEta] : jetRecoEnergyScale[iEta]);

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

      if (calcPtClosure)
        thisGraph->GetXaxis ()->SetRangeUser (50, 150);
      else
        thisGraph->GetXaxis ()->SetRangeUser (50, pbins[numpbins]);

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

    myText (0.20, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
    const int counts = (int)((iJet == 0 ? jetCalibRespCounts:jetRecoRespCounts)->Integral ());
    myText (0.20, 0.78, kBlack, Form ("%i jets", counts), 0.04);
    myText (0.20, 0.90, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);

    jesSummaryCanvas->SaveAs (Form ("%s/%s_JetEnergyScale.pdf", plotPath.Data (), iJet==0 ? "calib":"reco"));
    Delete1DArray (hArr, numetabins+1);

    hArr = Get1DArray <TGraphAsymmErrors*> (numetabins+1);
    for (int iEta = 0; iEta <= numetabins; iEta++) {
      thisHist = (iJet == 0 ? jetCalibEnergyRes[iEta] : jetRecoEnergyRes[iEta]);

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

      if (calcPtClosure)
        thisGraph->GetXaxis ()->SetRangeUser (50, 150);
      else
        thisGraph->GetXaxis ()->SetRangeUser (50, pbins[numpbins]);

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

    myText (0.20, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
    myText (0.20, 0.76, kBlack, Form ("%i jets", counts), 0.04);
    myText (0.20, 0.90, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);

    jesSummaryCanvas->SaveAs (Form ("%s/%s_JetEnergyResolution.pdf", plotPath.Data (), iJet==0 ? "calib":"reco"));
    Delete1DArray (hArr, numetabins+1);
  }

  if (jesSummaryCanvas) delete jesSummaryCanvas;

  TCanvas* ptSpectrumCanvas = new TCanvas ("ptSpectrumCanvas", "", 800, 600);
  gPad->SetLogx ();
  gPad->SetLogy ();
  jetPtSpectrum->Scale (1, "width");
  TGraphAsymmErrors** hArr = Get1DArray <TGraphAsymmErrors*> (numetabins);
  for (short iEta = 0; iEta < numetabins; iEta++) {
    TH1D* thisHist = jetPtSpectrum->ProjectionX (Form ("_px%i", iEta+1), iEta+1, iEta+1);

    thisGraph = make_graph (thisHist);
    hArr[iEta] = thisGraph;

    thisGraph->GetYaxis ()->SetRangeUser (0.01, 100);
    thisGraph->SetLineColor (colors[iEta]);
    thisGraph->SetMarkerColor (colors[iEta]);
    thisGraph->SetMarkerStyle (kFullDotLarge);
    thisGraph->SetMarkerSize (0.5);
    thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{Jet} #left[GeV#right]");
    thisGraph->GetXaxis ()->SetTitleOffset (1);
    thisGraph->GetYaxis ()->SetTitle ("d^{2}N / d#it{p}_{T} d#eta #left[GeV^{-1}#right]");
    thisGraph->GetYaxis ()->SetTitleOffset (1);

    if (iEta == 0)
      thisGraph->Draw ("ap");
    else
      thisGraph->Draw ("p");
    if (thisHist) delete thisHist;

    myMarkerText (0.62+ (iEta>= (numetabins/2))*0.18, 0.9-0.04*iEta+0.04* (numetabins/2)* (iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
  }
  const int counts = jetCalibRespCounts->Integral();
  myText (0.20, 0.26, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
  myText (0.20, 0.20, kBlack, Form ("%i jets", counts), 0.04);
  myText (0.20, 0.32, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
  ptSpectrumCanvas->SaveAs (Form ("%s/JetPtSpectra.pdf", plotPath.Data ()));
  Delete1DArray (hArr, numetabins);

  TCanvas* etaPhiCanvas = new TCanvas ("etaPhiCanvas", "", 800, 800);
  FormatTH2Canvas (etaPhiCanvas, true);
  gStyle->SetPalette (kRainBow);
  //jetEtaPhiCorr->RebinX (7);
  jetEtaPhiCorr->RebinY (2);
  jetEtaPhiCorr->GetXaxis ()->SetTitle ("#eta^{Jet}_{det}");
  jetEtaPhiCorr->GetXaxis ()->SetTitleOffset (1);
  jetEtaPhiCorr->GetYaxis ()->SetTitle ("#phi^{Jet}_{det}");
  jetEtaPhiCorr->GetYaxis ()->SetTitleOffset (1);
  jetEtaPhiCorr->GetZaxis ()->SetTitle ("Counts");
  jetEtaPhiCorr->GetZaxis ()->SetTitleOffset (1.2);
  jetEtaPhiCorr->Draw ("colz");

  TLine** lines = Get1DArray <TLine*> (4);
  const double etabounds[5] = {1.5, 1.5, 3.2, 3.2, 1.5};
  const double phibounds[5] = {-pi, -pi/2., -pi/2., -pi, -pi};
  for (short iLine = 0; iLine < 4; iLine++) {
    lines[iLine] = new TLine (etabounds[iLine], phibounds[iLine], etabounds[iLine+1], phibounds[iLine+1]);
    lines[iLine]->SetLineStyle (7);
    lines[iLine]->SetLineWidth (1);
    lines[iLine]->Draw ("same");
  }
  etaPhiCanvas->SaveAs (Form ("%s/JetEtaPhi.pdf", plotPath.Data ()));

  for (short iJS = 0; iJS < 28; iJS++) {

    TCanvas* jetSamplingCanvas = new TCanvas (Form ("jetSamplingCanvas_iJS%i", iJS), "", 800, 800);
    FormatTH2Canvas (jetSamplingCanvas, true);
    jetSamplingCanvas->cd ();
    
    proj2d = Project2D (TString (jetSamplingHist[iJS]->GetName ()) + "_pxy", jetSamplingHist[iJS], "x", "y", 1, jetSamplingHist[iJS]->GetNbinsZ ());
    //proj2d->Scale (1. / jetSamplingHist[iJS]->GetNbinsZ ());
    proj2d->GetXaxis ()->SetTitle ("#eta^{Jet}_{det}");
    proj2d->GetYaxis ()->SetTitle ("#phi^{Jet}_{det}");
    proj2d->GetZaxis ()->SetTitle ("Counts");
    proj2d->GetXaxis ()->SetTitleOffset (1);
    proj2d->GetYaxis ()->SetTitleOffset (1);
    proj2d->Draw ("colz");

    for (short iLine = 0; iLine < 4; iLine++) {
      lines[iLine]->Draw ("same");
    }

    jetSamplingCanvas->SaveAs (Form ("%s/Sampling/JS_%i.pdf", plotPath.Data (), iJS));

    if (jetSamplingCanvas) delete jetSamplingCanvas;
    if (proj2d) { delete proj2d; proj2d = NULL; }
  }

  Delete1DArray (jetCalibEnergyScale, numetabins+1);
  Delete1DArray (jetRecoEnergyScale, numetabins+1);
  Delete1DArray (jetCalibEnergyRes, numetabins+1);
  Delete1DArray (jetRecoEnergyRes, numetabins+1);

  if (jetCalibRespHist) delete jetCalibRespHist;
  if (jetRecoRespHist) delete jetRecoRespHist;
  if (jetCalibRespCounts) delete jetCalibRespCounts;
  if (jetRecoRespCounts) delete jetRecoRespCounts;

  if (jetPtSpectrum) delete jetPtSpectrum;
  if (jetEtaPhiCorr) delete jetEtaPhiCorr;

  return;
}

} // end namespace
