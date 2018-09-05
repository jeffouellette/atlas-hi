#include "EnergyScaleChecksHist.h"

#include <ArrayTemplates.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
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

void EnergyScaleChecksHist () {

  SetAtlasStyle();

  // Setup trigger vectors
  SetupDirectories("EnergyScaleChecks/", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<TString> gammaJetSampleIds(0);
  for (short i = 0; i < 6; i++) {
   gammaJetSampleIds.push_back(TString("Pbp") + (runValidation ? "_Valid":"_Overlay") + "_GammaJet_Slice" + to_string(i+1));
   gammaJetSampleIds.push_back(TString("pPb") + (runValidation ? "_Valid":"_Overlay") + "_GammaJet_Slice" + to_string(i+1));
  }

  vector<TString> zeeJetSampleIds(0);
  zeeJetSampleIds.push_back("Pbp_ZeeJet_Overlay");
  zeeJetSampleIds.push_back("pPb_ZeeJet_Overlay");

  vector<TString> zmumuJetSampleIds(0);
  zmumuJetSampleIds.push_back("Pbp_ZmumuJet");
  zmumuJetSampleIds.push_back("pPb_ZmumuJet");


  /**** Initialize histograms ****/
  TH1D* electronEnergyScale = new TH1D("electronEnergyScale", ";Electron #it{p}_{T}^{reco} / #it{p}_{T}^{truth};", 200, 0, 2.0);
  electronEnergyScale->Sumw2();

  TH2D* lowResponseEtaPhi[2];
  TH2D* highResponseEtaPhi[2];
  TH2D* lowResponseEtaPt[2];
  TH2D* highResponseEtaPt[2];
  TH2D* responsePt[2];
  TH1D* jetSpectrum[2];

  lowResponseEtaPhi[0] = new TH2D("akt4hi_lowResponseEtaPhi", ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
  lowResponseEtaPhi[1] = new TH2D("akt4emtopo_lowResponseEtaPhi", ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
  highResponseEtaPhi[0] = new TH2D("akt4hi_highResponseEtaPhi", ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
  highResponseEtaPhi[1] = new TH2D("akt4emtopo_highResponseEtaPhi", ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
  lowResponseEtaPt[0] = new TH2D("akt4hi_lowResponseEtaPt", ";#eta;#it{p}_{T}^{reco} #left[GeV#right];", 25, -2.5, 2.5, 66, 20, 350);
  lowResponseEtaPt[1] = new TH2D("akt4emtopo_lowResponseEtaPt", ";#eta;#it{p}_{T}^{reco} #left[GeV#right];", 25, -2.5, 2.5, 66, 20, 350);
  highResponseEtaPt[0] = new TH2D("akt4hi_highResponseEtaPt", ";#eta;#it{p}_{T}^{reco} #left[GeV#right];", 25, -2.5, 2.5, 66, 20, 350);
  highResponseEtaPt[1] = new TH2D("akt4emtopo_highResponseEtaPt", ";#eta;#it{p}_{T}^{reco} #left[GeV#right];", 25, -2.5, 2.5, 66, 20, 350);
  responsePt[0] = new TH2D("akt4hi_responsePt", ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};#it{p}_{T}^{reco} #left[GeV#right];", 200, 0, 4.0, 66, 20, 350);
  responsePt[1] = new TH2D("akt4emtopo_responsePt", ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};#it{p}_{T}^{reco} #left[GeV#right];", 200, 0, 4.0, 66, 20, 350);
  jetSpectrum[0] = new TH1D("akt4hi_jetSpectrum", ";#it{p}_{T}^{jet} #left[GeV#right];", 33, 20, 350);
  jetSpectrum[1] = new TH1D("akt4emtopo_jetSpectrum", ";#it{p}_{T}^{jet} #left[GeV#right];", 33, 20, 350);

  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   lowResponseEtaPhi[iAlgo]->Sumw2();
   highResponseEtaPhi[iAlgo]->Sumw2();
   lowResponseEtaPt[iAlgo]->Sumw2();
   highResponseEtaPt[iAlgo]->Sumw2();
   responsePt[iAlgo]->Sumw2();
   jetSpectrum[iAlgo]->Sumw2();
  }

  TH1D* jetEnergyResponseCalib[2][numpbins+1][numetabins+1];
  TH1D* jetEnergyResponseReco[2][numpbins+1][numetabins+1];
  TH1D* photonEnergyResponse[numpbins+1][numetabins+1];

  for (short iP = 0; iP <= numpbins; iP++) {
   for (short iEta = 0; iEta <= numetabins; iEta++) {
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     jetEnergyResponseCalib[iAlgo][iP][iEta] = new TH1D(Form("%s_jetEnergyResponseCalib_iP%i_iEta%i", algo.Data(), iP, iEta), "", 50, 0, 2);
     jetEnergyResponseCalib[iAlgo][iP][iEta]->Sumw2();
     jetEnergyResponseReco[iAlgo][iP][iEta] = new TH1D(Form("%s_jetEnergyResponseReco_iP%i_iEta%i", algo.Data(), iP, iEta), "", 50, 0, 2);
     jetEnergyResponseReco[iAlgo][iP][iEta]->Sumw2();
    }
    photonEnergyResponse[iP][iEta] = new TH1D(Form("photonEnergyResponse_iP%i_iEta%i", iP, iEta), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", 50, 0, 2);
    photonEnergyResponse[iP][iEta]->Sumw2();
   }
  }

  int *** nJet = Get3DArray <int> (2, numpbins+1, numetabins+1);
  int ** nGamma = Get2DArray <int> (numpbins+1, numetabins+1);

  {
   TSystemDirectory dir(rootPath.Data(), rootPath.Data());
   TList* sysfiles = dir.GetListOfFiles();
   if (!sysfiles) {
    cout << "Cannot get list of files! Exiting." << endl;
    return;
   }
   TSystemFile *sysfile;
   TString fname;
   TString histName;
   TIter next(sysfiles);
   TVectorD *nJetVec, *nGammaVec;
   int numFiles = 0;
   while ((sysfile=(TSystemFile*)next())) {
    fname = sysfile->GetName();
    if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
     if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data() << endl;

     // do this if gamma jet MC sample
     for (TString gammaJetSampleId : gammaJetSampleIds) { // check for gamma jet MC
      if (fname.Contains(gammaJetSampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       nJetVec = (TVectorD*)thisFile->Get(Form("nJetVec_%s", gammaJetSampleId.Data()));
       nGammaVec = (TVectorD*)thisFile->Get(Form("nGammaVec_%s", gammaJetSampleId.Data()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
        lowResponseEtaPhi[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_lowResponseEtaPhi_dataSet%s", algo.Data(), gammaJetSampleId.Data())));
        highResponseEtaPhi[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_highResponseEtaPhi_dataSet%s", algo.Data(), gammaJetSampleId.Data())));
        lowResponseEtaPt[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_lowResponseEtaPt_dataSet%s", algo.Data(), gammaJetSampleId.Data())));
        highResponseEtaPt[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_highResponseEtaPt_dataSet%s", algo.Data(), gammaJetSampleId.Data())));
        responsePt[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_responsePt_dataSet%s", algo.Data(), gammaJetSampleId.Data())));
        jetSpectrum[iAlgo]->Add((TH1D*)thisFile->Get(Form("%s_jetPt_dataSet%s", algo.Data(), gammaJetSampleId.Data())));
       }

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        for (short iP = 0; iP <= numpbins; iP++) {
         for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
          const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
          nJet[iAlgo][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo*(numpbins+1)*(numetabins+1)];
          jetEnergyResponseCalib[iAlgo][iP][iEta]->Add((TH1D*)thisFile->Get(Form("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data(), gammaJetSampleId.Data(), iP, iEta)));
          jetEnergyResponseReco[iAlgo][iP][iEta]->Add((TH1D*)thisFile->Get(Form("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data(), gammaJetSampleId.Data(), iP, iEta)));
         }
         nGamma[iP][iEta] += (*nGammaVec)[iP + (numpbins+1)*iEta];
         photonEnergyResponse[iP][iEta]->Add((TH1D*)thisFile->Get(Form("photonEnergyResponse_dataSet%s_iP%i_iEta%i", gammaJetSampleId.Data(), iP, iEta)));
        }
       }

       thisFile->Close();
       delete thisFile;
       break;
      }
     }
     // do this if Z->ee MC sample
     for (TString zeeJetSampleId : zeeJetSampleIds) { // check for Z->ee MC
      if (fname.Contains(zeeJetSampleId)) { // if Z->ee MC do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       nJetVec = (TVectorD*)thisFile->Get(Form("nJetVec_%s", zeeJetSampleId.Data()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
        lowResponseEtaPhi[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_lowResponseEtaPhi_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
        highResponseEtaPhi[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_highResponseEtaPhi_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
        lowResponseEtaPt[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_lowResponseEtaPt_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
        highResponseEtaPt[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_highResponseEtaPt_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
        responsePt[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_responsePt_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
        jetSpectrum[iAlgo]->Add((TH1D*)thisFile->Get(Form("%s_jetPt_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
       }

       electronEnergyScale->Add((TH1D*)thisFile->Get(Form("electronEnergyScale_dataSet%s", zeeJetSampleId.Data())));

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        for (short iP = 0; iP <= numpbins; iP++) {
         for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
          const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
          nJet[iAlgo][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo*(numpbins+1)*(numetabins+1)];
          jetEnergyResponseCalib[iAlgo][iP][iEta]->Add((TH1D*)thisFile->Get(Form("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data(), zeeJetSampleId.Data(), iP, iEta)));
          jetEnergyResponseReco[iAlgo][iP][iEta]->Add((TH1D*)thisFile->Get(Form("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data(), zeeJetSampleId.Data(), iP, iEta)));
         }
        }
       }

       thisFile->Close();
       delete thisFile;
       break;
      }
     }
     // do this if Z->mumu sample
     for (TString zmumuJetSampleId : zmumuJetSampleIds) { // check for Z->mumu MC
      if (fname.Contains(zmumuJetSampleId)) { // if Z->mumu sample do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       nJetVec = (TVectorD*)thisFile->Get(Form("nJetVec_%s", zmumuJetSampleId.Data()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
        lowResponseEtaPhi[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_lowResponseEtaPhi_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
        highResponseEtaPhi[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_highResponseEtaPhi_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
        lowResponseEtaPt[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_lowResponseEtaPt_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
        highResponseEtaPt[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_highResponseEtaPt_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
        responsePt[iAlgo]->Add((TH2D*)thisFile->Get(Form("%s_responsePt_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
        jetSpectrum[iAlgo]->Add((TH1D*)thisFile->Get(Form("%s_jetPt_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
       }

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        for (short iP = 0; iP <= numpbins; iP++) {
         for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
          const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
          nJet[iAlgo][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo*(numpbins+1)*(numetabins+1)];
          jetEnergyResponseCalib[iAlgo][iP][iEta]->Add((TH1D*)thisFile->Get(Form("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data(), zmumuJetSampleId.Data(), iP, iEta)));
          jetEnergyResponseReco[iAlgo][iP][iEta]->Add((TH1D*)thisFile->Get(Form("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data(), zmumuJetSampleId.Data(), iP, iEta)));
         }
        }
       }

       thisFile->Close();
       delete thisFile;
       break;
      }
     }
    }
   }
   cout << numFiles << " files read in." << endl;
  }
  /**** End loop over input files ****/


  /**** Plots the electron energy scale ****/
  TCanvas* energyScaleCanvas = new TCanvas("energyScaleCanvas", "", 800, 600);
  TH1D* thisHist = electronEnergyScale;
  TF1* gausFit = new TF1("gausFit", "gaus(0)", 0, 2.0);
  energyScaleCanvas->cd();
  const float ncounts = thisHist->Integral();
  thisHist->Scale(1./ncounts);
  thisHist->Fit(gausFit, "R", "L");
  const float m = gausFit->GetParameter(1);
  const float s = gausFit->GetParameter(2);
  if (gausFit) delete gausFit;
  gausFit = new TF1("gausFit2", "gaus(0)", m - 1.6*s, m + 1.6*s);
  thisHist->Fit(gausFit, "R", "L");

  thisHist->GetYaxis()->SetTitle ("Counts / Total");
  //thisHist->SetMarkerStyle(6);
  thisHist->Draw("e1 x0");
  myText(0.18, 0.85, kBlack, "Electron energy response");
  myText(0.18, 0.78, kBlack, Form("%i electrons", (int)ncounts));
  myText(0.18, 0.71, kBlack, Form("En. Scale = %.5f #pm %.5f", gausFit->GetParameter(1), gausFit->GetParError(1)));
  myText(0.18, 0.64, kBlack, Form("En. Res. = %.5f #pm %.5f", gausFit->GetParameter(2), gausFit->GetParError(2)));
  energyScaleCanvas->SaveAs(Form("%s/electronEnergyScale.pdf", plotPath.Data()));


  /**** Plots the jet energy response ****/
  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

   TH1D* jetEnergyScale[numetabins];
   TH1D* jetEnergyRes[numetabins];

   for (int iEta = 0; iEta <= numetabins; iEta++) {
    TCanvas* jesCanvas = new TCanvas (Form ("%s_jetEnergyScalePad_%i", algo.Data(), iEta), "", 800, 600);

    if (iEta < numetabins) {
     jetEnergyScale[iEta] = new TH1D (Form ("%s_jetEnergyScale_eta%i", algo.Data(), iEta), ";#it{p}_{T} #left[GeV#right];JES", numpbins, pbins);
     jetEnergyRes[iEta] = new TH1D (Form ("%s_jetEnergyRes_eta%i", algo.Data(), iEta), ";#it{p}_{T} #left[GeV#right];JER", numpbins, pbins);
    }

    jesCanvas->cd();
    jesCanvas->Divide(4, 3);

    TF1* recoFits[12];
    TF1* calibFits[12];
    for (int iP = 0; iP < numpbins; iP++) {
     if (iP <= 11) {
      jesCanvas->cd(iP+1);
     }
     TH1D* thisHist = jetEnergyResponseReco[iAlgo][iP][iEta];
     int nJets = nJet[iAlgo][iP][iEta];

     TF1* recoFit = new TF1 (Form ("recoFit_iP%i", iP), "gaus(0)", 0, 2.0);
     thisHist->Scale (1./thisHist->Integral());
     thisHist->Fit (recoFit, "RNL");
     double m = recoFit->GetParameter(1);
     double s = recoFit->GetParameter(2);
     if (recoFit) delete recoFit;

     recoFit = new TF1 (Form ("recoFit2_iP%i", iP), "gaus(0)", m - 1.3*s, m + 1.3*s);
     thisHist->Fit (recoFit, "RNL");
     //m = recoFit->GetParameter(1);
     //s = recoFit->GetParameter(2);
     //if (recoFit) delete recoFit;

     //recoFit = new TF1("recoFit3", "gaus(0)", m - 1.3*s, m + 1.3*s);
     //thisHist->Fit(recoFit, "RNL");

     thisHist->SetMarkerStyle(kFullDotMedium);
     thisHist->SetLineColor(kBlack);
     thisHist->SetMarkerColor(kBlack);
     thisHist->GetYaxis()->SetTitle ("Counts / Total");

     if (iP <= 11) {
      thisHist->Draw("e1 x0");
      recoFit->SetLineColor(kBlack);
      recoFit->Draw("same");
      recoFits[iP] = recoFit;
     }

     thisHist = jetEnergyResponseCalib[iAlgo][iP][iEta];
     thisHist->Scale (1./thisHist->Integral());

     TF1* calibFit = new TF1 (Form ("calibFit_iP%i", iP), "gaus(0)", 0, 2.0);
     thisHist->Fit(calibFit, "RNL");
     m = calibFit->GetParameter(1);
     s = calibFit->GetParameter(2);
     if (calibFit) delete calibFit;

     calibFit = new TF1 (Form ("calibFit2_iP%i", iP), "gaus(0)", m - 1.3*s, m + 1.3*s);
     thisHist->Fit(calibFit, "RNL");
     //m = calibFit->GetParameter(1);
     //s = calibFit->GetParameter(2);
     //if (calibFit) delete calibFit;

     //calibFit = new TF1("calibFit3", "gaus(0)", m - 1.3*s, m + 1.3*s);
     //thisHist->Fit(calibFit, "RNL");

     thisHist->SetMarkerStyle(kFullDotMedium);
     thisHist->SetLineColor(kBlue);
     thisHist->SetMarkerColor(kBlue);
     if (iP <= 11) {
      thisHist->Draw("same e1 x0");
      calibFit->SetLineColor(kBlue);
      calibFit->Draw("same");
      calibFits[iP] = calibFit;
     }

     if (iEta < numetabins) {
      jetEnergyScale[iEta]->SetBinContent (iP+1, calibFit->GetParameter(1));
      jetEnergyScale[iEta]->SetBinError (iP+1, calibFit->GetParError(1));
      //jetEnergyScale[iEta]->SetBinError (iP+1, 0.000001);
      jetEnergyRes[iEta]->SetBinContent (iP+1, calibFit->GetParameter(2));
      jetEnergyRes[iEta]->SetBinError (iP+1, calibFit->GetParError(2));
      //jetEnergyRes[iEta]->SetBinError (iP+1, 0.000001);
     }

     myText(0.18, 0.90, kBlack, Form("#mu = %s", FormatMeasurement (calibFit->GetParameter(1), calibFit->GetParError(1), 1)), 0.04 * 2);
     myText(0.18, 0.80, kBlack, Form("#sigma = %s", FormatMeasurement (calibFit->GetParameter(2), calibFit->GetParError(2), 1)), 0.04 * 2);

     //if (iEta < numetabins) {
     // myText (0.18, 0.62, kBlack, Form("%g < #eta_{Lab}^{Jet} < %g", etabins[iEta], etabins[iEta+1]));
     //}
     if (iP <= 11) {
      myText (0.18, 0.70, kBlack, Form("%g < #it{p}_{T}^{J} < %g", pbins[iP], pbins[iP+1]), 0.04 * 2);
     }

    }

    jesCanvas->cd(1);
    myText (0.2, 0.8, kBlack, Form("%g < #eta_{Lab} < %g", etabins[iEta], etabins[iEta+1]), 0.03*3);
    if (iAlgo == 0)
     myText (0.2, 0.65, kBlack, "HI jets");
    else
     myText (0.2, 0.65, kBlack, "EMTopo jets");

    jesCanvas->SaveAs(Form("%s/jetEnergyResponse/%s_iEta%i.pdf", plotPath.Data(), algo.Data(), iEta));
    if (jesCanvas) delete jesCanvas;

    for (short iP = 0; iP <= 11; iP++) {
     if (recoFits[iP]) delete recoFits[iP];
     if (calibFits[iP]) delete calibFits[iP];
    }
   }

   TCanvas* jesSummaryCanvas = new TCanvas("jesSummaryCanvas", "", 800, 600);
   jesSummaryCanvas->cd();
   gPad->SetLogx();
   Color_t colors[14] = {kAzure-2, kBlue, kViolet+7, kMagenta, kPink, kRed+2, kOrange+8, kSpring, kGreen+2, kCyan, kGray, kRed-9, kViolet-8, kBlue-9};
   for (int iEta = 0; iEta < numetabins; iEta++) {
     jetEnergyScale[iEta]->GetYaxis()->SetTitle ("Jet Energy Scale");
     jetEnergyScale[iEta]->GetYaxis()->SetRangeUser (0.85, 1.15);
     jetEnergyScale[iEta]->SetLineColor (colors[iEta]);
     jetEnergyScale[iEta]->SetMarkerColor (colors[iEta]);
     //jetEnergyScale[iEta]->SetMarkerStyle (kFullDotLarge);
     if (iEta == 0) jetEnergyScale[iEta]->Draw("hist ][");
     else jetEnergyScale[iEta]->Draw("same hist ][");
     myMarkerText (0.77, 0.88-0.04*iEta, colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
   }
   jesSummaryCanvas->SaveAs (Form ("%s/%s_JetEnergyScale.pdf", plotPath.Data(), algo.Data()));
   for (int iEta = 0; iEta < numetabins; iEta++) {
     jetEnergyRes[iEta]->GetYaxis()->SetTitle ("Jet Energy Resolution");
     jetEnergyRes[iEta]->GetYaxis()->SetRangeUser(0, 0.3);
     jetEnergyRes[iEta]->SetLineColor (colors[iEta]);
     jetEnergyRes[iEta]->SetMarkerColor (colors[iEta]);
     //jetEnergyRes[iEta]->SetMarkerStyle (kFullDotLarge);
     if (iEta == 0) jetEnergyRes[iEta]->Draw("hist ][");
     else jetEnergyRes[iEta]->Draw("same hist ][");
     myMarkerText (0.77, 0.88-0.04*iEta, colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
   }
   jesSummaryCanvas->SaveAs (Form ("%s/%s_JetEnergyResolution.pdf", plotPath.Data(), algo.Data()));
   for (int iEta = 0; iEta < numetabins; iEta++) {
     if (jetEnergyScale[iEta]) delete jetEnergyScale[iEta];
     if (jetEnergyRes[iEta]) delete jetEnergyRes[iEta];
   }
   if (jesSummaryCanvas) delete jesSummaryCanvas;
  }


  /**** Plots the photon energy response ****/
  //energyScaleCanvas->cd();
  //gPad->SetLogy(true);

  TH1D* photonEnergyScale[numetabins+1];
  TH1D* photonEnergyRes[numetabins+1];
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   TCanvas* pesCanvas = new TCanvas (Form ("photonEnergyScalePad_%i", iEta), "", 800, 600);

   //if (iEta < numetabins) {
   photonEnergyScale[iEta] = new TH1D (Form ("photonEnergyScale_eta%i", iEta), "", numpbins, pbins);
   photonEnergyRes[iEta] = new TH1D (Form ("photonEnergyRes_eta%i", iEta), "", numpbins, pbins);
   //}

   pesCanvas->cd();
   pesCanvas->Divide(4, 3);

   TF1* calibFits[12];

   for (short iP = 0; iP < numpbins; iP++) {
    if (iP <= 11)
     pesCanvas->cd(iP+1);

    gPad->SetLogy(true);
    TH1D* thisHist = photonEnergyResponse[iP][iEta];
    int nGammas = nGamma[iP][iEta];
    thisHist->Scale(1./thisHist->Integral());
    thisHist->SetAxisRange(3e-6, 2e1, "Y");

    TF1* calibFit = new TF1("calibFit", "gaus(0)", 0, 2.0);
    thisHist->Fit(calibFit, "RNL");
    float m = calibFit->GetParameter(1);
    float s = calibFit->GetParameter(2);
    if (calibFit) delete calibFit;

    calibFit = new TF1("calibFit2", "gaus(0)", m - 3.5*s, m + 3.5*s);
    thisHist->Fit(calibFit, "RNL");
    m = calibFit->GetParameter(1);
    s = calibFit->GetParameter(2);
    if (calibFit) delete calibFit;

    calibFit = new TF1("calibFit3", "gaus(0)", m - 3.5*s, m + 3.5*s);
    thisHist->Fit(calibFit, "RNL");

    thisHist->SetMarkerStyle(kFullDotMedium);
    thisHist->SetLineColor(kBlue);
    thisHist->SetMarkerColor(kBlue);
    thisHist->GetYaxis()->SetTitle ("Counts / Total");
    if (iP <= 11) {
     thisHist->Draw("e1 x0");
     calibFit->SetLineColor(kBlue);
     calibFit->Draw("same");
     calibFits[iP] = calibFit;
    }

    //if (iEta < numetabins) {
    photonEnergyScale[iEta]->SetBinContent (iP+1, calibFit->GetParameter(1));
    photonEnergyScale[iEta]->SetBinError (iP+1, calibFit->GetParError(1));
    //photonEnergyScale[iEta]->SetBinError (iP+1, 0.000001);
    photonEnergyRes[iEta]->SetBinContent (iP+1, calibFit->GetParameter(2));
    photonEnergyRes[iEta]->SetBinError (iP+1, calibFit->GetParError(2));
    //photonEnergyRes[iEta]->SetBinError (iP+1, 0.000001);
    //}

    //myText(0.18, 0.9, kBlack, "Photon energy response");
    //myText(0.18, 0.83, kBlack, Form("%i photons", nGammas));
    myText(0.18, 0.90, kBlack, Form("#mu = %s", FormatMeasurement (calibFit->GetParameter(1), calibFit->GetParError(1), 1)), 0.04 * 2);
    myText(0.18, 0.80, kBlack, Form("#sigma = %s", FormatMeasurement (calibFit->GetParameter(2), calibFit->GetParError(2), 1)), 0.04 * 2);

    //if (iEta < numetabins) {
    // myText (0.18, 0.62, kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]));
    //}
    if (iP <= 11) {
     myText (0.18, 0.70, kBlack, Form("%g < #it{p}_{T}^{reco} < %g", pbins[iP], pbins[iP+1]), 0.04 * 2);
    }

    //energyScaleCanvas->SaveAs(Form("%s/photonEnergyResponse/iP%i_iEta%i.pdf", plotPath.Data(), iP, iEta));

    //if (calibFit) delete calibFit;
   }
   pesCanvas->SaveAs(Form("%s/photonEnergyResponse/iEta%i.pdf", plotPath.Data(), iEta));
   if (pesCanvas) delete pesCanvas;

   for (short iP = 0; iP <= 11; iP++) if (calibFits[iP]) delete calibFits[iP];
  }

  TCanvas* pesSummaryCanvas = new TCanvas("pesSummaryCanvas", "", 800, 600);
  pesSummaryCanvas->cd();
  gPad->SetLogx();
  Color_t colors[14] = {kAzure-2, kBlue, kViolet+7, kMagenta, kPink, kRed+2, kOrange+8, kSpring, kGreen+2, kCyan, kGray, kRed-9, kViolet-8, kBlue-9};
  for (short iEta = 0; iEta < numetabins; iEta++) {
   photonEnergyScale[iEta]->GetYaxis()->SetTitle ("Photon Energy Scale");
   photonEnergyScale[iEta]->GetYaxis()->SetRangeUser (0.98, 1.02);
   photonEnergyScale[iEta]->SetLineColor (colors[iEta]);
   photonEnergyScale[iEta]->SetMarkerColor (colors[iEta]);
   //photonEnergyScale[iEta]->SetMarkerStyle (kFullDotLarge);
   if (iEta == 0) photonEnergyScale[iEta]->Draw("hist ][");
   else photonEnergyScale[iEta]->Draw("same hist ][");
   myMarkerText (0.77, 0.88-0.04*iEta, colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
  }
  pesSummaryCanvas->SaveAs (Form ("%s/PhotonEnergyScale.pdf", plotPath.Data()));
  for (short iEta = 0; iEta < numetabins; iEta++) {
   photonEnergyRes[iEta]->GetYaxis()->SetTitle ("Photon Energy Resolution");
   photonEnergyRes[iEta]->GetYaxis()->SetRangeUser(0, 0.04);
   photonEnergyRes[iEta]->SetLineColor (colors[iEta]);
   photonEnergyRes[iEta]->SetMarkerColor (colors[iEta]);
   //photonEnergyRes[iEta]->SetMarkerStyle (kFullDotLarge);
   if (iEta == 0) photonEnergyRes[iEta]->Draw("hist ][");
   else photonEnergyRes[iEta]->Draw("same hist ][");
   myMarkerText (0.77, 0.88-0.04*iEta, colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
  }
  pesSummaryCanvas->SaveAs (Form ("%s/PhotonEnergyResolution.pdf", plotPath.Data()));
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   if (photonEnergyScale[iEta]) delete photonEnergyScale[iEta];
   if (photonEnergyRes[iEta]) delete photonEnergyRes[iEta];
  }
  if (pesSummaryCanvas) delete pesSummaryCanvas;


// Miscellaneous plots
  //lowResponseEtaPhi->Draw("col");
  //energyScaleCanvas->SaveAs(Form("%s/lowJetResponseEtaPhi.pdf", plotPath.Data()));
  //highResponseEtaPhi->Draw("col");
  //energyScaleCanvas->SaveAs(Form("%s/highJetResponseEtaPhi.pdf", plotPath.Data()));
  //lowResponseEtaPt->Draw("col");
  //energyScaleCanvas->SaveAs(Form("%s/lowJetResponseEtaPt.pdf", plotPath.Data()));
  //highResponseEtaPt->Draw("col");
  //energyScaleCanvas->SaveAs(Form("%s/highJetResponseEtaPt.pdf", plotPath.Data()));
  //gPad->SetLogz(true);
  //TProfile* responsePtProf = responsePt->ProfileY("reponsePtProfY");
  //TGraphAsymmErrors* responsePtProfGraph = new TGraphAsymmErrors(responsePtProf->GetNbinsX());
  //for (int binx = 0; binx < responsePtProf->GetNbinsX(); binx++) {
  // responsePtProfGraph->SetPoint(binx, responsePtProf->GetBinContent(binx+1), responsePtProf->GetBinCenter(binx+1));
  // responsePtProfGraph->SetPointEYlow(binx, responsePtProf->GetBinLowEdge(binx+1)+responsePtProf->GetBinWidth(binx+1)-responsePtProf->GetBinCenter(binx+1));
  // responsePtProfGraph->SetPointEYhigh(binx, responsePtProf->GetBinCenter(binx+1)-responsePtProf->GetBinLowEdge(binx+1));
  // responsePtProfGraph->SetPointEXlow(binx, responsePtProf->GetBinError(binx+1));
  // responsePtProfGraph->SetPointEXhigh(binx, responsePtProf->GetBinError(binx+1));
  //}
  //responsePt->Draw("col");
  //responsePtProfGraph->SetMarkerStyle(kDot);
  //responsePtProfGraph->Draw("p2 same");
  //energyScaleCanvas->SaveAs(Form("%s/responsePt.pdf", plotPath.Data()));
  //responsePtProf->SetAxisRange(0.95, 1.10, "Y");
  //responsePtProf->Draw("e1");
  //energyScaleCanvas->SaveAs(Form("%s/responsePtProfY.pdf", plotPath.Data()));
  //gPad->SetLogy(true);
  //jetSpectrum->Scale(1e6, "width");
  //jetSpectrum->Draw("e1");
  //energyScaleCanvas->SaveAs(Form("%s/jetSpectrum.pdf", plotPath.Data()));
  

  return;
}

}
