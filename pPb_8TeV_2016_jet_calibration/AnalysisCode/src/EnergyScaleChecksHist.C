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
  vector<int> runNumbers(0);
  for (short i = 0; i < sizeof(full_run_list)/sizeof(full_run_list[0]); i++) runNumbers.push_back(full_run_list[i]);
  vector<TString> gammaJetOverlaySampleIds(0);
  for (short i = 0; i < 6; i++) {
   gammaJetOverlaySampleIds.push_back (TString ("Pbp_Overlay_GammaJet_Slice") + to_string (i+1));
   gammaJetOverlaySampleIds.push_back (TString ("pPb_Overlay_GammaJet_Slice") + to_string (i+1));
  }
  vector<TString> gammaJetSignalSampleIds(0);
  for (short i = 0; i < 6; i++) {
   gammaJetSignalSampleIds.push_back (TString ("pPb_Signal_GammaJet_Slice") + to_string (i+1));
  }
  vector<TString> zeeJetSampleIds(0);
  zeeJetSampleIds.push_back("Pbp_Overlay_ZeeJet");
  zeeJetSampleIds.push_back("pPb_Overlay_ZeeJet");

  vector<TString> zmumuJetSampleIds(0);
  zmumuJetSampleIds.push_back("Pbp_Overlay_ZmumuJet");
  zmumuJetSampleIds.push_back("pPb_Overlay_ZmumuJet");

  vector<TString> dijetSampleIds(0);
  dijetSampleIds.push_back("pPb_Signal_Dijet_Slice2");


  /**** Initialize histograms ****/
  TH1D* electronEnergyScale = new TH1D ("electronEnergyScale", ";Electron #it{p}_{T}^{reco} / #it{p}_{T}^{truth};", 200, 0, 2.0);
  electronEnergyScale->Sumw2();

  TH2D*** lowResponseEtaPhi = Get2DArray <TH2D*> (2, 2);
  TH2D*** highResponseEtaPhi = Get2DArray <TH2D*> (2, 2);
  TH2D*** lowResponseEtaPt = Get2DArray <TH2D*> (2, 2);
  TH2D*** highResponseEtaPt = Get2DArray <TH2D*> (2, 2);
  TH2D*** responsePt = Get2DArray <TH2D*> (2, 2);
  TH1D*** jetSpectrum = Get2DArray <TH1D*> (2, 2);

  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   const TString algo = (iAlgo == 0 ? "akt4hi":"akt4emtopo");
   for (short iMC = 0; iMC < 2; iMC++) {
    const TString mcType = (iMC == 0 ? "overlay":"signal");

    lowResponseEtaPhi[iAlgo][iMC] = new TH2D (Form ("%s_%s_lowResponseEtaPhi", algo.Data(), mcType.Data()), ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
    highResponseEtaPhi[iAlgo][iMC] = new TH2D (Form ("%s_%s_highResponseEtaPhi", algo.Data(), mcType.Data()), ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
    lowResponseEtaPt[iAlgo][iMC] = new TH2D (Form ("%s_%s_lowResponseEtaPt", algo.Data(), mcType.Data()), ";#eta;#it{p}_{T}^{reco} #left[GeV#right];", 25, -2.5, 2.5, 66, 20, 350);
    highResponseEtaPt[iAlgo][iMC] = new TH2D (Form ("%s_%s_highResponseEtaPt", algo.Data(), mcType.Data()), ";#eta;#it{p}_{T}^{reco} #left[GeV#right];", 25, -2.5, 2.5, 66, 20, 350);
    responsePt[iAlgo][iMC] = new TH2D (Form ("%s_%s_responsePt", algo.Data(), mcType.Data()), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};#it{p}_{T}^{reco} #left[GeV#right];", 200, 0, 4.0, 66, 20, 350);
    jetSpectrum[iAlgo][iMC] = new TH1D (Form ("%s_%s_jetSpectrum", algo.Data(), mcType.Data()), ";#it{p}_{T}^{jet} #left[GeV#right];", 33, 20, 350);
  
   }
  } 

  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   for (short iMC = 0; iMC < 2; iMC++) {
    lowResponseEtaPhi[iAlgo][iMC]->Sumw2();
    highResponseEtaPhi[iAlgo][iMC]->Sumw2();
    lowResponseEtaPt[iAlgo][iMC]->Sumw2();
    highResponseEtaPt[iAlgo][iMC]->Sumw2();
    responsePt[iAlgo][iMC]->Sumw2();
    jetSpectrum[iAlgo][iMC]->Sumw2();
   }
  }

  TH1D***** jetEnergyResponseCalib = Get4DArray <TH1D*> (2, 2, numpbins+1, numetabins+1);
  TH1D***** jetEnergyResponseReco = Get4DArray <TH1D*> (2, 2, numpbins+1, numetabins+1);
  TH1D**** photonEnergyResponse = Get3DArray <TH1D*> (2, numpbins+1, numetabins+1);

  for (short iP = 0; iP <= numpbins; iP++) {
   for (short iEta = 0; iEta <= numetabins; iEta++) {
    for (short iMC = 0; iMC < 2; iMC++) {
     const TString mcType = (iMC == 0 ? "overlay" : "signal");
     for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
      const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
      jetEnergyResponseCalib[iAlgo][iMC][iP][iEta] = new TH1D (Form ("%s_%s_jetEnergyResponseCalib_iP%i_iEta%i", algo.Data(), mcType.Data(), iP, iEta), "", 50, 0, 2);
      jetEnergyResponseCalib[iAlgo][iMC][iP][iEta]->Sumw2();
      jetEnergyResponseReco[iAlgo][iMC][iP][iEta] = new TH1D (Form ("%s_%s_jetEnergyResponseReco_iP%i_iEta%i", algo.Data(), mcType.Data(), iP, iEta), "", 50, 0, 2);
      jetEnergyResponseReco[iAlgo][iMC][iP][iEta]->Sumw2();
     }
     photonEnergyResponse[iMC][iP][iEta] = new TH1D (Form ("%s_photonEnergyResponse_iP%i_iEta%i", mcType.Data(), iP, iEta), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", 50, 0, 2);
     photonEnergyResponse[iMC][iP][iEta]->Sumw2();
    }
   }
  }

  int**** nJet = Get4DArray <int> (2, 2, numpbins+1, numetabins+1);
  int*** nGamma = Get3DArray <int> (2, numpbins+1, numetabins+1);

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

     // do this if gamma jet MC sample (OVERLAY)
     for (TString gammaJetOverlaySampleId : gammaJetOverlaySampleIds) { // check for gamma jet MC
      if (fname.Contains(gammaJetOverlaySampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", gammaJetOverlaySampleId.Data()));
       nGammaVec = (TVectorD*)thisFile->Get (Form ("nGammaVec_%s", gammaJetOverlaySampleId.Data()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
        lowResponseEtaPhi[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_lowResponseEtaPhi_dataSet%s", algo.Data(), gammaJetOverlaySampleId.Data())));
        highResponseEtaPhi[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_highResponseEtaPhi_dataSet%s", algo.Data(), gammaJetOverlaySampleId.Data())));
        lowResponseEtaPt[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_lowResponseEtaPt_dataSet%s", algo.Data(), gammaJetOverlaySampleId.Data())));
        highResponseEtaPt[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_highResponseEtaPt_dataSet%s", algo.Data(), gammaJetOverlaySampleId.Data())));
        responsePt[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_responsePt_dataSet%s", algo.Data(), gammaJetOverlaySampleId.Data())));
        jetSpectrum[iAlgo][0]->Add ((TH1D*)thisFile->Get (Form ("%s_jetPt_dataSet%s", algo.Data(), gammaJetOverlaySampleId.Data())));
       }

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        for (short iP = 0; iP <= numpbins; iP++) {
         for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
          const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
          nJet[iAlgo][0][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo*(numpbins+1)*(numetabins+1)];
          jetEnergyResponseCalib[iAlgo][0][iP][iEta]->Add ((TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data(), gammaJetOverlaySampleId.Data(), iP, iEta)));
          jetEnergyResponseReco[iAlgo][0][iP][iEta]->Add ((TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data(), gammaJetOverlaySampleId.Data(), iP, iEta)));
         }
         nGamma[0][iP][iEta] += (*nGammaVec)[iP + (numpbins+1)*iEta];
         photonEnergyResponse[0][iP][iEta]->Add ((TH1D*)thisFile->Get (Form ("photonEnergyResponse_dataSet%s_iP%i_iEta%i", gammaJetOverlaySampleId.Data(), iP, iEta)));
        }
       }

       thisFile->Close();
       delete thisFile;
       break;
      }
     }
     // do this if gamma jet MC sample (SIGNAL)
     for (TString gammaJetSignalSampleId : gammaJetSignalSampleIds) { // check for gamma jet MC
      if (fname.Contains(gammaJetSignalSampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", gammaJetSignalSampleId.Data()));
       nGammaVec = (TVectorD*)thisFile->Get (Form ("nGammaVec_%s", gammaJetSignalSampleId.Data()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
        lowResponseEtaPhi[iAlgo][1]->Add ((TH2D*)thisFile->Get (Form ("%s_lowResponseEtaPhi_dataSet%s", algo.Data(), gammaJetSignalSampleId.Data())));
        highResponseEtaPhi[iAlgo][1]->Add ((TH2D*)thisFile->Get (Form ("%s_highResponseEtaPhi_dataSet%s", algo.Data(), gammaJetSignalSampleId.Data())));
        lowResponseEtaPt[iAlgo][1]->Add ((TH2D*)thisFile->Get (Form ("%s_lowResponseEtaPt_dataSet%s", algo.Data(), gammaJetSignalSampleId.Data())));
        highResponseEtaPt[iAlgo][1]->Add ((TH2D*)thisFile->Get (Form ("%s_highResponseEtaPt_dataSet%s", algo.Data(), gammaJetSignalSampleId.Data())));
        responsePt[iAlgo][1]->Add ((TH2D*)thisFile->Get (Form ("%s_responsePt_dataSet%s", algo.Data(), gammaJetSignalSampleId.Data())));
        jetSpectrum[iAlgo][1]->Add ((TH1D*)thisFile->Get (Form ("%s_jetPt_dataSet%s", algo.Data(), gammaJetSignalSampleId.Data())));
       }

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        for (short iP = 0; iP <= numpbins; iP++) {
         for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
          const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
          nJet[iAlgo][1][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo*(numpbins+1)*(numetabins+1)];
          jetEnergyResponseCalib[iAlgo][1][iP][iEta]->Add ((TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data(), gammaJetSignalSampleId.Data(), iP, iEta)));
          jetEnergyResponseReco[iAlgo][1][iP][iEta]->Add ((TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data(), gammaJetSignalSampleId.Data(), iP, iEta)));
         }
         nGamma[1][iP][iEta] += (*nGammaVec)[iP + (numpbins+1)*iEta];
         photonEnergyResponse[1][iP][iEta]->Add ((TH1D*)thisFile->Get (Form ("photonEnergyResponse_dataSet%s_iP%i_iEta%i", gammaJetSignalSampleId.Data(), iP, iEta)));
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
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", zeeJetSampleId.Data()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
        lowResponseEtaPhi[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_lowResponseEtaPhi_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
        highResponseEtaPhi[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_highResponseEtaPhi_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
        lowResponseEtaPt[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_lowResponseEtaPt_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
        highResponseEtaPt[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_highResponseEtaPt_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
        responsePt[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_responsePt_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
        jetSpectrum[iAlgo][0]->Add ((TH1D*)thisFile->Get (Form ("%s_jetPt_dataSet%s", algo.Data(), zeeJetSampleId.Data())));
       }

       electronEnergyScale->Add ((TH1D*)thisFile->Get (Form ("electronEnergyScale_dataSet%s", zeeJetSampleId.Data())));

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        for (short iP = 0; iP <= numpbins; iP++) {
         for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
          const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
          nJet[iAlgo][0][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo*(numpbins+1)*(numetabins+1)];
          jetEnergyResponseCalib[iAlgo][0][iP][iEta]->Add ((TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data(), zeeJetSampleId.Data(), iP, iEta)));
          jetEnergyResponseReco[iAlgo][0][iP][iEta]->Add ((TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data(), zeeJetSampleId.Data(), iP, iEta)));
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
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", zmumuJetSampleId.Data()));

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
        lowResponseEtaPhi[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_lowResponseEtaPhi_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
        highResponseEtaPhi[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_highResponseEtaPhi_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
        lowResponseEtaPt[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_lowResponseEtaPt_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
        highResponseEtaPt[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_highResponseEtaPt_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
        responsePt[iAlgo][0]->Add ((TH2D*)thisFile->Get (Form ("%s_responsePt_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
        jetSpectrum[iAlgo][0]->Add ((TH1D*)thisFile->Get (Form ("%s_jetPt_dataSet%s", algo.Data(), zmumuJetSampleId.Data())));
       }

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        for (short iP = 0; iP <= numpbins; iP++) {
         for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
          const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
          nJet[iAlgo][0][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo*(numpbins+1)*(numetabins+1)];
          jetEnergyResponseCalib[iAlgo][0][iP][iEta]->Add ((TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data(), zmumuJetSampleId.Data(), iP, iEta)));
          jetEnergyResponseReco[iAlgo][0][iP][iEta]->Add ((TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data(), zmumuJetSampleId.Data(), iP, iEta)));
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
  thisHist->Draw ("e1 x0");
  myText(0.18, 0.85, kBlack, "Electron energy response");
  myText(0.18, 0.78, kBlack, Form ("%i electrons", (int)ncounts));
  myText(0.18, 0.71, kBlack, Form ("En. Scale = %.5f #pm %.5f", gausFit->GetParameter(1), gausFit->GetParError(1)));
  myText(0.18, 0.64, kBlack, Form ("En. Res. = %.5f #pm %.5f", gausFit->GetParameter(2), gausFit->GetParError(2)));
  energyScaleCanvas->SaveAs(Form ("%s/electronEnergyScale.pdf", plotPath.Data()));


  /**** Plots the jet energy response ****/
  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

   TH1D*** jetEnergyScale = Get2DArray <TH1D*> (2, numetabins+1);
   TH1D*** jetEnergyRes = Get2DArray <TH1D*> (2, numetabins+1);

   for (short iMC = 0; iMC < 2; iMC++) {
    for (int iEta = 0; iEta <= numetabins; iEta++) {
     TCanvas* jesCanvas = new TCanvas (Form ("%s_jetEnergyScalePad_%i", algo.Data(), iEta), "", 800, 600);

     jetEnergyScale[iMC][iEta] = new TH1D (Form ("%s_jetEnergyScale_eta%i", algo.Data(), iEta), "", numpbins, pbins);
     jetEnergyRes[iMC][iEta] = new TH1D (Form ("%s_jetEnergyRes_eta%i", algo.Data(), iEta), "", numpbins, pbins);

     jesCanvas->cd();
     jesCanvas->Divide(4, 3);

     TF1* recoFits[12];
     TF1* calibFits[12];
     for (int iP = 0; iP < numpbins; iP++) {
      if (iP <= 11) {
       jesCanvas->cd(iP+1);
      }
      TH1D* thisHist = jetEnergyResponseReco[iAlgo][iMC][iP][iEta];
      int nJets = nJet[iAlgo][iMC][iP][iEta];

      TF1* recoFit = new TF1 (Form ("recoFit_iP%i", iP), "gaus(0)", 0, 2.0);
      thisHist->Scale (1./thisHist->Integral());
      thisHist->Fit (recoFit, "RNL");
      double m = recoFit->GetParameter(1);
      double s = recoFit->GetParameter(2);
      if (recoFit) delete recoFit;

      recoFit = new TF1 (Form ("recoFit2_iP%i", iP), "gaus(0)", m - 1.5*s, m + 1.5*s);
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
       thisHist->Draw ("e1 x0");
       recoFit->SetLineColor(kBlack);
       recoFit->Draw ("same");
       recoFits[iP] = recoFit;
      }

      thisHist = jetEnergyResponseCalib[iAlgo][iMC][iP][iEta];
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
       thisHist->Draw ("same e1 x0");
       calibFit->SetLineColor(kBlue);
       calibFit->Draw ("same");
       calibFits[iP] = calibFit;
      }

      //jetEnergyScale[iMC][iEta]->SetBinContent (iP+1, calibFit->GetParameter(1));
      jetEnergyScale[iMC][iEta]->SetBinContent (iP+1, thisHist->GetMean());

      //jetEnergyScale[iMC][iEta]->SetBinError (iP+1, calibFit->GetParError(1));
      //jetEnergyScale[iMC][iEta]->SetBinError (iP+1, thisHist->GetMeanError());
      jetEnergyScale[iMC][iEta]->SetBinError (iP+1, 0.000001);

      //jetEnergyRes[iMC][iEta]->SetBinContent (iP+1, calibFit->GetParameter(2));
      jetEnergyRes[iMC][iEta]->SetBinContent (iP+1, thisHist->GetStdDev());

      //jetEnergyRes[iMC][iEta]->SetBinError (iP+1, calibFit->GetParError(2));
      //jetEnergyRes[iMC][iEta]->SetBinError (iP+1, thisHist->GetStdDevError());
      jetEnergyRes[iMC][iEta]->SetBinError (iP+1, 0.000001);

      myText(0.18, 0.90, kBlack, Form ("#mu = %s", FormatMeasurement (calibFit->GetParameter(1), calibFit->GetParError(1), 1)), 0.04 * 2);
      myText(0.18, 0.80, kBlack, Form ("#sigma = %s", FormatMeasurement (calibFit->GetParameter(2), calibFit->GetParError(2), 1)), 0.04 * 2);

      //if (iEta < numetabins) {
      // myText (0.18, 0.62, kBlack, Form ("%g < #eta_{Lab}^{Jet} < %g", etabins[iEta], etabins[iEta+1]));
      //}
      if (iP <= 11) {
       myText (0.18, 0.70, kBlack, Form ("%g < #it{p}_{T}^{J} < %g", pbins[iP], pbins[iP+1]), 0.04 * 2);
      }

     }
     

     jesCanvas->cd(1);
     myText (0.2, 0.8, kBlack, Form ("%g < #eta_{Lab} < %g", etabins[iEta], etabins[iEta+1]), 0.03*3);
     if (iAlgo == 0)
      myText (0.2, 0.65, kBlack, "HI jets", 0.03*3);
     else
      myText (0.2, 0.65, kBlack, "EMTopo jets");

     jesCanvas->SaveAs(Form ("%s/jetEnergyResponse/%s_iEta%i.pdf", plotPath.Data(), algo.Data(), iEta));
     if (jesCanvas) delete jesCanvas;

     for (short iP = 0; iP <= 11; iP++) {
      if (recoFits[iP]) delete recoFits[iP];
      if (calibFits[iP]) delete calibFits[iP];
     }
    }
   }

   TCanvas* jesSummaryCanvas = new TCanvas("jesSummaryCanvas", "", 800, 600);
   TH1D* thisHist = NULL;
   jesSummaryCanvas->cd();
   gPad->SetLogx();
   gStyle->SetErrorX(0.5);
   Color_t colors[14] = {kBlack, kAzure, kViolet, kMagenta, kPink+10, kPink, kOrange+10, kOrange, kSpring+10, kSpring, kTeal+10, kTeal, kAzure+10, kGray};
   for (int iEta = 0; iEta < numetabins; iEta++) {
     thisHist = jetEnergyScale[0][iEta];
     thisHist->GetXaxis()->SetTitle ("#it{p}_{T}^{jet} #left[GeV#right]");
     thisHist->GetYaxis()->SetTitle ("Jet Energy Scale");
     thisHist->GetYaxis()->SetRangeUser (0.85, 1.15);
     thisHist->SetLineColor (colors[iEta]);
     thisHist->SetMarkerColor (colors[iEta]);
     thisHist->SetMarkerStyle (kFullDotLarge);
     thisHist->SetMarkerSize (0.5);
     if (iEta == 0) thisHist->Draw ("L ][");
     else thisHist->Draw ("same L ][");
     myMarkerText (0.22+(iEta>=(numetabins/2))*0.18, 0.45-0.04*iEta+0.04*(numetabins/2)*(iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
   }
   myText (0.48, 0.9, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
   myText (0.60, 0.22, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
   jesSummaryCanvas->SaveAs (Form ("%s/%s_JetEnergyScale.pdf", plotPath.Data(), algo.Data()));

   for (short iMC = 0; iMC < 2; iMC++) {
    Color_t color = (iMC == 0 ? mcOverlayColor:mcSignalColor);
    thisHist = jetEnergyScale[iMC][numetabins];
    thisHist->GetXaxis()->SetTitle ("#it{p}_{T}^{jet} #left[GeV#right]");
    thisHist->GetYaxis()->SetTitle ("Jet Energy Scale");
    thisHist->GetYaxis()->SetRangeUser (0.85, 1.15);
    thisHist->SetLineColor (color);
    thisHist->SetMarkerColor (color);
    thisHist->SetMarkerStyle (kFullDotLarge);
    if (iMC == 0) thisHist->Draw ("L ][");
    else thisHist->Draw ("same L ][");
    myMarkerText (0.22, 0.28-iMC*0.06, color, kFullDotLarge, Form ("Pythia8 #it{pp} 8.16 TeV%s", (iMC == 0 ? " with #it{p}-Pb Overlay":"")), 1.25, 0.04);
   }
   myText (0.60, 0.82, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
   jesSummaryCanvas->SaveAs (Form ("%s/%s_JetEnergyScale_MCcomparison.pdf", plotPath.Data(), algo.Data()));

   for (int iEta = 0; iEta < numetabins; iEta++) {
     TH1D* thisHist = jetEnergyRes[0][iEta];
     thisHist->GetXaxis()->SetTitle ("#it{p}_{T}^{jet} #left[GeV#right]");
     thisHist->GetYaxis()->SetTitle ("Jet Energy Resolution");
     thisHist->GetYaxis()->SetRangeUser(0, 0.3);
     thisHist->SetLineColor (colors[iEta]);
     thisHist->SetMarkerColor (colors[iEta]);
     thisHist->SetMarkerStyle (kFullDotLarge);
     thisHist->SetMarkerSize (0.5);
     if (iEta == 0) thisHist->Draw ("L ][");
     else thisHist->Draw ("same L ][");
     myMarkerText (0.57+(iEta>=(numetabins/2))*0.18, 0.85-0.04*iEta+0.04*(numetabins/2)*(iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
   }
   myText (0.48, 0.9, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
   myText (0.20, 0.22, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
   jesSummaryCanvas->SaveAs (Form ("%s/%s_JetEnergyResolution.pdf", plotPath.Data(), algo.Data()));

   for (short iMC = 0; iMC < 2; iMC++) {
    Color_t color = (iMC == 0 ? mcOverlayColor:mcSignalColor);
    thisHist = jetEnergyRes[iMC][numetabins];
    thisHist->GetXaxis()->SetTitle ("#it{p}_{T}^{jet} #left[GeV#right]");
    thisHist->GetYaxis()->SetTitle ("Jet Energy Resolution");
    thisHist->GetYaxis()->SetRangeUser (0, 0.3);
    thisHist->SetLineColor (color);
    thisHist->SetMarkerColor (color);
    thisHist->SetMarkerStyle (kFullDotLarge);
    if (iMC == 0) thisHist->Draw ("L ][");
    else thisHist->Draw ("same L ][");
    myMarkerText (0.22, 0.28-iMC*0.06, color, kFullDotLarge, Form ("Pythia8 #it{pp} 8.16 TeV%s", (iMC == 0 ? " with #it{p}-Pb Overlay":"")), 1.25, 0.04);
   }
   myText (0.60, 0.82, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
   jesSummaryCanvas->SaveAs (Form ("%s/%s_JetEnergyResolution_MCcomparison.pdf", plotPath.Data(), algo.Data()));

   for (short iMC = 0; iMC < 2; iMC++) {
    for (int iEta = 0; iEta < numetabins; iEta++) {
     if (jetEnergyScale[iMC][iEta]) delete jetEnergyScale[iMC][iEta];
     if (jetEnergyRes[iMC][iEta]) delete jetEnergyRes[iMC][iEta];
    }
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
    TH1D* thisHist = photonEnergyResponse[0][iP][iEta];
    int nGammas = nGamma[0][iP][iEta];
    thisHist->Scale(1./thisHist->Integral());
    thisHist->SetAxisRange(3e-6, 2e1, "Y");

    TF1* calibFit = new TF1("calibFit", "gaus(0)", 0, 2.0);
    thisHist->Fit(calibFit, "RNL");
    float m = calibFit->GetParameter(1);
    float s = calibFit->GetParameter(2);
    if (calibFit) delete calibFit;

    calibFit = new TF1("calibFit2", "gaus(0)", m - 3.5*s, m + 3.5*s);
    thisHist->Fit(calibFit, "RNL");
    //m = calibFit->GetParameter(1);
    //s = calibFit->GetParameter(2);
    //if (calibFit) delete calibFit;

    //calibFit = new TF1("calibFit3", "gaus(0)", m - 3.5*s, m + 3.5*s);
    //thisHist->Fit(calibFit, "RNL");

    thisHist->SetMarkerStyle(kFullDotMedium);
    thisHist->SetLineColor(kBlue);
    thisHist->SetMarkerColor(kBlue);
    thisHist->GetYaxis()->SetTitle ("Counts / Total");
    if (iP <= 11) {
     thisHist->Draw ("e1 x0");
     calibFit->SetLineColor(kBlue);
     calibFit->Draw ("same");
     calibFits[iP] = calibFit;
    }

    //photonEnergyScale[iEta]->SetBinContent (iP+1, calibFit->GetParameter(1));
    photonEnergyScale[iEta]->SetBinContent (iP+1, thisHist->GetMean());

    //photonEnergyScale[iEta]->SetBinError (iP+1, calibFit->GetParError(1));
    //photonEnergyScale[iEta]->SetBinError (iP+1, thisHist->GetMeanError());
    photonEnergyScale[iEta]->SetBinError (iP+1, 0.000001);

    //photonEnergyRes[iEta]->SetBinContent (iP+1, calibFit->GetParameter(2));
    photonEnergyRes[iEta]->SetBinContent (iP+1, thisHist->GetStdDev());

    //photonEnergyRes[iEta]->SetBinError (iP+1, calibFit->GetParError(2));
    //photonEnergyRes[iEta]->SetBinError (iP+1, thisHist->GetStdDevError());
    photonEnergyRes[iEta]->SetBinError (iP+1, 0.000001);

    //myText(0.18, 0.9, kBlack, "Photon energy response");
    //myText(0.18, 0.83, kBlack, Form ("%i photons", nGammas));
    myText(0.18, 0.90, kBlack, Form ("#mu = %s", FormatMeasurement (calibFit->GetParameter(1), calibFit->GetParError(1), 1)), 0.04 * 2);
    myText(0.18, 0.80, kBlack, Form ("#sigma = %s", FormatMeasurement (calibFit->GetParameter(2), calibFit->GetParError(2), 1)), 0.04 * 2);

    //if (iEta < numetabins) {
    // myText (0.18, 0.62, kBlack, Form ("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]));
    //}
    if (iP <= 11) {
     myText (0.18, 0.70, kBlack, Form ("%g < #it{p}_{T}^{reco} < %g", pbins[iP], pbins[iP+1]), 0.04 * 2);
    }

    //energyScaleCanvas->SaveAs(Form ("%s/photonEnergyResponse/iP%i_iEta%i.pdf", plotPath.Data(), iP, iEta));

    //if (calibFit) delete calibFit;
   }
   pesCanvas->SaveAs(Form ("%s/photonEnergyResponse/iEta%i.pdf", plotPath.Data(), iEta));
   if (pesCanvas) delete pesCanvas;

   for (short iP = 0; iP <= 11; iP++) if (calibFits[iP]) delete calibFits[iP];
  }

  TCanvas* pesSummaryCanvas = new TCanvas("pesSummaryCanvas", "", 800, 600);
  pesSummaryCanvas->cd();
  gPad->SetLogx();
  Color_t colors[14] = {kBlack, kAzure, kViolet+10, kViolet, kPink+10, kPink, kOrange+10, kOrange, kSpring+10, kSpring, kTeal+10, kTeal, kAzure+10, kGray};
  for (short iEta = 0; iEta < numetabins; iEta++) {
   photonEnergyScale[iEta]->GetXaxis()->SetTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
   photonEnergyScale[iEta]->GetYaxis()->SetTitle ("Photon Energy Scale");
   photonEnergyScale[iEta]->GetYaxis()->SetRangeUser (0.985, 1.015);
   photonEnergyScale[iEta]->SetLineColor (colors[iEta]);
   photonEnergyScale[iEta]->SetMarkerColor (colors[iEta]);
   //photonEnergyScale[iEta]->SetMarkerStyle (kFullDotLarge);
   if (iEta == 0) photonEnergyScale[iEta]->Draw ("hist ][");
   else photonEnergyScale[iEta]->Draw ("same hist ][");
   myMarkerText (0.22+(iEta>=(numetabins/2))*0.18, 0.45-0.04*iEta+0.04*(numetabins/2)*(iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
  }
  myText (0.48, 0.9, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
  myText (0.60, 0.22, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
  pesSummaryCanvas->SaveAs (Form ("%s/PhotonEnergyScale.pdf", plotPath.Data()));
  for (short iEta = 0; iEta < numetabins; iEta++) {
   photonEnergyRes[iEta]->GetXaxis()->SetTitle ("#it{p}_{T}^{#gamma} #left[GeV#right]");
   photonEnergyRes[iEta]->GetYaxis()->SetTitle ("Photon Energy Resolution");
   photonEnergyRes[iEta]->GetYaxis()->SetRangeUser(0, 0.07);
   photonEnergyRes[iEta]->SetLineColor (colors[iEta]);
   photonEnergyRes[iEta]->SetMarkerColor (colors[iEta]);
   //photonEnergyRes[iEta]->SetMarkerStyle (kFullDotLarge);
   if (iEta == 0) photonEnergyRes[iEta]->Draw ("hist ][");
   else photonEnergyRes[iEta]->Draw ("same hist ][");
   myMarkerText (0.22+(iEta>=(numetabins/2))*0.18, 0.45-0.04*iEta+0.04*(numetabins/2)*(iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
  }
  myText (0.48, 0.9, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
  myText (0.20, 0.22, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
  pesSummaryCanvas->SaveAs (Form ("%s/PhotonEnergyResolution.pdf", plotPath.Data()));
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   if (photonEnergyScale[iEta]) delete photonEnergyScale[iEta];
   if (photonEnergyRes[iEta]) delete photonEnergyRes[iEta];
  }
  if (pesSummaryCanvas) delete pesSummaryCanvas;


// Miscellaneous plots
  //lowResponseEtaPhi->Draw ("col");
  //energyScaleCanvas->SaveAs(Form ("%s/lowJetResponseEtaPhi.pdf", plotPath.Data()));
  //highResponseEtaPhi->Draw ("col");
  //energyScaleCanvas->SaveAs(Form ("%s/highJetResponseEtaPhi.pdf", plotPath.Data()));
  //lowResponseEtaPt->Draw ("col");
  //energyScaleCanvas->SaveAs(Form ("%s/lowJetResponseEtaPt.pdf", plotPath.Data()));
  //highResponseEtaPt->Draw ("col");
  //energyScaleCanvas->SaveAs(Form ("%s/highJetResponseEtaPt.pdf", plotPath.Data()));
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
  //responsePt->Draw ("col");
  //responsePtProfGraph->SetMarkerStyle(kDot);
  //responsePtProfGraph->Draw ("p2 same");
  //energyScaleCanvas->SaveAs(Form ("%s/responsePt.pdf", plotPath.Data()));
  //responsePtProf->SetAxisRange(0.95, 1.10, "Y");
  //responsePtProf->Draw ("e1");
  //energyScaleCanvas->SaveAs(Form ("%s/responsePtProfY.pdf", plotPath.Data()));
  //gPad->SetLogy(true);
  //jetSpectrum->Scale(1e6, "width");
  //jetSpectrum->Draw ("e1");
  //energyScaleCanvas->SaveAs(Form ("%s/jetSpectrum.pdf", plotPath.Data()));
  

  return;
}

}
