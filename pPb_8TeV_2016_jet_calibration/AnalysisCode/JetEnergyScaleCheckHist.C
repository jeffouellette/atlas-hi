#include "../Params.C"
#include "../../Initialization.C"

void JetEnergyScaleCheckHist () {

  // Setup trigger vectors
  SetupDirectories("JetEnergyScaleCheck", "pPb_8TeV_2016_jet_calibration/");

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
  TH2D* lowResponseEtaPhi = new TH2D("lowResponseEtaPhi", ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
  TH2D* highResponseEtaPhi = new TH2D("highResponseEtaPhi", ";#eta;#phi;", 98, -4.9, 4.9, 100, -pi, pi);
  TH2D* lowResponseEtaPt = new TH2D("lowResponseEtaPt", ";#eta;#it{p}_{T}^{reco} #left[GeV#right];", 25, -2.5, 2.5, 66, 20, 350);
  TH2D* highResponseEtaPt = new TH2D("highResponseEtaPt", ";#eta;#it{p}_{T}^{reco} #left[GeV#right];", 25, -2.5, 2.5, 66, 20, 350);
  TH2D* responsePt = new TH2D("responsePt", ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};#it{p}_{T}^{reco} #left[GeV#right];", 200, 0, 4.0, 66, 20, 350);
  TH1D* jetSpectrum = new TH1D("jetSpectrum", ";#it{p}_{T}^{jet} #left[GeV#right];", 33, 20, 350);
  TH1D* jetEnergyResponseCalib[numpbins+1][numetabins+1];
  TH1D* jetEnergyResponseReco[numpbins+1][numetabins+1];
  TH1D* photonEnergyResponse[numpbins+1][numetabins+1];

  for (short pbin = 0; pbin <= numpbins; pbin++) {
   for (short etabin = 0; etabin <= numetabins; etabin++) {
    jetEnergyResponseCalib[pbin][etabin] = new TH1D(Form("jetEnergyResponseCalib_pbin%i_etabin%i", pbin, etabin), "", 50, 0, 2);
    jetEnergyResponseCalib[pbin][etabin]->Sumw2();
    jetEnergyResponseReco[pbin][etabin] = new TH1D(Form("jetEnergyResponseReco_pbin%i_etabin%i", pbin, etabin), "", 50, 0, 2);
    jetEnergyResponseReco[pbin][etabin]->Sumw2();
    photonEnergyResponse[pbin][etabin] = new TH1D(Form("photonEnergyResponse_pbin%i_etabin%i", pbin, etabin), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", 50, 0, 2);
    photonEnergyResponse[pbin][etabin]->Sumw2();
   }
  }

  int* nJet[numpbins+1] = {};
  int* nGamma[numpbins+1] = {};
  for (short pbin = 0; pbin <= numpbins; pbin++) {
   nJet[pbin] = new int[numetabins+1];
   nGamma[pbin] = new int[numetabins+1];
   for (short etabin = 0; etabin <= numetabins; etabin++) {
    nJet[pbin][etabin] = 0;
    nGamma[pbin][etabin] = 0;
   }
  }

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
   TVectorD *nJetVec, *nGammaVec, *nZeeMassVec, *nZeeJetVec, *nZmumuMassVec, *nZmumuJetVec, *nGammaJetVec;
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

       lowResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPhi_dataSet%s", gammaJetSampleId.Data())));
       highResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("highResponseEtaPhi_dataSet%s", gammaJetSampleId.Data())));
       lowResponseEtaPt->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPt_dataSet%s", gammaJetSampleId.Data())));
       highResponseEtaPt->Add((TH2D*)thisFile->Get(Form("highResponseEtaPt_dataSet%s", gammaJetSampleId.Data())));
       responsePt->Add((TH2D*)thisFile->Get(Form("responsePt_dataSet%s", gammaJetSampleId.Data())));
       jetSpectrum->Add((TH1D*)thisFile->Get(Form("jetPt_dataSet%s", gammaJetSampleId.Data())));

       for (short etabin = 0; etabin <= numetabins; etabin++) {
        for (short pbin = 0; pbin <= numpbins; pbin++) {
         nJet[pbin][etabin] += (*nJetVec)[pbin + (numpbins+1)*etabin];
         nGamma[pbin][etabin] += (*nGammaVec)[pbin + (numpbins+1)*etabin];
         jetEnergyResponseCalib[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("jetEnergyResponseCalib_dataSet%s_pbin%i_etabin%i", gammaJetSampleId.Data(), pbin, etabin)));
         jetEnergyResponseReco[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("jetEnergyResponseReco_dataSet%s_pbin%i_etabin%i", gammaJetSampleId.Data(), pbin, etabin)));
         photonEnergyResponse[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("photonEnergyResponse_dataSet%s_pbin%i_etabin%i", gammaJetSampleId.Data(), pbin, etabin)));
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
       nZeeMassVec = (TVectorD*)thisFile->Get(Form("nZeeMassVec_%s", zeeJetSampleId.Data()));
       nZeeJetVec = (TVectorD*)thisFile->Get(Form("nZeeJetVec_%s", zeeJetSampleId.Data()));

       lowResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPhi_dataSet%s", zeeJetSampleId.Data())));
       highResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("highResponseEtaPhi_dataSet%s", zeeJetSampleId.Data())));
       lowResponseEtaPt->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPt_dataSet%s", zeeJetSampleId.Data())));
       highResponseEtaPt->Add((TH2D*)thisFile->Get(Form("highResponseEtaPt_dataSet%s", zeeJetSampleId.Data())));
       responsePt->Add((TH2D*)thisFile->Get(Form("responsePt_dataSet%s", zeeJetSampleId.Data())));
       jetSpectrum->Add((TH1D*)thisFile->Get(Form("jetPt_dataSet%s", zeeJetSampleId.Data())));

       electronEnergyScale->Add((TH1D*)thisFile->Get(Form("electronEnergyScale_dataSet%s", zeeJetSampleId.Data())));
       for (short etabin = 0; etabin <= numetabins; etabin++) {
        for (short pbin = 0; pbin <= numpbins; pbin++) {
         nJet[pbin][etabin] += (*nJetVec)[pbin + (numpbins+1)*etabin];
         jetEnergyResponseCalib[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("jetEnergyResponseCalib_dataSet%s_pbin%i_etabin%i", zeeJetSampleId.Data(), pbin, etabin)));
         jetEnergyResponseReco[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("jetEnergyResponseReco_dataSet%s_pbin%i_etabin%i", zeeJetSampleId.Data(), pbin, etabin)));
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

       lowResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPhi_dataSet%s", zmumuJetSampleId.Data())));
       highResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("highResponseEtaPhi_dataSet%s", zmumuJetSampleId.Data())));
       lowResponseEtaPt->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPt_dataSet%s", zmumuJetSampleId.Data())));
       highResponseEtaPt->Add((TH2D*)thisFile->Get(Form("highResponseEtaPt_dataSet%s", zmumuJetSampleId.Data())));
       responsePt->Add((TH2D*)thisFile->Get(Form("responsePt_dataSet%s", zmumuJetSampleId.Data())));
       jetSpectrum->Add((TH1D*)thisFile->Get(Form("jetPt_dataSet%s", zmumuJetSampleId.Data())));

       for (short etabin = 0; etabin <= numetabins; etabin++) {
        for (short pbin = 0; pbin <= numpbins; pbin++) {
         nJet[pbin][etabin] += (*nJetVec)[pbin + (numpbins+1)*etabin];
         jetEnergyResponseCalib[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("jetEnergyResponseCalib_dataSet%s_pbin%i_etabin%i", zmumuJetSampleId.Data(), pbin, etabin)));
         jetEnergyResponseReco[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("jetEnergyResponseReco_dataSet%s_pbin%i_etabin%i", zmumuJetSampleId.Data(), pbin, etabin)));
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
  thisHist->Scale(1./ncounts, "width");
  thisHist->Fit(gausFit, "R", "L");
  const float m = gausFit->GetParameter(1);
  const float s = gausFit->GetParameter(2);
  if (gausFit) delete gausFit;
  gausFit = new TF1("gausFit2", "gaus(0)", m - 1.6*s, m + 1.6*s);
  thisHist->Fit(gausFit, "R", "L");

  thisHist->GetYaxis()->SetTitle("Fractional counts / 0.01 GeV");
  //thisHist->SetMarkerStyle(6);
  thisHist->Draw("e1 x0");
  myText(0.18, 0.85, kBlack, "Electron energy response");
  myText(0.18, 0.78, kBlack, Form("%i electrons", (int)ncounts));
  myText(0.18, 0.71, kBlack, Form("En. Scale = %.5f #pm %.5f", gausFit->GetParameter(1), gausFit->GetParError(1)));
  myText(0.18, 0.64, kBlack, Form("En. Res. = %.5f #pm %.5f", gausFit->GetParameter(2), gausFit->GetParError(2)));
  energyScaleCanvas->SaveAs(Form("%s/electronEnergyScale.pdf", plotPath.Data()));


  /**** Plots the jet energy response ****/
  for (int etabin = 0; etabin <= numetabins; etabin++) {
   TCanvas* jesCanvas = new TCanvas (Form ("jetEnergyScalePad_%i", etabin), "", 800, 600);
   jesCanvas->cd();
   jesCanvas->Divide(4, 3);
   for (int pbin = 0; pbin <= 11; pbin++) {
    jesCanvas->cd(pbin+1);
    TH1D* thisHist = jetEnergyResponseReco[pbin][etabin];
    int nJets = nJet[pbin][etabin];

    TF1* recoFit = new TF1("recoFit", "gaus(0)", 0, 2.0);
    thisHist->Scale(1./thisHist->Integral(), "width");
    thisHist->Fit(recoFit, "RNL");
    double m = recoFit->GetParameter(1);
    double s = recoFit->GetParameter(2);
    if (recoFit) delete recoFit;

    recoFit = new TF1("recoFit2", "gaus(0)", m - 1.3*s, m + 1.3*s);
    thisHist->Fit(recoFit, "RNL");
    m = recoFit->GetParameter(1);
    s = recoFit->GetParameter(2);
    if (recoFit) delete recoFit;

    recoFit = new TF1("recoFit3", "gaus(0)", m - 1.3*s, m + 1.3*s);
    thisHist->Fit(recoFit, "RNL");

    thisHist->SetLineColor(kBlack);
    thisHist->SetMarkerColor(kBlack);
    thisHist->GetYaxis()->SetTitle(Form("Fractional counts / %.2f GeV", thisHist->GetBinWidth(1)));
    thisHist->Draw("e1 x0");
    recoFit->SetLineColor(kBlack);
    recoFit->Draw("same");

    thisHist = jetEnergyResponseCalib[pbin][etabin];
    thisHist->Scale(1./thisHist->Integral(), "width");

    TF1* calibFit = new TF1("calibFit", "gaus(0)", 0, 2.0);
    thisHist->Fit(calibFit, "RNL");
    m = calibFit->GetParameter(1);
    s = calibFit->GetParameter(2);
    if (calibFit) delete calibFit;

    calibFit = new TF1("calibFit2", "gaus(0)", m - 1.3*s, m + 1.3*s);
    thisHist->Fit(calibFit, "RNL");
    m = calibFit->GetParameter(1);
    s = calibFit->GetParameter(2);
    if (calibFit) delete calibFit;

    calibFit = new TF1("calibFit3", "gaus(0)", m - 1.3*s, m + 1.3*s);
    thisHist->Fit(calibFit, "RNL");

    thisHist->SetLineColor(kBlue);
    thisHist->SetMarkerColor(kBlue);
    //thisHist->GetYaxis()->SetTitle(Form("Fractional counts / %.2f GeV", thisHist->GetBinWidth(1)));
    thisHist->Draw("same e1 x0");
    calibFit->SetLineColor(kBlue);
    calibFit->Draw("same");

    //myText(0.18, 0.9, kBlack, "Jet energy response");
    //myText(0.18, 0.83, kBlack, Form("%i jets", nJets));
    myText(0.18, 0.90, kBlack, Form("#mu = %s", FormatMeasurement (calibFit->GetParameter(1), calibFit->GetParError(1), 1)), 0.04 * 2);
    myText(0.18, 0.80, kBlack, Form("#sigma = %s", FormatMeasurement (calibFit->GetParameter(2), calibFit->GetParError(2), 1)), 0.04 * 2);
    //myMarkerText(0.58, 0.9, kBlack, kFullCircle, "Uncalibrated (EM scale)", 1.25, 0.04 * 4);
    //myMarkerText(0.58, 0.83, kBlue, kFullCircle, "Calibrated", 1.25, 0.04 * 4);

    //if (etabin < numetabins) {
    // myText (0.18, 0.62, kBlack, Form("%g < #eta_{Lab}^{Jet} < %g", etabins[etabin], etabins[etabin+1]));
    //}
    if (pbin < numpbins) {
     myText (0.18, 0.70, kBlack, Form("%g < #it{p}_{T}^{J} < %g", pbins[pbin], pbins[pbin+1]), 0.04 * 2);
    }

    //energyScaleCanvas->SaveAs(Form("%s/jetEnergyResponse/pbin%i_etabin%i.pdf", plotPath.Data(), pbin, etabin));

    if (calibFit) delete calibFit;
    if (recoFit) delete recoFit;
   }
   jesCanvas->SaveAs(Form("%s/jetEnergyResponse/etabin%i.pdf", plotPath.Data(), etabin));
   if (jesCanvas) delete jesCanvas;
  }


  /**** Plots the photon energy response ****/
  //energyScaleCanvas->cd();
  //gPad->SetLogy(true);
  for (short etabin = 0; etabin <= numetabins; etabin++) {
   TCanvas* jesCanvas = new TCanvas (Form ("photonEnergyScalePad_%i", etabin), "", 800, 600);
   jesCanvas->cd();
   jesCanvas->Divide(4, 3);
   for (short pbin = 0; pbin <= 11; pbin++) {
    jesCanvas->cd(pbin+1);
    gPad->SetLogy(true);
    TH1D* thisHist = photonEnergyResponse[pbin][etabin];
    int nGammas = nGamma[pbin][etabin];
    thisHist->Scale(1./thisHist->Integral(), "width");
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

    thisHist->SetLineColor(kBlue);
    thisHist->SetMarkerColor(kBlue);
    thisHist->GetYaxis()->SetTitle(Form("Fractional counts / %.2f GeV", thisHist->GetBinWidth(1)));
    thisHist->Draw("e1 x0");
    calibFit->SetLineColor(kBlue);
    calibFit->Draw("same");

    //myText(0.18, 0.9, kBlack, "Photon energy response");
    //myText(0.18, 0.83, kBlack, Form("%i photons", nGammas));
    myText(0.18, 0.90, kBlack, Form("#mu = %s", FormatMeasurement (calibFit->GetParameter(1), calibFit->GetParError(1), 1)), 0.04 * 2);
    myText(0.18, 0.80, kBlack, Form("#sigma = %s", FormatMeasurement (calibFit->GetParameter(2), calibFit->GetParError(2), 1)), 0.04 * 2);

    //if (etabin < numetabins) {
    // myText (0.18, 0.62, kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[etabin], etabins[etabin+1]));
    //}
    if (pbin < numpbins) {
     myText (0.18, 0.70, kBlack, Form("%g < #it{p}_{T}^{reco} < %g", pbins[pbin], pbins[pbin+1]), 0.04 * 2);
    }

    //energyScaleCanvas->SaveAs(Form("%s/photonEnergyResponse/pbin%i_etabin%i.pdf", plotPath.Data(), pbin, etabin));

    if (calibFit) delete calibFit;
   }
   jesCanvas->SaveAs(Form("%s/photonEnergyResponse/etabin%i.pdf", plotPath.Data(), etabin));
   if (jesCanvas) delete jesCanvas;
  }


  // Miscellaneous plots
  lowResponseEtaPhi->Draw("col");
  energyScaleCanvas->SaveAs(Form("%s/lowJetResponseEtaPhi.pdf", plotPath.Data()));
  highResponseEtaPhi->Draw("col");
  energyScaleCanvas->SaveAs(Form("%s/highJetResponseEtaPhi.pdf", plotPath.Data()));
  lowResponseEtaPt->Draw("col");
  energyScaleCanvas->SaveAs(Form("%s/lowJetResponseEtaPt.pdf", plotPath.Data()));
  highResponseEtaPt->Draw("col");
  energyScaleCanvas->SaveAs(Form("%s/highJetResponseEtaPt.pdf", plotPath.Data()));
  gPad->SetLogz(true);
  TProfile* responsePtProf = responsePt->ProfileY("reponsePtProfY");
  TGraphAsymmErrors* responsePtProfGraph = new TGraphAsymmErrors(responsePtProf->GetNbinsX());
  for (int binx = 0; binx < responsePtProf->GetNbinsX(); binx++) {
   responsePtProfGraph->SetPoint(binx, responsePtProf->GetBinContent(binx+1), responsePtProf->GetBinCenter(binx+1));
   responsePtProfGraph->SetPointEYlow(binx, responsePtProf->GetBinLowEdge(binx+1)+responsePtProf->GetBinWidth(binx+1)-responsePtProf->GetBinCenter(binx+1));
   responsePtProfGraph->SetPointEYhigh(binx, responsePtProf->GetBinCenter(binx+1)-responsePtProf->GetBinLowEdge(binx+1));
   responsePtProfGraph->SetPointEXlow(binx, responsePtProf->GetBinError(binx+1));
   responsePtProfGraph->SetPointEXhigh(binx, responsePtProf->GetBinError(binx+1));
  }
  responsePt->Draw("col");
  responsePtProfGraph->SetMarkerStyle(kDot);
  responsePtProfGraph->Draw("p2 same");
  energyScaleCanvas->SaveAs(Form("%s/responsePt.pdf", plotPath.Data()));
  responsePtProf->SetAxisRange(0.95, 1.10, "Y");
  responsePtProf->Draw("e1");
  energyScaleCanvas->SaveAs(Form("%s/responsePtProfY.pdf", plotPath.Data()));
  gPad->SetLogy(true);
  jetSpectrum->Scale(1e6, "width");
  jetSpectrum->Draw("e1");
  energyScaleCanvas->SaveAs(Form("%s/jetSpectrum.pdf", plotPath.Data()));
  

  return;
}
