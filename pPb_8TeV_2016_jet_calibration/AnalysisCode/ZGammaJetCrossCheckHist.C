#include "../Params.C"
#include "../../Initialization.C"
#include "../RooZfit.C"

using namespace RooFit;

TH1D* GetProfileX(const TString name, TH2D* hist, const int nbinsx, const double* xbins, const bool useFit) {
  TH1D* prof = new TH1D(name, "", nbinsx, xbins);
  for (int xbin = 1; xbin <= nbinsx; xbin++) {
   TH1D* projy = hist->ProjectionY("projy", xbin, xbin);
   projy->Rebin(rebinFactor);
   //projy->GetXaxis()->SetLimits(0, 2.0);
   double mean, mean_err;
   double chi_square = 0;
   int numNonzeroBins = 0;
   for (int xbin = 1; xbin <= projy->GetNbinsX(); xbin++)
    if (projy->GetBinContent(xbin) > 0) numNonzeroBins++;

   if (useFit && useGaussian && numNonzeroBins > 4) {
    TF1* gaus = new TF1("gaus", "gaus(0)", projy->GetXaxis()->GetBinLowEdge(1), projy->GetXaxis()->GetBinUpEdge(projy->GetNbinsX()));
    projy->Fit(gaus, "Q0R");
    mean = gaus->GetParameter(1);
    mean_err = gaus->GetParError(1);
    chi_square = gaus->GetChisquare() / (projy->GetNbinsX() - 3);
    if (gaus) delete gaus;
   }
   if (!useGaussian || !useFit || chi_square > 1.0 || numNonzeroBins <= 4) {
    mean = projy->GetMean();
    mean_err = projy->GetMeanError();
   }
   prof->SetBinContent(xbin, mean);
   prof->SetBinError(xbin, mean_err);
   if (projy) delete projy;
  }
  return prof;
}

TH1D* GetDataOverMC(const TString name, TH2D* data, TH2D* mc, const int numxbins, const double* xbins, const int numybins, const double* ybins, const bool useFit) {
  TH1D* dataOverMC = new TH1D(name, "", numxbins, xbins);
  for (int xbin = 1; xbin <= numxbins; xbin++) {
   TH1D* projy = data->ProjectionY(name + TString(Form("data_xbin%i", xbin)), xbin, xbin);
   projy->Rebin(rebinFactor);
   double dataAvg, dataErr, mcAvg, mcErr;
   double chi_square = 0;
   int numNonzeroBins = 0;
   for (int xbin = 1; xbin <= projy->GetNbinsX(); xbin++)
    if (projy->GetBinContent(xbin) > 0) numNonzeroBins++;

   if (useFit && useGaussian && numNonzeroBins > 4) {
    TF1* gaus = new TF1("gaus", "gaus(0)", projy->GetXaxis()->GetBinLowEdge(1), 2.0);//projy->GetXaxis()->GetBinUpEdge(projy->GetNbinsX()));
    projy->Fit(gaus, "Q0R");
    dataAvg = gaus->GetParameter(1);
    dataErr = gaus->GetParError(1);
    chi_square = gaus->GetChisquare() / (projy->GetNbinsX() - 3);
    if (gaus) delete gaus;
   }
   if (!useGaussian || !useFit || chi_square > 1.0 || numNonzeroBins <= 4) {
    dataAvg = projy->GetMean();
    dataErr = projy->GetMeanError();
   }
   if (projy) delete projy;

   projy = mc->ProjectionY(name + TString(Form("mc_xbin%i", xbin)), xbin, xbin);
   projy->Rebin(rebinFactor);
   chi_square = 0;
   numNonzeroBins = 0;
   for (int xbin = 1; xbin <= projy->GetNbinsX(); xbin++)
    if (projy->GetBinContent(xbin) > 0) numNonzeroBins++;

   if (useFit && useGaussian && numNonzeroBins > 4) {
    TF1* gaus = new TF1("gaus", "gaus(0)", projy->GetXaxis()->GetBinLowEdge(1), 2.0);//projy->GetXaxis()->GetBinUpEdge(projy->GetNbinsX()));
    projy->Fit(gaus, "Q0R");
    mcAvg = gaus->GetParameter(1);
    mcErr = gaus->GetParError(1);
    chi_square = gaus->GetChisquare() / (projy->GetNbinsX() - 3);
    if (gaus) delete gaus;
   }
   if (!useGaussian || !useFit || chi_square > 1.0 || numNonzeroBins <= 4) {
    mcAvg = projy->GetMean();
    mcErr = projy->GetMeanError();
   }
   if (projy) delete projy;

   const double dataOverMCavg = dataAvg/mcAvg;
   const double dataOverMCerr = dataOverMCavg * TMath::Sqrt(TMath::Power(dataErr/dataAvg, 2) + TMath::Power(mcErr/mcAvg, 2));
   if (!isnan(dataOverMCavg) && !isnan(dataOverMCerr)) {
    dataOverMC->SetBinContent(xbin, dataOverMCavg);
    dataOverMC->SetBinError(xbin, dataOverMCerr);
   }
  }
  dataOverMC->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());
  return dataOverMC;
}

void ZGammaJetCrossCheckHist () {

  // Setup trigger vectors
  SetupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers(0);
  for (short i = 0; i < sizeof(full_run_list)/sizeof(full_run_list[0]); i++) runNumbers.push_back(full_run_list[i]);
  vector<TString> gammaJetSampleIds(0);
  for (short i = 0; i < 6; i++) {
   gammaJetSampleIds.push_back(TString("Pbp") + (runValidation ? "_Valid":"_Overlay") + "_GammaJet_Slice" + to_string(i+1));
   gammaJetSampleIds.push_back(TString("pPb") + (runValidation ? "_Valid":"_Overlay") + "_GammaJet_Slice" + to_string(i+1));
  }
  vector<TString> zeeJetSampleIds(0);
  //for (short i = 0; i < 6; i++) {
  // zeeJetSampleIds.push_back(string("Pbp_ZeeJet") + to_string(i));
  // zeeJetSampleIds.push_back(string("pPb_ZeeJet") + to_string(i));
  //}
  zeeJetSampleIds.push_back("Pbp_ZeeJet_Overlay");
  zeeJetSampleIds.push_back("pPb_ZeeJet_Overlay");

  vector<TString> zmumuJetSampleIds(0);
  zmumuJetSampleIds.push_back("Pbp_ZmumuJet");
  zmumuJetSampleIds.push_back("pPb_ZmumuJet");

  TH2D* zeeJetHists[3][numetabins+1][2][3];
  //TH2D* zeeJetHistsSys[3][numetabins][2][3];
  TH2D* zmumuJetHists[3][numetabins+1][2][3];
  //TH2D* zmumuJetHistsSys[3][numetabins][2][3];
  TH2D* gJetHists[3][numetabins+1][2][3];
  TH2D* gJetHistsSys[3][numetabins+1][2][3];
  TH1D* zMassSpectra[2][2][numetabins+1];
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

  for (short etabin = 0; etabin <= numetabins; etabin++) {
   for (short dType = 0; dType < 2; dType++) { // dType is 0 for data, 1 for MC
    const TString dataType = (dType == 0 ? "data":"mc");

    for (short errType = 0; errType < 3; errType++) {
     TString error = "sys_lo";
     if (errType == 1) error = "stat";
     else if (errType == 2) error = "sys_hi";

     for (short pType = 0; pType < 3; pType++) {
      TString period = "periodA";
      if (pType == 1) period = "periodB";
      else if (pType == 2) period = "periodAB";

      zeeJetHists[pType][etabin][dType][errType] = new TH2D(Form("zeeJetPtRatio_etabin%i_%s_%s_%s", etabin, dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{Z} #left[GeV#right];#it{x}_{J}^{ref}", numpzbins, pzbins, numxjrefbins, xjrefbins);
      zeeJetHists[pType][etabin][dType][errType]->Sumw2();
      zmumuJetHists[pType][etabin][dType][errType] = new TH2D(Form("zmumuJetPtRatio_etabin%i_%s_%s_%s", etabin, dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{Z} #left[GeV#right];#it{x}_{J}^{ref}", numpzbins, pzbins, numxjrefbins, xjrefbins);
      zmumuJetHists[pType][etabin][dType][errType]->Sumw2();
      gJetHists[pType][etabin][dType][errType] = new TH2D(Form("gJetPtRatio_etabin%i_%s_%s_%s", etabin, dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{#gamma} #left[GeV#right];#it{x}_{J}^{ref}", numpbins, pbins, numxjrefbins, xjrefbins);
      gJetHists[pType][etabin][dType][errType]->Sumw2();

      gJetHistsSys[pType][etabin][dType][errType] = new TH2D(Form("gJetPtRatioSys_etabin%i_%s_%s_%s", etabin, dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{jet} #left[GeV#right];#Delta#it{x}_{J}^{ref}#it{p}_{T}^{ref}/#it{p}_{T}^{jet}", numpzbins, pzbins, numSigmaBins, -maxSigma, maxSigma);
      gJetHistsSys[pType][etabin][dType][errType]->Sumw2();
     }
    }

    for (short spcType = 0; spcType < 2; spcType++) {
     const TString species = (spcType == 0 ? "mumu":"ee");

     zMassSpectra[spcType][dType][etabin] = new TH1D(Form("z%sMassSpectrum_%s_etabin%i", species.Data(), dataType.Data(), etabin), "", 50, 60, 110);
     zMassSpectra[spcType][dType][etabin]->Sumw2();
    }
   }
  }
  //for (short etabin = 0; etabin <= numetabins; etabin++) {
  // for (short pType = 0; pType < 3; pType++) {
  //  TString period = "periodA";
  //  if (pType == 1) period = "periodB";
  //  else if (pType == 2) period = "periodAB";
  //  jetEnergyResponseCalib[pType][etabin] = new TH1D(Form("jetEnergyResponseCalib_etabin%i_%s", etabin, period.Data()), "Jet #it{p}_{T}^{reco} / #it{p}_{T}^{truth};", 50, 0, 2.0);
  //  jetEnergyResponseCalib[pType][etabin]->Sumw2();
  //  jetEnergyResponseReco[pType][etabin] = new TH1D(Form("jetEnergyResponseReco_etabin%i_%s", etabin, period.Data()), "Jet #it{p}_{T}^{reco} / #it{p}_{T}^{truth};", 50, 0, 2.0);
  //  jetEnergyResponseReco[pType][etabin]->Sumw2();
  // }
  //}

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
  int nZeeMass[3][2][numetabins+1] = {{{}, {}}, {{}, {}}, {{}, {}}};
  int nZeeJet[3][2][numetabins+1] = {{{}, {}}, {{}, {}}, {{}, {}}};
  int nZmumuMass[3][2][numetabins+1] = {{{}, {}}, {{}, {}}, {{}, {}}};
  int nZmumuJet[3][2][numetabins+1] = {{{}, {}}, {{}, {}}, {{}, {}}};
  int nGammaJet[3][2][numetabins+1] = {{{}, {}}, {{}, {}}, {{}, {}}};

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

     // do this if file is data
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains(to_string(runNumber))) { // if data, do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       const short pType = (runNumber < 313500 ? 0 : 1);
       //infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%i", runNumber));
       nZeeMassVec = (TVectorD*)thisFile->Get(Form("nZeeMassVec_%i", runNumber));
       nZeeJetVec = (TVectorD*)thisFile->Get(Form("nZeeJetVec_%i", runNumber));
       nZmumuMassVec = (TVectorD*)thisFile->Get(Form("nZmumuMassVec_%i", runNumber));
       nZmumuJetVec = (TVectorD*)thisFile->Get(Form("nZmumuJetVec_%i", runNumber));
       nGammaJetVec = (TVectorD*)thisFile->Get(Form("nGammaJetVec_%i", runNumber));

       for (short etabin = 0; etabin <= numetabins; etabin++) {
        const bool flipEta = runNumber < 313500 && etabin < numetabins;
        const short act_etabin = (flipEta ? (numetabins - etabin - 1) : etabin);

        nZeeMass[pType][0][etabin] += (*nZeeMassVec)[etabin];
        nZeeMass[2][0][etabin] += (*nZeeMassVec)[act_etabin];

        nZeeJet[pType][0][etabin] += (*nZeeJetVec)[etabin];
        nZeeJet[2][0][etabin] += (*nZeeJetVec)[act_etabin];

        nZmumuMass[pType][0][etabin] += (*nZmumuMassVec)[etabin];
        nZmumuMass[2][0][etabin] += (*nZmumuMassVec)[act_etabin];

        nZmumuJet[pType][0][etabin] += (*nZmumuJetVec)[etabin];
        nZmumuJet[2][0][etabin] += (*nZmumuJetVec)[act_etabin];

        nGammaJet[pType][0][etabin] += (*nGammaJetVec)[etabin];
        nGammaJet[2][0][etabin] += (*nGammaJetVec)[act_etabin];

        for (short errType = 0; errType < 3; errType++) {
         TString error = "sys_lo";
         if (errType == 1) error = "stat";
         else if (errType == 2) error = "sys_hi";

         TH2D* temp = (TH2D*)thisFile->Get(Form("zeeJetPtRatio_dataSet%i_etabin%i_data_%s", runNumber, etabin, error.Data()));
         zeeJetHists[pType][etabin][0][errType]->Add(temp);
         temp = (TH2D*)thisFile->Get(Form("zeeJetPtRatio_dataSet%i_etabin%i_data_%s", runNumber, act_etabin, error.Data()));
         zeeJetHists[2][etabin][0][errType]->Add(temp);

         temp = (TH2D*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%i_etabin%i_data_%s", runNumber, etabin, error.Data()));
         zmumuJetHists[pType][etabin][0][errType]->Add(temp);
         temp = (TH2D*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%i_etabin%i_data_%s", runNumber, act_etabin, error.Data()));
         zmumuJetHists[2][etabin][0][errType]->Add(temp);

         temp = (TH2D*)thisFile->Get(Form("gJetPtRatio_dataSet%i_etabin%i_data_%s", runNumber, etabin, error.Data()));
         gJetHists[pType][etabin][0][errType]->Add(temp);
         temp = (TH2D*)thisFile->Get(Form("gJetPtRatio_dataSet%i_etabin%i_data_%s", runNumber, act_etabin, error.Data()));
         gJetHists[2][etabin][0][errType]->Add(temp);

         if (errType == 1) 
          gJetHistsSys[pType][etabin][0][errType]->Add((TH2D*)thisFile->Get(Form("gJetPtRatioSys_dataSet%i_etabin%i_data_%s", runNumber, etabin, error.Data())));
         //gJetHistsSys[2][etabin][0][errType]->Add((TH2D*)thisFile->Get(Form("gJetPtRatioSys_dataSet%i_etabin%i_data_%s", runNumber, act_etabin, error.Data())));
        }
       }

       for (short spcType = 0; spcType < 2; spcType++) {
        const TString species = (spcType == 0 ? "mumu":"ee");

        for (short etabin = 0; etabin <= numetabins; etabin++) {
         zMassSpectra[spcType][0][etabin]->Add((TH1D*)thisFile->Get(Form("z%sMassSpectrum_dataSet%i_data_etabin%i", species.Data(), runNumber, etabin)));
        }
       }

       thisFile->Close();
       delete thisFile;
       break;
      }
     }
     // do this if gamma jet MC sample
     for (TString gammaJetSampleId : gammaJetSampleIds) { // check for gamma jet MC
      if (fname.Contains(gammaJetSampleId)) { // if gamma jet MC sample
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       const short pType = (gammaJetSampleId.Contains("pPb") ? 0 : 1);
       //infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%s", gammaJetSampleId.Data()));
       nJetVec = (TVectorD*)thisFile->Get(Form("nJetVec_%s", gammaJetSampleId.Data()));
       nGammaVec = (TVectorD*)thisFile->Get(Form("nGammaVec_%s", gammaJetSampleId.Data()));
       nGammaJetVec = (TVectorD*)thisFile->Get(Form("nGammaJetVec_%s", gammaJetSampleId.Data()));

       lowResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPhi_dataSet%s", gammaJetSampleId.Data())));
       highResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("highResponseEtaPhi_dataSet%s", gammaJetSampleId.Data())));
       lowResponseEtaPt->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPt_dataSet%s", gammaJetSampleId.Data())));
       highResponseEtaPt->Add((TH2D*)thisFile->Get(Form("highResponseEtaPt_dataSet%s", gammaJetSampleId.Data())));
       responsePt->Add((TH2D*)thisFile->Get(Form("responsePt_dataSet%s", gammaJetSampleId.Data())));
       jetSpectrum->Add((TH1D*)thisFile->Get(Form("jetPt_dataSet%s", gammaJetSampleId.Data())));

       for (short etabin = 0; etabin <= numetabins; etabin++) {
        const bool flipEta = gammaJetSampleId.Contains("pPb") && etabin < numetabins;
        const short act_etabin = (flipEta ? (numetabins - etabin - 1) : etabin); // period A condition

        nGammaJet[pType][1][etabin] += (*nGammaJetVec)[etabin];
        nGammaJet[2][1][etabin] += (*nGammaJetVec)[act_etabin];

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH2D* temp = (TH2D*)thisFile->Get(Form("gJetPtRatio_dataSet%s_etabin%i_mc_stat", gammaJetSampleId.Data(), etabin));
        gJetHists[pType][etabin][1][1]->Add(temp);
        temp = (TH2D*)thisFile->Get(Form("gJetPtRatio_dataSet%s_etabin%i_mc_stat", gammaJetSampleId.Data(), act_etabin));
        gJetHists[2][etabin][1][1]->Add(temp);

        gJetHistsSys[pType][etabin][1][1]->Add((TH2D*)thisFile->Get(Form("gJetPtRatioSys_dataSet%s_etabin%i_mc_stat", gammaJetSampleId.Data(), etabin)));
        //gJetHistsSys[2][etabin][1][1]->Add((TH2D*)thisFile->Get(Form("gJetPtRatioSys_dataSet%s_etabin%i_mc_stat", gammaJetSampleId.Data(), act_etabin)));
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
       const short pType = (zeeJetSampleId.Contains("pPb") ? 0 : 1);
       nJetVec = (TVectorD*)thisFile->Get(Form("nJetVec_%s", zeeJetSampleId.Data()));
       nZeeMassVec = (TVectorD*)thisFile->Get(Form("nZeeMassVec_%s", zeeJetSampleId.Data()));
       nZeeJetVec = (TVectorD*)thisFile->Get(Form("nZeeJetVec_%s", zeeJetSampleId.Data()));

       lowResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPhi_dataSet%s", zeeJetSampleId.Data())));
       highResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("highResponseEtaPhi_dataSet%s", zeeJetSampleId.Data())));
       lowResponseEtaPt->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPt_dataSet%s", zeeJetSampleId.Data())));
       highResponseEtaPt->Add((TH2D*)thisFile->Get(Form("highResponseEtaPt_dataSet%s", zeeJetSampleId.Data())));
       responsePt->Add((TH2D*)thisFile->Get(Form("responsePt_dataSet%s", zeeJetSampleId.Data())));
       //jetSpectrum->Add((TH1D*)thisFile->Get(Form("jetPt_dataSet%s", zeeJetSampleId.Data())));

       electronEnergyScale->Add((TH1D*)thisFile->Get(Form("electronEnergyScale_dataSet%s", zeeJetSampleId.Data())));
       for (short etabin = 0; etabin <= numetabins; etabin++) {
        const bool flipEta = zeeJetSampleId.Contains("pPb") && etabin < numetabins;
        const short act_etabin = (flipEta ? numetabins - etabin - 1 : etabin); // period A condition

        nZeeMass[pType][1][etabin] += (*nZeeMassVec)[etabin];
        nZeeMass[2][1][etabin] += (*nZeeMassVec)[act_etabin];
        nZeeJet[pType][1][etabin] += (*nZeeJetVec)[etabin];
        nZeeJet[2][1][etabin] += (*nZeeJetVec)[act_etabin];

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH2D* temp = (TH2D*)thisFile->Get(Form("zeeJetPtRatio_dataSet%s_etabin%i_mc_stat", zeeJetSampleId.Data(), etabin));
        zeeJetHists[pType][etabin][1][1]->Add(temp);
        temp = (TH2D*)thisFile->Get(Form("zeeJetPtRatio_dataSet%s_etabin%i_mc_stat", zeeJetSampleId.Data(), act_etabin));
        zeeJetHists[2][etabin][1][1]->Add(temp);

        for (short pbin = 0; pbin <= numpbins; pbin++) {
         nJet[pbin][etabin] += (*nJetVec)[pbin + (numpbins+1)*etabin];
         jetEnergyResponseCalib[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("jetEnergyResponseCalib_dataSet%s_pbin%i_etabin%i", zeeJetSampleId.Data(), pbin, etabin)));
         jetEnergyResponseReco[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("jetEnergyResponseReco_dataSet%s_pbin%i_etabin%i", zeeJetSampleId.Data(), pbin, etabin)));
        }
       }

       for (short etabin = 0; etabin <= numetabins; etabin++) {
        zMassSpectra[1][1][etabin]->Add((TH1D*)thisFile->Get(Form("zeeMassSpectrum_dataSet%s_mc_etabin%i", zeeJetSampleId.Data(), etabin)));
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
       const short pType = (zmumuJetSampleId.Contains("pPb") ? 0 : 1);
       nJetVec = (TVectorD*)thisFile->Get(Form("nJetVec_%s", zmumuJetSampleId.Data()));
       nZmumuMassVec = (TVectorD*)thisFile->Get(Form("nZmumuMassVec_%s", zmumuJetSampleId.Data()));
       nZmumuJetVec = (TVectorD*)thisFile->Get(Form("nZmumuJetVec_%s", zmumuJetSampleId.Data()));

       lowResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPhi_dataSet%s", zmumuJetSampleId.Data())));
       highResponseEtaPhi->Add((TH2D*)thisFile->Get(Form("highResponseEtaPhi_dataSet%s", zmumuJetSampleId.Data())));
       lowResponseEtaPt->Add((TH2D*)thisFile->Get(Form("lowResponseEtaPt_dataSet%s", zmumuJetSampleId.Data())));
       highResponseEtaPt->Add((TH2D*)thisFile->Get(Form("highResponseEtaPt_dataSet%s", zmumuJetSampleId.Data())));
       responsePt->Add((TH2D*)thisFile->Get(Form("responsePt_dataSet%s", zmumuJetSampleId.Data())));
       //jetSpectrum->Add((TH1D*)thisFile->Get(Form("jetPt_dataSet%s", zmumuJetSampleId.Data())));

       for (short etabin = 0; etabin <= numetabins; etabin++) {
        const bool flipEta = zmumuJetSampleId.Contains("pPb") && etabin < numetabins;
        const short act_etabin = (flipEta ? numetabins - etabin - 1 : etabin); // period A condition

        nZmumuMass[pType][1][etabin] += (*nZmumuMassVec)[etabin];
        nZmumuMass[2][1][etabin] += (*nZmumuMassVec)[act_etabin];
        nZmumuJet[pType][1][etabin] += (*nZmumuJetVec)[etabin];
        nZmumuJet[2][1][etabin] += (*nZmumuJetVec)[act_etabin];

        // Only add the statistical error plots for MC (don't need to
        // consider systematics)
        TH2D* temp = (TH2D*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%s_etabin%i_mc_stat", zmumuJetSampleId.Data(), etabin));
        zmumuJetHists[pType][etabin][1][1]->Add(temp);
        temp = (TH2D*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%s_etabin%i_mc_stat", zmumuJetSampleId.Data(), act_etabin));
        zmumuJetHists[2][etabin][1][1]->Add(temp);

        for (short pbin = 0; pbin <= numpbins; pbin++) {
         nJet[pbin][etabin] += (*nJetVec)[pbin + (numpbins+1)*etabin];
         jetEnergyResponseCalib[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("jetEnergyResponseCalib_dataSet%s_pbin%i_etabin%i", zmumuJetSampleId.Data(), pbin, etabin)));
         jetEnergyResponseReco[pbin][etabin]->Add((TH1D*)thisFile->Get(Form("jetEnergyResponseReco_dataSet%s_pbin%i_etabin%i", zmumuJetSampleId.Data(), pbin, etabin)));
        }
       }

       for (short etabin = 0; etabin <= numetabins; etabin++) {
        zMassSpectra[0][1][etabin]->Add((TH1D*)thisFile->Get(Form("zmumuMassSpectrum_dataSet%s_mc_etabin%i", zmumuJetSampleId.Data(), etabin)));
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

   zlines[i] = new TLine(pzbins[0], 1.0-2*dz+dz*i, pzbins[numpzbins], 1.0-2*dz+dz*i);
   glines[i] = new TLine(pbins[0], 1.0-1*dg+dg*i, pbins[numpbins], 1.0-1*dg+dg*i);
   xlines[i] = new TLine(xjrefbins[0], 1.0-2*dx+dx*i, xjrefbins[numxjrefbins], 1.0-2*dx+dx*i);
   dplines[i] = new TLine(dpbounds[i], 0.75, dpbounds[i], 2.15);
   dplines_bottom[i] = new TLine(dpbounds[i], 0.91, dpbounds[i], 1.09);

   if (1.0-2*dz+dz*i == 1) zlines[i]->SetLineStyle(1);
   else zlines[i]->SetLineStyle(3);
   if (1.0-1*dg+dg*i == 1) glines[i]->SetLineStyle(1);
   else glines[i]->SetLineStyle(3);
   if (1.0-2*dx+dx*i == 1) xlines[i]->SetLineStyle(1);
   else xlines[i]->SetLineStyle(3);
   dplines[i]->SetLineStyle(3);
   dplines_bottom[i]->SetLineStyle(3);
  }


  /**** Canvas definitions ****/
  TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
  const double padRatio = 1.5; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double dPadY = 1.0/(padRatio+1.0);
  const double uPadY = 1.0 - dPadY;
  TPad* topPad = new TPad("topPad", "", 0, dPadY, 1, 1);
  TPad* bottomPad = new TPad("bottomPad", "", 0, 0, 1, dPadY);
  topPad->SetBottomMargin(0);
  topPad->SetLeftMargin(-0.20);
  bottomPad->SetTopMargin(0);
  bottomPad->SetBottomMargin(0.30);
  bottomPad->SetLeftMargin(-0.20);
  topPad->Draw();
  bottomPad->Draw();

  /**** Define local histograms, graphs, etc. ****/
  TH1D *vJetHist, *vJetHist_mc, *vJetHist_lo, *vJetHist_hi, *vJetHist_rat, *vJetHist_rat_lo, *vJetHist_rat_hi;
  TGraphAsymmErrors *vJetGraph_sys, *vJetGraph_rat_sys;
  TH1D* gJetHistDifference[3][numetabins][3];

  TFile* outFile = new TFile(TString(rootPath) + "cc_difference.root", "recreate");

  for (short pType = 0; pType < 3; pType++) {

   TString period = "Period A";
   if (pType == 1) period = "Period B";
   else if (pType == 2) period = "Period A+B";

   for (short etabin = 0; etabin <= numetabins; etabin++) {

    /**** Plot ZmumuJet info ****/
    topPad->cd();
    topPad->SetLogx();
    vJetHist = GetProfileX("vJetHist", zmumuJetHists[pType][etabin][0][1], numpzbins, pzbins, false);
    vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
    vJetHist->SetYTitle("#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
    vJetHist->SetAxisRange(0.75, 2.15, "Y");
    vJetHist->SetMarkerColor(data_color);
    vJetHist->SetLineColor(data_color);
    vJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleOffset(uPadY);

    // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
    vJetHist_lo = GetProfileX("vJetHist_lo", zmumuJetHists[pType][etabin][0][0], numpzbins, pzbins, false);
    vJetHist_hi = GetProfileX("vJetHist_hi", zmumuJetHists[pType][etabin][0][2], numpzbins, pzbins, false);
    CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    if (vJetHist_lo) delete vJetHist_lo;
    if (vJetHist_hi) delete vJetHist_hi;
    vJetGraph_sys->SetFillColor(kBlack);
    vJetGraph_sys->SetFillStyle(3001);

    vJetHist_mc = GetProfileX("vJetHist_mc", zmumuJetHists[pType][etabin][1][1], numpzbins, pzbins, false);
    vJetHist_mc->SetMarkerColor(mc_color);
    vJetHist_mc->SetLineColor(mc_color);

    vJetHist->DrawCopy("E1 X0");
    vJetHist_mc->DrawCopy("SAME E1 X0");
    vJetGraph_sys->Draw("2");

    myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 #it{p}+Pb 8 TeV, with Insitu Corrections (%i events)", nZmumuJet[pType][0][etabin]), 1.25, 0.04/uPadY);
    myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV with Overlay (%i events)", nZmumuJet[pType][1][etabin]), 1.25, 0.04/uPadY);
    if (etabin < numetabins) {
     if (pType == 2) myText(0.155, 0.1,kBlack, Form("%g < #eta_{Proton}^{#mu#mu} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{Lab}^{#mu#mu} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
    }
    myText(0.155, 0.28, kBlack, "Z (#mu#mu) + Jet", 0.04/uPadY);
    myText(0.155, 0.19, kBlack, period.Data(), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx();

    vJetHist_rat = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_etabin%i", etabin)), zmumuJetHists[pType][etabin][0][1], zmumuJetHists[pType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
    vJetHist_rat_lo = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_lo_etabin%i", etabin)), zmumuJetHists[pType][etabin][0][0], zmumuJetHists[pType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetHist_rat_hi = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_hi_etabin%i", etabin)), zmumuJetHists[pType][etabin][0][2], zmumuJetHists[pType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    CalcSystematics(vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
    if (vJetHist_rat_lo) delete vJetHist_rat_lo;
    if (vJetHist_rat_hi) delete vJetHist_rat_hi;
    vJetGraph_rat_sys->SetFillColor(kBlack);
    vJetGraph_rat_sys->SetFillStyle(3001);

    vJetHist_rat->SetYTitle("Data / MC");
    vJetHist_rat->SetAxisRange(0.85, 1.15, "Y");
    //vJetHist_rat->SetAxisRange(0.75, 1.35, "Y");
    vJetHist_rat->GetYaxis()->SetNdivisions(405);
    vJetHist_rat->GetXaxis()->SetTitleSize(0.04/dPadY);
    vJetHist_rat->GetYaxis()->SetTitleSize(0.04/dPadY);
    vJetHist_rat->GetXaxis()->SetTitleOffset(1);
    vJetHist_rat->GetYaxis()->SetTitleOffset(dPadY);
    vJetHist_rat->GetYaxis()->CenterTitle(true);
    vJetHist_rat->GetXaxis()->SetLabelSize(0.04/dPadY);
    vJetHist_rat->GetYaxis()->SetLabelSize(0.04/dPadY);
    vJetHist_rat->GetXaxis()->SetTickLength(0.08);

    vJetHist_rat->Draw("E1 X0");
    vJetGraph_rat_sys->Draw("2");
    for (TLine* line : zlines) line->Draw("SAME");
    char* plotName;
    if (etabin < numetabins) plotName = Form("z_mumu_jet%i.pdf", etabin);
    else plotName = Form("z_mumu_jet_combined.pdf");

    switch (pType) {
     case 0:
      canvas->SaveAs(Form("%s/PeriodA/%s", plotPath.Data(), plotName));
      break;
     case 1:
      canvas->SaveAs(Form("%s/PeriodB/%s", plotPath.Data(), plotName));
      break;
     case 2:
      canvas->SaveAs(Form("%s/PeriodAB/%s", plotPath.Data(), plotName));
      break;
    }
    if (vJetHist) delete vJetHist;
    if (vJetHist_mc) delete vJetHist_mc;
    if (vJetGraph_sys) delete vJetGraph_sys;
    if (vJetHist_rat) delete vJetHist_rat;
    if (vJetGraph_rat_sys) delete vJetGraph_rat_sys;
    


    /**** Plots ZeeJet info ****/
    topPad->cd();
    topPad->SetLogx();
    vJetHist = GetProfileX("vJetHist", zeeJetHists[pType][etabin][0][1], numpzbins, pzbins, false);
    vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
    vJetHist->SetYTitle("#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
    vJetHist->SetAxisRange(0.75, 2.15, "Y");
    vJetHist->SetMarkerColor(data_color);
    vJetHist->SetLineColor(data_color);
    vJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleOffset(uPadY);

    // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
    vJetHist_lo = GetProfileX("vJetHist_lo", zeeJetHists[pType][etabin][0][0], numpzbins, pzbins, false);
    vJetHist_hi = GetProfileX("vJetHist_hi", zeeJetHists[pType][etabin][0][2], numpzbins, pzbins, false);
    CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    if (vJetHist_lo) delete vJetHist_lo;
    if (vJetHist_hi) delete vJetHist_hi;
    vJetGraph_sys->SetFillColor(kBlack);
    vJetGraph_sys->SetFillStyle(3001);

    vJetHist_mc = GetProfileX("vJetHist_mc", zeeJetHists[pType][etabin][1][1], numpzbins, pzbins, false);
    vJetHist_mc->SetMarkerColor(mc_color);
    vJetHist_mc->SetLineColor(mc_color);

    vJetHist->DrawCopy("E1 X0");
    vJetHist_mc->DrawCopy("SAME E1 X0");
    vJetGraph_sys->Draw("2");

    myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 #it{p}+Pb 8 TeV, with Insitu Corrections (%i events)", nZeeJet[pType][0][etabin]), 1.25, 0.04/uPadY);
    myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV with Overlay (%i events)", nZeeJet[pType][1][etabin]), 1.25, 0.04/uPadY);
    if (etabin < numetabins) {
     if (pType == 2) myText(0.155, 0.1,kBlack, Form("%g < #eta_{Proton}^{ee} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{Lab}^{ee} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
    }
    myText(0.155, 0.28, kBlack, "Z (ee) + Jet", 0.04/uPadY);
    myText(0.155, 0.19, kBlack, period.Data(), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx();

    vJetHist_rat = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_etabin%i", etabin)), zeeJetHists[pType][etabin][0][1], zeeJetHists[pType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
    vJetHist_rat_lo = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_lo_etabin%i", etabin)), zeeJetHists[pType][etabin][0][0], zeeJetHists[pType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetHist_rat_hi = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_hi_etabin%i", etabin)), zeeJetHists[pType][etabin][0][2], zeeJetHists[pType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    CalcSystematics(vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
    if (vJetHist_rat_lo) delete vJetHist_rat_lo;
    if (vJetHist_rat_hi) delete vJetHist_rat_hi;
    vJetGraph_rat_sys->SetFillColor(kBlack);
    vJetGraph_rat_sys->SetFillStyle(3001);

    vJetHist_rat->SetYTitle("Data / MC");
    vJetHist_rat->SetAxisRange(0.85, 1.15, "Y");
    //vJetHist_rat->SetAxisRange(0.75, 1.35, "Y");
    vJetHist_rat->GetYaxis()->SetNdivisions(405);
    vJetHist_rat->GetXaxis()->SetTitleSize(0.04/dPadY);
    vJetHist_rat->GetYaxis()->SetTitleSize(0.04/dPadY);
    vJetHist_rat->GetXaxis()->SetTitleOffset(1);
    vJetHist_rat->GetYaxis()->SetTitleOffset(dPadY);
    vJetHist_rat->GetYaxis()->CenterTitle(true);
    vJetHist_rat->GetXaxis()->SetLabelSize(0.04/dPadY);
    vJetHist_rat->GetYaxis()->SetLabelSize(0.04/dPadY);
    vJetHist_rat->GetXaxis()->SetTickLength(0.08);

    vJetHist_rat->Draw("E1 X0");
    vJetGraph_rat_sys->Draw("2");
    for (TLine* line : zlines) line->Draw("SAME");
    if (etabin < numetabins) plotName = Form("z_ee_jet%i.pdf", etabin);
    else plotName = Form("z_ee_jet_combined.pdf");
    switch (pType) {
     case 0:
      canvas->SaveAs(Form("%s/PeriodA/%s", plotPath.Data(), plotName));
      break;
     case 1:
      canvas->SaveAs(Form("%s/PeriodB/%s", plotPath.Data(), plotName));
      break;
     case 2:
      canvas->SaveAs(Form("%s/PeriodAB/%s", plotPath.Data(), plotName));
      break;
    }
    if (vJetHist) delete vJetHist;
    if (vJetHist_mc) delete vJetHist_mc;
    if (vJetGraph_sys) delete vJetGraph_sys;
    if (vJetHist_rat) delete vJetHist_rat;
    if (vJetGraph_rat_sys) delete vJetGraph_rat_sys;


    /**** Plots GammaJet info as a function of p_T^ref****/
    topPad->cd();
    topPad->SetLogx();
    vJetHist = GetProfileX("vJetHist", gJetHists[pType][etabin][0][1], numpbins, pbins, true);
    vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
    vJetHist->SetYTitle("#it{p}_{T}^{J} / #it{p}_{T}^{ref}>");
    vJetHist->SetAxisRange(0.75, 2.15, "Y");
    vJetHist->SetMarkerColor(data_color);
    vJetHist->SetLineColor(data_color);
    vJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleOffset(uPadY);

    // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
    vJetHist_lo = GetProfileX("vJetHist_lo", gJetHists[pType][etabin][0][0], numpbins, pbins, true);
    vJetHist_hi = GetProfileX("vJetHist_hi", gJetHists[pType][etabin][0][2], numpbins, pbins, true);
    CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    vJetGraph_sys->SetFillColor(kBlack);
    vJetGraph_sys->SetFillStyle(3001);

    vJetHist_mc = GetProfileX("vJetHist_mc", gJetHists[pType][etabin][1][1], numpbins, pbins, true);
    vJetHist_mc->SetMarkerColor(mc_color);
    vJetHist_mc->SetLineColor(mc_color);

    for (short errType = 0; errType < 3; errType++) {
     const TString error = (errType == 0 ? "syslo" : (errType == 1 ? "stat" : "syshi"));
     TString periodStr = "periodA";
     if (pType == 1) periodStr = "periodB";
     else if (pType == 2) periodStr = "periodAB";
     gJetHistDifference[pType][etabin][errType] = new TH1D(Form("gJetPtRatio_diff%i_%s_%s", etabin, error.Data(), periodStr.Data()), ";#it{p}_{T}^{ref} #left[GeV#right]", numpbins, pbins);
     for (short pbin = 1; pbin <= numpbins; pbin++) {
      double dataVal, dataErr;
      switch (errType) {
       case 0:
        dataVal = vJetHist_lo->GetBinContent(pbin);
        dataErr = vJetHist_lo->GetBinError(pbin);
        break;
       case 2:
        dataVal = vJetHist_hi->GetBinContent(pbin);
        dataErr = vJetHist_hi->GetBinError(pbin);
        break;
       default:
        dataVal = vJetHist->GetBinContent(pbin);
        dataErr = vJetHist->GetBinError(pbin);
      } 
      gJetHistDifference[pType][etabin][errType]->SetBinContent(pbin, dataVal - vJetHist_mc->GetBinContent(pbin));
      gJetHistDifference[pType][etabin][errType]->SetBinError(pbin, TMath::Sqrt(TMath::Power(dataErr,2) + TMath::Power(vJetHist_mc->GetBinError(pbin),2)));
     }
     gJetHistDifference[pType][etabin][errType]->Write();
    }
    if (vJetHist_lo) delete vJetHist_lo;
    if (vJetHist_hi) delete vJetHist_hi;

    vJetHist->DrawCopy("E1 X0");
    vJetHist_mc->DrawCopy("SAME E1 X0");
    vJetGraph_sys->Draw("2");
    for (TLine* line : dplines) line->Draw("same");

    myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 #it{p}+Pb 8 TeV, with Insitu Corrections (%i events)", nGammaJet[pType][0][etabin]), 1.25, 0.04/uPadY);
    myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV with Overlay (%i events)", nGammaJet[pType][1][etabin]), 1.25, 0.04/uPadY);
    if (etabin < numetabins) {
     if (pType == 2) myText(0.155, 0.1,kBlack, Form("%g < #eta_{Proton}^{#gamma} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
    }
    myText(0.155, 0.28, kBlack, "#gamma + Jet", 0.04/uPadY);
    myText(0.155, 0.19, kBlack, period.Data(), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx();
    vJetHist_rat = GetDataOverMC(TString(Form("gJetPtDataMCRatio_etabin%i", etabin)), gJetHists[pType][etabin][0][1], gJetHists[pType][etabin][1][1], numpbins, pbins, numxjrefbins, xjrefbins, true);
    vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
    vJetHist_rat_lo = GetDataOverMC(TString(Form("gJetPtDataMCRatio_lo_etabin%i", etabin)), gJetHists[pType][etabin][0][0], gJetHists[pType][etabin][1][1], numpbins, pbins, numxjrefbins, xjrefbins, true);
    vJetHist_rat_hi = GetDataOverMC(TString(Form("gJetPtDataMCRatio_hi_etabin%i", etabin)), gJetHists[pType][etabin][0][2], gJetHists[pType][etabin][1][1], numpbins, pbins, numxjrefbins, xjrefbins, true);
    CalcSystematics(vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
    if (vJetHist_rat_lo) delete vJetHist_rat_lo;
    if (vJetHist_rat_hi) delete vJetHist_rat_hi;
    vJetGraph_rat_sys->SetFillColor(kBlack);
    vJetGraph_rat_sys->SetFillStyle(3001);

    vJetHist_rat->SetYTitle("Data / MC");
    vJetHist_rat->SetAxisRange(0.91, 1.09, "Y");
    vJetHist_rat->GetYaxis()->SetNdivisions(405);
    vJetHist_rat->GetXaxis()->SetTitleSize(0.04/dPadY);
    vJetHist_rat->GetYaxis()->SetTitleSize(0.04/dPadY);
    vJetHist_rat->GetXaxis()->SetTitleOffset(1);
    vJetHist_rat->GetYaxis()->SetTitleOffset(dPadY);
    vJetHist_rat->GetYaxis()->CenterTitle(true);
    vJetHist_rat->GetXaxis()->SetLabelSize(0.04/dPadY);
    vJetHist_rat->GetYaxis()->SetLabelSize(0.04/dPadY);
    vJetHist_rat->GetXaxis()->SetTickLength(0.08);

    vJetHist_rat->Draw("e1 X0"); 
    vJetGraph_rat_sys->Draw("2");
    for (TLine* line : glines) line->Draw("same");
    for (TLine* line : dplines_bottom) line->Draw("same");

    if (etabin < numetabins) plotName = Form("gamma_jet%i.pdf", etabin);
    else plotName = Form("gamma_jet_combined.pdf");
    switch (pType) {
     case 0:
      canvas->SaveAs(Form("%s/PeriodA/%s", plotPath.Data(), plotName));
      break;
     case 1:
      canvas->SaveAs(Form("%s/PeriodB/%s", plotPath.Data(), plotName));
      break;
     case 2:
      canvas->SaveAs(Form("%s/PeriodAB/%s", plotPath.Data(), plotName));
      break;
    }
    if (vJetHist) delete vJetHist;
    if (vJetHist_mc) delete vJetHist_mc;
    if (vJetGraph_sys) delete vJetGraph_sys;
    if (vJetHist_rat) delete vJetHist_rat;
    if (vJetGraph_rat_sys) delete vJetGraph_rat_sys;
    for (short errType = 0; errType < 3; errType++)
     if (gJetHistDifference[pType][etabin][errType])
      delete gJetHistDifference[pType][etabin][errType];

    if (etabin == numetabins) continue;


    if (!plot_xjref) continue;
    /**** Plots xjref distributions, binned by ptref ****/
    for (int pbin = 0; pbin < numpbins; pbin++) {
     const double pref_lo = pbins[pbin];
     const double pref_hi =  pbins[pbin+1];
     topPad->cd();
     topPad->SetLogx(0);
     vJetHist = gJetHists[pType][etabin][0][1]->ProjectionY("vJetProjection", pbin, pbin);
     const float counts_data = vJetHist->Integral();
     const float total_data = gJetHists[pType][etabin][0][1]->Integral();
     vJetHist->Rebin(rebinFactor);
     vJetHist->Scale(1./vJetHist->Integral());
//     vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
     vJetHist->SetXTitle("#it{x}_{J}^{ref}");
     vJetHist->SetYTitle("Fractional counts");
     vJetHist->SetMarkerColor(data_color);
     vJetHist->SetLineColor(data_color);
     //vJetHist->GetXaxis()->SetRangeUser(0., 2.);
     vJetHist->GetYaxis()->SetRangeUser(0., 0.6);//vJetHist->GetYaxis()->GetXmax());
     
     vJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetTitleOffset(1.1*uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     vJetHist_lo = gJetHists[pType][etabin][0][0]->ProjectionY("vJetProjection_lo", pbin, pbin);
     vJetHist_lo->Rebin(rebinFactor);
     vJetHist_lo->Scale(1./vJetHist_lo->Integral()); 
     //vJetHist_lo->Scale(1./counts_data); 
     vJetHist_hi = gJetHists[pType][etabin][0][2]->ProjectionY("vJetProjection_hi", pbin, pbin);
     vJetHist_hi->Rebin(rebinFactor);
     vJetHist_hi->Scale(1./vJetHist_hi->Integral()); 
     //vJetHist_hi->Scale(1./counts_data); 
//     CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
//     vJetGraph_sys->SetFillColor(kBlack);
//     vJetGraph_sys->SetFillStyle(3001);

     vJetHist_lo->SetMarkerColor(data_color);
     vJetHist_lo->SetLineColor(data_color);
     vJetHist_hi->SetLineColor(data_color);
     vJetHist_hi->SetLineColor(data_color);

     vJetHist_mc = gJetHists[pType][etabin][1][1]->ProjectionY("vJetProjection_mc", pbin, pbin);
     const float counts_mc = vJetHist_mc->Integral();
     const float total_mc = gJetHists[pType][etabin][1][1]->Integral();
     vJetHist_mc->Rebin(rebinFactor);
     vJetHist_mc->Scale(1./vJetHist_mc->Integral()); 
     vJetHist_mc->SetMarkerColor(mc_color);
     vJetHist_mc->SetLineColor(mc_color);

     vJetHist->DrawCopy("E1 X0");
     vJetHist_mc->DrawCopy("SAME E1 X0");
//     vJetGraph_sys->Draw("2");

//     float n = 100. * counts_data / total_data;
     myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 #it{p}+Pb 8 TeV, with Insitu Corrections (%i events)", (int)counts_data), 1.25, 0.04/uPadY);

//     n = 100. * counts_mc / total_mc;
     myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV %s (%i events)", (runValidation ? "":"with Overlay"), (int)counts_mc), 1.25, 0.04/uPadY);

     float mean, mean_err, mean_mc, mean_mc_err, mean_lo, mean_hi;
   
     if (useGaussian) {
      TF1* gaus_data = new TF1("gaus_data", "gaus(0)", 0, 4.0);
      vJetHist->Fit(gaus_data, "Q0R");
      TF1* gaus_mc = new TF1("gaus_mc", "gaus(0)", 0, 4.0);
      vJetHist_mc->Fit(gaus_mc, "Q0R");
      TF1* gaus_data_lo = new TF1("gaus_data_lo", "gaus(0)", 0, 4.0);
      vJetHist_lo->Fit(gaus_data_lo, "Q0R");
      TF1* gaus_data_hi = new TF1("gaus_data_hi", "gaus(0)", 0, 4.0);
      vJetHist_hi->Fit(gaus_data_hi, "Q0R");
      mean = gaus_data->GetParameter(1);
      mean_err = gaus_data->GetParError(1);
      mean_mc = gaus_mc->GetParameter(1);
      mean_mc_err = gaus_mc->GetParError(1);
      mean_lo = gaus_data_lo->GetParameter(1);
      mean_hi = gaus_data_hi->GetParameter(1);
      if (gaus_data) delete gaus_data;
      if (gaus_mc) delete gaus_mc;
      if (gaus_data_lo) delete gaus_data_lo;
      if (gaus_data_hi) delete gaus_data_hi;
     }
     else {
      mean = vJetHist->GetMean();
      mean_err = vJetHist->GetMeanError();
      mean_mc = vJetHist_mc->GetMean();
      mean_mc_err = vJetHist_mc->GetMeanError();
      mean_lo = vJetHist_lo->GetMean();
      mean_hi = vJetHist_hi->GetMean();
     }

     const float sys_err = 0.5* (TMath::Abs(mean_lo-mean) + TMath::Abs(mean_hi-mean));

     myText(0.155, 0.73, kBlack, Form("<#it{x}_{J}^{ref}>^{data} = %.2f #pm %.2f #pm %.2f", mean, mean_err, sys_err), 0.04/uPadY);
     myText(0.155, 0.65, kBlack, Form("<#it{x}_{J}^{ref}>^{MC} = %.2f #pm %.2f", mean_mc, mean_mc_err), 0.04/uPadY);

     myText(0.155, 0.43, kBlack, "#gamma + Jet", 0.04/uPadY);
     myText(0.155, 0.34, kBlack, Form("%g < #it{p}_{T}^{ref} < %g", pref_lo, pref_hi), 0.04/uPadY);
     myText(0.155, 0.25, kBlack, period.Data(), 0.04/uPadY);
     if (pType == 2) myText(0.155, 0.16,kBlack, Form("%g < #eta_{Proton}^{#gamma} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.155, 0.16,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);

     bottomPad->cd();
     bottomPad->SetLogx(0);
     vJetHist->Divide(vJetHist_mc);
     vJetHist_lo->Divide(vJetHist_mc);
     vJetHist_hi->Divide(vJetHist_mc);

//     vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist);
//     CalcSystematics(vJetGraph_rat_sys, vJetHist, vJetHist_hi, vJetHist_lo);
//     vJetGraph_rat_sys->SetFillColor(kBlack);
//     vJetGraph_rat_sys->SetFillStyle(3001);

     vJetHist->SetYTitle("Data / MC");
     vJetHist->SetAxisRange(0.45, 1.65, "Y");
     vJetHist->GetYaxis()->SetNdivisions(605);
     vJetHist->GetXaxis()->SetTitleSize(0.04/dPadY);
     vJetHist->GetYaxis()->SetTitleSize(0.04/dPadY);
     vJetHist->GetXaxis()->SetTitleOffset(1);
     vJetHist->GetYaxis()->SetTitleOffset(1.1*dPadY);
     vJetHist->GetYaxis()->CenterTitle(true);
     vJetHist->GetXaxis()->SetLabelSize(0.04/dPadY);
     vJetHist->GetYaxis()->SetLabelSize(0.04/dPadY);
     vJetHist->GetXaxis()->SetTickLength(0.08);

     vJetHist->Draw("E1 X0"); 
//     vJetGraph_rat_sys->Draw("2");
     for (TLine* line : xlines) line->Draw("SAME");
     plotName = Form("pref_slices/gamma_jet%i_pbin%i.pdf", etabin, pbin);
     switch (pType) {
      case 0:
       canvas->SaveAs(Form("%s/PeriodA/%s", plotPath.Data(), plotName));
       break;
      case 1:
       canvas->SaveAs(Form("%s/PeriodB/%s", plotPath.Data(), plotName));
       break;
      case 2:
       canvas->SaveAs(Form("%s/PeriodAB/%s", plotPath.Data(), plotName));
       break;
     }
     if (vJetHist) delete vJetHist;
     if (vJetHist_mc) delete vJetHist_mc;
     if (vJetHist_lo) delete vJetHist_lo;
     if (vJetHist_hi) delete vJetHist_hi;
//     if (vJetGraph_sys) delete vJetGraph_sys;
//     if (vJetGraph_rat_sys) delete vJetGraph_rat_sys;
      
    }
   }
  }


  /**** Plots systematic errors vs jet pt ****/
  for (short pType = 0; pType < 2; pType++) {
   const TString period = (pType == 0 ? "Period A":"Period B");

   for (short etabin = 0; etabin < numetabins; etabin++) {
    topPad->cd();
    TH2D* thisHist = gJetHistsSys[pType][etabin][0][1];
    TH1D* rmsHist = new TH1D(Form("rms_etabin%i_%s", etabin, (pType==0?"pPb":"Pbp")), "", numpzbins, pzbins);
    for (short pzbin = 0; pzbin < numpzbins; pzbin++) {
     float rms = 0;
     float sumWeights = 0;
     for (short sigbin = 0; sigbin < numSigmaBins; sigbin++) {
      const float sig = thisHist->GetYaxis()->GetBinCenter(sigbin+1);
      const float weight = thisHist->GetBinContent(pzbin+1, sigbin+1);
      rms += pow(sig, 2) * weight;
      sumWeights += weight;
     }
     if (sumWeights > 0) rms = sqrt(rms) / sqrt(sumWeights);
     rmsHist->SetBinContent(pzbin+1, rms);
    }
    topPad->SetLogz();
    thisHist->Draw("col");
    thisHist->GetXaxis()->SetLabelSize(0.04/uPadY);
    thisHist->GetYaxis()->SetLabelSize(0.04/uPadY);
    thisHist->GetYaxis()->SetTitleSize(0.04/uPadY);
    thisHist->GetYaxis()->SetTitleOffset(1.1*uPadY);
    myText(0.72, 0.89, kBlack, period.Data(), 0.04/uPadY);
    myText(0.72, 0.8,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);

    bottomPad->cd();
    rmsHist->SetXTitle("#it{p}_{T}^{jet} #left[GeV#right]");
    //rmsHist->SetYTitle("RMS(#sigma)/#it{p}_{T}^{jet}");
    rmsHist->SetYTitle("RMS");
    rmsHist->SetAxisRange(0, 0.09, "Y");
    rmsHist->GetYaxis()->SetNdivisions(405);
    rmsHist->GetXaxis()->SetTitleSize(0.04/dPadY);
    rmsHist->GetYaxis()->SetTitleSize(0.04/dPadY);
    rmsHist->GetXaxis()->SetTitleOffset(1);
    rmsHist->GetYaxis()->SetTitleOffset(1.1*dPadY);
    rmsHist->GetYaxis()->CenterTitle(true);
    rmsHist->GetXaxis()->SetLabelSize(0.04/dPadY);
    rmsHist->GetYaxis()->SetLabelSize(0.04/dPadY);
    rmsHist->GetXaxis()->SetTickLength(0.08);
    rmsHist->Draw("hist");
    canvas->SaveAs(Form("%s/Period%s/jetSystematics_etabin%i.pdf", plotPath.Data(), (pType==0 ? "A":"B"), etabin));
    if (rmsHist) delete rmsHist;
   }
  }
  outFile->Write();
  if (outFile) delete outFile;


  /**** create new lines for Z mass spectra ****/
  TLine* lines[5] = {};
  for (short i = 0; i < 5; i++) {
   lines[i] = new TLine(60, 0.6+0.2*i, 110, 0.6+0.2*i);
   if (0.6+0.2*i == 1) lines[i]->SetLineStyle(1);
   else lines[i]->SetLineStyle(3);
  }


  /**** Plot mumu mass spectra ****/
  for (int etabin = 0; etabin <= numetabins; etabin++) {
   if (etabin != numetabins) continue;
   for (int species = 0; species < 2; species++) {
    topPad->cd();
    topPad->SetLogx(0);
    double mean[2] = {};
    double mean_err[2] = {};
    double sigma[2] = {};
    double sigma_err[2] = {};
    TF1* fits[2];
    //RooZfit* fits[2];
    for (short dType = 0; dType < 2; dType++) {
     TH1D* thisHist = zMassSpectra[species][dType][etabin];
     Color_t color = (dType==0 ? data_color : mc_color);
     thisHist->GetXaxis()->SetTitle("#font[12]{ll} Invariant Mass #left[GeV#right]");
     thisHist->GetYaxis()->SetTitle("Normalized Counts / 1 GeV");
     thisHist->GetYaxis()->SetTitleSize(0.04/uPadY);
     thisHist->GetYaxis()->SetTitleOffset(1.1*uPadY);
     thisHist->GetYaxis()->SetLabelSize(0.04/uPadY);
     thisHist->SetLineColor(color);
     thisHist->SetMarkerColor(color);
//     thisHist->GetXaxis()->SetNdivisions(50802, false);
     thisHist->Scale(1.0, "width");

     //TF1* gausFit = new TF1(Form("fit_guess_species%i_dType%i", species, dType), "gaus(0)", Z_mass - 10, Z_mass + 10);
     //thisHist->Fit(gausFit, "R", "L");
     //double m = gausFit->GetParameter(1);
     //double s = gausFit->GetParameter(2);
     //const double scale = 1.0 / thisHist->Integral(thisHist->FindBin(m-Z_mass_fitNsigma*s), thisHist->FindBin(m+Z_mass_fitNsigma*s));
     //thisHist->Scale(scale);
     //
     //RooZfit* fit = new RooZfit(thisHist, color);
     //fits[dType] = fit;

     //mean[dType] = fit->mean->getVal();
     //mean_err[dType] = fit->mean->getError();
     //sigma[dType] = fit->sigma->getVal();
     //sigma_err[dType] = fit->sigma->getError();

     TF1* fit = new TF1(Form("fit_guess_species%i_dType%i", species, dType), "gaus(0)", Z_mass - 10, Z_mass + 10);
     thisHist->Fit(fit, "R", "L");
     double m = fit->GetParameter(1);
     double s = fit->GetParameter(2);

     TF1* fit_better = new TF1(Form("fit_better_species%i_dType%i", species, dType), "gaus(0)", m-Z_mass_fitNsigma*s, m+Z_mass_fitNsigma*s);
     thisHist->Fit(fit_better, "R", "L");
     m = fit_better->GetParameter(1);
     s = fit_better->GetParameter(2);
     mean[dType] = m;
     mean_err[dType] = fit_better->GetParError(1);
     sigma[dType] = s;
     sigma_err[dType] = fit_better->GetParError(2);
     const double scale = 1.0 / thisHist->Integral(thisHist->FindBin(m-Z_mass_fitNsigma*s), thisHist->FindBin(m+Z_mass_fitNsigma*s));
     thisHist->Scale(scale);

     fits[dType] = fit_better;
     fit_better->SetLineColor(color);
     fit_better->SetParameter(0, scale*fit_better->GetParameter(0));

     thisHist->SetAxisRange(0, 0.20, "Y");
     thisHist->GetYaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    }
    for (short dType = 0; dType < 2; dType++) {
     TH1D* thisHist = zMassSpectra[species][dType][etabin];
     if (dType == 0) thisHist->Draw("p");
     else thisHist->Draw("hist same");
     //fits[dType]->plot->Draw("same");
     fits[dType]->Draw("same");
     if (etabin < numetabins) {
      if (species == 0)
       myText(0.175, 0.15,kBlack, Form("%g < #eta_{Lab}^{#mu#mu} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
      else if (species == 1)
       myText(0.175, 0.15,kBlack, Form("%g < #eta_{Lab}^{ee} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     }
    }
    if (species == 0) {
     myText(0.175, 0.88, kBlack, "Z (#mu#mu) + Jet", 0.04/uPadY);
     myMarkerText(0.175, 0.80, data_color, kFullCircle, Form("2016 #it{p}+{Pb} 8 TeV (%i events)", nZmumuMass[2][0][etabin]), 1.25, 0.04/uPadY);
     myMarkerText(0.175, 0.55, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV with Overlay (%i events)", nZmumuMass[2][1][etabin]), 1.25, 0.04/uPadY);
    }
    else if (species == 1) {
     myText(0.175, 0.88, kBlack, "Z (ee) + Jet", 0.04/uPadY);
     myMarkerText(0.175, 0.80, data_color, kFullCircle, Form("2016 #it{p}+{Pb} 8 TeV (%i events)", nZeeMass[2][0][etabin]), 1.25, 0.04/uPadY);
     myMarkerText(0.175, 0.55, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV with Overlay (%i events)", nZeeMass[2][1][etabin]), 1.25, 0.04/uPadY);
    }
    myText(0.175, 0.72, kBlack, Form("m_{Z}^{data} = %.2f #pm %.2f GeV", mean[0], mean_err[0]), 0.04/uPadY);
    myText(0.175, 0.64, kBlack, Form("#sigma_{Z}^{data} = %.2f #pm %.2f GeV", sigma[0], sigma_err[0]), 0.04/uPadY);
    myText(0.175, 0.47, kBlack, Form("m_{Z}^{mc} = %.2f #pm %.2f GeV", mean[1], mean_err[1]), 0.04/uPadY);
    myText(0.175, 0.39, kBlack, Form("#sigma_{Z}^{mc} = %.2f #pm %.2f GeV", sigma[1], sigma_err[1]), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx(0);
    TH1D* thisHist = (TH1D*)zMassSpectra[species][0][etabin]->Clone(Form("invMass_species%i_clone", species));
    thisHist->Divide(zMassSpectra[species][1][etabin]);
    thisHist->GetXaxis()->SetTitle("#font[12]{ll} Invariant Mass #left[GeV#right]");
    thisHist->GetYaxis()->SetTitle("Data / MC");
    thisHist->GetXaxis()->SetTitleSize(0.04/dPadY);
    thisHist->GetYaxis()->SetTitleSize(0.04/dPadY);
    thisHist->GetXaxis()->SetTitleOffset(1);
    thisHist->GetYaxis()->SetTitleOffset(1.1*dPadY);
    thisHist->GetYaxis()->CenterTitle(true);
    thisHist->GetXaxis()->SetLabelSize(0.04/dPadY);
    thisHist->GetYaxis()->SetLabelSize(0.04/dPadY);
    thisHist->SetLineColor(kBlack);
    thisHist->SetMarkerColor(kBlack);
    //thisHist->SetAxisRange(0.6, 1.6, "Y");
    thisHist->SetAxisRange(0.2, 2.2, "Y");

    thisHist->GetYaxis()->ChangeLabel(-1, -1, -1, -1, -1, -1, " ");
    thisHist->GetXaxis()->SetTickLength(0.08);
//    thisHist->GetXaxis()->SetNdivisions(50802, false);
    thisHist->GetYaxis()->SetNdivisions(405, false);
    thisHist->Draw("p");
    for (TLine* line : lines) line->Draw("same");

    if (etabin != numetabins) {
     if (species == 0) canvas->SaveAs(Form("%s/zmumu_mass_comparison_etabin%i.pdf", plotPath.Data(), etabin));
     else if (species == 1) canvas->SaveAs(Form("%s/zee_mass_comparison_etabin%i.pdf", plotPath.Data(), etabin));
    }
    else {
     if (species == 0) canvas->SaveAs(Form("%s/zmumu_mass_comparison.pdf", plotPath.Data()));
     else if (species == 1) canvas->SaveAs(Form("%s/zee_mass_comparison.pdf", plotPath.Data()));
    }

    for (short fit = 0; fit < 2; fit++) if (fits[fit]) delete fits[fit];
   }
  }


  /**** Plot 2d histograms ****/
  TCanvas* th2canvas = new TCanvas("th2canvas", "", 800, 600);
  const double padRatio_th2 = 1; // ratio of size of upper pad to lower pad. Used to scale plots and font sizes equally.
  const double rPadX = 1.0/(padRatio_th2+1.0);
  const double lPadX = 1.0 - rPadX;
  TPad* leftPad = new TPad("leftPad", "", 0, 0, lPadX, 1);
  TPad* rightPad = new TPad("rightPad", "", lPadX, 0, 1, 1);
  leftPad->SetRightMargin(0);
  rightPad->SetLeftMargin(0);
  leftPad->SetTopMargin(-0.05);
  rightPad->SetTopMargin(-0.05);
  leftPad->SetBottomMargin(-0.10);
  rightPad->SetBottomMargin(-0.10);
  leftPad->SetLeftMargin(-0.20);
  rightPad->SetRightMargin(-0.10);
  leftPad->Draw();
  rightPad->Draw();
  for (short pType = 0; pType < 3; pType++) {
   const TString period = (pType == 0 ? "Period A": (pType == 1 ? "Period B":"Period A+B"));

   for (short etabin = 0; etabin < numetabins; etabin++) {
    TH2D* dataHist = gJetHists[pType][etabin][0][1];
    TH2D* mcHist = gJetHists[pType][etabin][1][1];
    dataHist->Scale(1./dataHist->Integral());
    mcHist->Scale(1./mcHist->Integral());

    rightPad->cd();
    rightPad->SetLogx();
    rightPad->SetLogz();
    dataHist->GetXaxis()->SetTitleSize(0.02/rPadX);
    dataHist->GetYaxis()->SetTitleSize(0.02/rPadX);
    dataHist->GetXaxis()->SetTitleOffset(1);
    dataHist->GetYaxis()->SetTitleOffset(1);
    dataHist->GetXaxis()->SetLabelSize(0.02/rPadX);
    dataHist->GetYaxis()->SetLabelSize(0.02/rPadX);
    dataHist->Draw("col");
    myText(0.1, 0.15, kBlack, Form("2016 #it{p}+Pb 8 TeV, with Insitu Corrections (%i events)", nGammaJet[pType][0][etabin]), 0.02/rPadX);
    if (pType != 2) myText(0.6, 0.85,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[etabin], etabins[etabin+1]), 0.02/rPadX);
    else myText(0.6, 0.85,kBlack, Form("%g < #eta_{Proton}^{#gamma} < %g", etabins[etabin], etabins[etabin+1]), 0.02/rPadX);
    myText(0.6, 0.8,kBlack, period.Data(), 0.02/rPadX);

    leftPad->cd();
    leftPad->SetLogx();
    leftPad->SetLogz();
    mcHist->GetXaxis()->SetTitleSize(0.02/lPadX);
    mcHist->GetYaxis()->SetTitleSize(0.02/lPadX);
    mcHist->GetXaxis()->SetTitleOffset(1);
    mcHist->GetYaxis()->SetTitleOffset(1);
    mcHist->GetXaxis()->SetLabelSize(0.02/lPadX);
    mcHist->GetYaxis()->SetLabelSize(0.02/lPadX);
    mcHist->Draw("col");
    myText(0.2, 0.15, kBlack, Form("Pythia8 #it{pp} 8 TeV %s (%i events)", (runValidation ? "":"with Overlay"), nGammaJet[pType][1][etabin]), 0.02/lPadX);

    const char* plotName = Form("gamma_jet%i_th2.pdf", etabin);
    switch (pType) {
     case 0:
      th2canvas->SaveAs(Form("%s/PeriodA/%s", plotPath.Data(), plotName));
      break;
     case 1:
      th2canvas->SaveAs(Form("%s/PeriodB/%s", plotPath.Data(), plotName));
      break;
     case 2:
      th2canvas->SaveAs(Form("%s/PeriodAB/%s", plotPath.Data(), plotName));
      break;
    }
   }
  }


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
