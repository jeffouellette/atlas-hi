#include "../Params.C"
#include "../../Initialization.C"
#include "../RooZfit.C"

//using namespace RooFit;

TH1D* GetProfileX(const TString name, TH2F* hist, const int nbinsx, const double* xbins, const bool useFit) {
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

TH1D* GetDataOverMC(const TString name, TH2F* data, TH2F* mc, const int numxbins, const double* xbins, const int numybins, const double* ybins, const bool useFit) {
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
  setupDirectories("", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers(0);
  for (short i = 0; i < sizeof(full_run_list)/sizeof(full_run_list[0]); i++) runNumbers.push_back(full_run_list[i]);
  vector<TString> gammaJetSampleIds(0);
  for (short i = 0; i < 6; i++) {
   gammaJetSampleIds.push_back(TString("Pbp") + (runValidation ? "_Valid":"_Overlay") + "_GammaJet_Slice" + to_string(i+1));
   gammaJetSampleIds.push_back(TString("pPb") + (runValidation ? "_Valid":"_Overlay") + "_GammaJet_Slice" + to_string(i+1));
  }
  vector<TString> zeeJetSampleIds(0);
  for (short i = 0; i < 6; i++) {
   zeeJetSampleIds.push_back(string("Pbp_ZeeJet") + to_string(i));
   zeeJetSampleIds.push_back(string("pPb_ZeeJet") + to_string(i));
  }
  vector<TString> zmumuJetSampleIds(0);
  zmumuJetSampleIds.push_back("Pbp_ZmumuJet");
  zmumuJetSampleIds.push_back("pPb_ZmumuJet");

  TH2F* zeeJetHists[3][numetabins+1][2][3];
  //TH2F* zeeJetHistsSys[3][numetabins][2][3];
  TH2F* zmumuJetHists[3][numetabins+1][2][3];
  //TH2F* zmumuJetHistsSys[3][numetabins][2][3];
  TH2F* gJetHists[3][numetabins+1][2][3];
  TH2F* gJetHistsSys[3][numetabins+1][2][3];
  TH1F* zMassSpectra[2][2][numetabins+1];
  TH1F* zMassSpectra_AllSigns[2][2][numetabins+1];

  for (short etabin = 0; etabin <= numetabins; etabin++) {

   for (short dType = 0; dType < 2; dType++) { // dType is 0 for data, 1 for MC
    string dataType = "data";
    if (dType == 1) dataType = "mc";

    for (short errType = 0; errType < 3; errType++) {
     string error = "sys_lo";
     if (errType == 1) error = "stat";
     else if (errType == 2) error = "sys_hi";

     for (short periodType = 0; periodType < 3; periodType++) {
      string period = "periodA";
      if (periodType == 1) period = "periodB";
      else if (periodType == 2) period = "periodAB";

      zeeJetHists[periodType][etabin][dType][errType] = new TH2F(Form("zeeJetPtRatio_hist%i_%s_%s_%s", etabin, dataType.c_str(), error.c_str(), period.c_str()), ";#it{p}_{T}^{ref} #left[GeV#right];#it{p}_{T}^{j} / #it{p}_{T}^{ref}", numpzbins, pzbins, numxjrefbins, xjrefbins);
      zeeJetHists[periodType][etabin][dType][errType]->Sumw2();
      zmumuJetHists[periodType][etabin][dType][errType] = new TH2F(Form("zmumuJetPtRatio_hist%i_%s_%s_%s", etabin, dataType.c_str(), error.c_str(), period.c_str()), ";#it{p}_{T}^{ref} #left[GeV#right];#it{p}_{T}^{j} / #it{p}_{T}^{ref}", numpzbins, pzbins, numxjrefbins, xjrefbins);
      zmumuJetHists[periodType][etabin][dType][errType]->Sumw2();
      gJetHists[periodType][etabin][dType][errType] = new TH2F(Form("gJetPtRatio_hist%i_%s_%s_%s", etabin, dataType.c_str(), error.c_str(), period.c_str()), ";#it{p}_{T}^{ref} #left[GeV#right];#it{p}_{T}^{j} / #it{p}_{T}^{ref}", numpgammabins, pgammabins, numxjrefbins, xjrefbins);
      gJetHists[periodType][etabin][dType][errType]->Sumw2();

      gJetHistsSys[periodType][etabin][dType][errType] = new TH2F(Form("gJetPtRatioSys_hist%i_%s_%s_%s", etabin, dataType.c_str(), error.c_str(), period.c_str()), ";#it{p}_{T}^{jet} #left[GeV#right];#Delta#it{x}_{J}^{ref}#it{p}_{T}^{ref}/#it{p}_{T}^{jet}", numpzbins, pzbins, numSigmaBins, -maxSigma, maxSigma);
      gJetHistsSys[periodType][etabin][dType][errType]->Sumw2();
     }
    }

    for (short spcType = 0; spcType < 2; spcType++) {
     string species = "ee";
     if (spcType == 0) species = "mumu";

     zMassSpectra[spcType][dType][etabin] = new TH1F(Form("z%sMassSpectrum_%s_etabin%i", species.c_str(), dataType.c_str(), etabin), "", 50, 60, 110);
     zMassSpectra[spcType][dType][etabin]->Sumw2();
     zMassSpectra_AllSigns[spcType][dType][etabin] = new TH1F(Form("z%sMassSpectrum_AllSigns_%s_etabin%i", species.c_str(), dataType.c_str(), etabin), "", 50, 60, 110);
     zMassSpectra_AllSigns[spcType][dType][etabin]->Sumw2();
    }
   }
  }

  int Zee_n[3][2][numetabins+1] = {{{}, {}}, {{}, {}}, {{}, {}}};
  int Zmumu_n[3][2][numetabins+1] = {{{}, {}}, {{}, {}}, {{}, {}}};
  int g_n[3][2][numetabins+1] = {{{}, {}}, {{}, {}}, {{}, {}}};

  {
   TSystemDirectory dir(rootPath.c_str(), rootPath.c_str());
   TList* sysfiles = dir.GetListOfFiles();
   if (!sysfiles) {
    cout << "Cannot get list of files! Exiting." << endl;
    return;
   }
   TSystemFile *sysfile;
   TString fname;
   TString histName;
   TIter next(sysfiles);
   TVectorD* infoVec;
   int numFiles = 0;
   while ((sysfile=(TSystemFile*)next())) {
    fname = sysfile->GetName();
    if (!sysfile->IsDirectory() && fname.EndsWith(".root")) {
     if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data() << endl;
     for (int runNumber : runNumbers) { // check for data
      if (fname.Contains(to_string(runNumber))) {
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%i", runNumber));

       for (short etabin = 0; etabin < numetabins; etabin++) {
        const short periodType = (runNumber < 313500 ? 0 : 1);
        const short act_etabin = (runNumber < 313500 ? (numetabins - etabin - 1) : etabin);

        Zee_n[periodType][0][etabin] += (*infoVec)[2+etabin];
        Zee_n[2][0][act_etabin] += (*infoVec)[2+etabin];
        Zmumu_n[periodType][0][etabin] += (*infoVec)[2+numetabins+etabin];
        Zmumu_n[2][0][act_etabin] += (*infoVec)[2+numetabins+etabin];
        g_n[periodType][0][etabin] += (*infoVec)[2+2*numetabins+etabin];
        g_n[2][0][act_etabin] += (*infoVec)[2+2*numetabins+etabin];

        for (short errType = 0; errType < 3; errType++) {
         string error = "sys_lo";
         if (errType == 1) error = "stat";
         else if (errType == 2) error = "sys_hi";

         TH2F* temp = (TH2F*)thisFile->Get(Form("zeeJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, etabin, error.c_str()));
         zeeJetHists[periodType][etabin][0][errType]->Add(temp);
         zeeJetHists[periodType][numetabins][0][errType]->Add(temp);
         temp = (TH2F*)thisFile->Get(Form("zeeJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, act_etabin, error.c_str()));
         zeeJetHists[2][etabin][0][errType]->Add(temp);
         zeeJetHists[2][numetabins][0][errType]->Add(temp);

         temp = (TH2F*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, etabin, error.c_str()));
         zmumuJetHists[periodType][etabin][0][errType]->Add(temp);
         zmumuJetHists[periodType][numetabins][0][errType]->Add(temp);
         temp = (TH2F*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, act_etabin, error.c_str()));
         zmumuJetHists[2][etabin][0][errType]->Add(temp);
         zmumuJetHists[2][numetabins][0][errType]->Add(temp);

         temp = (TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, etabin, error.c_str()));
         gJetHists[periodType][etabin][0][errType]->Add(temp);
         gJetHists[periodType][numetabins][0][errType]->Add(temp);
         temp = (TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, act_etabin, error.c_str()));
         gJetHists[2][etabin][0][errType]->Add(temp);
         gJetHists[2][numetabins][0][errType]->Add(temp);

         //zmumuJetHists[periodType][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, etabin, error.c_str())));
         //zmumuJetHists[2][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, act_etabin, error.c_str())));
         //gJetHists[periodType][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, etabin, error.c_str())));
         //gJetHists[2][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, act_etabin, error.c_str())));
         if (errType == 1) 
          gJetHistsSys[periodType][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("gJetPtRatioSys_dataSet%i_hist%i_data_%s", runNumber, etabin, error.c_str())));
         //gJetHistsSys[2][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("gJetPtRatioSys_dataSet%i_hist%i_data_%s", runNumber, act_etabin, error.c_str())));
        }
       }

       for (short spcType = 0; spcType < 2; spcType++) {
        string species = "ee";
        if (spcType == 0) species = "mumu";

        for (short etabin = 0; etabin < numetabins; etabin++) {
         zMassSpectra[spcType][0][etabin]->Add((TH1F*)thisFile->Get(Form("z%sMassSpectrum_dataSet%i_data_etabin%i", species.c_str(), runNumber, etabin)));
         zMassSpectra_AllSigns[spcType][0][etabin]->Add((TH1F*)thisFile->Get(Form("z%sMassSpectrum_AllSigns_dataSet%i_data_etabin%i", species.c_str(), runNumber, etabin)));
        }
        zMassSpectra[spcType][0][numetabins]->Add((TH1F*)thisFile->Get(Form("z%sMassSpectrum_dataSet%i_data", species.c_str(), runNumber)));
        zMassSpectra_AllSigns[spcType][0][numetabins]->Add((TH1F*)thisFile->Get(Form("z%sMassSpectrum_AllSigns_dataSet%i_data", species.c_str(), runNumber)));
       }

       thisFile->Close();
       delete thisFile;
       break;
      }
     }
     for (TString gammaJetSampleId : gammaJetSampleIds) { // check for gamma jet MC
      if (fname.Contains(gammaJetSampleId)) {
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%s", gammaJetSampleId.Data()));

       for (short etabin = 0; etabin < numetabins; etabin++) {
        const short periodType = (gammaJetSampleId.Contains("pPb") ? 0 : 1);
        const short act_etabin = (gammaJetSampleId.Contains("pPb") ? (numetabins - etabin - 1) : etabin); // period A condition

        g_n[periodType][1][etabin] += (*infoVec)[2+2*numetabins+etabin];
        g_n[2][1][act_etabin] += (*infoVec)[2+2*numetabins+etabin];

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH2F* temp = (TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%s_hist%i_mc_stat", gammaJetSampleId.Data(), etabin));
        gJetHists[periodType][etabin][1][1]->Add(temp);
        gJetHists[periodType][numetabins][1][1]->Add(temp);
        temp = (TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%s_hist%i_mc_stat", gammaJetSampleId.Data(), act_etabin));
        gJetHists[2][etabin][1][1]->Add(temp);
        gJetHists[2][numetabins][1][1]->Add(temp);

        //gJetHists[periodType][etabin][1][1]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%s_hist%i_mc_stat", gammaJetSampleId.Data(), etabin)));
        //gJetHists[2][etabin][1][1]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%s_hist%i_mc_stat", gammaJetSampleId.Data(), act_etabin)));
        gJetHistsSys[periodType][etabin][1][1]->Add((TH2F*)thisFile->Get(Form("gJetPtRatioSys_dataSet%s_hist%i_mc_stat", gammaJetSampleId.Data(), etabin)));
        //gJetHistsSys[2][etabin][1][1]->Add((TH2F*)thisFile->Get(Form("gJetPtRatioSys_dataSet%s_hist%i_mc_stat", gammaJetSampleId.Data(), act_etabin)));
       }

       thisFile->Close();
       delete thisFile;
       break;
      }
     }
     for (TString zeeJetSampleId : zeeJetSampleIds) { // check for Z->ee MC
      if (fname.Contains(zeeJetSampleId)) {
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%s", zeeJetSampleId.Data()));

       for (short etabin = 0; etabin < numetabins; etabin++) {
        const short periodType = (zeeJetSampleId.Contains("pPb") ? 0 : 1);
        const short act_etabin = (zeeJetSampleId.Contains("pPb") ? numetabins - etabin - 1 : etabin); // period A condition

        Zee_n[periodType][1][etabin] += (*infoVec)[2+etabin];
        Zee_n[2][1][act_etabin] += (*infoVec)[2+etabin];

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH2F* temp = (TH2F*)thisFile->Get(Form("zeeJetPtRatio_dataSet%s_hist%i_mc_stat", zeeJetSampleId.Data(), etabin));
        zeeJetHists[periodType][etabin][1][1]->Add(temp);
        zeeJetHists[periodType][numetabins][1][1]->Add(temp);
        temp = (TH2F*)thisFile->Get(Form("zeeJetPtRatio_dataSet%s_hist%i_mc_stat", zeeJetSampleId.Data(), act_etabin));
        zeeJetHists[2][etabin][1][1]->Add(temp);
        zeeJetHists[2][numetabins][1][1]->Add(temp);
       }

       for (short etabin = 0; etabin < numetabins; etabin++) {
        zMassSpectra[1][1][etabin]->Add((TH1F*)thisFile->Get(Form("zeeMassSpectrum_dataSet%s_mc_etabin%i", zeeJetSampleId.Data(), etabin)));
        zMassSpectra_AllSigns[1][1][etabin]->Add((TH1F*)thisFile->Get(Form("zeeMassSpectrum_AllSigns_dataSet%s_mc_etabin%i", zeeJetSampleId.Data(), etabin)));
       }
       zMassSpectra[1][1][numetabins]->Add((TH1F*)thisFile->Get(Form("zeeMassSpectrum_dataSet%s_mc", zeeJetSampleId.Data())));
       zMassSpectra_AllSigns[1][1][numetabins]->Add((TH1F*)thisFile->Get(Form("zeeMassSpectrum_AllSigns_dataSet%s_mc", zeeJetSampleId.Data())));

       thisFile->Close();
       delete thisFile;
       break;
      }
     }
     for (TString zmumuJetSampleId : zmumuJetSampleIds) { // check for Z->mumu MC
      if (fname.Contains(zmumuJetSampleId)) {
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile(rootPath + fname, "READ");
       infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%s", zmumuJetSampleId.Data()));

       for (short etabin = 0; etabin < numetabins; etabin++) {
        const short periodType = (zmumuJetSampleId.Contains("pPb") ? 0 : 1);
        const short act_etabin = (zmumuJetSampleId.Contains("pPb") ? numetabins - etabin - 1 : etabin); // period A condition

        Zmumu_n[periodType][1][etabin] += (*infoVec)[2+numetabins+etabin];
        Zmumu_n[2][1][act_etabin] += (*infoVec)[2+numetabins+etabin];

        // Only add the statistical error plots for MC (don't need to
        // consider systematics)
        TH2F* temp = (TH2F*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%s_hist%i_mc_stat", zmumuJetSampleId.Data(), etabin));
        zmumuJetHists[periodType][etabin][1][1]->Add(temp);
        zmumuJetHists[periodType][numetabins][1][1]->Add(temp);
        temp = (TH2F*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%s_hist%i_mc_stat", zmumuJetSampleId.Data(), act_etabin));
        zmumuJetHists[2][etabin][1][1]->Add(temp);
        zmumuJetHists[2][numetabins][1][1]->Add(temp);
       }

       for (short etabin = 0; etabin < numetabins; etabin++) {
        zMassSpectra[0][1][etabin]->Add((TH1F*)thisFile->Get(Form("zmumuMassSpectrum_dataSet%s_mc_etabin%i", zmumuJetSampleId.Data(), etabin)));
        zMassSpectra_AllSigns[0][1][etabin]->Add((TH1F*)thisFile->Get(Form("zmumuMassSpectrum_AllSigns_dataSet%s_mc_etabin%i", zmumuJetSampleId.Data(), etabin)));
       }
       zMassSpectra[0][1][numetabins]->Add((TH1F*)thisFile->Get(Form("zmumuMassSpectrum_dataSet%s_mc", zmumuJetSampleId.Data())));
       zMassSpectra_AllSigns[0][1][numetabins]->Add((TH1F*)thisFile->Get(Form("zmumuMassSpectrum_AllSigns_dataSet%s_mc", zmumuJetSampleId.Data())));

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


  /**** Calculate total vector boson + jet counts ****/
  for (short periodType = 0; periodType < 3; periodType++) {
   for (short dType = 0; dType < 2; dType++) {
    Zee_n[periodType][dType][numetabins] = 0;
    Zmumu_n[periodType][dType][numetabins] = 0;
    g_n[periodType][dType][numetabins] = 0;
    for (short etabin = 0; etabin < numetabins; etabin++) {
     Zee_n[periodType][dType][numetabins] += Zee_n[periodType][dType][etabin];
     Zmumu_n[periodType][dType][numetabins] += Zmumu_n[periodType][dType][etabin];
     g_n[periodType][dType][numetabins] += g_n[periodType][dType][etabin];
    }
   }
  }

  //cout << Zee_n[0] << " Z(ee)+jet candidate events in data" << endl;
  //cout << Zmumu_n[0] << " Z(mumu)+jet candidate events in data" << endl;
  //cout << g_n[0] << " gamma+jet candidate events in data" << endl;
  //cout << Zee_n[1] << " Z(ee)+jet candidate events in MC" << endl;
  //cout << Zmumu_n[1] << " Z(mumu)+jet candidate events in MC" << endl;
  //cout << g_n[1] << " gamma+jet candidate events in MC" << endl;
 
  TLine* zlines[8] = {};
  TLine* glines[8] = {};
  TLine* xlines[8] = {};
  for (short i = 0; i < 8; i++) {
   const float dz = 0.1;
   const float dg = 0.1;
   const float dx = 0.2;
   zlines[i] = new TLine(pzbins[0], 0.8+dz*i, pzbins[numpzbins], 0.8+dz*i);
   glines[i] = new TLine(pgammabins[0], 0.8+dg*i, pgammabins[numpgammabins], 0.8+dg*i);
   xlines[i] = new TLine(xjrefbins[0], 0.8+dx*i, xjrefbins[numxjrefbins], 0.8+dx*i);
   if (0.8+dz*i == 1) zlines[i]->SetLineStyle(1);
   else zlines[i]->SetLineStyle(3);
   if (0.8+dg*i == 1) glines[i]->SetLineStyle(1);
   else glines[i]->SetLineStyle(3);
   if (0.8+dx*i == 1) xlines[i]->SetLineStyle(1);
   else xlines[i]->SetLineStyle(3);
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
  TH1F* gJetHistDifference[3][numetabins][3];

  TFile* outFile = new TFile(TString(rootPath) + "cc_difference.root", "recreate");

  for (short periodType = 0; periodType < 3; periodType++) {

   string period = "Period A";
   if (periodType == 1) period = "Period B";
   else if (periodType == 2) period = "Period A+B";

   for (short etabin = 0; etabin <= numetabins; etabin++) {

    /**** Plot ZmumuJet info ****/
    topPad->cd();
    topPad->SetLogx();
    vJetHist = GetProfileX("vJetHist", zmumuJetHists[periodType][etabin][0][1], numpzbins, pzbins, false);
    vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
    vJetHist->SetYTitle("#it{p}_{T}^{jet} / #it{p}_{T}^{ref}");
    vJetHist->SetAxisRange(0.65, 1.45, "Y");
    vJetHist->SetMarkerColor(data_color);
    vJetHist->SetLineColor(data_color);
    vJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleOffset(uPadY);

    // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
    vJetHist_lo = GetProfileX("vJetHist_lo", zmumuJetHists[periodType][etabin][0][0], numpzbins, pzbins, false);
    vJetHist_hi = GetProfileX("vJetHist_hi", zmumuJetHists[periodType][etabin][0][2], numpzbins, pzbins, false);
    CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    if (vJetHist_lo) delete vJetHist_lo;
    if (vJetHist_hi) delete vJetHist_hi;
    vJetGraph_sys->SetFillColor(kBlack);
    vJetGraph_sys->SetFillStyle(3001);

    vJetHist_mc = GetProfileX("vJetHist_mc", zmumuJetHists[periodType][etabin][1][1], numpzbins, pzbins, false);
    vJetHist_mc->SetMarkerColor(mc_color);
    vJetHist_mc->SetLineColor(mc_color);

    vJetHist->DrawCopy("E1 X0");
    vJetHist_mc->DrawCopy("SAME E1 X0");
    vJetGraph_sys->Draw("2");

    myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 Data, Cross-Calib Insitu (%i events)", Zmumu_n[periodType][0][etabin]), 1.25, 0.04/uPadY);
    myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("MC with Data Overlay (%i events)", Zmumu_n[periodType][1][etabin]), 1.25, 0.04/uPadY);
    if (etabin < numetabins) {
     if (periodType == 2) myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{Proton} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{Lab} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
    }
    myText(0.155, 0.28, kBlack, "Z (#mu#mu) + Jet", 0.04/uPadY);
    myText(0.155, 0.19, kBlack, period.c_str(), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx();

    vJetHist_rat = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_hist%i", etabin)), zmumuJetHists[periodType][etabin][0][1], zmumuJetHists[periodType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
    vJetHist_rat_lo = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_lo_hist%i", etabin)), zmumuJetHists[periodType][etabin][0][0], zmumuJetHists[periodType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetHist_rat_hi = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_hi_hist%i", etabin)), zmumuJetHists[periodType][etabin][0][2], zmumuJetHists[periodType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    CalcSystematics(vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
    if (vJetHist_rat_lo) delete vJetHist_rat_lo;
    if (vJetHist_rat_hi) delete vJetHist_rat_hi;
    vJetGraph_rat_sys->SetFillColor(kBlack);
    vJetGraph_rat_sys->SetFillStyle(3001);

    vJetHist_rat->SetYTitle("Data / MC");
    //vJetHist_rat->SetAxisRange(0.85, 1.15, "Y");
    vJetHist_rat->SetAxisRange(0.75, 1.35, "Y");
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

    switch (periodType) {
     case 0:
      canvas->SaveAs(Form("%s/PeriodA/%s", plotPath.c_str(), plotName));
      break;
     case 1:
      canvas->SaveAs(Form("%s/PeriodB/%s", plotPath.c_str(), plotName));
      break;
     case 2:
      canvas->SaveAs(Form("%s/PeriodAB/%s", plotPath.c_str(), plotName));
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
    vJetHist = GetProfileX("vJetHist", zeeJetHists[periodType][etabin][0][1], numpzbins, pzbins, false);
    vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
    vJetHist->SetYTitle("#it{p}_{T}^{jet} / #it{p}_{T}^{ref}");
    vJetHist->SetAxisRange(0.65, 1.45, "Y");
    vJetHist->SetMarkerColor(data_color);
    vJetHist->SetLineColor(data_color);
    vJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleOffset(uPadY);

    // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
    vJetHist_lo = GetProfileX("vJetHist_lo", zeeJetHists[periodType][etabin][0][0], numpzbins, pzbins, false);
    vJetHist_hi = GetProfileX("vJetHist_hi", zeeJetHists[periodType][etabin][0][2], numpzbins, pzbins, false);
    CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    if (vJetHist_lo) delete vJetHist_lo;
    if (vJetHist_hi) delete vJetHist_hi;
    vJetGraph_sys->SetFillColor(kBlack);
    vJetGraph_sys->SetFillStyle(3001);

    vJetHist_mc = GetProfileX("vJetHist_mc", zeeJetHists[periodType][etabin][1][1], numpzbins, pzbins, false);
    vJetHist_mc->SetMarkerColor(mc_color);
    vJetHist_mc->SetLineColor(mc_color);

    vJetHist->DrawCopy("E1 X0");
    vJetHist_mc->DrawCopy("SAME E1 X0");
    vJetGraph_sys->Draw("2");

    myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 Data, Cross-Calib Insitu (%i events)", Zee_n[periodType][0][etabin]), 1.25, 0.04/uPadY);
    myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("MC Signal Only (%i events)", Zee_n[periodType][1][etabin]), 1.25, 0.04/uPadY);
    if (etabin < numetabins) {
     if (periodType == 2) myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{Proton} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{Lab} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
    }
    myText(0.155, 0.28, kBlack, "Z (ee) + Jet", 0.04/uPadY);
    myText(0.155, 0.19, kBlack, period.c_str(), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx();

    vJetHist_rat = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_hist%i", etabin)), zeeJetHists[periodType][etabin][0][1], zeeJetHists[periodType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
    vJetHist_rat_lo = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_lo_hist%i", etabin)), zeeJetHists[periodType][etabin][0][0], zeeJetHists[periodType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetHist_rat_hi = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_hi_hist%i", etabin)), zeeJetHists[periodType][etabin][0][2], zeeJetHists[periodType][etabin][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    CalcSystematics(vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
    if (vJetHist_rat_lo) delete vJetHist_rat_lo;
    if (vJetHist_rat_hi) delete vJetHist_rat_hi;
    vJetGraph_rat_sys->SetFillColor(kBlack);
    vJetGraph_rat_sys->SetFillStyle(3001);

    vJetHist_rat->SetYTitle("Data / MC");
    //vJetHist_rat->SetAxisRange(0.85, 1.15, "Y");
    vJetHist_rat->SetAxisRange(0.75, 1.35, "Y");
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
    switch (periodType) {
     case 0:
      canvas->SaveAs(Form("%s/PeriodA/%s", plotPath.c_str(), plotName));
      break;
     case 1:
      canvas->SaveAs(Form("%s/PeriodB/%s", plotPath.c_str(), plotName));
      break;
     case 2:
      canvas->SaveAs(Form("%s/PeriodAB/%s", plotPath.c_str(), plotName));
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
    vJetHist = GetProfileX("vJetHist", gJetHists[periodType][etabin][0][1], numpgammabins, pgammabins, true);
    vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
    vJetHist->SetYTitle("#it{p}_{T}^{jet} / #it{p}_{T}^{ref}");
    vJetHist->SetAxisRange(0.65, 1.45, "Y");
    vJetHist->SetMarkerColor(data_color);
    vJetHist->SetLineColor(data_color);
    vJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
    vJetHist->GetYaxis()->SetTitleOffset(uPadY);

    // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
    vJetHist_lo = GetProfileX("vJetHist_lo", gJetHists[periodType][etabin][0][0], numpgammabins, pgammabins, true);
    vJetHist_hi = GetProfileX("vJetHist_hi", gJetHists[periodType][etabin][0][2], numpgammabins, pgammabins, true);
    CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    vJetGraph_sys->SetFillColor(kBlack);
    vJetGraph_sys->SetFillStyle(3001);

    vJetHist_mc = GetProfileX("vJetHist_mc", gJetHists[periodType][etabin][1][1], numpgammabins, pgammabins, true);
    vJetHist_mc->SetMarkerColor(mc_color);
    vJetHist_mc->SetLineColor(mc_color);

    for (short errType = 0; errType < 3; errType++) {
     const string error = (errType == 0 ? "syslo" : (errType == 1 ? "stat" : "syshi"));
     string periodStr = "periodA";
     if (periodType == 1) periodStr = "periodB";
     else if (periodType == 2) periodStr = "periodAB";
     gJetHistDifference[periodType][etabin][errType] = new TH1F(Form("gJetPtRatio_diff%i_%s_%s", etabin, error.c_str(), periodStr.c_str()), ";#it{p}_{T}^{ref} #left[GeV#right]", numpgammabins, pgammabins);
     for (short pgammabin = 1; pgammabin <= numpgammabins; pgammabin++) {
      double dataVal, dataErr;
      switch (errType) {
       case 0:
        dataVal = vJetHist_lo->GetBinContent(pgammabin);
        dataErr = vJetHist_lo->GetBinError(pgammabin);
        break;
       case 2:
        dataVal = vJetHist_hi->GetBinContent(pgammabin);
        dataErr = vJetHist_hi->GetBinError(pgammabin);
        break;
       default:
        dataVal = vJetHist->GetBinContent(pgammabin);
        dataErr = vJetHist->GetBinError(pgammabin);
      } 
      gJetHistDifference[periodType][etabin][errType]->SetBinContent(pgammabin, dataVal - vJetHist_mc->GetBinContent(pgammabin));
      gJetHistDifference[periodType][etabin][errType]->SetBinError(pgammabin, TMath::Sqrt(TMath::Power(dataErr,2) + TMath::Power(vJetHist_mc->GetBinError(pgammabin),2)));
     }
     gJetHistDifference[periodType][etabin][errType]->Write();
    }
    if (vJetHist_lo) delete vJetHist_lo;
    if (vJetHist_hi) delete vJetHist_hi;

    vJetHist->DrawCopy("E1 X0");
    vJetHist_mc->DrawCopy("SAME E1 X0");
    vJetGraph_sys->Draw("2");
    myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 Data, Cross-Calib Insitu (%i events)", g_n[periodType][0][etabin]), 1.25, 0.04/uPadY);
    myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("MC %s (%i events)", (runValidation ? "Signal Only":"with Data Overlay"), g_n[periodType][1][etabin]), 1.25, 0.04/uPadY);
    if (etabin < numetabins) {
     if (periodType == 2) myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{Proton} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{Lab} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
    }
    myText(0.155, 0.28, kBlack, "#gamma + Jet", 0.04/uPadY);
    myText(0.155, 0.19, kBlack, period.c_str(), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx();
    vJetHist_rat = GetDataOverMC(TString(Form("gJetPtDataMCRatio_hist%i", etabin)), gJetHists[periodType][etabin][0][1], gJetHists[periodType][etabin][1][1], numpgammabins, pgammabins, numxjrefbins, xjrefbins, true);
    vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
    vJetHist_rat_lo = GetDataOverMC(TString(Form("gJetPtDataMCRatio_lo_hist%i", etabin)), gJetHists[periodType][etabin][0][0], gJetHists[periodType][etabin][1][1], numpgammabins, pgammabins, numxjrefbins, xjrefbins, true);
    vJetHist_rat_hi = GetDataOverMC(TString(Form("gJetPtDataMCRatio_hi_hist%i", etabin)), gJetHists[periodType][etabin][0][2], gJetHists[periodType][etabin][1][1], numpgammabins, pgammabins, numxjrefbins, xjrefbins, true);
    CalcSystematics(vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
    if (vJetHist_rat_lo) delete vJetHist_rat_lo;
    if (vJetHist_rat_hi) delete vJetHist_rat_hi;
    vJetGraph_rat_sys->SetFillColor(kBlack);
    vJetGraph_rat_sys->SetFillStyle(3001);

    vJetHist_rat->SetYTitle("Data / MC");
    vJetHist_rat->SetAxisRange(0.85, 1.15, "Y");
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
    for (TLine* line : glines) line->Draw("SAME");
    if (etabin < numetabins) plotName = Form("gamma_jet%i.pdf", etabin);
    else plotName = Form("gamma_jet_combined.pdf");
    switch (periodType) {
     case 0:
      canvas->SaveAs(Form("%s/PeriodA/%s", plotPath.c_str(), plotName));
      break;
     case 1:
      canvas->SaveAs(Form("%s/PeriodB/%s", plotPath.c_str(), plotName));
      break;
     case 2:
      canvas->SaveAs(Form("%s/PeriodAB/%s", plotPath.c_str(), plotName));
      break;
    }
    if (vJetHist) delete vJetHist;
    if (vJetHist_mc) delete vJetHist_mc;
    if (vJetGraph_sys) delete vJetGraph_sys;
    if (vJetHist_rat) delete vJetHist_rat;
    if (vJetGraph_rat_sys) delete vJetGraph_rat_sys;
    for (short errType = 0; errType < 3; errType++)
     if (gJetHistDifference[periodType][etabin][errType])
      delete gJetHistDifference[periodType][etabin][errType];

    if (etabin == numetabins) continue;


    /**** Plots xjref distributions, binned by ptref ****/
    for (int pgammabin = 0; pgammabin < numpgammabins; pgammabin++) {
     const double pref_lo = pgammabins[pgammabin];
     const double pref_hi =  pgammabins[pgammabin+1];
     topPad->cd();
     topPad->SetLogx(0);
     vJetHist = gJetHists[periodType][etabin][0][1]->ProjectionY("vJetProjection", pgammabin, pgammabin);
     const float counts_data = vJetHist->Integral();
     const float total_data = gJetHists[periodType][etabin][0][1]->Integral();
     vJetHist->Rebin(rebinFactor);
     vJetHist->Scale(1./vJetHist->Integral());
//     vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
     vJetHist->SetXTitle("#it{p}_{T}^{jet} / #it{p}_{T}^{ref}");
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
     vJetHist_lo = gJetHists[periodType][etabin][0][0]->ProjectionY("vJetProjection_lo", pgammabin, pgammabin);
     vJetHist_lo->Rebin(rebinFactor);
     vJetHist_lo->Scale(1./vJetHist_lo->Integral()); 
     //vJetHist_lo->Scale(1./counts_data); 
     vJetHist_hi = gJetHists[periodType][etabin][0][2]->ProjectionY("vJetProjection_hi", pgammabin, pgammabin);
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

     vJetHist_mc = gJetHists[periodType][etabin][1][1]->ProjectionY("vJetProjection_mc", pgammabin, pgammabin);
     const float counts_mc = vJetHist_mc->Integral();
     const float total_mc = gJetHists[periodType][etabin][1][1]->Integral();
     vJetHist_mc->Rebin(rebinFactor);
     vJetHist_mc->Scale(1./vJetHist_mc->Integral()); 
     vJetHist_mc->SetMarkerColor(mc_color);
     vJetHist_mc->SetLineColor(mc_color);

     vJetHist->DrawCopy("E1 X0");
     vJetHist_mc->DrawCopy("SAME E1 X0");
//     vJetGraph_sys->Draw("2");

     float n = 100. * counts_data / total_data;
     myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 Data, Cross-Calib Insitu (%.1f%% of events)", n), 1.25, 0.04/uPadY);

     n = 100. * counts_mc / total_mc;
     myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("MC %s (%.1f%% of events)", (runValidation ? "Signal Only":"with Data Overlay"), n), 1.25, 0.04/uPadY);

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

     myText(0.655, 0.43, kBlack, "#gamma + Jet", 0.04/uPadY);
     myText(0.655, 0.34, kBlack, Form("%g < #it{p}_{T}^{ref} < %g", pref_lo, pref_hi), 0.04/uPadY);
     myText(0.655, 0.25, kBlack, period.c_str(), 0.04/uPadY);
     if (periodType == 2) myText(0.655, 0.16,kBlack, Form("%g < #eta_{jet}^{Proton} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.655, 0.16,kBlack, Form("%g < #eta_{jet}^{Lab} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);

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
     plotName = Form("pref_slices/gamma_jet%i_pbin%i.pdf", etabin, pgammabin);
     switch (periodType) {
      case 0:
       canvas->SaveAs(Form("%s/PeriodA/%s", plotPath.c_str(), plotName));
       break;
      case 1:
       canvas->SaveAs(Form("%s/PeriodB/%s", plotPath.c_str(), plotName));
       break;
      case 2:
       canvas->SaveAs(Form("%s/PeriodAB/%s", plotPath.c_str(), plotName));
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
  for (short periodType = 0; periodType < 2; periodType++) {
   const string period = (periodType == 0 ? "Period A":"Period B");

   for (short etabin = 0; etabin < numetabins; etabin++) {
    topPad->cd();
    TH2F* thisHist = gJetHistsSys[periodType][etabin][0][1];
    TH1F* rmsHist = new TH1F(Form("rms_etabin%i_%s", etabin, period.c_str()), "", numpzbins, pzbins);
    for (short pzbin = 0; pzbin < numpzbins; pzbin++) {
     //const float pt = 0.5*(pzbins[pzbin+1]+pzbins[pzbin]);
     float rms = 0;
     float sumWeights = 0;
     for (short sigbin = 0; sigbin < numSigmaBins; sigbin++) {
      const float sig = thisHist->GetYaxis()->GetBinCenter(sigbin+1);
      //const float sig = -maxSigma + 2*maxSigma*(sigbin+0.5) / numSigmaBins;
      const float weight = thisHist->GetBinContent(pzbin+1, sigbin+1);
      rms += pow(sig, 2) * weight;
      sumWeights += weight;
     }
     rms = sqrt(rms) / sqrt(sumWeights);
     rmsHist->SetBinContent(pzbin+1, rms);
    }
    topPad->SetLogz();
    thisHist->Draw("col");
    thisHist->GetXaxis()->SetLabelSize(0.04/uPadY);
    thisHist->GetYaxis()->SetLabelSize(0.04/uPadY);
    thisHist->GetYaxis()->SetTitleSize(0.04/uPadY);
    thisHist->GetYaxis()->SetTitleOffset(1.1*uPadY);
    myText(0.72, 0.89, kBlack, period.c_str(), 0.04/uPadY);
    myText(0.72, 0.8,kBlack, Form("%g < #eta_{jet}^{Lab} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);

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
    canvas->SaveAs(Form("%s/Period%s/jetSystematics_etabin%i.pdf", plotPath.c_str(), (periodType==0 ? "A":"B"), etabin));
    if (rmsHist) delete rmsHist;
   }
  }
  outFile->Write();
  if (outFile) delete outFile;


  /**** delete objects ****/
  TLine* lines[8] = {};
  for (short i = 0; i < 8; i++) {
   //if (glines[i]) delete glines[i];
   //if (zlines[i]) delete zlines[i];
   //if (xlines[i]) delete xlines[i];
   lines[i] = new TLine(60, 0.6+0.2*i, 110, 0.6+0.2*i);
   if (0.6+0.2*i == 1) lines[i]->SetLineStyle(1);
   else lines[i]->SetLineStyle(3);
  }


  /**** Plot mumu mass spectra ****/
  for (int etabin = 0; etabin <= numetabins; etabin++) {
   for (int species = 0; species < 2; species++) {
    topPad->cd();
    topPad->SetLogx();
    double mean[2] = {};
    double mean_err[2] = {};
    double sigma[2] = {};
    double sigma_err[2] = {};
    TF1* fits[2];
    //RooZfit* fits[2];
    for (short dType = 0; dType < 2; dType++) {
     TH1F* thisHist = zMassSpectra[species][dType][etabin];
     Color_t color = (dType==0 ? data_color : mc_color);
     thisHist->GetXaxis()->SetTitle("#font[12]{ll} Invariant Mass #left[GeV#right]");
     thisHist->GetYaxis()->SetTitle("Normalized Counts / 1 GeV");
     thisHist->GetYaxis()->SetTitleSize(0.04/uPadY);
     thisHist->GetYaxis()->SetTitleOffset(1.1*uPadY);
     thisHist->GetYaxis()->SetLabelSize(0.04/uPadY);
     thisHist->SetLineColor(color);
     thisHist->SetMarkerColor(color);
     thisHist->GetXaxis()->SetNdivisions(50802, false);
     thisHist->Scale(1.0, "width");

     //TF1* gausFit = new TF1(Form("fit_guess_species%i_dType%i", species, dType), "gaus(0)", Z_mass - 10, Z_mass + 10);
     //thisHist->Fit(gausFit, "R", "L");
     //double m = gausFit->GetParameter(1);
     //double s = gausFit->GetParameter(2);
     //const double scale = 1.0 / thisHist->Integral(thisHist->FindBin(m-Z_mass_fitNsigma*s), thisHist->FindBin(m+Z_mass_fitNsigma*s));
     //thisHist->Scale(scale);
     
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
     TH1F* thisHist = zMassSpectra[species][dType][etabin];
     if (dType == 0) thisHist->Draw("hist");
     else thisHist->Draw("hist same");
     fits[dType]->Draw("same");
     if (etabin < numetabins)
      myText(0.175, 0.15,kBlack, Form("%g < #eta_{jet} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
    }
    if (species == 0) {
     myText(0.175, 0.88, kBlack, "Z (#mu#mu) + Jet", 0.04/uPadY);
     myMarkerText(0.175, 0.80, data_color, kFullCircle, Form("2016 Data (%i events)", Zmumu_n[2][0][etabin]), 1.25, 0.04/uPadY);
     myMarkerText(0.175, 0.55, mc_color, kFullCircle, Form("MC with Data Overlay (%i events)", Zmumu_n[2][1][etabin]), 1.25, 0.04/uPadY);
    }
    else if (species == 1) {
     myText(0.175, 0.88, kBlack, "Z (ee) + Jet", 0.04/uPadY);
     myMarkerText(0.175, 0.80, data_color, kFullCircle, Form("2016 Data (%i events)", Zee_n[2][0][etabin]), 1.25, 0.04/uPadY);
     myMarkerText(0.175, 0.55, mc_color, kFullCircle, Form("MC Signal Only (%i events)", Zee_n[2][1][etabin]), 1.25, 0.04/uPadY);
    }
    myText(0.175, 0.72, kBlack, Form("m_{Z}^{data} = %.2f #pm %.2f GeV", mean[0], mean_err[0]), 0.04/uPadY);
    myText(0.175, 0.64, kBlack, Form("#sigma_{Z}^{data} = %.2f #pm %.2f GeV", sigma[0], sigma_err[0]), 0.04/uPadY);
    myText(0.175, 0.47, kBlack, Form("m_{Z}^{mc} = %.2f #pm %.2f GeV", mean[1], mean_err[1]), 0.04/uPadY);
    myText(0.175, 0.39, kBlack, Form("#sigma_{Z}^{mc} = %.2f #pm %.2f GeV", sigma[1], sigma_err[1]), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx();
    TH1F* thisHist = (TH1F*)zMassSpectra[species][0][etabin]->Clone(Form("invMass_species%i_clone", species));
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
    thisHist->SetAxisRange(0.6, 1.6, "Y");
    thisHist->GetYaxis()->ChangeLabel(-1, -1, -1, -1, -1, -1, " ");
    thisHist->GetXaxis()->SetTickLength(0.08);
    thisHist->GetXaxis()->SetNdivisions(50802, false);
    thisHist->GetYaxis()->SetNdivisions(405, false);
    thisHist->Draw("p");
    for (TLine* line : lines) line->Draw("same");

    if (etabin != numetabins) {
     if (species == 0) canvas->SaveAs(Form("%s/zmumu_mass_comparison_etabin%i.pdf", plotPath.c_str(), etabin));
     else if (species == 1) canvas->SaveAs(Form("%s/zee_mass_comparison_etabin%i.pdf", plotPath.c_str(), etabin));
    }
    else {
     if (species == 0) canvas->SaveAs(Form("%s/zmumu_mass_comparison.pdf", plotPath.c_str()));
     else if (species == 1) canvas->SaveAs(Form("%s/zee_mass_comparison.pdf", plotPath.c_str()));
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
  for (short periodType = 0; periodType < 3; periodType++) {
   const string period = (periodType == 0 ? "Period A": (periodType == 1 ? "Period B":"Period A+B"));

   for (short etabin = 0; etabin < numetabins; etabin++) {
    TH2F* dataHist = gJetHists[periodType][etabin][0][1];
    TH2F* mcHist = gJetHists[periodType][etabin][1][1];
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
    myText(0.1, 0.15, kBlack, Form("2016 Data, Cross-Calib Insitu (%i events)", g_n[periodType][0][etabin]), 0.02/rPadX);
    myText(0.6, 0.85,kBlack, Form("%g < #eta_{jet} < %g", etabins[etabin], etabins[etabin+1]), 0.02/rPadX);
    myText(0.6, 0.8,kBlack, period.c_str(), 0.02/rPadX);

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
    myText(0.2, 0.15, kBlack, Form("MC %s (%i events)", (runValidation ? "Signal Only":"with Data Overlay"), g_n[periodType][1][etabin]), 0.02/lPadX);

    const char* plotName = Form("gamma_jet%i_th2.pdf", etabin);
    switch (periodType) {
     case 0:
      th2canvas->SaveAs(Form("%s/PeriodA/%s", plotPath.c_str(), plotName));
      break;
     case 1:
      th2canvas->SaveAs(Form("%s/PeriodB/%s", plotPath.c_str(), plotName));
      break;
     case 2:
      th2canvas->SaveAs(Form("%s/PeriodAB/%s", plotPath.c_str(), plotName));
      break;
    }
   }
  }

  return;
}
