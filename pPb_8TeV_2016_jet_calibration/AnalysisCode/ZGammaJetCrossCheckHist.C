#include "../Params.C"
#include "../../Initialization.C"

TH1D* GetDataOverMC(const TString name, TH2F* data, TH2F* mc, const int numxbins, const double* xbins, const int numybins, const double* ybins, const int ylo, const int yhi) {
  TH1D* dataOverMC = new TH1D(name, "", numxbins, xbins);
  for (int xbin = 1; xbin <= numxbins; xbin++) {
   TH1D* temp = data->ProjectionY(name + TString(Form("data_xbin%i", xbin)), xbin, xbin);
   // if we only want a subrange in xjref (y-axis), set other bins & errors to 0
   if (yhi != -1) {
    for (int ybin = 1; ybin <= numybins; ybin++) {
     if (ylo <= ybin && ybin <= yhi) continue;
     temp->SetBinContent(ybin, 0);
     temp->SetBinError(ybin, 0);
    }
   }
   const double dataAvg = temp->GetMean();
   const double dataErr = temp->GetMeanError();
   delete temp;
   temp = mc->ProjectionY(name + TString(Form("mc_xbin%i", xbin)), xbin, xbin);
   if (yhi != -1) {
    for (int ybin = 1; ybin <= numybins; ybin++) {
     if (ylo <= ybin && ybin <= yhi) continue;
     temp->SetBinContent(ybin, 0);
     temp->SetBinError(ybin, 0);
    }
   }
   const double mcAvg = temp->GetMean();
   const double mcErr = temp->GetMeanError();
   delete temp;
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

  TH2F* zeeJetHists[3][numetabins][2][3];
  //TH2F* zeeJetHistsSys[3][numetabins][2][3];
  TH2F* zmumuJetHists[3][numetabins][2][3];
  //TH2F* zmumuJetHistsSys[3][numetabins][2][3];
  TH2F* gJetHists[3][numetabins][2][3];
  TH2F* gJetHistsSys[3][numetabins][2][3];
  TH1F* zMassSpectra[2][2][numetabins+1];
  TH1F* zMassSpectra_AllSigns[2][2][numetabins+1];

  for (short etabin = 0; etabin < numetabins; etabin++) {

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

      zeeJetHists[periodType][etabin][dType][errType] = new TH2F(Form("zeeJetPtRatio_hist%i_%s_%s_%s", etabin, dataType.c_str(), error.c_str(), period.c_str()), ";#it{p}_{T}^{ref} #left[GeV#right];#it{p}_{T}^{j} / #it{p}_{T}^{ref}", numpbins, pbins, numxjrefbins, xjrefbins);
      zeeJetHists[periodType][etabin][dType][errType]->Sumw2();
      zmumuJetHists[periodType][etabin][dType][errType] = new TH2F(Form("zmumuJetPtRatio_hist%i_%s_%s_%s", etabin, dataType.c_str(), error.c_str(), period.c_str()), ";#it{p}_{T}^{ref} #left[GeV#right];#it{p}_{T}^{j} / #it{p}_{T}^{ref}", numpbins, pbins, numxjrefbins, xjrefbins);
      zmumuJetHists[periodType][etabin][dType][errType]->Sumw2();
      gJetHists[periodType][etabin][dType][errType] = new TH2F(Form("gJetPtRatio_hist%i_%s_%s_%s", etabin, dataType.c_str(), error.c_str(), period.c_str()), ";#it{p}_{T}^{ref} #left[GeV#right];#it{p}_{T}^{j} / #it{p}_{T}^{ref}", numprefbins, prefbins, numxjrefbins, xjrefbins);
      gJetHists[periodType][etabin][dType][errType]->Sumw2();

      gJetHistsSys[periodType][etabin][dType][errType] = new TH2F(Form("gJetPtRatioSys_hist%i_%s_%s_%s", etabin, dataType.c_str(), error.c_str(), period.c_str()), ";#it{p}_{T}^{jet} #left[GeV#right];#sigma(#it{p}_{T}^{ref})/#it{p}_{T}^{jet}", numpbins, pbins, numSigmaBins, -maxSigma, maxSigma);
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
     if (etabin == 0) {
      zMassSpectra[spcType][dType][numetabins] = new TH1F(Form("z%sMassSpectrum_%s", species.c_str(), dataType.c_str()), "", 50, 60, 110);
      zMassSpectra[spcType][dType][numetabins]->Sumw2();
      zMassSpectra_AllSigns[spcType][dType][numetabins] = new TH1F(Form("z%sMassSpectrum_AllSigns_%s", species.c_str(), dataType.c_str()), "", 50, 60, 110);
      zMassSpectra_AllSigns[spcType][dType][numetabins]->Sumw2();
     }
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

         zeeJetHists[periodType][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("zeeJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, etabin, error.c_str())));
         zeeJetHists[2][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("zeeJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, act_etabin, error.c_str())));
         zmumuJetHists[periodType][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, etabin, error.c_str())));
         zmumuJetHists[2][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, act_etabin, error.c_str())));
         gJetHists[periodType][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, etabin, error.c_str())));
         gJetHists[2][etabin][0][errType]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%i_hist%i_data_%s", runNumber, act_etabin, error.c_str())));
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
        gJetHists[periodType][etabin][1][1]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%s_hist%i_mc_stat", gammaJetSampleId.Data(), etabin)));
        gJetHists[2][etabin][1][1]->Add((TH2F*)thisFile->Get(Form("gJetPtRatio_dataSet%s_hist%i_mc_stat", gammaJetSampleId.Data(), act_etabin)));
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
        zeeJetHists[periodType][etabin][1][1]->Add((TH2F*)thisFile->Get(Form("zeeJetPtRatio_dataSet%s_hist%i_mc_stat", zeeJetSampleId.Data(), etabin)));
        zeeJetHists[2][etabin][1][1]->Add((TH2F*)thisFile->Get(Form("zeeJetPtRatio_dataSet%s_hist%i_mc_stat", zeeJetSampleId.Data(), act_etabin)));
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
        zmumuJetHists[periodType][etabin][1][1]->Add((TH2F*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%s_hist%i_mc_stat", zmumuJetSampleId.Data(), etabin)));
        zmumuJetHists[2][etabin][1][1]->Add((TH2F*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%s_hist%i_mc_stat", zmumuJetSampleId.Data(), act_etabin)));
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
  for (short i = 0; i < 8; i++) {
   zlines[i] = new TLine(pbins[0], 0.8+0.2*i, pbins[numpbins], 0.8+0.2*i);
   glines[i] = new TLine(prefbins[0], 0.8+0.2*i, prefbins[numprefbins], 0.8+0.2*i);
   if (0.8+0.2*i == 1) {
    zlines[i]->SetLineStyle(1);
    glines[i]->SetLineStyle(1);
   }
   else {
    zlines[i]->SetLineStyle(3);
    glines[i]->SetLineStyle(3);
   }
  }

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

  TProfile *zJetHist, *zJetHist_mc, *zJetHist_lo, *zJetHist_hi, *gJetHist, *gJetHist_mc, *gJetHist_lo, *gJetHist_hi;
  TH1D *zJetHist_rat, *zJetHist_rat_lo, *zJetHist_rat_hi, *gJetHist_rat, *gJetHist_rat_lo, *gJetHist_rat_hi;
  TGraphAsymmErrors *zJetGraph_sys, *zJetGraph_rat_sys, *gJetGraph_sys, *gJetGraph_rat_sys;
  TH1F* gJetHistDifference[3][numetabins][3];
  TFile* outFile = new TFile(TString(rootPath) + "cc_difference.root", "recreate");
  const short xjrefmax = (runallxjrefs ? numxjrefbins/xjrefRebinFactor : 0);
  for (short periodType = 0; periodType < 3; periodType++) {
   string period = "Period A";
   if (periodType == 1) period = "Period B";
   else if (periodType == 2) period = "Period A+B";
   for (short etabin = 0; etabin < numetabins; etabin++) {

    /**** Plot ZmumuJet info ****/

    for (int xjrefbin = 0; xjrefbin <= xjrefmax; xjrefbin++) {
     const int xjrefbin_lo = (xjrefbin==0 ? 1 : (xjrefbin-1)*xjrefRebinFactor+1);
     const int xjrefbin_hi = (xjrefbin==0 ? -1 : xjrefbin*xjrefRebinFactor);
     const double xjref_lo = (xjrefbin==0 ? 0.65 : xjrefbins[xjrefbin_lo - 1]);
     const double xjref_hi = (xjrefbin==0 ? 1.45 : xjrefbins[xjrefbin_hi]);
     topPad->cd();
     gPad->SetLogx();
     zJetHist = zmumuJetHists[periodType][etabin][0][1]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     zJetGraph_sys = new TGraphAsymmErrors((TH1F*)zJetHist); // for plotting systematics
     //TProfile* zJetHist_sys = (TProfile*)zJetHist->Clone(Form("%s_systematics", zJetHist->GetName()));
     zJetHist->SetYTitle("#it{p}_{T}^{jet} / #it{p}_{T}^{ref}");
     zJetHist->SetAxisRange(xjref_lo, xjref_hi, "Y");
     zJetHist->SetMarkerColor(data_color);
     zJetHist->SetLineColor(data_color);
     zJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
     zJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
     zJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
     zJetHist->GetYaxis()->SetTitleOffset(uPadY);

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     zJetHist_lo = zmumuJetHists[periodType][etabin][0][0]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     zJetHist_hi = zmumuJetHists[periodType][etabin][0][2]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     for (short pbin = 1; pbin <= numpbins; pbin++) {
      zJetGraph_sys->SetPoint(pbin-1, 0.5*(pbins[pbin-1]+pbins[pbin]), zJetHist->GetBinContent(pbin));
      zJetGraph_sys->SetPointEYhigh(pbin-1, TMath::Abs(zJetHist_hi->GetBinContent(pbin) - zJetHist->GetBinContent(pbin))); // set high systematics
      zJetGraph_sys->SetPointEYlow(pbin-1, TMath::Abs(zJetHist->GetBinContent(pbin) - zJetHist_lo->GetBinContent(pbin))); // set low systematics
     }
     if (zJetHist_lo) delete zJetHist_lo;
     if (zJetHist_hi) delete zJetHist_hi;
     zJetGraph_sys->SetFillColor(kBlack);
     zJetGraph_sys->SetFillStyle(3001);

     zJetHist_mc = zmumuJetHists[periodType][etabin][1][1]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     zJetHist_mc->SetMarkerColor(mc_color);
     zJetHist_mc->SetLineColor(mc_color);

     zJetHist->DrawCopy("E1 X0");
     zJetHist_mc->DrawCopy("SAME E1 X0");
     zJetGraph_sys->Draw("2");
     if (xjrefbin == 0) {
      myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 Data, Cross-Calib Insitu (%i events)", Zmumu_n[periodType][0][etabin]), 1.25, 0.04/uPadY);
      myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("MC with Data Overlay (%i events)", Zmumu_n[periodType][1][etabin]), 1.25, 0.04/uPadY);
     }
     else {
      const float counts_data = zmumuJetHists[periodType][etabin][0][1]->Integral(1, numpbins, xjrefbin_lo, xjrefbin_hi);
      const float total_data = zmumuJetHists[periodType][etabin][0][1]->Integral();
      float n = 100. * counts_data / total_data; 
      myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 Data, Cross-Calib Insitu (%.1f%% of events)", n), 1.25, 0.04/uPadY);

      const float counts_mc = zmumuJetHists[periodType][etabin][1][1]->Integral(1, numpbins, xjrefbin_lo, xjrefbin_hi);
      const float total_mc = zmumuJetHists[periodType][etabin][1][1]->Integral();
      n = 100. * counts_mc / total_mc;
      myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("MC with Data Overlay (%.1f%% events)", n), 1.25, 0.04/uPadY);
     }
     if (periodType == 0) myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{A} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{B} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     if (xjrefbin != 0) {
      myText(0.155, 0.37, kBlack, "Z (#mu#mu) + Jet", 0.04/uPadY);
      myText(0.155, 0.28, kBlack, Form("%g < #it{x}_{J}^{ref} < %g", xjrefbins[xjrefbin_lo-1], xjrefbins[xjrefbin_hi]), 0.04/uPadY);
     }
     else myText(0.155, 0.28, kBlack, "Z (#mu#mu) + Jet", 0.04/uPadY);
     myText(0.155, 0.19, kBlack, period.c_str(), 0.04/uPadY);

     bottomPad->cd();
     gPad->SetLogx();

     zJetHist_rat = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_hist%i", etabin)), zmumuJetHists[periodType][etabin][0][1], zmumuJetHists[periodType][etabin][1][1], numpbins, pbins, numxjrefbins, xjrefbins, xjrefbin_lo, xjrefbin_hi);
     zJetGraph_rat_sys = new TGraphAsymmErrors(zJetHist_rat);
     zJetHist_rat_lo = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_lo_hist%i", etabin)), zmumuJetHists[periodType][etabin][0][0], zmumuJetHists[periodType][etabin][1][1], numpbins, pbins, numxjrefbins, xjrefbins, xjrefbin_lo, xjrefbin_hi);
     zJetHist_rat_hi = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_hi_hist%i", etabin)), zmumuJetHists[periodType][etabin][0][2], zmumuJetHists[periodType][etabin][1][1], numpbins, pbins, numxjrefbins, xjrefbins, xjrefbin_lo, xjrefbin_hi);
     for (short pbin = 1; pbin <= numpbins; pbin++) {
      zJetGraph_rat_sys->SetPoint(pbin-1, 0.5*(pbins[pbin-1]+pbins[pbin]), zJetHist_rat->GetBinContent(pbin));
      zJetGraph_rat_sys->SetPointEYlow(pbin-1, TMath::Abs(zJetHist_rat->GetBinContent(pbin) - zJetHist_rat_lo->GetBinContent(pbin)));
      zJetGraph_rat_sys->SetPointEYhigh(pbin-1, TMath::Abs(zJetHist_rat_hi->GetBinContent(pbin) - zJetHist_rat->GetBinContent(pbin)));
     }
     if (zJetHist_rat_lo) delete zJetHist_rat_lo;
     if (zJetHist_rat_hi) delete zJetHist_rat_hi;
     zJetGraph_rat_sys->SetFillColor(kBlack);
     zJetGraph_rat_sys->SetFillStyle(3001);

     zJetHist_rat->SetYTitle("Data / MC");
     zJetHist_rat->SetAxisRange(0.65, 1.45, "Y");
     zJetHist_rat->GetYaxis()->SetNdivisions(405);
     zJetHist_rat->GetXaxis()->SetTitleSize(0.04/dPadY);
     zJetHist_rat->GetYaxis()->SetTitleSize(0.04/dPadY);
     zJetHist_rat->GetXaxis()->SetTitleOffset(1);
     zJetHist_rat->GetYaxis()->SetTitleOffset(dPadY);
     zJetHist_rat->GetXaxis()->SetLabelSize(0.04/dPadY);
     zJetHist_rat->GetYaxis()->SetLabelSize(0.04/dPadY);
     zJetHist_rat->GetXaxis()->SetTickLength(0.08);

     zJetHist_rat->Draw("E1 X0");
     zJetGraph_rat_sys->Draw("2");
     for (TLine* line : zlines) line->Draw("SAME");
     char* plotName;
     if (xjrefbin == 0) plotName = Form("z_mumu_jet%i.pdf", etabin);
     else plotName = Form("xjref/z_mumu_jet%i_xjrefbin%i.pdf", etabin, xjrefbin);

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
     if (zJetHist_rat) delete zJetHist_rat;
     if (zJetGraph_rat_sys) delete zJetGraph_rat_sys;
    }


    /**** Plots ZeeJet info ****/
    for (int xjrefbin = 0; xjrefbin <= xjrefmax; xjrefbin++) {
     const int xjrefbin_lo = (xjrefbin==0 ? 1 : (xjrefbin-1)*xjrefRebinFactor+1);
     const int xjrefbin_hi = (xjrefbin==0 ? -1 : xjrefbin*xjrefRebinFactor);
     const double xjref_lo = (xjrefbin==0 ? 0.65 : xjrefbins[xjrefbin_lo - 1]);
     const double xjref_hi = (xjrefbin==0 ? 1.45 : xjrefbins[xjrefbin_hi]); 
     topPad->cd();
     gPad->SetLogx();
     zJetHist = zeeJetHists[periodType][etabin][0][1]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     zJetGraph_sys = new TGraphAsymmErrors((TH1F*)zJetHist); // for plotting systematics
     zJetHist->SetYTitle("#it{p}_{T}^{jet} / #it{p}_{T}^{ref}");
     zJetHist->SetAxisRange(xjref_lo, xjref_hi, "Y");
     zJetHist->SetMarkerColor(data_color);
     zJetHist->SetLineColor(data_color);
     zJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
     zJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
     zJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
     zJetHist->GetYaxis()->SetTitleOffset(uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     zJetHist_lo = zeeJetHists[periodType][etabin][0][0]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     zJetHist_hi = zeeJetHists[periodType][etabin][0][2]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     for (short pbin = 1; pbin <= numpbins; pbin++) {
      zJetGraph_sys->SetPoint(pbin-1, 0.5*(pbins[pbin-1]+pbins[pbin]), zJetHist->GetBinContent(pbin));
      zJetGraph_sys->SetPointEYhigh(pbin-1, TMath::Abs(zJetHist_hi->GetBinContent(pbin) - zJetHist->GetBinContent(pbin))); // set high systematics
      zJetGraph_sys->SetPointEYlow(pbin-1, TMath::Abs(zJetHist->GetBinContent(pbin) - zJetHist_lo->GetBinContent(pbin))); // set low systematics
     }
     if (zJetHist_lo) delete zJetHist_lo;
     if (zJetHist_hi) delete zJetHist_hi;
     zJetGraph_sys->SetFillColor(kBlack);
     zJetGraph_sys->SetFillStyle(3001);

     zJetHist_mc = zeeJetHists[periodType][etabin][1][1]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     zJetHist_mc->SetMarkerColor(mc_color);
     zJetHist_mc->SetLineColor(mc_color);

     zJetHist->DrawCopy("E1 X0");
     zJetHist_mc->DrawCopy("SAME E1 X0");
     zJetGraph_sys->Draw("2");
     if (xjrefbin == 0) {
      myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 Data, Cross-Calib Insitu (%i events)", Zee_n[periodType][0][etabin]), 1.25, 0.04/uPadY);
      myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("MC Signal Only (%i events)", Zee_n[periodType][1][etabin]), 1.25, 0.04/uPadY);
     }
     else {
      const float counts_data = zeeJetHists[periodType][etabin][0][1]->Integral(1, numpbins, xjrefbin_lo, xjrefbin_hi);
      const float total_data = zeeJetHists[periodType][etabin][0][1]->Integral();
      float n = 100. * counts_data / total_data; 
      myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 Data, Cross-Calib Insitu (%.1f%% of events)", n), 1.25, 0.04/uPadY);

      const float counts_mc = zeeJetHists[periodType][etabin][1][1]->Integral(1, numpbins, xjrefbin_lo, xjrefbin_hi);
      const float total_mc = zeeJetHists[periodType][etabin][1][1]->Integral();
      n = 100. * counts_mc / total_mc;
      myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("MC Signal Only (%.1f%% events)", n), 1.25, 0.04/uPadY);
     }
     if (periodType == 0) myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{A} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{B} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     if (xjrefbin != 0) {
      myText(0.155, 0.37, kBlack, "Z (ee) + Jet", 0.04/uPadY);
      myText(0.155, 0.28, kBlack, Form("%g < #it{x}_{J}^{ref} < %g", xjrefbins[xjrefbin_lo-1], xjrefbins[xjrefbin_hi]), 0.04/uPadY);
     }
     else myText(0.155, 0.28, kBlack, "Z (ee) + Jet", 0.04/uPadY);
     myText(0.155, 0.19, kBlack, period.c_str(), 0.04/uPadY);

     bottomPad->cd();
     gPad->SetLogx();

     zJetHist_rat = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_hist%i", etabin)), zeeJetHists[periodType][etabin][0][1], zeeJetHists[periodType][etabin][1][1], numpbins, pbins, numxjrefbins, xjrefbins, xjrefbin_lo, xjrefbin_hi);
     zJetGraph_rat_sys = new TGraphAsymmErrors(zJetHist_rat);
     zJetHist_rat_lo = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_lo_hist%i", etabin)), zeeJetHists[periodType][etabin][0][0], zeeJetHists[periodType][etabin][1][1], numpbins, pbins, numxjrefbins, xjrefbins, xjrefbin_lo, xjrefbin_hi);
     zJetHist_rat_hi = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_hi_hist%i", etabin)), zeeJetHists[periodType][etabin][0][2], zeeJetHists[periodType][etabin][1][1], numpbins, pbins, numxjrefbins, xjrefbins, xjrefbin_lo, xjrefbin_hi);
     for (short pbin = 1; pbin <= numpbins; pbin++) {
      zJetGraph_rat_sys->SetPoint(pbin-1, 0.5*(pbins[pbin-1]+pbins[pbin]), zJetHist_rat->GetBinContent(pbin));
      zJetGraph_rat_sys->SetPointEYlow(pbin-1, TMath::Abs(zJetHist_rat->GetBinContent(pbin) - zJetHist_rat_lo->GetBinContent(pbin)));
      zJetGraph_rat_sys->SetPointEYhigh(pbin-1, TMath::Abs(zJetHist_rat_hi->GetBinContent(pbin) - zJetHist_rat->GetBinContent(pbin)));
     }
     if (zJetHist_rat_lo) delete zJetHist_rat_lo;
     if (zJetHist_rat_hi) delete zJetHist_rat_hi;
     zJetGraph_rat_sys->SetFillColor(kBlack);
     zJetGraph_rat_sys->SetFillStyle(3001);

     zJetHist_rat->SetYTitle("Data / MC");
     zJetHist_rat->SetAxisRange(0.65, 1.45, "Y");
     zJetHist_rat->GetYaxis()->SetNdivisions(405);
     zJetHist_rat->GetXaxis()->SetTitleSize(0.04/dPadY);
     zJetHist_rat->GetYaxis()->SetTitleSize(0.04/dPadY);
     zJetHist_rat->GetXaxis()->SetTitleOffset(1);
     zJetHist_rat->GetYaxis()->SetTitleOffset(dPadY);
     zJetHist_rat->GetXaxis()->SetLabelSize(0.04/dPadY);
     zJetHist_rat->GetYaxis()->SetLabelSize(0.04/dPadY);
     zJetHist_rat->GetXaxis()->SetTickLength(0.08);

     zJetHist_rat->Draw("E1 X0");
     zJetGraph_rat_sys->Draw("2");
     for (TLine* line : zlines) line->Draw("SAME");
     char* plotName;
     if (xjrefbin == 0) plotName = Form("z_ee_jet%i.pdf", etabin);
     else plotName = Form("xjref/z_ee_jet%i_xjrefbin%i.pdf", etabin, xjrefbin);
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
     if (zJetHist_rat) delete zJetHist_rat;
     if (zJetGraph_rat_sys) delete zJetGraph_rat_sys;
    }


    /**** Plots GammaJet info as a function of p_T^ref****/
    for (int xjrefbin = 0; xjrefbin <= xjrefmax; xjrefbin++) {
     const int xjrefbin_lo = (xjrefbin==0 ? 1 : (xjrefbin-1)*xjrefRebinFactor+1);
     const int xjrefbin_hi = (xjrefbin==0 ? -1 : xjrefbin*xjrefRebinFactor);
     const double xjref_lo = (xjrefbin==0 ? 0.65 : xjrefbins[xjrefbin_lo - 1]);
     const double xjref_hi = (xjrefbin==0 ? 1.45 : xjrefbins[xjrefbin_hi]); 
     topPad->cd();
     gPad->SetLogx();
     gJetHist = gJetHists[periodType][etabin][0][1]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     gJetGraph_sys = new TGraphAsymmErrors((TH1F*)gJetHist); // for plotting systematics
     gJetHist->SetYTitle("#it{p}_{T}^{jet} / #it{p}_{T}^{ref}");
     gJetHist->SetAxisRange(xjref_lo, xjref_hi, "Y");
     gJetHist->SetMarkerColor(data_color);
     gJetHist->SetLineColor(data_color);
     gJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
     gJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
     gJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
     gJetHist->GetYaxis()->SetTitleOffset(uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     gJetHist_lo = gJetHists[periodType][etabin][0][0]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     gJetHist_hi = gJetHists[periodType][etabin][0][2]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     for (short prefbin = 1; prefbin <= numprefbins; prefbin++) {
      gJetGraph_sys->SetPoint(prefbin-1, 0.5*(prefbins[prefbin-1]+prefbins[prefbin]), gJetHist->GetBinContent(prefbin));
      gJetGraph_sys->SetPointEYhigh(prefbin-1, TMath::Abs(gJetHist_hi->GetBinContent(prefbin) - gJetHist->GetBinContent(prefbin))); // set high systematics
      gJetGraph_sys->SetPointEYlow(prefbin-1, TMath::Abs(gJetHist->GetBinContent(prefbin) - gJetHist_lo->GetBinContent(prefbin))); // set low systematics
     }
     gJetGraph_sys->SetFillColor(kBlack);
     gJetGraph_sys->SetFillStyle(3001);

     gJetHist_mc = gJetHists[periodType][etabin][1][1]->ProfileX("_pfx", xjrefbin_lo, xjrefbin_hi);
     gJetHist_mc->SetMarkerColor(mc_color);
     gJetHist_mc->SetLineColor(mc_color);

     for (short errType = 0; errType < 3; errType++) {
      const string error = (errType == 0 ? "syslo" : (errType == 1 ? "stat" : "syshi"));
      string periodStr = "periodA";
      if (periodType == 1) periodStr = "periodB";
      else if (periodType == 2) periodStr = "periodAB";
      gJetHistDifference[periodType][etabin][errType] = new TH1F(Form("gJetPtRatio_diff%i_%s_%s", etabin, error.c_str(), periodStr.c_str()), ";#it{p}_{T}^{ref} #left[GeV#right]", numprefbins, prefbins);
      for (short prefbin = 1; prefbin <= numprefbins; prefbin++) {
       double dataVal, dataErr;
       switch (errType) {
        case 0:
         dataVal = gJetHist_lo->GetBinContent(prefbin);
         dataErr = gJetHist_lo->GetBinError(prefbin);
         break;
        case 2:
         dataVal = gJetHist_hi->GetBinContent(prefbin);
         dataErr = gJetHist_hi->GetBinError(prefbin);
         break;
        default:
         dataVal = gJetHist->GetBinContent(prefbin);
         dataErr = gJetHist->GetBinError(prefbin);
       } 
       gJetHistDifference[periodType][etabin][errType]->SetBinContent(prefbin, dataVal - gJetHist_mc->GetBinContent(prefbin));
       gJetHistDifference[periodType][etabin][errType]->SetBinError(prefbin, TMath::Sqrt(TMath::Power(dataErr,2) + TMath::Power(gJetHist_mc->GetBinError(prefbin),2)));
      }
      gJetHistDifference[periodType][etabin][errType]->Write();
     }
     if (gJetHist_lo) delete gJetHist_lo;
     if (gJetHist_hi) delete gJetHist_hi;

     gJetHist->DrawCopy("E1 X0");
     gJetHist_mc->DrawCopy("SAME E1 X0");
     gJetGraph_sys->Draw("2");
     if (xjrefbin == 0) {
      myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 Data, Cross-Calib Insitu (%i events)", g_n[periodType][0][etabin]), 1.25, 0.04/uPadY);
      myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("MC %s (%i events)", (runValidation ? "Signal Only":"with Data Overlay"), g_n[periodType][1][etabin]), 1.25, 0.04/uPadY);
     }
     else {
      const float counts_data = gJetHists[periodType][etabin][0][1]->Integral(1, numprefbins, xjrefbin_lo, xjrefbin_hi);
      const float total_data = gJetHists[periodType][etabin][0][1]->Integral();
      float n = 100. * counts_data / total_data;
      myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 Data, Cross-Calib Insitu (%.1f%% of events)", n), 1.25, 0.04/uPadY);

      const float counts_mc = gJetHists[periodType][etabin][1][1]->Integral(1, numprefbins, xjrefbin_lo, xjrefbin_hi);
      const float total_mc = gJetHists[periodType][etabin][1][1]->Integral();
      n = 100. * counts_mc / total_mc;
      myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("MC %s (%.1f%% of events)", (runValidation ? "Signal Only":"with Data Overlay"), n), 1.25, 0.04/uPadY);
     }
     if (periodType == 0) myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{A} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{jet}^{B} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
     if (xjrefbin != 0) {
      myText(0.155, 0.37, kBlack, "#gamma + Jet", 0.04/uPadY);
      myText(0.155, 0.28, kBlack, Form("%g < #it{x}_{J}^{ref} < %g", xjrefbins[xjrefbin_lo-1], xjrefbins[xjrefbin_hi]), 0.04/uPadY);
     }
     else myText(0.155, 0.28, kBlack, "#gamma + Jet", 0.04/uPadY);
     myText(0.155, 0.19, kBlack, period.c_str(), 0.04/uPadY);

     bottomPad->cd();
     gPad->SetLogx();
     gJetHist_rat = GetDataOverMC(TString(Form("gJetPtDataMCRatio_hist%i", etabin)), gJetHists[periodType][etabin][0][1], gJetHists[periodType][etabin][1][1], numprefbins, prefbins, numxjrefbins, xjrefbins, xjrefbin_lo, xjrefbin_hi);
     gJetGraph_rat_sys = new TGraphAsymmErrors(gJetHist_rat);
     gJetHist_rat_lo = GetDataOverMC(TString(Form("gJetPtDataMCRatio_lo_hist%i", etabin)), gJetHists[periodType][etabin][0][0], gJetHists[periodType][etabin][1][1], numprefbins, prefbins, numxjrefbins, xjrefbins, xjrefbin_lo, xjrefbin_hi);
     gJetHist_rat_hi = GetDataOverMC(TString(Form("gJetPtDataMCRatio_hi_hist%i", etabin)), gJetHists[periodType][etabin][0][2], gJetHists[periodType][etabin][1][1], numprefbins, prefbins, numxjrefbins, xjrefbins, xjrefbin_lo, xjrefbin_hi);
     for (short prefbin = 1; prefbin <= numprefbins; prefbin++) {
      gJetGraph_rat_sys->SetPoint(prefbin-1, 0.5*(prefbins[prefbin-1]+prefbins[prefbin]), gJetHist_rat->GetBinContent(prefbin));
      gJetGraph_rat_sys->SetPointEYlow(prefbin-1, TMath::Abs(gJetHist_rat->GetBinContent(prefbin) - gJetHist_rat_lo->GetBinContent(prefbin)));
      gJetGraph_rat_sys->SetPointEYhigh(prefbin-1, TMath::Abs(gJetHist_rat_hi->GetBinContent(prefbin) - gJetHist_rat->GetBinContent(prefbin)));
     }
     if (gJetHist_rat_lo) delete gJetHist_rat_lo;
     if (gJetHist_rat_hi) delete gJetHist_rat_hi;
     gJetGraph_rat_sys->SetFillColor(kBlack);
     gJetGraph_rat_sys->SetFillStyle(3001);

     gJetHist_rat->SetYTitle("Data / MC");
     gJetHist_rat->SetAxisRange(0.65, 1.45, "Y");
     gJetHist_rat->GetYaxis()->SetNdivisions(405);
     gJetHist_rat->GetXaxis()->SetTitleSize(0.04/dPadY);
     gJetHist_rat->GetYaxis()->SetTitleSize(0.04/dPadY);
     gJetHist_rat->GetXaxis()->SetTitleOffset(1);
     gJetHist_rat->GetYaxis()->SetTitleOffset(dPadY);
     gJetHist_rat->GetXaxis()->SetLabelSize(0.04/dPadY);
     gJetHist_rat->GetYaxis()->SetLabelSize(0.04/dPadY);
     gJetHist_rat->GetXaxis()->SetTickLength(0.08);

     gJetHist_rat->Draw("e1 X0"); 
     gJetGraph_rat_sys->Draw("2");
     for (TLine* line : glines) line->Draw("SAME");
     char* plotName;
     if (xjrefbin == 0) plotName = Form("gamma_jet%i.pdf", etabin);
     else plotName = Form("xjref/gamma_jet%i_xjrefbin%i.pdf", etabin, xjrefbin);
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
     if (gJetHist_rat) delete gJetHist_rat;
     for (short errType = 0; errType < 3; errType++)
      if (gJetHistDifference[periodType][etabin][errType])
       delete gJetHistDifference[periodType][etabin][errType];
    }
   }
  }

  for (short periodType = 0; periodType < 2; periodType++) {
   const string period = (periodType == 0 ? "Period A":"Period B");

   for (short etabin = 0; etabin < numetabins; etabin++) {
    topPad->cd();
    TH2F* thisHist = gJetHistsSys[periodType][etabin][0][1];
    TH1F* rmsHist = new TH1F(Form("rms_etabin%i_%s", etabin, period.c_str()), "", numpbins, pbins);
    for (short pbin = 0; pbin < numpbins; pbin++) {
     //const float pt = 0.5*(pbins[pbin+1]+pbins[pbin]);
     float rms = 0;
     float sumWeights = 0;
     for (short sigbin = 0; sigbin < numSigmaBins; sigbin++) {
      const float sig = thisHist->GetYaxis()->GetBinCenter(sigbin+1);
      //const float sig = -maxSigma + 2*maxSigma*(sigbin+0.5) / numSigmaBins;
      const float weight = thisHist->GetBinContent(pbin+1, sigbin+1);
      rms += pow(sig, 2) * weight;
      sumWeights += weight;
     }
     rms = sqrt(rms) / sqrt(sumWeights);
     rmsHist->SetBinContent(pbin+1, rms);
    }
    //gPad->SetLogz();
    thisHist->Draw("col");
    thisHist->GetXaxis()->SetLabelSize(0.04/uPadY);
    thisHist->GetYaxis()->SetLabelSize(0.04/uPadY);
    thisHist->GetYaxis()->SetTitleSize(0.04/uPadY);
    thisHist->GetYaxis()->SetTitleOffset(1.1*uPadY);
    myText(0.72, 0.89, kBlack, period.c_str(), 0.04/uPadY);
    myText(0.72, 0.8,kBlack, Form("%g < #eta_{jet}^{Lab} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);

    bottomPad->cd();
    rmsHist->SetXTitle("#it{p}_{T}^{jet} #left[GeV#right]");
    rmsHist->SetYTitle("RMS(#sigma)/#it{p}_{T}^{jet}");
    rmsHist->SetAxisRange(0, 0.09, "Y");
    rmsHist->GetYaxis()->SetNdivisions(405);
    rmsHist->GetXaxis()->SetTitleSize(0.04/dPadY);
    rmsHist->GetYaxis()->SetTitleSize(0.04/dPadY);
    rmsHist->GetXaxis()->SetTitleOffset(1);
    rmsHist->GetYaxis()->SetTitleOffset(1.1*dPadY);
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

  TLine* lines[8];
  for (short i = 0; i < 8; i++) {
   if (glines[i]) delete glines[i];
   if (zlines[i]) delete zlines[i];
   lines[i] = new TLine(60, 0.6+0.2*i, 110, 0.6+0.2*i);
   if (0.6+0.2*i == 1) lines[i]->SetLineStyle(1);
   else lines[i]->SetLineStyle(3);
  }

  // Plot mumu mass spectra
  for (int etabin = 0; etabin < numetabins+1; etabin++) {
   for (int species = 0; species < 2; species++) {
    topPad->cd();
    gPad->SetLogx();
    double mean[2] = {};
    double mean_err[2] = {};
    double sigma[2] = {};
    double sigma_err[2] = {};
    TF1* fits[2];
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
     if (etabin != numetabins) myText(0.175, 0.15,kBlack, Form("%g < #eta_{jet} < %g", etabins[etabin], etabins[etabin+1]), 0.04/uPadY);
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
    gPad->SetLogx();
    TH1F* thisHist = (TH1F*)zMassSpectra[species][0][etabin]->Clone(Form("invMass_species%i_clone", species));
    thisHist->Divide(zMassSpectra[species][1][etabin]);
    thisHist->GetXaxis()->SetTitle("#font[12]{ll} Invariant Mass #left[GeV#right]");
    thisHist->GetYaxis()->SetTitle("Data / MC");
    thisHist->GetXaxis()->SetTitleSize(0.04/dPadY);
    thisHist->GetYaxis()->SetTitleSize(0.04/dPadY);
    thisHist->GetXaxis()->SetTitleOffset(1);
    thisHist->GetYaxis()->SetTitleOffset(1.1*dPadY);
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
    gPad->SetLogx();
    gPad->SetLogz();
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
    gPad->SetLogx();
    gPad->SetLogz();
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


}
