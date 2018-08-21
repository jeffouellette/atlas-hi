#include "../Params.h"
#include <Initialization.h>
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

  for (short iEta = 0; iEta <= numetabins; iEta++) {
   for (short iData = 0; iData < 2; iData++) { // iData is 0 for data, 1 for MC
    const TString dataType = (iData == 0 ? "data":"mc");

    for (short iErr = 0; iErr < 3; iErr++) {
     TString error = "sys_lo";
     if (iErr == 1) error = "stat";
     else if (iErr == 2) error = "sys_hi";

     for (short iPer = 0; iPer < 3; iPer++) {
      TString period = "periodA";
      if (iPer == 1) period = "periodB";
      else if (iPer == 2) period = "periodAB";

      zeeJetHists[iPer][iEta][iData][iErr] = new TH2D(Form("zeeJetPtRatio_iEta%i_%s_%s_%s", iEta, dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{Z} #left[GeV#right];#it{x}_{J}^{ref}", numpzbins, pzbins, numxjrefbins, xjrefbins);
      zeeJetHists[iPer][iEta][iData][iErr]->Sumw2();
      zmumuJetHists[iPer][iEta][iData][iErr] = new TH2D(Form("zmumuJetPtRatio_iEta%i_%s_%s_%s", iEta, dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{Z} #left[GeV#right];#it{x}_{J}^{ref}", numpzbins, pzbins, numxjrefbins, xjrefbins);
      zmumuJetHists[iPer][iEta][iData][iErr]->Sumw2();
      gJetHists[iPer][iEta][iData][iErr] = new TH2D(Form("gJetPtRatio_iEta%i_%s_%s_%s", iEta, dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{#gamma} #left[GeV#right];#it{x}_{J}^{ref}", numpbins, pbins, numxjrefbins, xjrefbins);
      gJetHists[iPer][iEta][iData][iErr]->Sumw2();

      gJetHistsSys[iPer][iEta][iData][iErr] = new TH2D(Form("gJetPtRatioSys_iEta%i_%s_%s_%s", iEta, dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{jet} #left[GeV#right];#Delta#it{x}_{J}^{ref}#it{p}_{T}^{ref}/#it{p}_{T}^{jet}", numpzbins, pzbins, numSigmaBins, -maxSigma, maxSigma);
      gJetHistsSys[iPer][iEta][iData][iErr]->Sumw2();
     }
    }

    for (short iSpc = 0; iSpc < 2; iSpc++) {
     const TString species = (iSpc == 0 ? "mumu":"ee");

     zMassSpectra[iSpc][iData][iEta] = new TH1D(Form("z%sMassSpectrum_%s_iEta%i", species.Data(), dataType.Data(), iEta), "", 50, 60, 110);
     zMassSpectra[iSpc][iData][iEta]->Sumw2();
    }
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
   TVectorD *nZeeMassVec, *nZeeJetVec, *nZmumuMassVec, *nZmumuJetVec, *nGammaJetVec;
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
       const short iPer = (runNumber < 313500 ? 0 : 1);
       //infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%i", runNumber));
       nZeeMassVec = (TVectorD*)thisFile->Get(Form("nZeeMassVec_%i", runNumber));
       nZeeJetVec = (TVectorD*)thisFile->Get(Form("nZeeJetVec_%i", runNumber));
       nZmumuMassVec = (TVectorD*)thisFile->Get(Form("nZmumuMassVec_%i", runNumber));
       nZmumuJetVec = (TVectorD*)thisFile->Get(Form("nZmumuJetVec_%i", runNumber));
       nGammaJetVec = (TVectorD*)thisFile->Get(Form("nGammaJetVec_%i", runNumber));

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        //const bool flipEta = runNumber < 313500 && iEta < numetabins;
        const bool flipEta = false;
        const short act_iEta = (flipEta ? (numetabins - iEta - 1) : iEta);

        nZeeMass[iPer][0][iEta] += (*nZeeMassVec)[iEta];
        nZeeMass[2][0][iEta] += (*nZeeMassVec)[act_iEta];

        nZeeJet[iPer][0][iEta] += (*nZeeJetVec)[iEta];
        nZeeJet[2][0][iEta] += (*nZeeJetVec)[act_iEta];

        nZmumuMass[iPer][0][iEta] += (*nZmumuMassVec)[iEta];
        nZmumuMass[2][0][iEta] += (*nZmumuMassVec)[act_iEta];

        nZmumuJet[iPer][0][iEta] += (*nZmumuJetVec)[iEta];
        nZmumuJet[2][0][iEta] += (*nZmumuJetVec)[act_iEta];

        nGammaJet[iPer][0][iEta] += (*nGammaJetVec)[iEta];
        nGammaJet[2][0][iEta] += (*nGammaJetVec)[act_iEta];

        for (short iErr = 0; iErr < 3; iErr++) {
         TString error = "sys_lo";
         if (iErr == 1) error = "stat";
         else if (iErr == 2) error = "sys_hi";

         TH2D* temp = (TH2D*)thisFile->Get(Form("zeeJetPtRatio_dataSet%i_iEta%i_data_%s", runNumber, iEta, error.Data()));
         zeeJetHists[iPer][iEta][0][iErr]->Add(temp);
         temp = (TH2D*)thisFile->Get(Form("zeeJetPtRatio_dataSet%i_iEta%i_data_%s", runNumber, act_iEta, error.Data()));
         zeeJetHists[2][iEta][0][iErr]->Add(temp);

         temp = (TH2D*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%i_iEta%i_data_%s", runNumber, iEta, error.Data()));
         zmumuJetHists[iPer][iEta][0][iErr]->Add(temp);
         temp = (TH2D*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%i_iEta%i_data_%s", runNumber, act_iEta, error.Data()));
         zmumuJetHists[2][iEta][0][iErr]->Add(temp);

         temp = (TH2D*)thisFile->Get(Form("gJetPtRatio_dataSet%i_iEta%i_data_%s", runNumber, iEta, error.Data()));
         gJetHists[iPer][iEta][0][iErr]->Add(temp);
         temp = (TH2D*)thisFile->Get(Form("gJetPtRatio_dataSet%i_iEta%i_data_%s", runNumber, act_iEta, error.Data()));
         gJetHists[2][iEta][0][iErr]->Add(temp);

         if (iErr == 1) 
          gJetHistsSys[iPer][iEta][0][iErr]->Add((TH2D*)thisFile->Get(Form("gJetPtRatioSys_dataSet%i_iEta%i_data_%s", runNumber, iEta, error.Data())));
         //gJetHistsSys[2][iEta][0][iErr]->Add((TH2D*)thisFile->Get(Form("gJetPtRatioSys_dataSet%i_iEta%i_data_%s", runNumber, act_iEta, error.Data())));
        }
       }

       for (short iSpc = 0; iSpc < 2; iSpc++) {
        const TString species = (iSpc == 0 ? "mumu":"ee");

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         zMassSpectra[iSpc][0][iEta]->Add((TH1D*)thisFile->Get(Form("z%sMassSpectrum_dataSet%i_data_iEta%i", species.Data(), runNumber, iEta)));
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
       const short iPer = (gammaJetSampleId.Contains("pPb") ? 0 : 1);
       //infoVec = (TVectorD*)thisFile->Get(Form("infoVec_%s", gammaJetSampleId.Data()));
       nGammaJetVec = (TVectorD*)thisFile->Get(Form("nGammaJetVec_%s", gammaJetSampleId.Data()));

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        //const bool flipEta = gammaJetSampleId.Contains("pPb") && iEta < numetabins;
        const bool flipEta = false;
        const short act_iEta = (flipEta ? (numetabins - iEta - 1) : iEta); // period A condition

        nGammaJet[iPer][1][iEta] += (*nGammaJetVec)[iEta];
        nGammaJet[2][1][iEta] += (*nGammaJetVec)[act_iEta];

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH2D* temp = (TH2D*)thisFile->Get(Form("gJetPtRatio_dataSet%s_iEta%i_mc_stat", gammaJetSampleId.Data(), iEta));
        gJetHists[iPer][iEta][1][1]->Add(temp);
        temp = (TH2D*)thisFile->Get(Form("gJetPtRatio_dataSet%s_iEta%i_mc_stat", gammaJetSampleId.Data(), act_iEta));
        gJetHists[2][iEta][1][1]->Add(temp);

        gJetHistsSys[iPer][iEta][1][1]->Add((TH2D*)thisFile->Get(Form("gJetPtRatioSys_dataSet%s_iEta%i_mc_stat", gammaJetSampleId.Data(), iEta)));
        //gJetHistsSys[2][iEta][1][1]->Add((TH2D*)thisFile->Get(Form("gJetPtRatioSys_dataSet%s_iEta%i_mc_stat", gammaJetSampleId.Data(), act_iEta)));
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
       const short iPer = (zeeJetSampleId.Contains("pPb") ? 0 : 1);
       nZeeMassVec = (TVectorD*)thisFile->Get(Form("nZeeMassVec_%s", zeeJetSampleId.Data()));
       nZeeJetVec = (TVectorD*)thisFile->Get(Form("nZeeJetVec_%s", zeeJetSampleId.Data()));

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        //const bool flipEta = zeeJetSampleId.Contains("pPb") && iEta < numetabins;
        const bool flipEta = false;
        const short act_iEta = (flipEta ? numetabins - iEta - 1 : iEta); // period A condition

        nZeeMass[iPer][1][iEta] += (*nZeeMassVec)[iEta];
        nZeeMass[2][1][iEta] += (*nZeeMassVec)[act_iEta];
        nZeeJet[iPer][1][iEta] += (*nZeeJetVec)[iEta];
        nZeeJet[2][1][iEta] += (*nZeeJetVec)[act_iEta];

        // Only add the statistical error plots for MC (don't need to consider systematics)
        TH2D* temp = (TH2D*)thisFile->Get(Form("zeeJetPtRatio_dataSet%s_iEta%i_mc_stat", zeeJetSampleId.Data(), iEta));
        zeeJetHists[iPer][iEta][1][1]->Add(temp);
        temp = (TH2D*)thisFile->Get(Form("zeeJetPtRatio_dataSet%s_iEta%i_mc_stat", zeeJetSampleId.Data(), act_iEta));
        zeeJetHists[2][iEta][1][1]->Add(temp);

       }

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        zMassSpectra[1][1][iEta]->Add((TH1D*)thisFile->Get(Form("zeeMassSpectrum_dataSet%s_mc_iEta%i", zeeJetSampleId.Data(), iEta)));
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
       const short iPer = (zmumuJetSampleId.Contains("pPb") ? 0 : 1);
       nZmumuMassVec = (TVectorD*)thisFile->Get(Form("nZmumuMassVec_%s", zmumuJetSampleId.Data()));
       nZmumuJetVec = (TVectorD*)thisFile->Get(Form("nZmumuJetVec_%s", zmumuJetSampleId.Data()));

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        //const bool flipEta = zmumuJetSampleId.Contains("pPb") && iEta < numetabins;
        const bool flipEta = false;
        const short act_iEta = (flipEta ? numetabins - iEta - 1 : iEta); // period A condition

        nZmumuMass[iPer][1][iEta] += (*nZmumuMassVec)[iEta];
        nZmumuMass[2][1][iEta] += (*nZmumuMassVec)[act_iEta];
        nZmumuJet[iPer][1][iEta] += (*nZmumuJetVec)[iEta];
        nZmumuJet[2][1][iEta] += (*nZmumuJetVec)[act_iEta];

        // Only add the statistical error plots for MC (don't need to
        // consider systematics)
        TH2D* temp = (TH2D*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%s_iEta%i_mc_stat", zmumuJetSampleId.Data(), iEta));
        zmumuJetHists[iPer][iEta][1][1]->Add(temp);
        temp = (TH2D*)thisFile->Get(Form("zmumuJetPtRatio_dataSet%s_iEta%i_mc_stat", zmumuJetSampleId.Data(), act_iEta));
        zmumuJetHists[2][iEta][1][1]->Add(temp);

       }

       for (short iEta = 0; iEta <= numetabins; iEta++) {
        zMassSpectra[0][1][iEta]->Add((TH1D*)thisFile->Get(Form("zmumuMassSpectrum_dataSet%s_mc_iEta%i", zmumuJetSampleId.Data(), iEta)));
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

  for (short iPer = 0; iPer < 3; iPer++) {

   TString period = "Period A";
   if (iPer == 1) period = "Period B";
   else if (iPer == 2) period = "Period A+B";

   for (short iEta = 0; iEta <= numetabins; iEta++) {

    /**** Plot ZmumuJet info ****/
    topPad->cd();
    topPad->SetLogx();
    vJetHist = GetProfileX("vJetHist", zmumuJetHists[iPer][iEta][0][1], numpzbins, pzbins, false);
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
    vJetHist_lo = GetProfileX("vJetHist_lo", zmumuJetHists[iPer][iEta][0][0], numpzbins, pzbins, false);
    vJetHist_hi = GetProfileX("vJetHist_hi", zmumuJetHists[iPer][iEta][0][2], numpzbins, pzbins, false);
    CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    if (vJetHist_lo) delete vJetHist_lo;
    if (vJetHist_hi) delete vJetHist_hi;
    vJetGraph_sys->SetFillColor(kBlack);
    vJetGraph_sys->SetFillStyle(3001);

    vJetHist_mc = GetProfileX("vJetHist_mc", zmumuJetHists[iPer][iEta][1][1], numpzbins, pzbins, false);
    vJetHist_mc->SetMarkerColor(mc_color);
    vJetHist_mc->SetLineColor(mc_color);

    vJetHist->DrawCopy("E1 X0");
    vJetHist_mc->DrawCopy("SAME E1 X0");
    vJetGraph_sys->Draw("2");

    myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 #it{p}+Pb 8 TeV, with Insitu Corrections (%i events)", nZmumuJet[iPer][0][iEta]), 1.25, 0.04/uPadY);
    myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV with Overlay (%i events)", nZmumuJet[iPer][1][iEta]), 1.25, 0.04/uPadY);
    if (iEta < numetabins) {
     if (iPer == 2) myText(0.155, 0.1,kBlack, Form("%g < #eta_{Lab}^{#mu#mu} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{Lab}^{#mu#mu} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
    }
    myText(0.155, 0.28, kBlack, "Z (#mu#mu) + Jet", 0.04/uPadY);
    myText(0.155, 0.19, kBlack, period.Data(), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx();

    vJetHist_rat = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_iEta%i", iEta)), zmumuJetHists[iPer][iEta][0][1], zmumuJetHists[iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
    vJetHist_rat_lo = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_lo_iEta%i", iEta)), zmumuJetHists[iPer][iEta][0][0], zmumuJetHists[iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetHist_rat_hi = GetDataOverMC(TString(Form("zmumuJetPtDataMCRatio_hi_iEta%i", iEta)), zmumuJetHists[iPer][iEta][0][2], zmumuJetHists[iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
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
    if (iEta < numetabins) plotName = Form("z_mumu_jet%i.pdf", iEta);
    else plotName = Form("z_mumu_jet_combined.pdf");

    switch (iPer) {
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
    vJetHist = GetProfileX("vJetHist", zeeJetHists[iPer][iEta][0][1], numpzbins, pzbins, false);
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
    vJetHist_lo = GetProfileX("vJetHist_lo", zeeJetHists[iPer][iEta][0][0], numpzbins, pzbins, false);
    vJetHist_hi = GetProfileX("vJetHist_hi", zeeJetHists[iPer][iEta][0][2], numpzbins, pzbins, false);
    CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    if (vJetHist_lo) delete vJetHist_lo;
    if (vJetHist_hi) delete vJetHist_hi;
    vJetGraph_sys->SetFillColor(kBlack);
    vJetGraph_sys->SetFillStyle(3001);

    vJetHist_mc = GetProfileX("vJetHist_mc", zeeJetHists[iPer][iEta][1][1], numpzbins, pzbins, false);
    vJetHist_mc->SetMarkerColor(mc_color);
    vJetHist_mc->SetLineColor(mc_color);

    vJetHist->DrawCopy("E1 X0");
    vJetHist_mc->DrawCopy("SAME E1 X0");
    vJetGraph_sys->Draw("2");

    myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 #it{p}+Pb 8 TeV, with Insitu Corrections (%i events)", nZeeJet[iPer][0][iEta]), 1.25, 0.04/uPadY);
    myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV with Overlay (%i events)", nZeeJet[iPer][1][iEta]), 1.25, 0.04/uPadY);
    if (iEta < numetabins) {
     if (iPer == 2) myText(0.155, 0.1,kBlack, Form("%g < #eta_{Lab}^{ee} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{Lab}^{ee} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
    }
    myText(0.155, 0.28, kBlack, "Z (ee) + Jet", 0.04/uPadY);
    myText(0.155, 0.19, kBlack, period.Data(), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx();

    vJetHist_rat = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_iEta%i", iEta)), zeeJetHists[iPer][iEta][0][1], zeeJetHists[iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
    vJetHist_rat_lo = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_lo_iEta%i", iEta)), zeeJetHists[iPer][iEta][0][0], zeeJetHists[iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
    vJetHist_rat_hi = GetDataOverMC(TString(Form("zeeJetPtDataMCRatio_hi_iEta%i", iEta)), zeeJetHists[iPer][iEta][0][2], zeeJetHists[iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
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
    if (iEta < numetabins) plotName = Form("z_ee_jet%i.pdf", iEta);
    else plotName = Form("z_ee_jet_combined.pdf");
    switch (iPer) {
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
    vJetHist = GetProfileX("vJetHist", gJetHists[iPer][iEta][0][1], numpbins, pbins, true);
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
    vJetHist_lo = GetProfileX("vJetHist_lo", gJetHists[iPer][iEta][0][0], numpbins, pbins, true);
    vJetHist_hi = GetProfileX("vJetHist_hi", gJetHists[iPer][iEta][0][2], numpbins, pbins, true);
    CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
    vJetGraph_sys->SetFillColor(kBlack);
    vJetGraph_sys->SetFillStyle(3001);

    vJetHist_mc = GetProfileX("vJetHist_mc", gJetHists[iPer][iEta][1][1], numpbins, pbins, true);
    vJetHist_mc->SetMarkerColor(mc_color);
    vJetHist_mc->SetLineColor(mc_color);

    for (short iErr = 0; iErr < 3; iErr++) {
     const TString error = (iErr == 0 ? "syslo" : (iErr == 1 ? "stat" : "syshi"));
     TString periodStr = "periodA";
     if (iPer == 1) periodStr = "periodB";
     else if (iPer == 2) periodStr = "periodAB";
     gJetHistDifference[iPer][iEta][iErr] = new TH1D(Form("gJetPtRatio_diff%i_%s_%s", iEta, error.Data(), periodStr.Data()), ";#it{p}_{T}^{ref} #left[GeV#right]", numpbins, pbins);
     for (short iP = 1; iP <= numpbins; iP++) {
      double dataVal, dataErr;
      switch (iErr) {
       case 0:
        dataVal = vJetHist_lo->GetBinContent(iP);
        dataErr = vJetHist_lo->GetBinError(iP);
        break;
       case 2:
        dataVal = vJetHist_hi->GetBinContent(iP);
        dataErr = vJetHist_hi->GetBinError(iP);
        break;
       default:
        dataVal = vJetHist->GetBinContent(iP);
        dataErr = vJetHist->GetBinError(iP);
      } 
      gJetHistDifference[iPer][iEta][iErr]->SetBinContent(iP, dataVal - vJetHist_mc->GetBinContent(iP));
      gJetHistDifference[iPer][iEta][iErr]->SetBinError(iP, TMath::Sqrt(TMath::Power(dataErr,2) + TMath::Power(vJetHist_mc->GetBinError(iP),2)));
     }
     gJetHistDifference[iPer][iEta][iErr]->Write();
    }
    if (vJetHist_lo) delete vJetHist_lo;
    if (vJetHist_hi) delete vJetHist_hi;

    vJetHist->DrawCopy("E1 X0");
    vJetHist_mc->DrawCopy("SAME E1 X0");
    vJetGraph_sys->Draw("2");
    for (TLine* line : dplines) line->Draw("same");

    myMarkerText(0.175, 0.88, data_color, kFullCircle, Form("2016 #it{p}+Pb 8 TeV, with Insitu Corrections (%i events)", nGammaJet[iPer][0][iEta]), 1.25, 0.04/uPadY);
    myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV with Overlay (%i events)", nGammaJet[iPer][1][iEta]), 1.25, 0.04/uPadY);
    if (iEta < numetabins) {
     if (iPer == 2) myText(0.155, 0.1,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
     else myText(0.155, 0.1,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
    }
    myText(0.155, 0.28, kBlack, "#gamma + Jet", 0.04/uPadY);
    myText(0.155, 0.19, kBlack, period.Data(), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx();
    vJetHist_rat = GetDataOverMC(TString(Form("gJetPtDataMCRatio_iEta%i", iEta)), gJetHists[iPer][iEta][0][1], gJetHists[iPer][iEta][1][1], numpbins, pbins, numxjrefbins, xjrefbins, true);
    vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
    vJetHist_rat_lo = GetDataOverMC(TString(Form("gJetPtDataMCRatio_lo_iEta%i", iEta)), gJetHists[iPer][iEta][0][0], gJetHists[iPer][iEta][1][1], numpbins, pbins, numxjrefbins, xjrefbins, true);
    vJetHist_rat_hi = GetDataOverMC(TString(Form("gJetPtDataMCRatio_hi_iEta%i", iEta)), gJetHists[iPer][iEta][0][2], gJetHists[iPer][iEta][1][1], numpbins, pbins, numxjrefbins, xjrefbins, true);
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

    if (iEta < numetabins) plotName = Form("gamma_jet%i.pdf", iEta);
    else plotName = Form("gamma_jet_combined.pdf");
    switch (iPer) {
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
    for (short iErr = 0; iErr < 3; iErr++)
     if (gJetHistDifference[iPer][iEta][iErr])
      delete gJetHistDifference[iPer][iEta][iErr];

    if (iEta == numetabins) continue;


    if (!plot_xjref) continue;
    /**** Plots xjref distributions, binned by ptref ****/
    for (short iP = 0; iP < numpbins; iP++) {
     const double pref_lo = pbins[iP];
     const double pref_hi =  pbins[iP+1];
     topPad->cd();
     topPad->SetLogx(0);
     vJetHist = gJetHists[iPer][iEta][0][1]->ProjectionY("vJetProjection", iP, iP);
     const float counts_data = vJetHist->Integral();
     const float total_data = gJetHists[iPer][iEta][0][1]->Integral();
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
     vJetHist_lo = gJetHists[iPer][iEta][0][0]->ProjectionY("vJetProjection_lo", iP, iP);
     vJetHist_lo->Rebin(rebinFactor);
     vJetHist_lo->Scale(1./vJetHist_lo->Integral()); 
     //vJetHist_lo->Scale(1./counts_data); 
     vJetHist_hi = gJetHists[iPer][iEta][0][2]->ProjectionY("vJetProjection_hi", iP, iP);
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

     vJetHist_mc = gJetHists[iPer][iEta][1][1]->ProjectionY("vJetProjection_mc", iP, iP);
     const float counts_mc = vJetHist_mc->Integral();
     const float total_mc = gJetHists[iPer][iEta][1][1]->Integral();
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
     if (iPer == 2) myText(0.155, 0.16,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
     else myText(0.155, 0.16,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);

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
     plotName = Form("pref_slices/gamma_jet%i_iP%i.pdf", iEta, iP);
     switch (iPer) {
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
  for (short iPer = 0; iPer < 2; iPer++) {
   const TString period = (iPer == 0 ? "Period A":"Period B");

   for (short iEta = 0; iEta < numetabins; iEta++) {
    topPad->cd();
    TH2D* thisHist = gJetHistsSys[iPer][iEta][0][1];
    TH1D* rmsHist = new TH1D(Form("rms_iEta%i_%s", iEta, (iPer==0?"pPb":"Pbp")), "", numpzbins, pzbins);
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
    myText(0.72, 0.8,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);

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
    canvas->SaveAs(Form("%s/Period%s/jetSystematics_iEta%i.pdf", plotPath.Data(), (iPer==0 ? "A":"B"), iEta));
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
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   if (iEta != numetabins) continue;
   for (short iSpc = 0; iSpc < 2; iSpc++) {
    topPad->cd();
    topPad->SetLogx(0);
    double mean[2] = {};
    double mean_err[2] = {};
    double sigma[2] = {};
    double sigma_err[2] = {};
    TF1* fits[2];
    //RooZfit* fits[2];
    for (short iData = 0; iData < 2; iData++) {
     TH1D* thisHist = zMassSpectra[iSpc][iData][iEta];
     Color_t color = (iData==0 ? data_color : mc_color);
     thisHist->GetXaxis()->SetTitle("#font[12]{ll} Invariant Mass #left[GeV#right]");
     thisHist->GetYaxis()->SetTitle("Normalized Counts / 1 GeV");
     thisHist->GetYaxis()->SetTitleSize(0.04/uPadY);
     thisHist->GetYaxis()->SetTitleOffset(1.1*uPadY);
     thisHist->GetYaxis()->SetLabelSize(0.04/uPadY);
     thisHist->SetLineColor(color);
     thisHist->SetMarkerColor(color);
//     thisHist->GetXaxis()->SetNdivisions(50802, false);
     thisHist->Scale(1.0, "width");

     //TF1* gausFit = new TF1(Form("fit_guess_iSpc%i_iData%i", iSpc, iData), "gaus(0)", Z_mass - 10, Z_mass + 10);
     //thisHist->Fit(gausFit, "R", "L");
     //double m = gausFit->GetParameter(1);
     //double s = gausFit->GetParameter(2);
     //const double scale = 1.0 / thisHist->Integral(thisHist->FindBin(m-Z_mass_fitNsigma*s), thisHist->FindBin(m+Z_mass_fitNsigma*s));
     //thisHist->Scale(scale);
     //
     //RooZfit* fit = new RooZfit(thisHist, color);
     //fits[iData] = fit;

     //mean[iData] = fit->mean->getVal();
     //mean_err[iData] = fit->mean->getError();
     //sigma[iData] = fit->sigma->getVal();
     //sigma_err[iData] = fit->sigma->getError();

     TF1* fit = new TF1(Form("fit_guess_iSpc%i_iData%i", iSpc, iData), "gaus(0)", Z_mass - 10, Z_mass + 10);
     thisHist->Fit(fit, "R", "L");
     double m = fit->GetParameter(1);
     double s = fit->GetParameter(2);

     TF1* fit_better = new TF1(Form("fit_better_iSpc%i_iData%i", iSpc, iData), "gaus(0)", m-Z_mass_fitNsigma*s, m+Z_mass_fitNsigma*s);
     thisHist->Fit(fit_better, "R", "L");
     m = fit_better->GetParameter(1);
     s = fit_better->GetParameter(2);
     mean[iData] = m;
     mean_err[iData] = fit_better->GetParError(1);
     sigma[iData] = s;
     sigma_err[iData] = fit_better->GetParError(2);
     const double scale = 1.0 / thisHist->Integral(thisHist->FindBin(m-Z_mass_fitNsigma*s), thisHist->FindBin(m+Z_mass_fitNsigma*s));
     thisHist->Scale(scale);

     fits[iData] = fit_better;
     fit_better->SetLineColor(color);
     fit_better->SetParameter(0, scale*fit_better->GetParameter(0));

     thisHist->SetAxisRange(0, 0.20, "Y");
     thisHist->GetYaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
    }
    for (short iData = 0; iData < 2; iData++) {
     TH1D* thisHist = zMassSpectra[iSpc][iData][iEta];
     if (iData == 0) thisHist->Draw("p");
     else thisHist->Draw("hist same");
     //fits[iData]->plot->Draw("same");
     fits[iData]->Draw("same");
     if (iEta < numetabins) {
      if (iSpc == 0)
       myText(0.175, 0.15,kBlack, Form("%g < #eta_{Lab}^{#mu#mu} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
      else if (iSpc == 1)
       myText(0.175, 0.15,kBlack, Form("%g < #eta_{Lab}^{ee} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
     }
    }
    if (iSpc == 0) {
     myText(0.175, 0.88, kBlack, "Z (#mu#mu) + Jet", 0.04/uPadY);
     myMarkerText(0.175, 0.80, data_color, kFullCircle, Form("2016 #it{p}+{Pb} 8 TeV (%i events)", nZmumuMass[2][0][iEta]), 1.25, 0.04/uPadY);
     myMarkerText(0.175, 0.55, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV with Overlay (%i events)", nZmumuMass[2][1][iEta]), 1.25, 0.04/uPadY);
    }
    else if (iSpc == 1) {
     myText(0.175, 0.88, kBlack, "Z (ee) + Jet", 0.04/uPadY);
     myMarkerText(0.175, 0.80, data_color, kFullCircle, Form("2016 #it{p}+{Pb} 8 TeV (%i events)", nZeeMass[2][0][iEta]), 1.25, 0.04/uPadY);
     myMarkerText(0.175, 0.55, mc_color, kFullCircle, Form("Pythia8 #it{pp} 8 TeV with Overlay (%i events)", nZeeMass[2][1][iEta]), 1.25, 0.04/uPadY);
    }
    myText(0.175, 0.72, kBlack, Form("m_{Z}^{data} = %.2f #pm %.2f GeV", mean[0], mean_err[0]), 0.04/uPadY);
    myText(0.175, 0.64, kBlack, Form("#sigma_{Z}^{data} = %.2f #pm %.2f GeV", sigma[0], sigma_err[0]), 0.04/uPadY);
    myText(0.175, 0.47, kBlack, Form("m_{Z}^{mc} = %.2f #pm %.2f GeV", mean[1], mean_err[1]), 0.04/uPadY);
    myText(0.175, 0.39, kBlack, Form("#sigma_{Z}^{mc} = %.2f #pm %.2f GeV", sigma[1], sigma_err[1]), 0.04/uPadY);

    bottomPad->cd();
    bottomPad->SetLogx(0);
    TH1D* thisHist = (TH1D*)zMassSpectra[iSpc][0][iEta]->Clone(Form("invMass_iSpc%i_clone", iSpc));
    thisHist->Divide(zMassSpectra[iSpc][1][iEta]);
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

    if (iEta != numetabins) {
     if (iSpc == 0) canvas->SaveAs(Form("%s/zmumu_mass_comparison_iEta%i.pdf", plotPath.Data(), iEta));
     else if (iSpc == 1) canvas->SaveAs(Form("%s/zee_mass_comparison_iEta%i.pdf", plotPath.Data(), iEta));
    }
    else {
     if (iSpc == 0) canvas->SaveAs(Form("%s/zmumu_mass_comparison.pdf", plotPath.Data()));
     else if (iSpc == 1) canvas->SaveAs(Form("%s/zee_mass_comparison.pdf", plotPath.Data()));
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
  for (short iPer = 0; iPer < 3; iPer++) {
   const TString period = (iPer == 0 ? "Period A": (iPer == 1 ? "Period B":"Period A+B"));

   for (short iEta = 0; iEta < numetabins; iEta++) {
    TH2D* dataHist = gJetHists[iPer][iEta][0][1];
    TH2D* mcHist = gJetHists[iPer][iEta][1][1];
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
    myText(0.1, 0.15, kBlack, Form("2016 #it{p}+Pb 8 TeV, with Insitu Corrections (%i events)", nGammaJet[iPer][0][iEta]), 0.02/rPadX);
    if (iPer != 2) myText(0.6, 0.85,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.02/rPadX);
    else myText(0.6, 0.85,kBlack, Form("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.02/rPadX);
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
    myText(0.2, 0.15, kBlack, Form("Pythia8 #it{pp} 8 TeV %s (%i events)", (runValidation ? "":"with Overlay"), nGammaJet[iPer][1][iEta]), 0.02/lPadX);

    const char* plotName = Form("gamma_jet%i_th2.pdf", iEta);
    switch (iPer) {
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


  return;
}
