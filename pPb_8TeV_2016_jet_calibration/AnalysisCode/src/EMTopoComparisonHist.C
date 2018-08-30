#include "EMTopoComparisonHist.h"

#include <GlobalParams.h>

#include <TF1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TVectorT.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include <AtlasUtils.h>
#include <AtlasStyle.h>

using namespace atlashi;

namespace pPb8TeV2016JetCalibration {

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
   TH1D* projy = data->ProjectionY(name + TString(Form ("data_xbin%i", xbin)), xbin, xbin);
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

   projy = mc->ProjectionY(name + TString(Form ("mc_xbin%i", xbin)), xbin, xbin);
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

void EMTopoComparisonHist () {

  SetAtlasStyle();

  // Setup trigger vectors
  SetupDirectories("EMTopoComparison/", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers(0);
  for (short i = 0; i < sizeof(full_run_list)/sizeof(full_run_list[0]); i++) runNumbers.push_back(full_run_list[i]);
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

  TH2D****** zeeJetHists = Get5DArray <TH2D*> (2, 3, numetabins+1, 2, 3);
  TH2D****** zmumuJetHists = Get5DArray <TH2D*> (2, 3, numetabins+1, 2, 3);
  TH2D****** gJetHists = Get5DArray <TH2D*> (2, 3, numetabins+1, 2, 3);

  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
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

       zeeJetHists[iAlgo][iPer][iEta][iData][iErr] = new TH2D(Form ("zeeJetPtRatio_iEta%i_%s_%s_%s_%s", iEta, algo.Data(), dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{Z} #left[GeV#right];#it{x}_{J}^{ref}", numpzbins, pzbins, numxjrefbins, xjrefbins);
       zeeJetHists[iAlgo][iPer][iEta][iData][iErr]->Sumw2();
       zmumuJetHists[iAlgo][iPer][iEta][iData][iErr] = new TH2D(Form ("zmumuJetPtRatio_iEta%i_%s_%s_%s_%s", iEta, algo.Data(), dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{Z} #left[GeV#right];#it{x}_{J}^{ref}", numpzbins, pzbins, numxjrefbins, xjrefbins);
       zmumuJetHists[iAlgo][iPer][iEta][iData][iErr]->Sumw2();
       gJetHists[iAlgo][iPer][iEta][iData][iErr] = new TH2D(Form ("gJetPtRatio_iEta%i_%s_%s_%s_%s", iEta, algo.Data(), dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{#gamma} #left[GeV#right];#it{x}_{J}^{ref}", numpbins, pbins, numxjrefbins, xjrefbins);
       gJetHists[iAlgo][iPer][iEta][iData][iErr]->Sumw2();

       //gJetHistsSys[iAlgo][iPer][iEta][iData][iErr] = new TH2D(Form ("gJetPtRatioSys_iEta%i_%s_%s_%s_%s", iEta, algo.Data(), dataType.Data(), error.Data(), period.Data()), ";#it{p}_{T}^{jet} #left[GeV#right];#Delta#it{x}_{J}^{ref}#it{p}_{T}^{ref}/#it{p}_{T}^{jet}", numpzbins, pzbins, numSigmaBins, -maxSigma, maxSigma);
       //gJetHistsSys[iAlgo][iPer][iEta][iData][iErr]->Sumw2();
      }
     }
    }
   }
  }

  int**** nZeeJet = Get4DArray <int> (2, 3, 2, numetabins+1);
  int**** nZmumuJet = Get4DArray <int> (2, 3, 2, numetabins+1);
  int**** nGammaJet = Get4DArray <int> (2, 3, 2, numetabins+1);

  int nTotalJets[2] = {};
  int nCleanJets[2] = {};

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
   TVectorD *nZeeJetVec, *nZmumuJetVec, *nGammaJetVec, *jetCleaningVec;
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
       //infoVec = (TVectorD*) thisFile->Get (Form ("infoVec_%i", runNumber));
       nZeeJetVec = (TVectorD*) thisFile->Get (Form ("nZeeJetVec_%i", runNumber));
       nZmumuJetVec = (TVectorD*) thisFile->Get (Form ("nZmumuJetVec_%i", runNumber));
       nGammaJetVec = (TVectorD*) thisFile->Get (Form ("nGammaJetVec_%i", runNumber));

       jetCleaningVec = (TVectorD*) thisFile->Get (Form ("jetCleaningVec_%i", runNumber));
       if (jetCleaningVec) {
        nTotalJets[0] += (int) (*jetCleaningVec)[0];
        nCleanJets[0] += (int) (*jetCleaningVec)[1];
       }

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         //const bool flipEta = runNumber < 313500 && iEta < numetabins;
         const bool flipEta = false;
         const short act_iEta = (flipEta ? (numetabins - iEta - 1) : iEta);

         nZeeJet[iAlgo][iPer][0][iEta] += (*nZeeJetVec)[iEta + iAlgo*(numetabins+1)];
         nZeeJet[iAlgo][2][0][iEta] += (*nZeeJetVec)[act_iEta + iAlgo*(numetabins+1)];

         nZmumuJet[iAlgo][iPer][0][iEta] += (*nZmumuJetVec)[iEta + iAlgo*(numetabins+1)];
         nZmumuJet[iAlgo][2][0][iEta] += (*nZmumuJetVec)[act_iEta + iAlgo*(numetabins+1)];

         nGammaJet[iAlgo][iPer][0][iEta] += (*nGammaJetVec)[iEta + iAlgo*(numetabins+1)];
         nGammaJet[iAlgo][2][0][iEta] += (*nGammaJetVec)[act_iEta + iAlgo*(numetabins+1)];

         for (short iErr = 0; iErr < 3; iErr++) {
          TString error = "sys_lo";
          if (iErr == 1) error = "stat";
          else if (iErr == 2) error = "sys_hi";

          TH2D* temp = (TH2D*) thisFile->Get (Form ("zeeJetPtRatio_dataSet%i_%s_iEta%i_data_%s", runNumber, algo.Data(), iEta, error.Data()));
          zeeJetHists[iAlgo][iPer][iEta][0][iErr]->Add(temp);
          temp = (TH2D*) thisFile->Get (Form ("zeeJetPtRatio_dataSet%i_%s_iEta%i_data_%s", runNumber, algo.Data(), act_iEta, error.Data()));
          zeeJetHists[iAlgo][2][iEta][0][iErr]->Add(temp);

          temp = (TH2D*) thisFile->Get (Form ("zmumuJetPtRatio_dataSet%i_%s_iEta%i_data_%s", runNumber, algo.Data(), iEta, error.Data()));
          zmumuJetHists[iAlgo][iPer][iEta][0][iErr]->Add(temp);
          temp = (TH2D*) thisFile->Get (Form ("zmumuJetPtRatio_dataSet%i_%s_iEta%i_data_%s", runNumber, algo.Data(), act_iEta, error.Data()));
          zmumuJetHists[iAlgo][2][iEta][0][iErr]->Add(temp);

          temp = (TH2D*) thisFile->Get (Form ("gJetPtRatio_dataSet%i_%s_iEta%i_data_%s", runNumber, algo.Data(), iEta, error.Data()));
          gJetHists[iAlgo][iPer][iEta][0][iErr]->Add(temp);
          temp = (TH2D*) thisFile->Get (Form ("gJetPtRatio_dataSet%i_%s_iEta%i_data_%s", runNumber, algo.Data(), act_iEta, error.Data()));
          gJetHists[iAlgo][2][iEta][0][iErr]->Add(temp);

          //if (iErr == 1)  {
          // gJetHistsSys[iAlgo][iPer][iEta][0][iErr]->Add((TH2D*) thisFile->Get (Form ("gJetPtRatioSys_dataSet%i_iEta%i_data_%s", runNumber, iEta, error.Data())));
          ////gJetHistsSys[2][iEta][0][iErr]->Add((TH2D*) thisFile->Get (Form ("gJetPtRatioSys_dataSet%i_iEta%i_data_%s", runNumber, act_iEta, error.Data())));
          //}
         }
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
       nGammaJetVec = (TVectorD*) thisFile->Get (Form ("nGammaJetVec_%s", gammaJetSampleId.Data()));

       jetCleaningVec = (TVectorD*) thisFile->Get (Form ("jetCleaningVec_%s", gammaJetSampleId.Data()));
       if (jetCleaningVec) {
        nTotalJets[1] += (int) (*jetCleaningVec)[0];
        nCleanJets[1] += (int) (*jetCleaningVec)[1];
       }

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         //const bool flipEta = gammaJetSampleId.Contains("pPb") && iEta < numetabins;
         const bool flipEta = false;
         const short act_iEta = (flipEta ? (numetabins - iEta - 1) : iEta); // period A condition

         nGammaJet[iAlgo][iPer][1][iEta] += (*nGammaJetVec)[iEta + iAlgo*(numetabins+1)];
         nGammaJet[iAlgo][2][1][iEta] += (*nGammaJetVec)[act_iEta + iAlgo*(numetabins+1)];

         // Only add the statistical error plots for MC (don't need to consider systematics)
         TH2D* temp = (TH2D*) thisFile->Get (Form ("gJetPtRatio_dataSet%s_%s_iEta%i_mc_stat", gammaJetSampleId.Data(), algo.Data(), iEta));
         gJetHists[iAlgo][iPer][iEta][1][1]->Add(temp);
         temp = (TH2D*) thisFile->Get (Form ("gJetPtRatio_dataSet%s_%s_iEta%i_mc_stat", gammaJetSampleId.Data(), algo.Data(), act_iEta));
         gJetHists[iAlgo][2][iEta][1][1]->Add(temp);

         //gJetHistsSys[iAlgo][iPer][iEta][1][1]->Add((TH2D*) thisFile->Get (Form ("gJetPtRatioSys_dataSet%s_iEta%i_mc_stat", gammaJetSampleId.Data(), iEta)));
         ////gJetHistsSys[iAlgo][2][iEta][1][1]->Add((TH2D*) thisFile->Get (Form ("gJetPtRatioSys_dataSet%s_iEta%i_mc_stat", gammaJetSampleId.Data(), act_iEta)));
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
       const short iPer = (zeeJetSampleId.Contains("pPb") ? 0 : 1);
       nZeeJetVec = (TVectorD*) thisFile->Get (Form ("nZeeJetVec_%s", zeeJetSampleId.Data()));

       jetCleaningVec = (TVectorD*) thisFile->Get (Form ("jetCleaningVec_%s", zeeJetSampleId.Data()));
       if (jetCleaningVec) {
        nTotalJets[1] += (int) (*jetCleaningVec)[0];
        nCleanJets[1] += (int) (*jetCleaningVec)[1];
       }

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         //const bool flipEta = zeeJetSampleId.Contains("pPb") && iEta < numetabins;
         const bool flipEta = false;
         const short act_iEta = (flipEta ? numetabins - iEta - 1 : iEta); // period A condition

         nZeeJet[iAlgo][iPer][1][iEta] += (*nZeeJetVec)[iEta + iAlgo*(numetabins+1)];
         nZeeJet[iAlgo][2][1][iEta] += (*nZeeJetVec)[act_iEta + iAlgo*(numetabins+1)];

         // Only add the statistical error plots for MC (don't need to consider systematics)
         TH2D* temp = (TH2D*) thisFile->Get (Form ("zeeJetPtRatio_dataSet%s_%s_iEta%i_mc_stat", zeeJetSampleId.Data(), algo.Data(), iEta));
         zeeJetHists[iAlgo][iPer][iEta][1][1]->Add(temp);
         temp = (TH2D*) thisFile->Get (Form ("zeeJetPtRatio_dataSet%s_%s_iEta%i_mc_stat", zeeJetSampleId.Data(), algo.Data(), act_iEta));
         zeeJetHists[iAlgo][2][iEta][1][1]->Add(temp);
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
       const short iPer = (zmumuJetSampleId.Contains("pPb") ? 0 : 1);
       nZmumuJetVec = (TVectorD*) thisFile->Get (Form ("nZmumuJetVec_%s", zmumuJetSampleId.Data()));

       jetCleaningVec = (TVectorD*) thisFile->Get (Form ("jetCleaningVec_%s", zmumuJetSampleId.Data()));
       if (jetCleaningVec) {
        nTotalJets[1] += (int) (*jetCleaningVec)[0];
        nCleanJets[1] += (int) (*jetCleaningVec)[1];
       }

       for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
        const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         //const bool flipEta = zmumuJetSampleId.Contains("pPb") && iEta < numetabins;
         const bool flipEta = false;
         const short act_iEta = (flipEta ? numetabins - iEta - 1 : iEta); // period A condition

         nZmumuJet[iAlgo][iPer][1][iEta] += (*nZmumuJetVec)[iEta + iAlgo*(numetabins+1)];
         nZmumuJet[iAlgo][2][1][iEta] += (*nZmumuJetVec)[act_iEta + iAlgo*(numetabins+1)];

         // Only add the statistical error plots for MC (don't need to
         // consider systematics)
         TH2D* temp = (TH2D*) thisFile->Get (Form ("zmumuJetPtRatio_dataSet%s_%s_iEta%i_mc_stat", zmumuJetSampleId.Data(), algo.Data(), iEta));
         zmumuJetHists[iAlgo][iPer][iEta][1][1]->Add(temp);
         temp = (TH2D*) thisFile->Get (Form ("zmumuJetPtRatio_dataSet%s_%s_iEta%i_mc_stat", zmumuJetSampleId.Data(), algo.Data(), act_iEta));
         zmumuJetHists[iAlgo][2][iEta][1][1]->Add(temp);
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

  for (short iPer = 0; iPer < 3; iPer++) { // loop over period configurations
   TString period = "Period A";
   if (iPer == 1) period = "Period B";
   else if (iPer == 2) period = "Period A+B";

   TString perType = "pA";
   if (iPer == 1) perType = "pB";
   else if (iPer == 2) perType = "pAB";

   for (short iEta = 0; iEta <= numetabins; iEta++) { // loop over bins in eta

    /**** Plot ZmumuJet info ****/
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) { // loop over jet algorithms
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     const Style_t markerStyle = (iAlgo == 0 ? 20 : 24);
     topPad->cd();
     topPad->SetLogx();
     vJetHist = GetProfileX(Form ("zmumuJetHist_%s_%s_iEta%i", algo.Data(), perType.Data(), iEta), zmumuJetHists[iAlgo][iPer][iEta][0][1], numpzbins, pzbins, false);
     vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
     vJetHist->SetYTitle("<#it{x}_{J}^{ref}>");
     vJetHist->SetAxisRange(0.75, 2.15, "Y");
     vJetHist->SetMarkerStyle(markerStyle);
     vJetHist->SetMarkerColor(data_color);
     vJetHist->SetLineColor(data_color);
     vJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetTitleOffset(uPadY);

     // Now calculate systematics by taking the TProfile, then set as the errors to the TGraphAsymmErrors object
     vJetHist_lo = GetProfileX("vJetHist_lo", zmumuJetHists[iAlgo][iPer][iEta][0][0], numpzbins, pzbins, false);
     vJetHist_hi = GetProfileX("vJetHist_hi", zmumuJetHists[iAlgo][iPer][iEta][0][2], numpzbins, pzbins, false);
     CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) delete vJetHist_lo;
     if (vJetHist_hi) delete vJetHist_hi;
     vJetGraph_sys->SetFillColor(kBlack);
     vJetGraph_sys->SetFillStyle(3001);

     vJetHist_mc = GetProfileX(Form ("zmumuJetHist_mc_%s_%s_iEta%i", algo.Data(), perType.Data(), iEta), zmumuJetHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, false);
     vJetHist_mc->SetMarkerStyle(markerStyle);
     vJetHist_mc->SetMarkerColor(mc_color);
     vJetHist_mc->SetLineColor(mc_color);

     if (iAlgo == 0) vJetHist->DrawCopy("e1 x0");
     else vJetHist->DrawCopy("same e1 x0");
     vJetHist_mc->DrawCopy("same e1 x0");
     vJetGraph_sys->Draw("2");

     if (iAlgo == 0) {
      myMarkerText(0.175, 0.88, data_color, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nZmumuJet[iAlgo][iPer][0][iEta]), 1.25, 0.04/uPadY);
      myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with Overlay (%i events)", nZmumuJet[iAlgo][iPer][1][iEta]), 1.25, 0.04/uPadY);
      if (iEta < numetabins) {
       if (iPer == 2) myText(0.155, 0.1,kBlack, Form ("%g < #eta_{Lab}^{#mu#mu} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
       else myText(0.155, 0.1,kBlack, Form ("%g < #eta_{Lab}^{#mu#mu} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
      }
      myText(0.155, 0.28, kBlack, "Z (#mu#mu) + Jet", 0.04/uPadY);
      myText(0.155, 0.19, kBlack, period.Data(), 0.04/uPadY);
     }

     bottomPad->cd();
     bottomPad->SetLogx();

     vJetHist_rat = GetDataOverMC(TString(Form ("zmumuJetPtDataMCRatio_%s_%s_iEta%i", algo.Data(), perType.Data(), iEta)), zmumuJetHists[iAlgo][iPer][iEta][0][1], zmumuJetHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
     vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
     vJetHist_rat_lo = GetDataOverMC(TString(Form ("zmumuJetPtDataMCRatio_lo_%s_iEta%i", algo.Data(), iEta)), zmumuJetHists[iAlgo][iPer][iEta][0][0], zmumuJetHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
     vJetHist_rat_hi = GetDataOverMC(TString(Form ("zmumuJetPtDataMCRatio_hi_%s_iEta%i", algo.Data(), iEta)), zmumuJetHists[iAlgo][iPer][iEta][0][2], zmumuJetHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
     CalcSystematics(vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) delete vJetHist_rat_lo;
     if (vJetHist_rat_hi) delete vJetHist_rat_hi;
     vJetGraph_rat_sys->SetFillColor(kBlack);
     vJetGraph_rat_sys->SetFillStyle(3001);

     vJetHist_rat->SetYTitle("Data / MC");
     vJetHist_rat->SetAxisRange(0.85, 1.15, "Y");
     vJetHist_rat->SetMarkerStyle(markerStyle);
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

     if (iAlgo == 0) vJetHist_rat->Draw("e1 x0");
     else vJetHist_rat->Draw("same e1 x0");
     vJetGraph_rat_sys->Draw("2");
     for (TLine* line : zlines) line->Draw("same");
    }

    char* plotName;
    if (iEta < numetabins) plotName = Form ("z_mumu_jet%i.pdf", iEta);
    else plotName = Form ("z_mumu_jet_combined.pdf");

    switch (iPer) {
     case 0:
      canvas->SaveAs(Form ("%s/PeriodA/%s", plotPath.Data(), plotName));
      break;
     case 1:
      canvas->SaveAs(Form ("%s/PeriodB/%s", plotPath.Data(), plotName));
      break;
     case 2:
      canvas->SaveAs(Form ("%s/PeriodAB/%s", plotPath.Data(), plotName));
      break;
    }
    if (vJetHist) delete vJetHist;
    if (vJetHist_mc) delete vJetHist_mc;
    if (vJetGraph_sys) delete vJetGraph_sys;
    if (vJetHist_rat) delete vJetHist_rat;
    if (vJetGraph_rat_sys) delete vJetGraph_rat_sys;


    /**** Plots ZeeJet info ****/
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     const Style_t markerStyle = (iAlgo == 0 ? 20 : 24);
     topPad->cd();
     topPad->SetLogx();
     vJetHist = GetProfileX(Form ("zeeJetHist_%s_%s_iEta%i", algo.Data(), perType.Data(), iEta), zeeJetHists[iAlgo][iPer][iEta][0][1], numpzbins, pzbins, false);
     vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
     vJetHist->SetYTitle("<#it{x}_{J}^{ref}>");
     vJetHist->SetAxisRange(0.75, 2.15, "Y");
     vJetHist->SetMarkerStyle(markerStyle);
     vJetHist->SetMarkerColor(data_color);
     vJetHist->SetLineColor(data_color);
     vJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetTitleOffset(uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     vJetHist_lo = GetProfileX("vJetHist_lo", zeeJetHists[iAlgo][iPer][iEta][0][0], numpzbins, pzbins, false);
     vJetHist_hi = GetProfileX("vJetHist_hi", zeeJetHists[iAlgo][iPer][iEta][0][2], numpzbins, pzbins, false);
     CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) delete vJetHist_lo;
     if (vJetHist_hi) delete vJetHist_hi;
     vJetGraph_sys->SetFillColor(kBlack);
     vJetGraph_sys->SetFillStyle(3001);

     vJetHist_mc = GetProfileX(Form ("zeeJetHist_mc_%s_%s_iEta%i", algo.Data(), perType.Data(), iEta), zeeJetHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, false);
     vJetHist_mc->SetMarkerStyle(markerStyle);
     vJetHist_mc->SetMarkerColor(mc_color);
     vJetHist_mc->SetLineColor(mc_color);

     if (iAlgo == 0) vJetHist->DrawCopy("e1 x0");
     else vJetHist->DrawCopy("same e1 x0");
     vJetHist_mc->DrawCopy("same e1 x0");
     vJetGraph_sys->Draw("2");

     if (iAlgo == 0) {
      myMarkerText(0.175, 0.88, data_color, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nZeeJet[iAlgo][iPer][0][iEta]), 1.25, 0.04/uPadY);
      myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with Overlay (%i events)", nZeeJet[iAlgo][iPer][1][iEta]), 1.25, 0.04/uPadY);
      if (iEta < numetabins) {
       if (iPer == 2) myText(0.155, 0.1,kBlack, Form ("%g < #eta_{Lab}^{ee} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
       else myText(0.155, 0.1,kBlack, Form ("%g < #eta_{Lab}^{ee} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
      }
      myText(0.155, 0.28, kBlack, "Z (ee) + Jet", 0.04/uPadY);
      myText(0.155, 0.19, kBlack, period.Data(), 0.04/uPadY);
     }

     bottomPad->cd();
     bottomPad->SetLogx();

     vJetHist_rat = GetDataOverMC(TString(Form ("zeeJetPtDataMCRatio_%s_%s_iEta%i", algo.Data(), perType.Data(), iEta)), zeeJetHists[iAlgo][iPer][iEta][0][1], zeeJetHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
     vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
     vJetHist_rat_lo = GetDataOverMC(TString(Form ("zeeJetPtDataMCRatio_lo_%s_iEta%i", algo.Data(), iEta)), zeeJetHists[iAlgo][iPer][iEta][0][0], zeeJetHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
     vJetHist_rat_hi = GetDataOverMC(TString(Form ("zeeJetPtDataMCRatio_hi_%s_iEta%i", algo.Data(), iEta)), zeeJetHists[iAlgo][iPer][iEta][0][2], zeeJetHists[iAlgo][iPer][iEta][1][1], numpzbins, pzbins, numxjrefbins, xjrefbins, false);
     CalcSystematics(vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) delete vJetHist_rat_lo;
     if (vJetHist_rat_hi) delete vJetHist_rat_hi;
     vJetGraph_rat_sys->SetFillColor(kBlack);
     vJetGraph_rat_sys->SetFillStyle(3001);

     vJetHist_rat->SetYTitle("Data / MC");
     vJetHist_rat->SetAxisRange(0.85, 1.15, "Y");
     vJetHist_rat->SetMarkerStyle(markerStyle);
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

     if (iAlgo == 0) vJetHist_rat->Draw("e1 x0");
     else vJetHist_rat->Draw("same e1 x0");
     vJetGraph_rat_sys->Draw("2");
     for (TLine* line : zlines) line->Draw("same");
    }

    if (iEta < numetabins) plotName = Form ("z_ee_jet%i.pdf", iEta);
    else plotName = Form ("z_ee_jet_combined.pdf");
    switch (iPer) {
     case 0:
      canvas->SaveAs(Form ("%s/PeriodA/%s", plotPath.Data(), plotName));
      break;
     case 1:
      canvas->SaveAs(Form ("%s/PeriodB/%s", plotPath.Data(), plotName));
      break;
     case 2:
      canvas->SaveAs(Form ("%s/PeriodAB/%s", plotPath.Data(), plotName));
      break;
    }
    if (vJetHist) delete vJetHist;
    if (vJetHist_mc) delete vJetHist_mc;
    if (vJetGraph_sys) delete vJetGraph_sys;
    if (vJetHist_rat) delete vJetHist_rat;
    if (vJetGraph_rat_sys) delete vJetGraph_rat_sys;


    /**** Plots GammaJet info as a function of p_T^ref****/
    for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     const Style_t markerStyle = (iAlgo == 0 ? 20 : 24);
     topPad->cd();
     topPad->SetLogx();
     vJetHist = GetProfileX(Form ("gammaJetHist_%s_%s_iEta%i", algo.Data(), perType.Data(), iEta), gJetHists[iAlgo][iPer][iEta][0][1], numpbins, pbins, true);
     vJetGraph_sys = new TGraphAsymmErrors(vJetHist); // for plotting systematics
     vJetHist->SetYTitle("<#it{x}_{J}^{ref}>");
     vJetHist->SetAxisRange(0.75, 2.15, "Y");
     vJetHist->SetMarkerStyle(markerStyle);
     vJetHist->SetMarkerColor(data_color);
     vJetHist->SetLineColor(data_color);
     vJetHist->GetXaxis()->SetLabelSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetLabelSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetTitleSize(0.04/uPadY);
     vJetHist->GetYaxis()->SetTitleOffset(uPadY);

     // Now calculate systematics by taking the TProfile of the pt+err and pt-err samples, then set as the errors to the TGraphAsymmErrors object
     vJetHist_lo = GetProfileX("vJetHist_lo", gJetHists[iAlgo][iPer][iEta][0][0], numpbins, pbins, true);
     vJetHist_hi = GetProfileX("vJetHist_hi", gJetHists[iAlgo][iPer][iEta][0][2], numpbins, pbins, true);
     CalcSystematics(vJetGraph_sys, vJetHist, vJetHist_hi, vJetHist_lo);
     if (vJetHist_lo) delete vJetHist_lo;
     if (vJetHist_hi) delete vJetHist_hi;
     vJetGraph_sys->SetFillColor(kBlack);
     vJetGraph_sys->SetFillStyle(3001);

     vJetHist_mc = GetProfileX(Form ("gammaJetHist_mc_%s_%s_iEta%i", algo.Data(), perType.Data(), iEta), gJetHists[iAlgo][iPer][iEta][1][1], numpbins, pbins, true);
     vJetHist_mc->SetMarkerStyle(markerStyle);
     vJetHist_mc->SetMarkerColor(mc_color);
     vJetHist_mc->SetLineColor(mc_color);

     if (iAlgo == 0) vJetHist->DrawCopy("e1 x0");
     else vJetHist->DrawCopy("same e1 x0");
     vJetHist_mc->DrawCopy("same e1 x0");
     vJetGraph_sys->Draw("2");
     for (TLine* line : dplines) line->Draw("same");

     if (iAlgo == 0) {
      myMarkerText(0.175, 0.88, data_color, kFullCircle, Form ("2016 #it{p}+Pb 8.16 TeV, with Insitu Corrections (%i events)", nGammaJet[iAlgo][iPer][0][iEta]), 1.25, 0.04/uPadY);
      myMarkerText(0.175, 0.81, mc_color, kFullCircle, Form ("Pythia8 #it{pp} 8.16 TeV with Overlay (%i events)", nGammaJet[iAlgo][iPer][1][iEta]), 1.25, 0.04/uPadY);
      if (iEta < numetabins) {
       if (iPer == 2) myText(0.155, 0.1,kBlack, Form ("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
       else myText(0.155, 0.1,kBlack, Form ("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);
      }
      myText(0.155, 0.28, kBlack, "#gamma + Jet", 0.04/uPadY);
      myText(0.155, 0.19, kBlack, period.Data(), 0.04/uPadY);
     }

     bottomPad->cd();
     bottomPad->SetLogx();
     vJetHist_rat = GetDataOverMC(TString(Form ("gJetPtDataMCRatio_%s_%s_iEta%i", algo.Data(), perType.Data(), iEta)), gJetHists[iAlgo][iPer][iEta][0][1], gJetHists[iAlgo][iPer][iEta][1][1], numpbins, pbins, numxjrefbins, xjrefbins, true);
     vJetGraph_rat_sys = new TGraphAsymmErrors(vJetHist_rat);
     vJetHist_rat_lo = GetDataOverMC(TString(Form ("gJetPtDataMCRatio_lo_%s_iEta%i", algo.Data(), iEta)), gJetHists[iAlgo][iPer][iEta][0][0], gJetHists[iAlgo][iPer][iEta][1][1], numpbins, pbins, numxjrefbins, xjrefbins, true);
     vJetHist_rat_hi = GetDataOverMC(TString(Form ("gJetPtDataMCRatio_hi_iEta%i", iEta)), gJetHists[iAlgo][iPer][iEta][0][2], gJetHists[iAlgo][iPer][iEta][1][1], numpbins, pbins, numxjrefbins, xjrefbins, true);
     CalcSystematics(vJetGraph_rat_sys, vJetHist_rat, vJetHist_rat_hi, vJetHist_rat_lo);
     if (vJetHist_rat_lo) delete vJetHist_rat_lo;
     if (vJetHist_rat_hi) delete vJetHist_rat_hi;
     vJetGraph_rat_sys->SetFillColor(kBlack);
     vJetGraph_rat_sys->SetFillStyle(3001);

     vJetHist_rat->SetYTitle("Data / MC");
     vJetHist_rat->SetAxisRange(0.91, 1.09, "Y");
     vJetHist_rat->SetMarkerStyle(markerStyle);
     vJetHist_rat->GetYaxis()->SetNdivisions(405);
     vJetHist_rat->GetXaxis()->SetTitleSize(0.04/dPadY);
     vJetHist_rat->GetYaxis()->SetTitleSize(0.04/dPadY);
     vJetHist_rat->GetXaxis()->SetTitleOffset(1);
     vJetHist_rat->GetYaxis()->SetTitleOffset(dPadY);
     vJetHist_rat->GetYaxis()->CenterTitle(true);
     vJetHist_rat->GetXaxis()->SetLabelSize(0.04/dPadY);
     vJetHist_rat->GetYaxis()->SetLabelSize(0.04/dPadY);
     vJetHist_rat->GetXaxis()->SetTickLength(0.08);

     if (iAlgo == 0) vJetHist_rat->Draw("e1 x0"); 
     else vJetHist_rat->Draw("same e1 x0");
     vJetGraph_rat_sys->Draw("2");
     for (TLine* line : glines) line->Draw("same");
     for (TLine* line : dplines_bottom) line->Draw("same");
    }

    if (iEta < numetabins) plotName = Form ("gamma_jet%i.pdf", iEta);
    else plotName = Form ("gamma_jet_combined.pdf");
    switch (iPer) {
     case 0:
      canvas->SaveAs(Form ("%s/PeriodA/%s", plotPath.Data(), plotName));
      break;
     case 1:
      canvas->SaveAs(Form ("%s/PeriodB/%s", plotPath.Data(), plotName));
      break;
     case 2:
      canvas->SaveAs(Form ("%s/PeriodAB/%s", plotPath.Data(), plotName));
      break;
    }
    if (vJetHist) delete vJetHist;
    if (vJetHist_mc) delete vJetHist_mc;
    if (vJetGraph_sys) delete vJetGraph_sys;
    if (vJetHist_rat) delete vJetHist_rat;
    if (vJetGraph_rat_sys) delete vJetGraph_rat_sys;
   }
  }


  /**** Plots systematic errors vs jet pt ****/
  //for (short iPer = 0; iPer < 2; iPer++) {
  // const TString period = (iPer == 0 ? "Period A":"Period B");

  // for (short iEta = 0; iEta < numetabins; iEta++) {
  //  topPad->cd();
  //  TH2D* thisHist = gJetHistsSys[iPer][iEta][0][1];
  //  TH1D* rmsHist = new TH1D(Form ("rms_iEta%i_%s", iEta, (iPer==0?"pPb":"Pbp")), "", numpzbins, pzbins);
  //  for (short pzbin = 0; pzbin < numpzbins; pzbin++) {
  //   float rms = 0;
  //   float sumWeights = 0;
  //   for (short sigbin = 0; sigbin < numSigmaBins; sigbin++) {
  //    const float sig = thisHist->GetYaxis()->GetBinCenter(sigbin+1);
  //    const float weight = thisHist->GetBinContent(pzbin+1, sigbin+1);
  //    rms += pow(sig, 2) * weight;
  //    sumWeights += weight;
  //   }
  //   if (sumWeights > 0) rms = sqrt(rms) / sqrt(sumWeights);
  //   rmsHist->SetBinContent(pzbin+1, rms);
  //  }
  //  topPad->SetLogz();
  //  thisHist->Draw("col");
  //  thisHist->GetXaxis()->SetLabelSize(0.04/uPadY);
  //  thisHist->GetYaxis()->SetLabelSize(0.04/uPadY);
  //  thisHist->GetYaxis()->SetTitleSize(0.04/uPadY);
  //  thisHist->GetYaxis()->SetTitleOffset(1.1*uPadY);
  //  myText(0.72, 0.89, kBlack, period.Data(), 0.04/uPadY);
  //  myText(0.72, 0.8,kBlack, Form ("%g < #eta_{Lab}^{#gamma} < %g", etabins[iEta], etabins[iEta+1]), 0.04/uPadY);

  //  bottomPad->cd();
  //  rmsHist->SetXTitle("#it{p}_{T}^{jet} #left[GeV#right]");
  //  //rmsHist->SetYTitle("RMS(#sigma)/#it{p}_{T}^{jet}");
  //  rmsHist->SetYTitle("RMS");
  //  rmsHist->SetAxisRange(0, 0.09, "Y");
  //  rmsHist->GetYaxis()->SetNdivisions(405);
  //  rmsHist->GetXaxis()->SetTitleSize(0.04/dPadY);
  //  rmsHist->GetYaxis()->SetTitleSize(0.04/dPadY);
  //  rmsHist->GetXaxis()->SetTitleOffset(1);
  //  rmsHist->GetYaxis()->SetTitleOffset(1.1*dPadY);
  //  rmsHist->GetYaxis()->CenterTitle(true);
  //  rmsHist->GetXaxis()->SetLabelSize(0.04/dPadY);
  //  rmsHist->GetYaxis()->SetLabelSize(0.04/dPadY);
  //  rmsHist->GetXaxis()->SetTickLength(0.08);
  //  rmsHist->Draw("hist");
  //  canvas->SaveAs(Form ("%s/Period%s/jetSystematics_iEta%i.pdf", plotPath.Data(), (iPer==0 ? "A":"B"), iEta));
  //  if (rmsHist) delete rmsHist;
  // }
  //}

  Delete4DArray (nZeeJet, 2, 3, 2, numetabins+1);
  Delete4DArray (nZmumuJet, 2, 3, 2, numetabins+1);
  Delete4DArray (nGammaJet, 2, 3, 2, numetabins+1);
  Delete5DArray (zeeJetHists, 2, 3, numetabins+1, 2, 3);
  Delete5DArray (zmumuJetHists, 2, 3, numetabins+1, 2, 3);
  Delete5DArray (gJetHists, 2, 3, numetabins+1, 2, 3);

  double num = (double)nCleanJets[0];
  double den = (double)nTotalJets[0];
  double frac = num / den;
  double fracErr = sqrt ( ((num+1.)*(num+2.)) / ((den+2.)*(den+3.)) - ((num+1.)*(num+1.))/((den+2.)*(den+2.)) );

  cout << "Clean / total jets in data = " << nCleanJets[0]
       << " / " << nTotalJets[0]
       << " = " << frac
       << " +/- " << fracErr
       << endl;

  num = (double)nCleanJets[1];
  den = (double)nTotalJets[1];
  frac = num / den;
  fracErr = sqrt ( ((num+1.)*(num+2.)) / ((den+2.)*(den+3.)) - ((num+1.)*(num+1.))/((den+2.)*(den+2.)) );

  cout << "Clean / total jets in MC   = " << nCleanJets[1]
       << " / " << nTotalJets[1]
       << " = " << frac
       << " +/- " << fracErr
       << endl;

  return;
}

} // end namespace
