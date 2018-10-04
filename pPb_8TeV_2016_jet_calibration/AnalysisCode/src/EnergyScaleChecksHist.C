#include "EnergyScaleChecksHist.h"
#include "Params.h"
#include "Utils.h"

#include <ArrayTemplates.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
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

const short plotMC = 1; // 0 plots overlay


TString GetMCType (const short iMC) {
  if (iMC == 0) return "overlay";
  else if (iMC == 1) return "signal";
  else if (iMC == 2) return "valid_signal";
  else return "";
}


void EnergyScaleChecksHist () {

  SetAtlasStyle ();

  // Setup trigger vectors
  SetupDirectories ("EnergyScaleChecks/", "pPb_8TeV_2016_jet_calibration/");

  // Setup list of data and lists of MC samples
  vector<int> runNumbers (0);
  for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);
  vector<TString> gammaJetOverlaySampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetOverlaySampleIds.push_back (TString ("Pbp_Overlay_GammaJet_Slice") + to_string (i+1));
   gammaJetOverlaySampleIds.push_back (TString ("pPb_Overlay_GammaJet_Slice") + to_string (i+1));
  }
  vector<TString> gammaJetSignalSampleIds (0);
  for (short i = 0; i < 6; i++) {
   gammaJetSignalSampleIds.push_back (TString ("pPb_Signal_GammaJet_Slice") + to_string (i+1));
  }
  vector<TString> zeeJetSampleIds (0);
  zeeJetSampleIds.push_back ("Pbp_Overlay_ZeeJet");
  zeeJetSampleIds.push_back ("pPb_Overlay_ZeeJet");

  vector<TString> zmumuJetSampleIds (0);
  zmumuJetSampleIds.push_back ("Pbp_Overlay_ZmumuJet");
  zmumuJetSampleIds.push_back ("pPb_Overlay_ZmumuJet");

  vector<TString> dijetSampleIds (0);
  dijetSampleIds.push_back ("pPb_Signal_Dijet_Slice2");

  vector<TString> dijetValidSampleIds (0);
  dijetValidSampleIds.push_back ("pPb_Signal_Valid_Dijet_Slice2");
  dijetValidSampleIds.push_back ("Pbp_Signal_Valid_Dijet_Slice2");


  /**** Initialize histograms ****/
  TH1D* electronEnergyScale = new TH1D ("electronEnergyScale", ";Electron #it{p}_{T}^{reco} / #it{p}_{T}^{truth};", 200, 0, 2.0);
  electronEnergyScale->Sumw2 ();

  TH1D****** jetEnergyResponseCalib = Get5DArray <TH1D*> (numdpbins, 2, 3, numpbins+1, numetabins+1); // DP slice, jet algo, MC type, pt bin, eta bin
  TH1D****** jetEnergyResponseReco = Get5DArray <TH1D*> (numdpbins, 2, 3, numpbins+1, numetabins+1);
  TH1D***** photonEnergyResponse = Get4DArray <TH1D*> (numdpbins, 3, numpbins+1, numetabins+1);

  for (short iDP = 0; iDP < numdpbins; iDP++) {
   for (short iP = 0; iP <= numpbins; iP++) {
    for (short iEta = 0; iEta <= numetabins; iEta++) {
     for (short iMC = 0; iMC < 3; iMC++) {
      const TString mcType = GetMCType (iMC);
      for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
       const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
       jetEnergyResponseCalib[iDP][iAlgo][iMC][iP][iEta] = new TH1D (Form ("%s_%s_jetEnergyResponseCalib_iP%i_iEta%i_Slice%i", algo.Data (), mcType.Data (), iP, iEta, iDP+1), "", 100, 0, 2);
       jetEnergyResponseCalib[iDP][iAlgo][iMC][iP][iEta]->Sumw2 ();
       jetEnergyResponseReco[iDP][iAlgo][iMC][iP][iEta] = new TH1D (Form ("%s_%s_jetEnergyResponseReco_iP%i_iEta%i_Slice%i", algo.Data (), mcType.Data (), iP, iEta, iDP+1), "", 100, 0, 2);
       jetEnergyResponseReco[iDP][iAlgo][iMC][iP][iEta]->Sumw2 ();
      }
      photonEnergyResponse[iDP][iMC][iP][iEta] = new TH1D (Form ("%s_photonEnergyResponse_iP%i_iEta%i_Slice%i", mcType.Data (), iP, iEta, iDP+1), ";#it{p}_{T}^{reco} / #it{p}_{T}^{truth};Counts", 100, 0, 2);
      photonEnergyResponse[iDP][iMC][iP][iEta]->Sumw2 ();
     }
    }
   }
  }

  int***** nJet = Get5DArray <int> (numdpbins, 2, 3, numpbins+1, numetabins+1);
  int**** nGamma = Get4DArray <int> (numdpbins, 3, numpbins+1, numetabins+1);

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
   TVectorD *nJetVec, *nGammaVec;
   int numFiles = 0;
   while ( (sysfile= (TSystemFile*)next ())) {
    fname = sysfile->GetName ();
    if (!sysfile->IsDirectory () && fname.EndsWith (".root")) {
     if (debugStatements) cout << "Status: In triggers_pt_counts.C (breakpoint I): Found " << fname.Data () << endl;

     // do this if gamma jet MC sample (OVERLAY)
     if (plotMC == 0) {
      for (TString gammaJetOverlaySampleId : gammaJetOverlaySampleIds) { // check for gamma jet MC
       if (fname.Contains (gammaJetOverlaySampleId)) { // if gamma jet MC sample
        numFiles++;
        cout << "Reading in " << rootPath+fname << endl;
        TFile* thisFile = new TFile (rootPath + fname, "READ");
        nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", gammaJetOverlaySampleId.Data ()));
        nGammaVec = (TVectorD*)thisFile->Get (Form ("nGammaVec_%s", gammaJetOverlaySampleId.Data ()));

        int iDP = -1;
        while (iDP < numdpbins && !gammaJetOverlaySampleId.Contains (Form ("Slice%i", iDP+1))) iDP++;
        if (iDP == -1 || 6 <= iDP) continue;

        //const int iDP = 0;

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         for (short iP = 0; iP <= numpbins; iP++) {
          for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
            const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
            nJet[iDP][iAlgo][0][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo* (numpbins+1)* (numetabins+1)];
            jetEnergyResponseCalib[iDP][iAlgo][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data (), gammaJetOverlaySampleId.Data (), iP, iEta)));
            jetEnergyResponseReco[iDP][iAlgo][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data (), gammaJetOverlaySampleId.Data (), iP, iEta)));
          }
          nGamma[iDP][0][iP][iEta] += (*nGammaVec)[iP + (numpbins+1)*iEta];
          photonEnergyResponse[iDP][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("photonEnergyResponse_dataSet%s_iP%i_iEta%i", gammaJetOverlaySampleId.Data (), iP, iEta)));
         }
        }

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
        nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", gammaJetSignalSampleId.Data ()));
        nGammaVec = (TVectorD*)thisFile->Get (Form ("nGammaVec_%s", gammaJetSignalSampleId.Data ()));

        int iDP = -1;
        while (iDP < numdpbins && !gammaJetSignalSampleId.Contains (Form ("Slice%i", iDP+1))) iDP++;
        if (iDP == -1 || 6 <= iDP) continue;

        //const int iDP = 0;

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         for (short iP = 0; iP <= numpbins; iP++) {
          for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
            const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
            nJet[iDP][iAlgo][1][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo* (numpbins+1)* (numetabins+1)];
            jetEnergyResponseCalib[iDP][iAlgo][1][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data (), gammaJetSignalSampleId.Data (), iP, iEta)));
            jetEnergyResponseReco[iDP][iAlgo][1][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data (), gammaJetSignalSampleId.Data (), iP, iEta)));
          }
          nGamma[iDP][1][iP][iEta] += (*nGammaVec)[iP + (numpbins+1)*iEta];
          photonEnergyResponse[iDP][1][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("photonEnergyResponse_dataSet%s_iP%i_iEta%i", gammaJetSignalSampleId.Data (), iP, iEta)));
         }
        }

        thisFile->Close ();
        delete thisFile;
        break;
       }
      }
     }
     // do this if Z->ee MC sample
     for (TString zeeJetSampleId : zeeJetSampleIds) { // check for Z->ee MC
      if (fname.Contains (zeeJetSampleId)) { // if Z->ee MC do this
       numFiles++;
       cout << "Reading in " << rootPath+fname << endl;
       TFile* thisFile = new TFile (rootPath + fname, "READ");
       nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", zeeJetSampleId.Data ()));

       electronEnergyScale->Add ( (TH1D*)thisFile->Get (Form ("electronEnergyScale_dataSet%s", zeeJetSampleId.Data ())));

       //for (short iEta = 0; iEta <= numetabins; iEta++) {
       // for (short iP = 0; iP <= numpbins; iP++) {
       //  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
       //   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
       //   nJet[iAlgo][0][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo* (numpbins+1)* (numetabins+1)];
       //   jetEnergyResponseCalib[iAlgo][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data (), zeeJetSampleId.Data (), iP, iEta)));
       //   jetEnergyResponseReco[iAlgo][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data (), zeeJetSampleId.Data (), iP, iEta)));
       //  }
       // }
       //}

       thisFile->Close ();
       delete thisFile;
       break;
      }
     }
     //// do this if Z->mumu sample
     //for (TString zmumuJetSampleId : zmumuJetSampleIds) { // check for Z->mumu MC
     // if (fname.Contains (zmumuJetSampleId)) { // if Z->mumu sample do this
     //  numFiles++;
     //  cout << "Reading in " << rootPath+fname << endl;
     //  TFile* thisFile = new TFile (rootPath + fname, "READ");
     //  nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", zmumuJetSampleId.Data ()));

     //  //for (short iEta = 0; iEta <= numetabins; iEta++) {
     //  // for (short iP = 0; iP <= numpbins; iP++) {
     //  //  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
     //  //   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
     //  //   nJet[iAlgo][0][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo* (numpbins+1)* (numetabins+1)];
     //  //   jetEnergyResponseCalib[iAlgo][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data (), zmumuJetSampleId.Data (), iP, iEta)));
     //  //   jetEnergyResponseReco[iAlgo][0][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data (), zmumuJetSampleId.Data (), iP, iEta)));
     //  //  }
     //  // }
     //  //}

     //  thisFile->Close ();
     //  delete thisFile;
     //  break;
     // }
     //}

     // do this if dijet MC sample (SIGNAL)
     if (plotMC == 1) {
      for (TString dijetSampleId : dijetSampleIds) { // check for gamma jet MC
       if (fname.Contains (dijetSampleId)) { // if gamma jet MC sample
        numFiles++;
        cout << "Reading in " << rootPath+fname << endl;
        TFile* thisFile = new TFile (rootPath + fname, "READ");
        nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", dijetSampleId.Data ()));
        nGammaVec = (TVectorD*)thisFile->Get (Form ("nGammaVec_%s", dijetSampleId.Data ()));

        const short iDP = 0;

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         for (short iP = 0; iP <= numpbins; iP++) {
          for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
            const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
            nJet[iDP][iAlgo][1][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo* (numpbins+1)* (numetabins+1)];
            jetEnergyResponseCalib[iDP][iAlgo][1][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data (), dijetSampleId.Data (), iP, iEta)));
            jetEnergyResponseReco[iDP][iAlgo][1][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data (), dijetSampleId.Data (), iP, iEta)));
          }
          //nGamma[iDP][1][iP][iEta] += (*nGammaVec)[iP + (numpbins+1)*iEta];
          //photonEnergyResponse[iDP][1][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("photonEnergyResponse_dataSet%s_iP%i_iEta%i", dijetSampleId.Data (), iP, iEta)));
         }
        }

        thisFile->Close ();
        delete thisFile;
        break;
       }
      }

      // do this if dijet MC sample (SIGNAL)
      for (TString dijetValidSampleId : dijetValidSampleIds) { // check for gamma jet MC
       if (fname.Contains (dijetValidSampleId)) { // if gamma jet MC sample
        numFiles++;
        cout << "Reading in " << rootPath+fname << endl;
        TFile* thisFile = new TFile (rootPath + fname, "READ");
        nJetVec = (TVectorD*)thisFile->Get (Form ("nJetVec_%s", dijetValidSampleId.Data ()));
        nGammaVec = (TVectorD*)thisFile->Get (Form ("nGammaVec_%s", dijetValidSampleId.Data ()));

        const short iDP = 0;

        for (short iEta = 0; iEta <= numetabins; iEta++) {
         for (short iP = 0; iP <= numpbins; iP++) {
          for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
            const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");
            nJet[iDP][iAlgo][2][iP][iEta] += (*nJetVec)[iP + (numpbins+1)*iEta + iAlgo* (numpbins+1)* (numetabins+1)];
            jetEnergyResponseCalib[iDP][iAlgo][2][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseCalib_dataSet%s_iP%i_iEta%i", algo.Data (), dijetValidSampleId.Data (), iP, iEta)));
            jetEnergyResponseReco[iDP][iAlgo][2][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("%s_jetEnergyResponseReco_dataSet%s_iP%i_iEta%i", algo.Data (), dijetValidSampleId.Data (), iP, iEta)));
          }
          //nGamma[iDP][2][iP][iEta] += (*nGammaVec)[iP + (numpbins+1)*iEta];
          //photonEnergyResponse[iDP][2][iP][iEta]->Add ( (TH1D*)thisFile->Get (Form ("photonEnergyResponse_dataSet%s_iP%i_iEta%i", dijetValidSampleId.Data (), iP, iEta)));
         }
        }

        thisFile->Close ();
        delete thisFile;
        break;
       }
      }
     }
    }
   }
   cout << numFiles << " files read in." << endl;
  }
  /**** End loop over input files ****/


  /////**** Plots the electron energy scale ****/
  TH1D* thisHist = NULL;
  TGraphAsymmErrors* thisGraph = NULL;
  Color_t colors[15] = {kGray+2, kAzure, kViolet, kMagenta, kPink+10, kPink, kOrange+10, kOrange, kSpring+10, kSpring, kTeal+10, kTeal, kAzure+10, kGray, kBlack};

  TCanvas* energyScaleCanvas = new TCanvas ("energyScaleCanvas", "", 800, 600);
  thisHist = electronEnergyScale;
  TF1* gausFit = new TF1 ("gausFit", "gaus (0)", 0, 2.0);
  energyScaleCanvas->cd ();
  const float ncounts = thisHist->Integral ();
  thisHist->Scale (1./ncounts);
  thisHist->Fit (gausFit, "R");
  const float m = gausFit->GetParameter (1);
  const float s = gausFit->GetParameter (2);
  if (gausFit) delete gausFit;
  gausFit = new TF1 ("gausFit2", "gaus (0)", m - 1.6*s, m + 1.6*s);
  thisHist->Fit (gausFit, "R");

  thisHist->GetYaxis ()->SetTitle ("Counts / Total");
  //thisHist->SetMarkerStyle (6);
  thisHist->Draw ("e1 x0");
  myText (0.18, 0.85, kBlack, "Electron energy response");
  myText (0.18, 0.78, kBlack, Form ("%i electrons", (int)ncounts));
  myText (0.18, 0.71, kBlack, Form ("En. Scale = %.5f #pm %.5f", gausFit->GetParameter (1), gausFit->GetParError (1)));
  myText (0.18, 0.64, kBlack, Form ("En. Res. = %.5f #pm %.5f", gausFit->GetParameter (2), gausFit->GetParError (2)));
  energyScaleCanvas->SaveAs (Form ("%s/electronEnergyScale.pdf", plotPath.Data ()));


  int** jetCounter = Get2DArray <int> (2, 3);
  /**** Plots the jet energy response ****/
  for (short iAlgo = 0; iAlgo < 2; iAlgo++) {
   const TString algo = (iAlgo == 0 ? "akt4hi" : "akt4emtopo");

   TH1D*** jetEnergyScale = Get2DArray <TH1D*> (3, numetabins+1);
   TH1D*** jetEnergyRes = Get2DArray <TH1D*> (3, numetabins+1);

   for (short iMC = plotMC; iMC < plotMC+2; iMC++) {
    const TString mcType = GetMCType (iMC);
    for (int iEta = 0; iEta <= numetabins; iEta++) {
     TCanvas* jesCanvas = new TCanvas (Form ("%s_jetEnergyScalePad_%i", algo.Data (), iEta), "", 800, 600);

     jetEnergyScale[iMC][iEta] = new TH1D (Form ("%s_%s_jetEnergyScale_eta%i", algo.Data (), mcType.Data (), iEta), "", numpbins, pbins);
     jetEnergyRes[iMC][iEta] = new TH1D (Form ("%s_%s_jetEnergyRes_eta%i", algo.Data (), mcType.Data (), iEta), "", numpbins, pbins);

     jesCanvas->cd ();
     jesCanvas->Divide (4, 3);

     TF1** recoFits = Get1DArray<TF1*> (12);
     TF1** calibFits = Get1DArray<TF1*> (12);
     for (short iDP = 0; iDP < numdpbins; iDP++) {
      for (int iP = 0; iP < numpbins; iP++) {

       if (! (dpbins[iDP] <= pbins[iP]) || (numdpbins > 1 && ! (pbins[iP+1] <= dpbins[iDP+1])))
        continue;

       if (iP <= 11)
        jesCanvas->cd (iP+1);

       thisHist = jetEnergyResponseReco[iDP][iAlgo][iMC][iP][iEta];
       if (thisHist->Integral () != 0)
        thisHist->Scale (1./thisHist->Integral ());

       int nJets = nJet[iDP][iAlgo][iMC][iP][iEta];
       jetCounter [iAlgo][iMC] += nJets;

       TF1* recoFit = new TF1 (Form ("recoFit_iP%i", iP), "gaus (0)", 0, 2.0);
       if (nJets > 0) {
        recoFit->SetParameter (1, thisHist->GetMean ());
        recoFit->SetParameter (2, thisHist->GetStdDev ());
        thisHist->Fit (recoFit, "RN");
        double m = recoFit->GetParameter (1);
        double s = recoFit->GetParameter (2);
        if (recoFit) delete recoFit;

        recoFit = new TF1 (Form ("recoFit2_iP%i", iP), "gaus (0)", m - 2.0*s, m + 2.0*s);
        recoFit->SetParameter (1, thisHist->GetMean ());
        recoFit->SetParameter (2, thisHist->GetStdDev ());
        thisHist->Fit (recoFit, "RN");
       }

       if (iMC == plotMC && iP <= 11) {
        thisHist->SetMarkerStyle (kFullDotMedium);
        thisHist->SetLineColor (kBlack);
        thisHist->SetMarkerColor (kBlack);
        thisHist->GetYaxis ()->SetTitle ("Counts / Total");
        thisHist->Draw ("e1 x0");

        recoFit->SetLineColor (kBlack);
        recoFit->Draw ("same");
        recoFits[iP] = recoFit;
       }
       

       thisHist = jetEnergyResponseCalib[iDP][iAlgo][iMC][iP][iEta];
       if (thisHist->Integral () != 0) {
        thisHist->Scale (1./thisHist->Integral ());
       }

       TF1* calibFit = new TF1 (Form ("calibFit_iP%i", iP), "gaus (0)", 0, 2.0);
       if (nJets > 0) {
        calibFit->SetParameter (1, thisHist->GetMean ());
        calibFit->SetParameter (2, thisHist->GetStdDev ());
        thisHist->Fit (calibFit, "RN");
        double m = calibFit->GetParameter (1);
        double s = calibFit->GetParameter (2);
        if (calibFit) delete calibFit;

        calibFit = new TF1 (Form ("calibFit2_iP%i", iP), "gaus (0)", m - 2.0*s, m + 2.0*s);
        calibFit->SetParameter (1, thisHist->GetMean ());
        calibFit->SetParameter (2, thisHist->GetStdDev ());
        thisHist->Fit (calibFit, "RN");
       }

       if (iMC == plotMC && iP <= 11) {
        thisHist->SetMarkerStyle (kFullDotMedium);
        thisHist->SetLineColor (kBlue);
        thisHist->SetMarkerColor (kBlue);
        thisHist->Draw ("same e1 x0");

        calibFit->SetLineColor (kBlue);
        calibFit->Draw ("same");
        calibFits[iP] = calibFit;
       }

       jetEnergyScale[iMC][iEta]->SetBinContent (iP+1, calibFit->GetParameter (1));
       //jetEnergyScale[iMC][iEta]->SetBinContent (iP+1, thisHist->GetMean ());

       jetEnergyScale[iMC][iEta]->SetBinError (iP+1, calibFit->GetParError (1));
       //jetEnergyScale[iMC][iEta]->SetBinError (iP+1, thisHist->GetMeanError ());
       //jetEnergyScale[iMC][iEta]->SetBinError (iP+1, 0.000001);

       jetEnergyRes[iMC][iEta]->SetBinContent (iP+1, calibFit->GetParameter (2));
       //jetEnergyRes[iMC][iEta]->SetBinContent (iP+1, thisHist->GetStdDev ());

       jetEnergyRes[iMC][iEta]->SetBinError (iP+1, calibFit->GetParError (2));
       //jetEnergyRes[iMC][iEta]->SetBinError (iP+1, thisHist->GetStdDevError ());
       //jetEnergyRes[iMC][iEta]->SetBinError (iP+1, 0.000001);

       if (iMC == plotMC && iP <= 11) {
        myText (0.18, 0.90, kBlack, Form ("#mu = %s", FormatMeasurement (calibFit->GetParameter (1), calibFit->GetParError (1), 1)), 0.04 * 2);
        myText (0.18, 0.80, kBlack, Form ("#sigma = %s", FormatMeasurement (calibFit->GetParameter (2), calibFit->GetParError (2), 1)), 0.04 * 2);
        if (calcPtClosure) myText (0.18, 0.70, kBlack, Form ("%g < #it{p}_{T}^{true} < %g", pbins[iP], pbins[iP+1]), 0.04 * 2);
        else myText (0.18, 0.70, kBlack, Form ("%g < #it{E}_{true} < %g", pbins[iP], pbins[iP+1]), 0.04 * 2);
       }
       else {
        if (recoFit) delete recoFit;
        if (calibFit) delete calibFit;
       }
      }
     }

     if (iMC == plotMC) {
      jesCanvas->cd (1);
      myText (0.2, 0.8, kBlack, Form ("%g < #eta_{Lab} < %g", etabins[iEta], etabins[iEta+1]), 0.03*3);
      if (iAlgo == 0)
       myText (0.2, 0.65, kBlack, "HI jets", 0.03*3);
      else
       myText (0.2, 0.65, kBlack, "EMTopo jets");

      jesCanvas->SaveAs (Form ("%s/jetEnergyResponse/%s_iEta%i.pdf", plotPath.Data (), algo.Data (), iEta));

      Delete1DArray (recoFits, 12);
      Delete1DArray (calibFits, 12);
     }
     if (jesCanvas) delete jesCanvas;
    }
   }

   TCanvas* jesSummaryCanvas = new TCanvas ("jesSummaryCanvas", "", 800, 600);

   jesSummaryCanvas->cd ();
   gPad->SetLogx ();
   gStyle->SetErrorX (0.5);
   Color_t colors[15] = {kGray+2, kAzure, kViolet, kMagenta, kPink+10, kPink, kOrange+10, kOrange, kSpring+10, kSpring, kTeal+10, kTeal, kAzure+10, kGray, kBlack};

   double delta = 0.005;
   for (short iMC = plotMC; iMC < plotMC+2; iMC++) {
    const TString mcType = GetMCType (iMC);

    TGraphAsymmErrors** hArr = Get1DArray <TGraphAsymmErrors*> (numetabins+1);
    for (int iEta = 0; iEta <= numetabins; iEta++) {
     thisGraph = make_graph (jetEnergyScale[iMC][iEta]);
     double t_delta = 1.;
     if (iEta < numetabins) t_delta += delta * (iEta - (numetabins/2) + (int) (iEta>= (numetabins/2)));
     deltaize (thisGraph, t_delta, true);

     if (calcPtClosure) {
      thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{true} #left[GeV#right]");
      thisGraph->GetYaxis ()->SetTitle ("<#it{p}_{T}^{calib} / #it{p}_{T}^{true}>");
     }
     else {
      thisGraph->GetXaxis ()->SetTitle ("#it{E}_{true} #left[GeV#right]");
      thisGraph->GetYaxis ()->SetTitle ("<#it{E}_{calib} / #it{E}_{true}>");
     }

     if (plotMC == 1) {
      if (calcPtClosure) thisGraph->GetXaxis ()->SetRangeUser (50, 150);
      else thisGraph->GetXaxis ()->SetRangeUser (50, pbins[numpbins]);
     }
     else {
      thisGraph->GetXaxis ()->SetRangeUser (pbins[0], pbins[numpbins]);
     }

     if (plotMC == 1) thisGraph->GetYaxis ()->SetRangeUser (0.85, 1.10);
     else thisGraph->GetYaxis ()->SetRangeUser (0.85, 1.20);

     thisGraph->SetLineColor (colors[iEta]);
     thisGraph->SetMarkerColor (colors[iEta]);
     thisGraph->SetMarkerStyle (kFullDotLarge);
     if (iEta < numetabins) thisGraph->SetMarkerSize (0.5);
     else thisGraph->SetLineWidth (2);

     if (iEta == 0) thisGraph->Draw ("AP");
     else thisGraph->Draw ("P");
     hArr[iEta] = thisGraph;

     if (iEta == numetabins) myMarkerText (0.22, 0.49, colors[iEta], kFullCircle, "#eta-integrated", 1.25, 0.03);
     else myMarkerText (0.22+ (iEta>= (numetabins/2))*0.18, 0.45-0.04*iEta+0.04* (numetabins/2)* (iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
    }

    if (iMC == 0) myText (0.20, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
    else myText (0.20, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV", 0.04);
    myText (0.20, 0.78, kBlack, Form ("%i jets", jetCounter[iAlgo][iMC]), 0.04);
    myText (0.20, 0.90, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);

    jesSummaryCanvas->SaveAs (Form ("%s/%s_%s_JetEnergyScale.pdf", plotPath.Data (), algo.Data (), mcType.Data ()));
cout<<"e";
    Delete1DArray (hArr, numetabins+1);
   }

   delta = 0.014;
   for (short iMC = plotMC; iMC < plotMC+2; iMC++) {
    Color_t color = (iMC == plotMC ? mcOverlayColor:mcSignalColor);
    thisGraph = make_graph (jetEnergyScale[iMC][numetabins]);
    deltaize (thisGraph, 1 + (iMC-plotMC)*delta, true);

    if (calcPtClosure) {
     thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{true} #left[GeV#right]");
     thisGraph->GetYaxis ()->SetTitle ("<#it{p}_{T}^{calib} / #it{p}_{T}^{true}>");
    }
    else {
     thisGraph->GetXaxis ()->SetTitle ("#it{E}_{true} #left[GeV#right]");
     thisGraph->GetYaxis ()->SetTitle ("<#it{E}_{calib} / #it{E}_{true}>");
    }

    if (plotMC == 1) {
     if (calcPtClosure) thisGraph->GetXaxis ()->SetRangeUser (50, 150);
     else thisGraph->GetXaxis ()->SetRangeUser (50, pbins[numpbins]);
    }
    else {
     thisGraph->GetXaxis ()->SetRangeUser (pbins[0], pbins[numpbins]);
    }

    if (plotMC == 1) thisGraph->GetYaxis ()->SetRangeUser (0.98, 1.01);
    else thisGraph->GetYaxis ()->SetRangeUser (0.94, 1.05);

    thisGraph->SetLineWidth (2);
    thisGraph->SetLineColor (color);
    thisGraph->SetMarkerColor (color);
    thisGraph->SetMarkerStyle (kFullDotLarge);

    if (iMC == plotMC) thisGraph->Draw ("AP");
    else thisGraph->Draw ("P");

    if (plotMC == 1) myMarkerText (0.22, 0.28- (iMC-1)*0.06, color, kFullDotLarge, Form ("B.S. = %i ns", (iMC==1?25:100)), 1.25, 0.04);
    else myMarkerText (0.22, 0.28-iMC*0.06, color, kFullDotLarge, Form ("Pythia8 #it{pp} 8.16 TeV%s", (iMC==plotMC?" with #it{p}-Pb Data Overlay":"")), 1.25, 0.04);
   }

   jetEnergyScale[plotMC+1][numetabins]->Divide (jetEnergyScale[plotMC][numetabins]);
   thisGraph = make_graph (jetEnergyScale[plotMC+1][numetabins]);
   deltaize (thisGraph, 1-delta, true);
   thisGraph->SetLineWidth (2);
   thisGraph->SetLineColor (kBlack);
   thisGraph->SetMarkerColor (kBlack);
   thisGraph->Draw ("P");

   if (plotMC == 1) myMarkerText (0.22, 0.34, kBlack, kFullDotLarge, "Pythia8 #it{pp} 8.16 TeV, New / Old", 1.25, 0.04);
   else myMarkerText (0.22, 0.34, kBlack, kFullDotLarge, "Signal / Overlay", 1.25, 0.04);
   myText (0.20, 0.40, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
   jesSummaryCanvas->SaveAs (Form ("%s/%s_JetEnergyScale_MCcomparison.pdf", plotPath.Data (), algo.Data ()));

   delta = 0.005;
   for (short iMC = plotMC; iMC < plotMC+2; iMC++) {
    const TString mcType = GetMCType (iMC);
    TGraphAsymmErrors** hArr = Get1DArray <TGraphAsymmErrors*> (numetabins+1);
    for (int iEta = 0; iEta <= numetabins; iEta++) {
     thisGraph = make_graph (jetEnergyRes[iMC][iEta]);
     double t_delta = 1;
     if (iEta < numetabins) t_delta += delta * (iEta - (numetabins/2) + (int) (iEta>= (numetabins/2)));
     deltaize (thisGraph, t_delta, true);

     if (calcPtClosure) {
      thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{true} #left[GeV#right]");
      thisGraph->GetYaxis ()->SetTitle ("#sigma#left[<#it{p}_{T}^{calib} / #it{p}_{T}^{true}>#right]");
     }
     else {
      thisGraph->GetXaxis ()->SetTitle ("#it{E}_{true} #left[GeV#right]");
      thisGraph->GetYaxis ()->SetTitle ("#sigma#left[<#it{E}_{calib} / #it{E}_{true}>#right]");
     }

     if (plotMC == 1) {
      if (calcPtClosure) thisGraph->GetXaxis ()->SetRangeUser (50, 150);
      else thisGraph->GetXaxis ()->SetRangeUser (50, pbins[numpbins]);
     }
     else {
      thisGraph->GetXaxis ()->SetRangeUser (pbins[0], pbins[numpbins]);
     }

     thisGraph->GetYaxis ()->SetRangeUser (0, 0.3);

     thisGraph->SetLineColor (colors[iEta]);
     thisGraph->SetMarkerColor (colors[iEta]);
     thisGraph->SetMarkerStyle (kFullDotLarge);
     if (iEta < numetabins) thisGraph->SetMarkerSize (0.5);
     else thisGraph->SetLineWidth (2);

     if (iEta == 0) thisGraph->Draw ("AP");
     else thisGraph->Draw ("P");
     hArr[iEta] = thisGraph;

     if (iEta == numetabins) myMarkerText (0.57, 0.89, colors[iEta], kFullCircle, "#eta-integrated", 1.25, 0.03);
     else myMarkerText (0.57+ (iEta>= (numetabins/2))*0.18, 0.85-0.04*iEta+0.04* (numetabins/2)* (iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
    }

    if (iMC == 0) myText (0.20, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
    else myText (0.20, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV", 0.04);
    myText (0.20, 0.76, kBlack, Form ("%i jets", jetCounter[iAlgo][iMC]), 0.04);
    myText (0.20, 0.90, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);

    jesSummaryCanvas->SaveAs (Form ("%s/%s_%s_JetEnergyResolution.pdf", plotPath.Data (), algo.Data (), mcType.Data ()));
    Delete1DArray (hArr, numetabins+1);
   }

   //for (short iMC = 0; iMC < 3; iMC++) {
   // Color_t color = (iMC == 0 ? mcOverlayColor:mcSignalColor);
   // thisGraph = jetEnergyRes[iMC][numetabins];
   //  //thisGraph->GetXaxis ()->SetTitle ("#it{p}_{T}^{true} #left[GeV#right]");
   // thisGraph->GetXaxis ()->SetTitle ("#it{E}_{true} #left[GeV#right]");
   // thisGraph->GetYaxis ()->SetTitle ("Jet Energy Resolution");
   // thisGraph->GetYaxis ()->SetRangeUser (0, 0.3);
   // thisGraph->SetLineColor (color);
   // thisGraph->SetMarkerColor (color);
   // thisGraph->SetMarkerStyle (kFullDotLarge);
   // if (iMC == 0) thisGraph->Draw ("L ][");
   // else thisGraph->Draw ("same L ][");
   // myMarkerText (0.22, 0.28-iMC*0.06, color, kFullDotLarge, Form ("Pythia8 #it{pp} 8.16 TeV%s", (iMC == 0 ? " with #it{p}-Pb Overlay":"")), 1.25, 0.04);
   //}
   //myText (0.60, 0.82, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
   //jesSummaryCanvas->SaveAs (Form ("%s/%s_JetEnergyResolution_MCcomparison.pdf", plotPath.Data (), algo.Data ()));

   for (short iMC = plotMC; iMC < plotMC+2; iMC++) {
    for (int iEta = 0; iEta <= numetabins; iEta++) {
     if (jetEnergyScale[iMC][iEta]) delete jetEnergyScale[iMC][iEta];
     if (jetEnergyRes[iMC][iEta]) delete jetEnergyRes[iMC][iEta];
    }
   }
   if (jesSummaryCanvas) delete jesSummaryCanvas;
  }


  /**** Plots the photon energy response ****/
  TH1D** photonEnergyScale = Get1DArray <TH1D*> (numetabins+1);
  TH1D** photonEnergyRes = Get1DArray <TH1D*> (numetabins+1);
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   TCanvas* pesCanvas = new TCanvas (Form ("photonEnergyScalePad_%i", iEta), "", 800, 600);

   //if (iEta < numetabins) {
   photonEnergyScale[iEta] = new TH1D (Form ("photonEnergyScale_eta%i", iEta), "", numpbins, pbins);
   photonEnergyRes[iEta] = new TH1D (Form ("photonEnergyRes_eta%i", iEta), "", numpbins, pbins);
   //}

   pesCanvas->cd ();
   pesCanvas->Divide (4, 3);

   TF1* calibFits[12];

   for (short iDP = 0; iDP < numdpbins; iDP++) {
    for (short iP = 0; iP < numpbins; iP++) {

     if (! (dpbins[iDP] <= pbins[iP]) || (numdpbins > 1 && ! (pbins[iP+1] <= dpbins[iDP+1])))
      continue;

     if (iP <= 11)
      pesCanvas->cd (iP+1);

     gPad->SetLogy (true);

     thisHist = photonEnergyResponse[iDP][0][iP][iEta];
     if (thisHist->Integral () != 0)
      thisHist->Scale (1./thisHist->Integral ());

     int nGammas = nGamma[iDP][0][iP][iEta];

     TF1* calibFit = new TF1 ("calibFit", "gaus (0)", 0, 2.0);
     if (nGammas > 0) {
      calibFit->SetParameter (1, thisHist->GetMean ());
      calibFit->SetParameter (2, thisHist->GetStdDev ());
      thisHist->Fit (calibFit, "RN");
      float m = calibFit->GetParameter (1);
      float s = calibFit->GetParameter (2);
      if (calibFit) delete calibFit;

      calibFit = new TF1 ("calibFit2", "gaus (0)", m - 3.5*s, m + 3.5*s);
      thisHist->Fit (calibFit, "RN");
      //m = calibFit->GetParameter (1);
      //s = calibFit->GetParameter (2);
      //if (calibFit) delete calibFit;

      //calibFit = new TF1 ("calibFit3", "gaus (0)", m - 3.5*s, m + 3.5*s);
      //thisHist->Fit (calibFit, "RNL");
     }

     if (iP <= 11) {
      thisHist->SetAxisRange (3e-6, 2e1, "Y");
      thisHist->SetMarkerStyle (kFullDotMedium);
      thisHist->SetLineColor (kBlue);
      thisHist->SetMarkerColor (kBlue);
      thisHist->GetYaxis ()->SetTitle ("Counts / Total");
      thisHist->Draw ("e1 x0");

      calibFit->SetLineColor (kBlue);
      calibFit->Draw ("same");
      calibFits[iP] = calibFit;
     }

     photonEnergyScale[iEta]->SetBinContent (iP+1, calibFit->GetParameter (1));
     //photonEnergyScale[iEta]->SetBinContent (iP+1, thisHist->GetMean ());

     photonEnergyScale[iEta]->SetBinError (iP+1, calibFit->GetParError (1));
     //photonEnergyScale[iEta]->SetBinError (iP+1, thisHist->GetMeanError ());
     //photonEnergyScale[iEta]->SetBinError (iP+1, 0.000001);

     photonEnergyRes[iEta]->SetBinContent (iP+1, calibFit->GetParameter (2));
     //photonEnergyRes[iEta]->SetBinContent (iP+1, thisHist->GetStdDev ());

     photonEnergyRes[iEta]->SetBinError (iP+1, calibFit->GetParError (2));
     //photonEnergyRes[iEta]->SetBinError (iP+1, thisHist->GetStdDevError ());
     //photonEnergyRes[iEta]->SetBinError (iP+1, 0.000001);

     //myText (0.18, 0.9, kBlack, "Photon energy response");
     //myText (0.18, 0.83, kBlack, Form ("%i photons", nGammas));
     if (iP <= 11) {
      myText (0.18, 0.90, kBlack, Form ("#mu = %s", FormatMeasurement (calibFit->GetParameter (1), calibFit->GetParError (1), 1)), 0.04 * 2);
      myText (0.18, 0.80, kBlack, Form ("#sigma = %s", FormatMeasurement (calibFit->GetParameter (2), calibFit->GetParError (2), 1)), 0.04 * 2);
      myText (0.18, 0.70, kBlack, Form ("%g < #it{p}_{T}^{reco} < %g", pbins[iP], pbins[iP+1]), 0.04 * 2);
     }
     else {
      if (calibFit) delete calibFit;
     }
    }
   }
   pesCanvas->SaveAs (Form ("%s/photonEnergyResponse/iEta%i.pdf", plotPath.Data (), iEta));
   if (pesCanvas) delete pesCanvas;

   for (short iP = 0; iP <= 11; iP++)
    if (calibFits[iP]) delete calibFits[iP];
  }

  TCanvas* pesSummaryCanvas = new TCanvas ("pesSummaryCanvas", "", 800, 600);
  pesSummaryCanvas->cd ();
  gPad->SetLogx ();
  double delta = 0.007;
  TGraphAsymmErrors** hArr = Get1DArray <TGraphAsymmErrors*> (numetabins+1);
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   thisGraph = make_graph (photonEnergyScale[iEta]);
   double t_delta = 1.;
   if (iEta < numetabins) t_delta += delta * (iEta - (numetabins/2) + (int) (iEta>= (numetabins/2)));
   deltaize (thisGraph, t_delta, true);

   thisGraph->GetXaxis ()->SetTitle ("#it{E}_{T}^{#gamma} #left[GeV#right]");
   thisGraph->GetYaxis ()->SetTitle ("<#it{E}_{T}^{calib} / #it{E}_{T}^{true}>");
   thisGraph->GetYaxis ()->SetRangeUser (0.985, 1.015);
   thisGraph->SetLineColor (colors[iEta]);
   thisGraph->SetMarkerColor (colors[iEta]);
   thisGraph->SetMarkerStyle (kFullDotLarge);
   if (iEta < numetabins) thisGraph->SetMarkerSize (0.5);
   else thisGraph->SetLineWidth (2);

   if (iEta == 0) thisGraph->Draw ("AP");
   else thisGraph->Draw ("P");
   hArr[iEta] = thisGraph;

   if (iEta == numetabins) myMarkerText (0.22, 0.49, colors[iEta], kFullCircle, "#eta-integrated", 1.25, 0.03);
   else myMarkerText (0.22+ (iEta>= (numetabins/2))*0.18, 0.45-0.04*iEta+0.04* (numetabins/2)* (iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
  }
  myText (0.60, 0.90, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
  myText (0.48, 0.84, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
  pesSummaryCanvas->SaveAs (Form ("%s/PhotonEnergyScale.pdf", plotPath.Data ()));
  Delete1DArray (hArr, numetabins+1);

  hArr = Get1DArray <TGraphAsymmErrors*> (numetabins+1);
  for (short iEta = 0; iEta <= numetabins; iEta++) {
   thisGraph = make_graph (photonEnergyRes[iEta]);
   double t_delta = 1.;
   if (iEta < numetabins) t_delta += delta * (iEta - (numetabins/2) + (int) (iEta>= (numetabins/2)));
   deltaize (thisGraph, t_delta, true);

   thisGraph->GetXaxis ()->SetTitle ("#it{E}_{T}^{#gamma} #left[GeV#right]");
   thisGraph->GetYaxis ()->SetTitle ("#sigma#left[<#it{E}_{T}^{calib} / #it{E}_{T}^{true}>#right]");
   thisGraph->GetYaxis ()->SetRangeUser (0, 0.07);
   thisGraph->SetLineColor (colors[iEta]);
   thisGraph->SetMarkerColor (colors[iEta]);
   thisGraph->SetMarkerStyle (kFullDotLarge);
   if (iEta < numetabins) thisGraph->SetMarkerSize (0.5);
   else thisGraph->SetLineWidth (2);

   if (iEta == 0) thisGraph->Draw ("AP");
   else thisGraph->Draw ("P");
   hArr[iEta] = thisGraph;

   if (iEta == numetabins) myMarkerText (0.22, 0.90, colors[iEta], kFullCircle, "#eta-integrated", 1.25, 0.03);
   else myMarkerText (0.22+ (iEta>= (numetabins/2))*0.18, 0.86-0.04*iEta+0.04* (numetabins/2)* (iEta>=numetabins/2), colors[iEta], kFullCircle, Form ("%g < #eta < %g", etabins[iEta], etabins[iEta+1]), 1.25, 0.03);
  }
  myText (0.20, 0.28, kBlack, "#bf{#it{ATLAS}} Simulation Internal", 0.04);
  myText (0.20, 0.22, kBlack, "Pythia8 #it{pp} 8.16 TeV with #it{p}-Pb Overlay", 0.04);
  pesSummaryCanvas->SaveAs (Form ("%s/PhotonEnergyResolution.pdf", plotPath.Data ()));

  Delete1DArray (photonEnergyScale, numetabins+1);
  Delete1DArray (photonEnergyRes, numetabins+1);
  if (pesSummaryCanvas) delete pesSummaryCanvas;
  Delete1DArray (hArr, numetabins+1);


  return;
}

}
