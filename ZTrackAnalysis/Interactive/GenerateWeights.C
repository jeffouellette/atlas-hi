#ifndef __GenerateWeights_C__
#define __GenerateWeights_C__

#include "Params.h"
#include "ZTrackUtilities.h"

#include <ArrayTemplates.h>

#include <TEfficiency.h>
#include <TTree.h>

#include <iostream>

void SafeWrite (TObject* tobj) {
  if (tobj)
    tobj->Write ();
}

const int gw_nFCalBins = 100;
const double* gw_fcalBins = linspace (0, 5200, gw_nFCalBins);
const int gw_nQ2Bins = 20;
const double* gw_q2Bins = linspace (0, 0.3, gw_nQ2Bins);

const int gw_nNchBins = 160;
const double* gw_nchBins = linspace (-0.5, 160.5, gw_nNchBins);

//const double gw_fcalBins[11]   = {0, 63.719, 144.14, 289.595, 525.092, 875.41, 1368.75, 2046.51, 2989.31, 3618.44, 5200};
//const int gw_nFCalBins = sizeof (gw_fcalBins) / sizeof (gw_fcalBins[0]) - 1;
//const int gw_nQ2Bins = 10;
//const double* gw_q2Bins = linspace (0, 1, gw_nQ2Bins);

TH1D* referenceFCalDist = nullptr;
TH1D* referenceQ2Dist[numFinerCentBins];

TH1D* h_PbPbFCalDist = nullptr;
TH1D* h_PbPbFCal_weights = nullptr;
TH1D* h_PbPbQ2Dist[numFinerCentBins];
TH1D* h_PbPbQ2_weights[numFinerCentBins];

TH1D* h_ppNchDist = nullptr;
TH1D* h_ppNch_weights = nullptr;


void GenerateWeights (const TString name, const TString inFileName = "outFile.root", const TString outFileName = "eventWeightsFile.root") {

  if (name == "data")
    SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
  else if (name == "mc")
    SetupDirectories ("MCAnalysis/", "ZTrackAnalysis/");
  else if (name == "minbias")
    SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
  else if (name == "hijing")
    SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
  else if (name == "truth")
    SetupDirectories ("TruthAnalysis/", "ZTrackAnalysis/");
  else {
    cout << "Unrecognized analysis name! Quitting." << endl;
    return;
  }

  TFile* inFile = new TFile (Form ("%s/Nominal/%s", rootPath.Data (), inFileName.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  h_PbPbFCalDist = new TH1D (Form ("h_PbPbFCalDist_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins);
  h_PbPbFCal_weights = new TH1D (Form ("h_PbPbFCal_weights_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins);
  h_PbPbFCalDist->Sumw2 ();
  h_PbPbFCal_weights->Sumw2 ();
  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    h_PbPbQ2Dist[iCent] = new TH1D (Form ("h_PbPbQ2Dist_iCent%i_%s", iCent, name.Data ()), "", gw_nQ2Bins, gw_q2Bins);
    h_PbPbQ2_weights[iCent] = new TH1D (Form ("h_PbPbQ2_weights_iCent%i_%s", iCent, name.Data ()), "", gw_nQ2Bins, gw_q2Bins);
    h_PbPbQ2Dist[iCent]->Sumw2 ();
    h_PbPbQ2_weights[iCent]->Sumw2 ();
  }
  h_ppNchDist = new TH1D (Form ("h_ppNchDist_%s", name.Data ()), "", gw_nNchBins, gw_nchBins);
  h_ppNchDist->Sumw2 ();
  h_ppNch_weights = new TH1D (Form ("h_ppNch_weights_%s", name.Data ()), "", gw_nNchBins, gw_nchBins);
  h_ppNch_weights->Sumw2 ();


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    float event_weight = 1, fcal_et = 0, q2 = 0, zpt = 0;

    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    //PbPbTree->SetBranchAddress ("z_pt",         &zpt);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      //if (zpt < 25)
      //  continue; // candidate events

      short iCent = 0;
      while (iCent < numFinerCentBins) {
        if (fcal_et < finerCentBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numFinerCentBins-1)
        continue;

      h_PbPbFCalDist->Fill (fcal_et, event_weight);
      h_PbPbQ2Dist[iCent]->Fill (q2, event_weight);
    }
    cout << "Done 1st Pb+Pb loop." << endl;
    if (h_PbPbFCalDist->Integral () > 0)
      h_PbPbFCalDist->Scale (1./h_PbPbFCalDist->Integral ());
    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
      if (h_PbPbQ2Dist[iCent]->Integral () > 0)
        h_PbPbQ2Dist[iCent]->Scale (1./h_PbPbQ2Dist[iCent]->Integral ());
    }

    if (name == "data") {
      for (int ix = 1; ix <= h_PbPbFCal_weights->GetNbinsX (); ix++) {
        h_PbPbFCal_weights->SetBinContent (ix, 1);
      }
      for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
        for (int ix = 1; ix <= h_PbPbQ2_weights[iCent]->GetNbinsX (); ix++) {
          h_PbPbQ2_weights[iCent]->SetBinContent (ix, 1);
        }
      }
    }
    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceFCalDist = (TH1D*)ztrackFile->Get ("h_PbPbFCalDist_data");

      TH1D* referenceQ2Dist[numFinerCentBins];
      for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
        referenceQ2Dist[iCent] = (TH1D*)ztrackFile->Get (Form ("h_PbPbQ2Dist_iCent%i_data", iCent));
      }

      for (int ix = 1; ix <= h_PbPbFCal_weights->GetNbinsX (); ix++) {
        const double fcal_weight = (h_PbPbFCalDist->GetBinContent (ix) != 0 ? referenceFCalDist->GetBinContent (ix) / h_PbPbFCalDist->GetBinContent (ix) : 0);
        h_PbPbFCal_weights->SetBinContent (ix, fcal_weight);
      }

      for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
        for (int ix = 1; ix <= h_PbPbQ2Dist[iCent]->GetNbinsX (); ix++) {
          h_PbPbQ2Dist[iCent]->SetBinContent (ix, 0);
          h_PbPbQ2Dist[iCent]->SetBinError (ix, 0);
        }
      }

      for (int iEvt = 0; iEvt < nEvts; iEvt++) {
        if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
          cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
        PbPbTree->GetEntry (iEvt);

        //if (zpt < 25)
        //  continue; // candidate events

        short iCent = 0;
        while (iCent < numFinerCentBins) {
          if (fcal_et < finerCentBins[iCent])
            break;
          else
            iCent++;
        }
        if (iCent < 1 || iCent > numFinerCentBins-1)
          continue;

        event_weight *= h_PbPbFCal_weights->GetBinContent (h_PbPbFCal_weights->FindBin (fcal_et));

        h_PbPbQ2Dist[iCent]->Fill (q2, event_weight);
      }
      cout << "Done 2nd Pb+Pb loop." << endl;

      for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
        if (h_PbPbQ2Dist[iCent]->Integral () > 0)
          h_PbPbQ2Dist[iCent]->Scale (1./h_PbPbQ2Dist[iCent]->Integral ());
      }

      for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
        for (int ix = 1; ix <= h_PbPbQ2_weights[iCent]->GetNbinsX (); ix++) {
          const double q2_weight = (h_PbPbQ2Dist[iCent]->GetBinContent (ix) != 0 ? referenceQ2Dist[iCent]->GetBinContent (ix) / h_PbPbQ2Dist[iCent]->GetBinContent (ix) : 0);
          h_PbPbQ2_weights[iCent]->SetBinContent (ix, q2_weight);
        }
      }
      ztrackFile->Close ();

      if (name == "data")
        SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      else if (name == "mc")
        SetupDirectories ("MCAnalysis/", "ZTrackAnalysis/");
      else if (name == "minbias")
        SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
      else if (name == "hijing")
        SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
      else if (name == "truth")
        SetupDirectories ("TruthAnalysis/", "ZTrackAnalysis/");
    }
  }

  if (ppTree) {
    float event_weight = 0;
    //float zpt = 0;
    int ntrk = 0;
    ppTree->SetBranchAddress ("event_weight", &event_weight);
    //ppTree->SetBranchAddress ("z_pt", &zpt);
    ppTree->SetBranchAddress ("ntrk", &ntrk);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      //if (zpt < 25)
      //  continue; // candidate events

      h_ppNchDist->Fill (ntrk, event_weight);
    }
    cout << "Done pp loop." << endl;

    if (name == "data") {
      for (int ix = 1; ix <= h_ppNch_weights->GetNbinsX (); ix++) {
        h_ppNch_weights->SetBinContent (ix, 1);
      }
    } else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceNchDist = (TH1D*)ztrackFile->Get ("h_ppNchDist_data");

      for (int ix = 1; ix <= h_ppNch_weights->GetNbinsX (); ix++) {
        const double nch_weight = (h_ppNchDist->GetBinContent (ix) != 0 ? referenceNchDist->GetBinContent (ix) / h_ppNchDist->GetBinContent (ix) : 0);
        h_ppNch_weights->SetBinContent (ix, nch_weight);
      }

      ztrackFile->Close ();

      if (name == "data")
        SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      else if (name == "mc")
        SetupDirectories ("MCAnalysis/", "ZTrackAnalysis/");
      else if (name == "minbias")
        SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
      else if (name == "hijing")
        SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
      else if (name == "truth")
        SetupDirectories ("TruthAnalysis/", "ZTrackAnalysis/");
    }
  }


  TFile* eventWeightsFile = new TFile (Form ("%s/%s", rootPath.Data (), outFileName.Data ()), "recreate");

  SafeWrite (h_PbPbFCalDist);
  SafeWrite (h_PbPbFCal_weights);
  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    SafeWrite (h_PbPbQ2Dist[iCent]);
    SafeWrite (h_PbPbQ2_weights[iCent]);
  }
  SafeWrite (h_ppNchDist);
  SafeWrite (h_ppNch_weights);

  eventWeightsFile->Close ();
}


void GenerateDataWeights () {
  GenerateWeights ("data");
}

void GenerateMCWeights () {
  GenerateWeights ("mc", "eventWeightsTree.root");
}

void GenerateTruthWeights () {
  GenerateWeights ("truth");
}

void Generate2015HijingWeights () {
  GenerateWeights ("hijing", "PbPb_Hijing_15.root", "PbPb_Hijing_15_eventWeights.root");
}

void Generate2018HijingWeights () {
  GenerateWeights ("hijing", "PbPb_Hijing_18.root", "PbPb_Hijing_18_eventWeights.root");
}

#endif
