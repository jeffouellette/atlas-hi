#ifndef __GenerateWeights_C__
#define __GenerateWeights_C__

#include "Params.h"
#include "ZTrackUtilities.h"

#include <ArrayTemplates.h>

#include <iostream>

void SafeWrite (TObject* tobj) {
  if (tobj)
    tobj->Write ();
}

TH3D** PbPbEventInfoDist  = Get1DArray <TH3D*> (2); // iData
TH3D** PbPbEventReweights = Get1DArray <TH3D*> (2); // iData
TH1D** ppEventInfoDist    = Get1DArray <TH1D*> (2); // iData
TH1D** ppEventReweights   = Get1DArray <TH1D*> (2); // iData

const double gw_fcalBins[11]   = {0, 63.719, 144.14, 289.595, 525.092, 875.41, 1368.75, 2046.51, 2989.31, 3618.44, 5200};
const int gw_nFCalBins = sizeof (gw_fcalBins) / sizeof (gw_fcalBins[0]) - 1;
const int gw_nPsi2Bins = 10;
const double* gw_psi2Bins = linspace (-pi/2, pi/2, nPsi2Bins);
const int gw_nVertZBins = 10;
const double* gw_vertZBins = linspace (-150, 150, nVertZBins);


void GenerateWeights (const TString name) {

  if (name == "data")
    SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
  else if (name == "mc")
    SetupDirectories ("MCAnalysis/", "ZTrackAnalysis/");
  else if (name == "minbias")
    SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
  else if (name == "truth")
    SetupDirectories ("TruthAnalysis/", "ZTrackAnalysis/");
  else {
    cout << "Unrecognized analysis name! Quitting." << endl;
    return;
  }

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  TH3D* PbPbEventInfoDist = new TH3D (Form ("PbPbEventInfoDist_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins, gw_nPsi2Bins, gw_psi2Bins, gw_nVertZBins, gw_vertZBins);
  PbPbEventInfoDist->Sumw2 ();
  TH3D* PbPbEventReweights = new TH3D (Form ("PbPbEventReweights_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins, gw_nPsi2Bins, gw_psi2Bins, gw_nVertZBins, gw_vertZBins);
  PbPbEventReweights->Sumw2 ();

  TH1D* ppEventInfoDist = new TH1D (Form ("ppEventInfoDist_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  ppEventInfoDist->Sumw2 ();
  TH1D* ppEventReweights = new TH1D (Form ("ppEventReweights_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  ppEventReweights->Sumw2 ();


  float event_weight = 0, fcal_et = 0, q2 = 0, psi2 = 0, vz = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("vz",           &vz);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      PbPbEventInfoDist->Fill (fcal_et, psi2, vz, event_weight);
    }
    cout << endl;
    if (PbPbEventInfoDist->Integral () > 0)
      PbPbEventInfoDist->Scale (1./PbPbEventInfoDist->Integral (), "width");

    if (name == "data") {
      PbPbEventReweights->Add (PbPbEventInfoDist);
      PbPbEventReweights->Divide (PbPbEventInfoDist);
    }
    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      PbPbEventReweights->Add ((TH3D*)ztrackFile->Get ("PbPbEventInfoDist_data"));
      PbPbEventReweights->Divide (PbPbEventInfoDist);
      ztrackFile->Close ();
    }

    if (name == "data")
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
    else if (name == "mc")
      SetupDirectories ("MCAnalysis/", "ZTrackAnalysis/");
    else if (name == "minbias")
      SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
    else if (name == "truth")
      SetupDirectories ("TruthAnalysis/", "ZTrackAnalysis/");
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("event_weight", &event_weight);
    ppTree->SetBranchAddress ("vz",           &vz);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      ppEventInfoDist->Fill (vz, event_weight);
    }
    cout << endl;
    if (ppEventInfoDist->Integral () > 0)
      ppEventInfoDist->Scale (1./ppEventInfoDist->Integral ());

    if (name == "data") {
      ppEventReweights->Add (ppEventInfoDist);
      ppEventReweights->Divide (ppEventInfoDist);
    }
    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      ppEventReweights->Add ((TH3D*)ztrackFile->Get ("ppEventInfoDist_data"));
      ppEventReweights->Divide (ppEventInfoDist);
      ztrackFile->Close ();
    }

    if (name == "data")
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
    else if (name == "mc")
      SetupDirectories ("MCAnalysis/", "ZTrackAnalysis/");
    else if (name == "minbias")
      SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
    else if (name == "truth")
      SetupDirectories ("TruthAnalysis/", "ZTrackAnalysis/");
  }

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "recreate");

  SafeWrite (PbPbEventInfoDist);
  SafeWrite (PbPbEventReweights);
  SafeWrite (ppEventInfoDist);
  SafeWrite (ppEventReweights);

  eventWeightsFile->Close ();
}

#endif
