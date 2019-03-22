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

const int gw_nFCalBins = 100;
const double* gw_fcalBins = linspace (0, 5200, gw_nFCalBins);
const int gw_nQ2Bins = 25;
const double* gw_q2Bins = linspace (0, 1, gw_nQ2Bins);
const int gw_nVertZBins = 50;
const double* gw_vertZBins = linspace (-200, 200, gw_nVertZBins);

//const double gw_fcalBins[11]   = {0, 63.719, 144.14, 289.595, 525.092, 875.41, 1368.75, 2046.51, 2989.31, 3618.44, 5200};
//const int gw_nFCalBins = sizeof (gw_fcalBins) / sizeof (gw_fcalBins[0]) - 1;
//const int gw_nQ2Bins = 10;
//const double* gw_q2Bins = linspace (0, 1, gw_nQ2Bins);
//const int gw_nVertZBins = 18;
//const double gw_vertZBins[19] = {-200, -175, -150, -125, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 125, 150, 175, 200};

TH1D* PbPbFCalDist = nullptr;
TH1D* PbPbVZDist = nullptr;
TH1D* PbPbQ2Dist = nullptr;
TH3D* PbPbEventReweights = nullptr;

TH1D* ppVZDist = nullptr;
TH1D* ppEventReweights = nullptr;


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

  PbPbFCalDist = new TH1D (Form ("PbPbFCalDist_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins);
  PbPbQ2Dist = new TH1D (Form ("PbPbQ2Dist_%s", name.Data ()), "", gw_nQ2Bins, gw_q2Bins);
  PbPbVZDist = new TH1D (Form ("PbPbVZDist_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  PbPbFCalDist->Sumw2 ();
  PbPbQ2Dist->Sumw2 ();
  PbPbVZDist->Sumw2 ();
  PbPbEventReweights = new TH3D (Form ("PbPbEventReweights_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins, gw_nQ2Bins, gw_q2Bins, gw_nVertZBins, gw_vertZBins);
  PbPbEventReweights->Sumw2 ();

  ppVZDist = new TH1D (Form ("ppVZDist_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  ppVZDist->Sumw2 ();
  ppEventReweights = new TH1D (Form ("ppEventReweights_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  ppEventReweights->Sumw2 ();


  float event_weight = 0, fcal_et = 0, q2 = 0, vz = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("vz",           &vz);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      PbPbFCalDist->Fill (fcal_et, event_weight);
      PbPbQ2Dist->Fill (q2, event_weight);
      PbPbVZDist->Fill (vz, event_weight);
    }
    cout << endl;
    if (PbPbFCalDist->Integral () > 0)
      PbPbFCalDist->Scale (1./PbPbFCalDist->Integral (), "width");
    if (PbPbQ2Dist->Integral () > 0)
      PbPbQ2Dist->Scale (1./PbPbQ2Dist->Integral (), "width");
    if (PbPbVZDist->Integral () > 0)
      PbPbVZDist->Scale (1./PbPbVZDist->Integral (), "width");

    if (name == "data") {
      for (int ix = 1; ix <= PbPbEventReweights->GetNbinsX (); ix++)
        for (int iy = 1; iy <= PbPbEventReweights->GetNbinsY (); iy++)
          for (int iz = 1; iz <= PbPbEventReweights->GetNbinsZ (); iz++)
            PbPbEventReweights->SetBinContent (ix, iy, iz, 1); // should be 1 in the reference
    }
    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceFCalDist = (TH1D*)ztrackFile->Get ("PbPbFCalDist_data");
      TH1D* referenceVZDist = (TH1D*)ztrackFile->Get ("PbPbVZDist_data");
      TH1D* referenceQ2Dist = (TH1D*)ztrackFile->Get ("PbPbQ2Dist_data");
      for (int ix = 1; ix <= PbPbEventReweights->GetNbinsX (); ix++) {
        const double fcalReweight = (PbPbFCalDist->GetBinContent (ix) != 0 ? referenceFCalDist->GetBinContent (ix) / PbPbFCalDist->GetBinContent (ix) : 0);
        for (int iy = 1; iy <= PbPbEventReweights->GetNbinsY (); iy++) {
          const double q2Reweight = (PbPbQ2Dist->GetBinContent (iy) != 0 ? referenceQ2Dist->GetBinContent (iy) / PbPbQ2Dist->GetBinContent (iy) : 0);
          for (int iz = 1; iz <= PbPbEventReweights->GetNbinsZ (); iz++) {
            const double vzReweight = (PbPbVZDist->GetBinContent (iz) != 0 ? referenceVZDist->GetBinContent (iz) / PbPbVZDist->GetBinContent (iz) : 0);
            PbPbEventReweights->SetBinContent (ix, iy, iz, fcalReweight * vzReweight * q2Reweight);
          }
        }
      }
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

      ppVZDist->Fill (vz, event_weight);
    }
    cout << endl;
    if (ppVZDist->Integral () > 0)
      ppVZDist->Scale (1./ppVZDist->Integral ());

    if (name == "data") {
      for (int ix = 1; ix <= ppEventReweights->GetNbinsX (); ix++)
        ppEventReweights->SetBinContent (ix, 1);
    }
    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceVZDist = (TH1D*)ztrackFile->Get ("ppVZDist_data");
      for (int ix = 1; ix <= ppEventReweights->GetNbinsX (); ix++)
        ppEventReweights->SetBinContent (ix, ppVZDist->GetBinContent (ix) != 0. ? referenceVZDist->GetBinContent (ix) / ppVZDist->GetBinContent (ix) : 0);
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

  SafeWrite (PbPbFCalDist);
  SafeWrite (PbPbVZDist);
  SafeWrite (PbPbQ2Dist);
  SafeWrite (PbPbEventReweights);
  SafeWrite (ppVZDist);
  SafeWrite (ppEventReweights);

  eventWeightsFile->Close ();
}

#endif
