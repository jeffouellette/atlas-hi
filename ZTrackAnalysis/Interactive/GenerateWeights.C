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

//const int gw_nFCalBins = 100;
//const double* gw_fcalBins = linspace (0, 5200, gw_nFCalBins);
//const int gw_nQ2Bins = 25;
//const double* gw_q2Bins = linspace (0, 1, gw_nQ2Bins);
//const int gw_nVertZBins = 50;
//const double* gw_vertZBins = linspace (-200, 200, gw_nVertZBins);

const double gw_fcalBins[11]   = {0, 63.719, 144.14, 289.595, 525.092, 875.41, 1368.75, 2046.51, 2989.31, 3618.44, 5200};
const int gw_nFCalBins = sizeof (gw_fcalBins) / sizeof (gw_fcalBins[0]) - 1;
const int gw_nQ2Bins = 10;
const double* gw_q2Bins = linspace (0, 1, gw_nQ2Bins);
const int gw_nVertZBins = 18;
const double gw_vertZBins[19] = {-200, -175, -150, -125, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 125, 150, 175, 200};

TH1D* h_PbPbFCalDist = nullptr;
TH1D* h_PbPbVZDist = nullptr;
TH1D* h_PbPbQ2Dist = nullptr;
TH3D* h_PbPbEventReweights = nullptr;

TH1D* h_ppVZDist = nullptr;
TH1D* h_ppEventReweights = nullptr;


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

  h_PbPbFCalDist = new TH1D (Form ("h_PbPbFCalDist_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins);
  h_PbPbQ2Dist = new TH1D (Form ("h_PbPbQ2Dist_%s", name.Data ()), "", gw_nQ2Bins, gw_q2Bins);
  h_PbPbVZDist = new TH1D (Form ("h_PbPbVZDist_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  h_PbPbFCalDist->Sumw2 ();
  h_PbPbQ2Dist->Sumw2 ();
  h_PbPbVZDist->Sumw2 ();
  h_PbPbEventReweights = new TH3D (Form ("h_PbPbEventReweights_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins, gw_nQ2Bins, gw_q2Bins, gw_nVertZBins, gw_vertZBins);
  h_PbPbEventReweights->Sumw2 ();

  h_ppVZDist = new TH1D (Form ("h_ppVZDist_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  h_ppVZDist->Sumw2 ();
  h_ppEventReweights = new TH1D (Form ("h_ppEventReweights_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  h_ppEventReweights->Sumw2 ();


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

      h_PbPbFCalDist->Fill (fcal_et, event_weight);
      h_PbPbQ2Dist->Fill (q2, event_weight);
      h_PbPbVZDist->Fill (vz, event_weight);
    }
    cout << endl;
    if (h_PbPbFCalDist->Integral () > 0)
      h_PbPbFCalDist->Scale (1./h_PbPbFCalDist->Integral (), "width");
    if (h_PbPbQ2Dist->Integral () > 0)
      h_PbPbQ2Dist->Scale (1./h_PbPbQ2Dist->Integral (), "width");
    if (h_PbPbVZDist->Integral () > 0)
      h_PbPbVZDist->Scale (1./h_PbPbVZDist->Integral (), "width");

    if (name == "data") {
      for (int ix = 1; ix <= h_PbPbEventReweights->GetNbinsX (); ix++)
        for (int iy = 1; iy <= h_PbPbEventReweights->GetNbinsY (); iy++)
          for (int iz = 1; iz <= h_PbPbEventReweights->GetNbinsZ (); iz++)
            h_PbPbEventReweights->SetBinContent (ix, iy, iz, 1); // should be 1 in the reference
    }
    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceFCalDist = (TH1D*)ztrackFile->Get ("h_PbPbFCalDist_data");
      TH1D* referenceVZDist = (TH1D*)ztrackFile->Get ("h_PbPbVZDist_data");
      TH1D* referenceQ2Dist = (TH1D*)ztrackFile->Get ("h_PbPbQ2Dist_data");
      for (int ix = 1; ix <= h_PbPbEventReweights->GetNbinsX (); ix++) {
        const double fcalReweight = (h_PbPbFCalDist->GetBinContent (ix) != 0 ? referenceFCalDist->GetBinContent (ix) / h_PbPbFCalDist->GetBinContent (ix) : 0);
        for (int iy = 1; iy <= h_PbPbEventReweights->GetNbinsY (); iy++) {
          const double q2Reweight = (h_PbPbQ2Dist->GetBinContent (iy) != 0 ? referenceQ2Dist->GetBinContent (iy) / h_PbPbQ2Dist->GetBinContent (iy) : 0);
          for (int iz = 1; iz <= h_PbPbEventReweights->GetNbinsZ (); iz++) {
            const double vzReweight = (h_PbPbVZDist->GetBinContent (iz) != 0 ? referenceVZDist->GetBinContent (iz) / h_PbPbVZDist->GetBinContent (iz) : 0);
            h_PbPbEventReweights->SetBinContent (ix, iy, iz, fcalReweight * vzReweight * q2Reweight);
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

      h_ppVZDist->Fill (vz, event_weight);
    }
    cout << endl;
    if (h_ppVZDist->Integral () > 0)
      h_ppVZDist->Scale (1./h_ppVZDist->Integral ());

    if (name == "data") {
      for (int ix = 1; ix <= h_ppEventReweights->GetNbinsX (); ix++)
        h_ppEventReweights->SetBinContent (ix, 1);
    }
    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceVZDist = (TH1D*)ztrackFile->Get ("h_ppVZDist_data");
      for (int ix = 1; ix <= h_ppEventReweights->GetNbinsX (); ix++)
        h_ppEventReweights->SetBinContent (ix, h_ppVZDist->GetBinContent (ix) != 0. ? referenceVZDist->GetBinContent (ix) / h_ppVZDist->GetBinContent (ix) : 0);
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

  SafeWrite (h_PbPbFCalDist);
  SafeWrite (h_PbPbVZDist);
  SafeWrite (h_PbPbQ2Dist);
  SafeWrite (h_PbPbEventReweights);
  SafeWrite (h_ppVZDist);
  SafeWrite (h_ppEventReweights);

  eventWeightsFile->Close ();
}

#endif
