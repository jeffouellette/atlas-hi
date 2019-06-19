#ifndef __GenerateWeights_C__
#define __GenerateWeights_C__

#include "Params.h"
#include "ZTrackUtilities.h"

#include <ArrayTemplates.h>

#include <TEfficiency.h>

#include <iostream>

void SafeWrite (TObject* tobj) {
  if (tobj)
    tobj->Write ();
}

const int gw_nFCalBins = 100;
const double* gw_fcalBins = linspace (0, 5200, gw_nFCalBins);
const int gw_nQ2Bins = 50;
const double* gw_q2Bins = linspace (0, 1, gw_nQ2Bins);
const int gw_nVertZBins = 50;
const double* gw_vertZBins = linspace (-200, 200, gw_nVertZBins);

//const double gw_fcalBins[11]   = {0, 63.719, 144.14, 289.595, 525.092, 875.41, 1368.75, 2046.51, 2989.31, 3618.44, 5200};
//const int gw_nFCalBins = sizeof (gw_fcalBins) / sizeof (gw_fcalBins[0]) - 1;
//const int gw_nQ2Bins = 10;
//const double* gw_q2Bins = linspace (0, 1, gw_nQ2Bins);
//const int gw_nVertZBins = 18;
//const double gw_vertZBins[19] = {-200, -175, -150, -125, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 125, 150, 175, 200};

TH1D* h_PbPbFCalDist = nullptr;
TH1D* h_PbPbFCal_weights = nullptr;
TH1D* h_PbPbVZDist = nullptr;
TH1D* h_PbPbVZ_weights = nullptr;
TH1D* h_PbPbQ2Dist = nullptr;
TH1D* h_PbPbQ2_weights = nullptr;
//TH3D* h_PbPbEventReweights = nullptr;

TH1D* h_ppVZDist = nullptr;
TH1D* h_ppVZ_weights = nullptr;
//TH1D* h_ppEventReweights = nullptr;


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

  TFile* inFile = new TFile (Form ("%s/Nominal/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  h_PbPbFCalDist = new TH1D (Form ("h_PbPbFCalDist_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins);
  h_PbPbFCal_weights = new TH1D (Form ("h_PbPbFCal_weights_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins);
  h_PbPbQ2Dist = new TH1D (Form ("h_PbPbQ2Dist_%s", name.Data ()), "", gw_nQ2Bins, gw_q2Bins);
  h_PbPbQ2_weights = new TH1D (Form ("h_PbPbQ2_weights_%s", name.Data ()), "", gw_nQ2Bins, gw_q2Bins);
  h_PbPbVZDist = new TH1D (Form ("h_PbPbVZDist_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  h_PbPbVZ_weights = new TH1D (Form ("h_PbPbVZ_weights_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  h_PbPbFCalDist->Sumw2 ();
  h_PbPbFCal_weights->Sumw2 ();
  h_PbPbQ2Dist->Sumw2 ();
  h_PbPbQ2_weights->Sumw2 ();
  h_PbPbVZDist->Sumw2 ();
  h_PbPbVZ_weights->Sumw2 ();
  //h_PbPbEventReweights = new TH3D (Form ("h_PbPbEventReweights_%s", name.Data ()), "", gw_nFCalBins, gw_fcalBins, gw_nQ2Bins, gw_q2Bins, gw_nVertZBins, gw_vertZBins);
  //h_PbPbEventReweights->Sumw2 ();

  h_ppVZDist = new TH1D (Form ("h_ppVZDist_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  h_ppVZ_weights = new TH1D (Form ("h_ppVZ_weights_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  h_ppVZDist->Sumw2 ();
  h_ppVZ_weights->Sumw2 ();
  //h_ppEventReweights = new TH1D (Form ("h_ppEventReweights_%s", name.Data ()), "", gw_nVertZBins, gw_vertZBins);
  //h_ppEventReweights->Sumw2 ();


  float event_weight = 0, fcal_et = 0, q2 = 0, vz = 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("fcalA_et",      &fcalA_et);
    PbPbTree->SetBranchAddress ("fcalC_et",      &fcalC_et);
    PbPbTree->SetBranchAddress ("fcalA_et_Cos",  &fcalA_et_Cos);
    PbPbTree->SetBranchAddress ("fcalC_et_Cos",  &fcalC_et_Cos);
    PbPbTree->SetBranchAddress ("fcalA_et_Sin",  &fcalA_et_Sin);
    PbPbTree->SetBranchAddress ("fcalC_et_Sin",  &fcalC_et_Sin);
    PbPbTree->SetBranchAddress ("vz",            &vz);

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
      h_PbPbFCalDist->Scale (1./h_PbPbFCalDist->Integral ());
    if (h_PbPbQ2Dist->Integral () > 0)
      h_PbPbQ2Dist->Scale (1./h_PbPbQ2Dist->Integral ());
    if (h_PbPbVZDist->Integral () > 0)
      h_PbPbVZDist->Scale (1./h_PbPbVZDist->Integral ());

    if (name == "data") {
      for (int ix = 1; ix <= h_PbPbFCal_weights->GetNbinsX (); ix++) {
        h_PbPbFCal_weights->SetBinContent (ix, 1);
      }
      for (int ix = 1; ix <= h_PbPbQ2_weights->GetNbinsX (); ix++) {
        h_PbPbQ2_weights->SetBinContent (ix, 1);
      }
      for (int ix = 1; ix <= h_PbPbVZ_weights->GetNbinsX (); ix++) {
        h_PbPbVZ_weights->SetBinContent (ix, 1);
      }
      //for (int ix = 1; ix <= h_PbPbEventReweights->GetNbinsX (); ix++)
      //  for (int iy = 1; iy <= h_PbPbEventReweights->GetNbinsY (); iy++)
      //    for (int iz = 1; iz <= h_PbPbEventReweights->GetNbinsZ (); iz++)
      //      h_PbPbEventReweights->SetBinContent (ix, iy, iz, 1); // should be 1 in the reference
    }
    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceFCalDist = (TH1D*)ztrackFile->Get ("h_PbPbFCalDist_data");
      TH1D* referenceVZDist = (TH1D*)ztrackFile->Get ("h_PbPbVZDist_data");
      TH1D* referenceQ2Dist = (TH1D*)ztrackFile->Get ("h_PbPbQ2Dist_data");

      for (int ix = 1; ix <= h_PbPbFCal_weights->GetNbinsX (); ix++) {
        const double fcal_weight = (h_PbPbFCalDist->GetBinContent (ix) != 0 ? referenceFCalDist->GetBinContent (ix) / h_PbPbFCalDist->GetBinContent (ix) : 0);
        h_PbPbFCal_weights->SetBinContent (ix, fcal_weight);
      }
      for (int ix = 1; ix <= h_PbPbQ2_weights->GetNbinsX (); ix++) {
        const double q2_weight = (h_PbPbQ2Dist->GetBinContent (ix) != 0 ? referenceQ2Dist->GetBinContent (ix) / h_PbPbQ2Dist->GetBinContent (ix) : 0);
        h_PbPbQ2_weights->SetBinContent (ix, q2_weight);
      }
      for (int ix = 1; ix <= h_PbPbVZ_weights->GetNbinsX (); ix++) {
        const double vz_weight = (h_PbPbVZDist->GetBinContent (ix) != 0 ? referenceVZDist->GetBinContent (ix) / h_PbPbVZDist->GetBinContent (ix) : 0);
        h_PbPbVZ_weights->SetBinContent (ix, vz_weight);
      }
      //for (int ix = 1; ix <= h_PbPbEventReweights->GetNbinsX (); ix++) {
      //  const double fcalReweight = (h_PbPbFCalDist->GetBinContent (ix) != 0 ? referenceFCalDist->GetBinContent (ix) / h_PbPbFCalDist->GetBinContent (ix) : 0);
      //  for (int iy = 1; iy <= h_PbPbEventReweights->GetNbinsY (); iy++) {
      //    const double q2Reweight = (h_PbPbQ2Dist->GetBinContent (iy) != 0 ? referenceQ2Dist->GetBinContent (iy) / h_PbPbQ2Dist->GetBinContent (iy) : 0);
      //    for (int iz = 1; iz <= h_PbPbEventReweights->GetNbinsZ (); iz++) {
      //      const double vzReweight = (h_PbPbVZDist->GetBinContent (iz) != 0 ? referenceVZDist->GetBinContent (iz) / h_PbPbVZDist->GetBinContent (iz) : 0);
      //      h_PbPbEventReweights->SetBinContent (ix, iy, iz, fcalReweight * vzReweight * q2Reweight);
      //    }
      //  }
      //}
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
      for (int ix = 1; ix <= h_ppVZ_weights->GetNbinsX (); ix++) {
        h_ppVZ_weights->SetBinContent (ix, 1);
      }
      //for (int ix = 1; ix <= h_ppEventReweights->GetNbinsX (); ix++)
      //  h_ppEventReweights->SetBinContent (ix, 1);
    }
    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceVZDist = (TH1D*)ztrackFile->Get ("h_ppVZDist_data");
      for (int ix = 1; ix <= h_ppVZ_weights->GetNbinsX (); ix++) {
        const double vz_weight = (h_ppVZDist->GetBinContent (ix) != 0. ? referenceVZDist->GetBinContent (ix) / h_ppVZDist->GetBinContent (ix) : 0);
        h_ppVZ_weights->SetBinContent (ix, vz_weight);
      }
      //for (int ix = 1; ix <= h_ppEventReweights->GetNbinsX (); ix++)
      //  h_ppEventReweights->SetBinContent (ix, h_ppVZDist->GetBinContent (ix) != 0. ? referenceVZDist->GetBinContent (ix) / h_ppVZDist->GetBinContent (ix) : 0);
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
  SafeWrite (h_PbPbFCal_weights);
  SafeWrite (h_PbPbVZ_weights);
  SafeWrite (h_PbPbQ2_weights);
  //SafeWrite (h_PbPbEventReweights);
  SafeWrite (h_ppVZDist);
  SafeWrite (h_ppVZ_weights);
  //SafeWrite (h_ppEventReweights);

  eventWeightsFile->Close ();
}

#endif
