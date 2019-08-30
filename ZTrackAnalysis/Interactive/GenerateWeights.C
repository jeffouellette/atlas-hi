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

const int gw_nFCalBins = 40;
const double* gw_fcalBins = linspace (0, 5000, gw_nFCalBins);
const int gw_nQ2Bins = 20;
const double* gw_q2Bins = linspace (0, 0.3, gw_nQ2Bins);
const int gw_nPsi2Bins = 20;
const double* gw_psi2Bins = linspace (-pi/2, pi/2, gw_nQ2Bins);

const int gw_nNchBins = 160;
const double* gw_nchBins = linspace (-0.5, 160.5, gw_nNchBins);

//const double gw_fcalBins[11]   = {0, 63.719, 144.14, 289.595, 525.092, 875.41, 1368.75, 2046.51, 2989.31, 3618.44, 5200};
//const int gw_nFCalBins = sizeof (gw_fcalBins) / sizeof (gw_fcalBins[0]) - 1;
//const int gw_nQ2Bins = 10;
//const double* gw_q2Bins = linspace (0, 1, gw_nQ2Bins);

TH1D* referenceFCalDist[nPtZBins+1];
TH1D* referenceQ2Dist[numFinerCentBins][nPtZBins+1];
TH1D* referencePsi2Dist[numFinerCentBins][nPtZBins+1];

TH1D* h_PbPbFCalDist[nPtZBins+1];
TH1D* h_PbPbFCal_weights[nPtZBins+1];
TH1D* h_PbPbQ2Dist[numFinerCentBins][nPtZBins+1];
TH1D* h_PbPbQ2_weights[numFinerCentBins][nPtZBins+1];
TH1D* h_PbPbPsi2Dist[numFinerCentBins][nPtZBins+1];
TH1D* h_PbPbPsi2_weights[numFinerCentBins][nPtZBins+1];

TH1D* h_PbPbNchDist = nullptr;
TH1D* h_PbPbNch_weights = nullptr;

TH1D* h_ppNchDist = nullptr;
TH1D* h_ppNch_weights = nullptr;


void GenerateWeights (const TString name, const TString inFileName = "outFile.root", const TString outFileName = "eventWeightsFile.root") {

  if (name == "data")
    SetupDirectories ("DataAnalysis/Nominal/", "ZTrackAnalysis/");
  else if (name == "mc")
    SetupDirectories ("MCAnalysis/Nominal/", "ZTrackAnalysis/");
  else if (name == "minbias")
    SetupDirectories ("MinbiasAnalysis/Nominal/", "ZTrackAnalysis/");
  else if (name == "hijing")
    SetupDirectories ("MinbiasAnalysis/Nominal/", "ZTrackAnalysis/");
  else if (name == "truth")
    SetupDirectories ("TruthAnalysis/Nominal/", "ZTrackAnalysis/");
  else if (name == "minbias_runvar")
    SetupDirectories ("MinbiasAnalysis/Variations/RunVariation/", "ZTrackAnalysis/");
  else {
    cout << "Unrecognized analysis name! Quitting." << endl;
    return;
  }

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  for (int iPtZ = 0; iPtZ < nPtZBins+1; iPtZ++) {
    h_PbPbFCalDist[iPtZ] = new TH1D (Form ("h_PbPbFCalDist_iPtZ%i_%s", iPtZ, name.Data ()), "", gw_nFCalBins, gw_fcalBins);
    h_PbPbFCal_weights[iPtZ] = new TH1D (Form ("h_PbPbFCal_weights_iPtZ%i_%s", iPtZ, name.Data ()), "", gw_nFCalBins, gw_fcalBins);
    h_PbPbFCalDist[iPtZ]->Sumw2 ();
    h_PbPbFCal_weights[iPtZ]->Sumw2 ();

    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
      h_PbPbQ2Dist[iCent][iPtZ] = new TH1D (Form ("h_PbPbQ2Dist_iCent%i_iPtZ%i_%s", iCent, iPtZ, name.Data ()), "", gw_nQ2Bins, gw_q2Bins);
      h_PbPbQ2_weights[iCent][iPtZ] = new TH1D (Form ("h_PbPbQ2_weights_iCent%i_iPtZ%i_%s", iCent, iPtZ, name.Data ()), "", gw_nQ2Bins, gw_q2Bins);
      h_PbPbQ2Dist[iCent][iPtZ]->Sumw2 ();
      h_PbPbQ2_weights[iCent][iPtZ]->Sumw2 ();

      h_PbPbPsi2Dist[iCent][iPtZ] = new TH1D (Form ("h_PbPbPsi2Dist_iCent%i_iPtZ%i_%s", iCent, iPtZ, name.Data ()), "", gw_nPsi2Bins, gw_psi2Bins);
      h_PbPbPsi2_weights[iCent][iPtZ] = new TH1D (Form ("h_PbPbPsi2_weights_iCent%i_iPtZ%i_%s", iCent, iPtZ, name.Data ()), "", gw_nPsi2Bins, gw_psi2Bins);
      h_PbPbPsi2Dist[iCent][iPtZ]->Sumw2 ();
      h_PbPbPsi2_weights[iCent][iPtZ]->Sumw2 ();
    }
  }

  h_PbPbNchDist = new TH1D (Form ("h_PbPbNchDist_%s", name.Data ()), "", 80, -0.5, 3999.5);
  h_PbPbNchDist->Sumw2 ();
  h_PbPbNch_weights = new TH1D (Form ("h_PbPbNch_weights_%s", name.Data ()), "", 80, -0.5, 3999.5);
  h_PbPbNch_weights->Sumw2 ();

  h_ppNchDist = new TH1D (Form ("h_ppNchDist_%s", name.Data ()), "", gw_nNchBins, gw_nchBins);
  h_ppNchDist->Sumw2 ();
  h_ppNch_weights = new TH1D (Form ("h_ppNch_weights_%s", name.Data ()), "", gw_nNchBins, gw_nchBins);
  h_ppNch_weights->Sumw2 ();


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    int ntrk = 0;
    float b_event_weight = 1, event_weight = 1, fcal_et = 0, q2 = 0, psi2 = 0, z_pt = 0;

    if (name != "minbias" && name != "minbias_runvar")
      PbPbTree->SetBranchAddress ("event_weight", &b_event_weight);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("ntrk_all",     &ntrk);
    if (name == "data" || name == "mc")
      PbPbTree->SetBranchAddress ("z_pt",       &z_pt);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      short iCent = 0;
      while (iCent < numFinerCentBins) {
        if (fcal_et < finerCentBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numFinerCentBins-1)
        continue;

      short iPtZ = 0;
      if (name == "data" || name == "mc") {
        while (iPtZ < nPtZBins) {
          if (z_pt < zPtBins[iPtZ+1])
            break;
          else
            iPtZ++;
        }
        h_PbPbFCalDist[iPtZ]->Fill (fcal_et, b_event_weight);
      }
      else {
        for (iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          h_PbPbFCalDist[iPtZ]->Fill (fcal_et, b_event_weight);
        }
      }
      h_PbPbNchDist->Fill (ntrk, event_weight);

      h_PbPbFCalDist[nPtZBins]->Fill (fcal_et, b_event_weight);

      if (name == "data") {
        h_PbPbQ2Dist[iCent][iPtZ]->Fill (q2, b_event_weight);
        h_PbPbPsi2Dist[iCent][iPtZ]->Fill (psi2, b_event_weight);
        h_PbPbQ2Dist[iCent][nPtZBins]->Fill (q2, b_event_weight);
        h_PbPbPsi2Dist[iCent][nPtZBins]->Fill (psi2, b_event_weight);
      }
    }
    cout << "Done 1st Pb+Pb loop." << endl;



    // Normalize histograms
    for (int iPtZ = 0; iPtZ < nPtZBins+1; iPtZ++) {
      if (h_PbPbFCalDist[iPtZ]->Integral () > 0)
        h_PbPbFCalDist[iPtZ]->Scale (1./h_PbPbFCalDist[iPtZ]->Integral ());

      for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
        if (h_PbPbQ2Dist[iCent][iPtZ]->Integral () > 0)
          h_PbPbQ2Dist[iCent][iPtZ]->Scale (1./h_PbPbQ2Dist[iCent][iPtZ]->Integral ());
        if (h_PbPbPsi2Dist[iCent][iPtZ]->Integral () > 0)
          h_PbPbPsi2Dist[iCent][iPtZ]->Scale (1./h_PbPbPsi2Dist[iCent][iPtZ]->Integral ());
      }
    }
    if (h_PbPbNchDist->Integral () > 0)
      h_PbPbNchDist->Scale (1./h_PbPbNchDist->Integral ());


    // Set weights to 1 in data
    if (name == "data") {
      for (int iPtZ = 0; iPtZ < nPtZBins+1; iPtZ++) {
        for (int ix = 1; ix <= h_PbPbFCal_weights[iPtZ]->GetNbinsX (); ix++)
          h_PbPbFCal_weights[iPtZ]->SetBinContent (ix, 1);
        for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
          for (int ix = 1; ix <= h_PbPbQ2_weights[iCent][iPtZ]->GetNbinsX (); ix++) {
            h_PbPbQ2_weights[iCent][iPtZ]->SetBinContent (ix, 1);
            h_PbPbPsi2_weights[iCent][iPtZ]->SetBinContent (ix, 1);
          }
        }
      }

      for (int ix = 1; ix <= h_PbPbNch_weights->GetNbinsX (); ix++)
        h_PbPbNch_weights->SetBinContent (ix, 1);
    }
    // Calculate reweighting factors otherwise (requires obtaining raw distributions in data)
    else {
      // get reference (data) distributions
      SetupDirectories ("DataAnalysis/Nominal/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceNchDist = (TH1D*)ztrackFile->Get ("h_PbPbNchDist_data");

      //TH1D* referenceFCalDist[nPtZBins];
      //TH1D* referencePsi2Dist[numFinerCentBins][nPtZBins];
      //TH1D* referenceQ2Dist[numFinerCentBins][nPtZBins];
      for (int iPtZ = 0; iPtZ < nPtZBins+1; iPtZ++) {
        referenceFCalDist[iPtZ] = (TH1D*) ztrackFile->Get (Form ("h_PbPbFCalDist_iPtZ%i_data", iPtZ));
        for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
          referenceQ2Dist[iCent][iPtZ] = (TH1D*)ztrackFile->Get (Form ("h_PbPbQ2Dist_iCent%i_iPtZ%i_data", iCent, iPtZ));
          referencePsi2Dist[iCent][iPtZ] = (TH1D*)ztrackFile->Get (Form ("h_PbPbPsi2Dist_iCent%i_iPtZ%i_data", iCent, iPtZ));
        }
      }

      // set FCal distribution weights, then calculate potential additional weights
      for (int iPtZ = 0; iPtZ < nPtZBins+1; iPtZ++) {
        for (int ix = 1; ix <= h_PbPbFCal_weights[iPtZ]->GetNbinsX (); ix++) {
          const double fcal_weight = (h_PbPbFCalDist[iPtZ]->GetBinContent (ix) != 0 ? referenceFCalDist[iPtZ]->GetBinContent (ix) / h_PbPbFCalDist[iPtZ]->GetBinContent (ix) : 0);
          h_PbPbFCal_weights[iPtZ]->SetBinContent (ix, fcal_weight);
        }
      }
      for (int ix = 1; ix <= h_PbPbNch_weights->GetNbinsX (); ix++) {
        const double fcal_weight = (h_PbPbNchDist->GetBinContent (ix) != 0 ? referenceNchDist->GetBinContent (ix) / h_PbPbNchDist->GetBinContent (ix) : 0);
        h_PbPbNch_weights->SetBinContent (ix, fcal_weight);
      }

      // restore other histograms (so we don't fill it twice)
      for (int iPtZ = 0; iPtZ < nPtZBins+1; iPtZ++) {
        for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
          for (int ix = 1; ix <= h_PbPbQ2Dist[iCent][iPtZ]->GetNbinsX (); ix++) {
            h_PbPbQ2Dist[iCent][iPtZ]->SetBinContent (ix, 0);
            h_PbPbQ2Dist[iCent][iPtZ]->SetBinError (ix, 0);
          }
          for (int ix = 1; ix <= h_PbPbPsi2Dist[iCent][iPtZ]->GetNbinsX (); ix++) {
            h_PbPbPsi2Dist[iCent][iPtZ]->SetBinContent (ix, 0);
            h_PbPbPsi2Dist[iCent][iPtZ]->SetBinError (ix, 0);
          }
        }
      }

      for (int iEvt = 0; iEvt < nEvts; iEvt++) {
        if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
          cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
        PbPbTree->GetEntry (iEvt);

        short iCent = 0;
        while (iCent < numFinerCentBins) {
          if (fcal_et < finerCentBins[iCent])
            break;
          else
            iCent++;
        }
        if (iCent < 1 || iCent > numFinerCentBins-1)
          continue;

        short iPtZ = 0;
        if (name == "mc") {
          while (iPtZ < nPtZBins) {
            if (z_pt < zPtBins[iPtZ+1])
              break;
            else
              iPtZ++;
          }
          event_weight = b_event_weight * h_PbPbFCal_weights[iPtZ]->GetBinContent (h_PbPbFCal_weights[iPtZ]->FindBin (fcal_et));
          h_PbPbQ2Dist[iCent][iPtZ]->Fill (q2, event_weight);
          h_PbPbPsi2Dist[iCent][iPtZ]->Fill (psi2, event_weight);
        }
        else {
          for (iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
            event_weight = b_event_weight * h_PbPbFCal_weights[iPtZ]->GetBinContent (h_PbPbFCal_weights[iPtZ]->FindBin (fcal_et));
            h_PbPbQ2Dist[iCent][iPtZ]->Fill (q2, event_weight);
            h_PbPbPsi2Dist[iCent][iPtZ]->Fill (psi2, event_weight);
          }
        }
        event_weight = b_event_weight * h_PbPbFCal_weights[nPtZBins]->GetBinContent (h_PbPbFCal_weights[nPtZBins]->FindBin (fcal_et));
        h_PbPbQ2Dist[iCent][nPtZBins]->Fill (q2, event_weight);
        h_PbPbPsi2Dist[iCent][nPtZBins]->Fill (psi2, event_weight);
      }
      cout << "Done 2nd Pb+Pb loop." << endl;

      // now set other weighting histograms
      for (int iPtZ = 0; iPtZ < nPtZBins+1; iPtZ++) {
        for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
          if (h_PbPbQ2Dist[iCent][iPtZ]->Integral () > 0)
            h_PbPbQ2Dist[iCent][iPtZ]->Scale (1./h_PbPbQ2Dist[iCent][iPtZ]->Integral ());
          if (h_PbPbPsi2Dist[iCent][iPtZ]->Integral () > 0)
            h_PbPbPsi2Dist[iCent][iPtZ]->Scale (1./h_PbPbPsi2Dist[iCent][iPtZ]->Integral ());
        }
      }
      for (int iPtZ = 0; iPtZ < nPtZBins+1; iPtZ++) {
        for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
          for (int ix = 1; ix <= h_PbPbQ2_weights[iCent][iPtZ]->GetNbinsX (); ix++) {
            const double q2_weight = (h_PbPbQ2Dist[iCent][iPtZ]->GetBinContent (ix) != 0 ? referenceQ2Dist[iCent][iPtZ]->GetBinContent (ix) / h_PbPbQ2Dist[iCent][iPtZ]->GetBinContent (ix) : 0);
            const double q2_weight_err = (h_PbPbQ2Dist[iCent][iPtZ]->GetBinContent (ix) != 0 ? q2_weight * sqrt (pow (h_PbPbQ2Dist[iCent][iPtZ]->GetBinError (ix) / h_PbPbQ2Dist[iCent][iPtZ]->GetBinContent (ix), 2) + pow (referenceQ2Dist[iCent][iPtZ]->GetBinError (ix) / referenceQ2Dist[iCent][iPtZ]->GetBinContent (ix), 2)) : 0);
            h_PbPbQ2_weights[iCent][iPtZ]->SetBinContent (ix, q2_weight);
            h_PbPbQ2_weights[iCent][iPtZ]->SetBinError (ix, q2_weight_err);
          }
          for (int ix = 1; ix <= h_PbPbPsi2_weights[iCent][iPtZ]->GetNbinsX (); ix++) {
            const double psi2_weight = (h_PbPbPsi2Dist[iCent][iPtZ]->GetBinContent (ix) != 0 ? referencePsi2Dist[iCent][iPtZ]->GetBinContent (ix) / h_PbPbPsi2Dist[iCent][iPtZ]->GetBinContent (ix) : 0);
            const double psi2_weight_err = (h_PbPbPsi2Dist[iCent][iPtZ]->GetBinContent (ix) != 0 ? psi2_weight * sqrt (pow (h_PbPbPsi2Dist[iCent][iPtZ]->GetBinError (ix) / h_PbPbPsi2Dist[iCent][iPtZ]->GetBinContent (ix), 2) + pow (referencePsi2Dist[iCent][iPtZ]->GetBinError (ix) / referencePsi2Dist[iCent][iPtZ]->GetBinContent (ix), 2)) : 0);
            h_PbPbPsi2_weights[iCent][iPtZ]->SetBinContent (ix, psi2_weight);
            h_PbPbPsi2_weights[iCent][iPtZ]->SetBinError (ix, psi2_weight_err);
          }
        }
      }
      ztrackFile->Close ();

      if (name == "data")
        SetupDirectories ("DataAnalysis/Nominal/", "ZTrackAnalysis/");
      else if (name == "mc")
        SetupDirectories ("MCAnalysis/Nominal/", "ZTrackAnalysis/");
      else if (name == "minbias")
        SetupDirectories ("MinbiasAnalysis/Nominal/", "ZTrackAnalysis/");
      else if (name == "hijing")
        SetupDirectories ("MinbiasAnalysis/Nominal/", "ZTrackAnalysis/");
      else if (name == "truth")
        SetupDirectories ("TruthAnalysis/Nominal/", "ZTrackAnalysis/");
      else if (name == "minbias_runvar")
        SetupDirectories ("MinbiasAnalysis/Variations/RunVariation/", "ZTrackAnalysis/");
    }
  }




  // Now calculate possible weights for pp collisions.
  if (ppTree) {
    float event_weight = 0;
    int ntrk = 0;
    ppTree->SetBranchAddress ("event_weight", &event_weight);
    ppTree->SetBranchAddress ("ntrk_all", &ntrk);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      h_ppNchDist->Fill (ntrk, event_weight);
    }
    cout << "Done pp loop." << endl;

    if (name == "data") {
      for (int ix = 1; ix <= h_ppNch_weights->GetNbinsX (); ix++) {
        h_ppNch_weights->SetBinContent (ix, 1);
      }
    } else {
      SetupDirectories ("DataAnalysis/Nominal/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceNchDist = (TH1D*)ztrackFile->Get ("h_ppNchDist_data");

      for (int ix = 1; ix <= h_ppNch_weights->GetNbinsX (); ix++) {
        const double nch_weight = (h_ppNchDist->GetBinContent (ix) != 0 ? referenceNchDist->GetBinContent (ix) / h_ppNchDist->GetBinContent (ix) : 0);
        h_ppNch_weights->SetBinContent (ix, nch_weight);
      }

      ztrackFile->Close ();

      if (name == "data")
        SetupDirectories ("DataAnalysis/Nominal/", "ZTrackAnalysis/");
      else if (name == "mc")
        SetupDirectories ("MCAnalysis/Nominal/", "ZTrackAnalysis/");
      else if (name == "minbias")
        SetupDirectories ("MinbiasAnalysis/Nominal/", "ZTrackAnalysis/");
      else if (name == "hijing")
        SetupDirectories ("MinbiasAnalysis/Nominal/", "ZTrackAnalysis/");
      else if (name == "truth")
        SetupDirectories ("TruthAnalysis/Nominal/", "ZTrackAnalysis/");
      else if (name == "minbias_runvar")
        SetupDirectories ("MinbiasAnalysis/Variations/RunVariation/", "ZTrackAnalysis/");
    }
  }


  TFile* eventWeightsFile = new TFile (Form ("%s/%s", rootPath.Data (), outFileName.Data ()), "recreate");

  for (int iPtZ = 0; iPtZ < nPtZBins+1; iPtZ++) {
    SafeWrite (h_PbPbFCalDist[iPtZ]);
    SafeWrite (h_PbPbFCal_weights[iPtZ]);

    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
      SafeWrite (h_PbPbQ2Dist[iCent][iPtZ]);
      SafeWrite (h_PbPbQ2_weights[iCent][iPtZ]);

      SafeWrite (h_PbPbPsi2Dist[iCent][iPtZ]);
      SafeWrite (h_PbPbPsi2_weights[iCent][iPtZ]);
    }
  }

  SafeWrite (h_PbPbNchDist);
  SafeWrite (h_PbPbNch_weights);

  SafeWrite (h_ppNchDist);
  SafeWrite (h_ppNch_weights);

  eventWeightsFile->Close ();
}


void GenerateDataWeights () {
  GenerateWeights ("data");
}

void GenerateMinbiasWeights () {
  GenerateWeights ("minbias", "eventWeightsTree.root");
}

void GenerateMinbiasRunVarWeights () {
  GenerateWeights ("minbias_runvar", "eventWeightsTree.root");
}

void GenerateMCWeights () {
  GenerateWeights ("mc", "eventWeightsTree.root");
}

void GenerateTruthWeights () {
  GenerateWeights ("truth", "eventWeightsTree.root");
}

void Generate2015HijingWeights () {
  GenerateWeights ("hijing", "PbPb_Hijing_15.root", "PbPb_Hijing_15_eventWeights.root");
}

void Generate2018HijingWeights () {
  GenerateWeights ("hijing", "PbPb_Hijing_18.root", "PbPb_Hijing_18_eventWeights.root");
}

#endif
