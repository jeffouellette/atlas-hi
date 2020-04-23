#ifndef __GenerateWeights_C__
#define __GenerateWeights_C__

#include "Params.h"
#include "Utilities.h"

#include <ArrayTemplates.h>

#include <TEfficiency.h>
#include <TTree.h>

#include <iostream>

void SafeWrite (TObject* tobj) {
  if (tobj)
    tobj->Write ();
}

//const int gw_nFCalBins = 40;
//const double* gw_fcalBins = linspace (0, 5000, gw_nFCalBins);
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

TH1D* referenceFCalDist;
TH1D* referenceQ2Dist[numFineCentBins];
TH1D* referencePsi2Dist[numFineCentBins];
TH1D* referenceNchDist;

TH1D* h_PbPbFCalDist;
TH1D* h_PbPbFCal_weights;
TH1D* h_PbPbQ2Dist[numFineCentBins];
TH1D* h_PbPbQ2_weights[numFineCentBins];
TH1D* h_PbPbPsi2Dist[numFineCentBins];
TH1D* h_PbPbPsi2_weights[numFineCentBins];

TH1D* h_PbPbNchDist = nullptr;
TH1D* h_PbPbNch_weights = nullptr;

TH1D* h_ppNchDist = nullptr;
TH1D* h_ppNch_weights = nullptr;


void GenerateWeights (const char* name, const char* inFileName = "outFile.root", const char* outFileName = "eventWeightsFile.root") {

  TFile* inFile = new TFile (inFileName, "read");

  TTree* PbPbTree = (TTree*)inFile->Get (strcmp (name, "minbias") != 0 ? "PbPbZTrackTree" : "PbPbMixedTree");
  TTree* ppTree = (TTree*)inFile->Get (strcmp (name, "minbias") != 0 ? "ppZTrackTree" : "ppMixedTree");

  h_PbPbFCalDist = new TH1D (Form ("h_PbPbFCalDist_%s", name), "", numSuperFineCentBins-1, superFineCentBins);
  h_PbPbFCal_weights = new TH1D (Form ("h_PbPbFCal_weights_%s", name), "", numSuperFineCentBins-1, superFineCentBins);
  h_PbPbFCalDist->Sumw2 ();
  h_PbPbFCal_weights->Sumw2 ();

  for (int iCent = 0; iCent < numFineCentBins; iCent++) {
    h_PbPbQ2Dist[iCent] = new TH1D (Form ("h_PbPbQ2Dist_iCent%i_%s", iCent, name), "", gw_nQ2Bins, gw_q2Bins);
    h_PbPbQ2_weights[iCent] = new TH1D (Form ("h_PbPbQ2_weights_iCent%i_%s", iCent, name), "", gw_nQ2Bins, gw_q2Bins);
    h_PbPbQ2Dist[iCent]->Sumw2 ();
    h_PbPbQ2_weights[iCent]->Sumw2 ();

    h_PbPbPsi2Dist[iCent] = new TH1D (Form ("h_PbPbPsi2Dist_iCent%i_%s", iCent, name), "", gw_nPsi2Bins, gw_psi2Bins);
    h_PbPbPsi2_weights[iCent] = new TH1D (Form ("h_PbPbPsi2_weights_iCent%i_%s", iCent, name), "", gw_nPsi2Bins, gw_psi2Bins);
    h_PbPbPsi2Dist[iCent]->Sumw2 ();
    h_PbPbPsi2_weights[iCent]->Sumw2 ();
  }
  

  h_PbPbNchDist = new TH1D (Form ("h_PbPbNchDist_%s", name), "", 80, -0.5, 3999.5);
  h_PbPbNchDist->Sumw2 ();
  h_PbPbNch_weights = new TH1D (Form ("h_PbPbNch_weights_%s", name), "", 80, -0.5, 3999.5);
  h_PbPbNch_weights->Sumw2 ();

  h_ppNchDist = new TH1D (Form ("h_ppNchDist_%s", name), "", gw_nNchBins, gw_nchBins);
  h_ppNchDist->Sumw2 ();
  h_ppNch_weights = new TH1D (Form ("h_ppNch_weights_%s", name), "", gw_nNchBins, gw_nchBins);
  h_ppNch_weights->Sumw2 ();


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    int ntrk = 0;
    float b_event_weight = 1, event_weight = 1, fcal_et = 0, q2 = 0, psi2 = 0; 

    PbPbTree->SetBranchAddress ("event_weight", &b_event_weight);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("ntrk_all",     &ntrk);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      short iCent = 0;
      while (iCent < numFineCentBins) {
        if (fcal_et < fineCentBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numFineCentBins-1)
        continue;

      h_PbPbFCalDist->Fill (fcal_et, b_event_weight);
      h_PbPbNchDist->Fill (ntrk, event_weight);

      if (strcmp (name, "data") == 0) {
        h_PbPbQ2Dist[iCent]->Fill (q2, b_event_weight);
        h_PbPbPsi2Dist[iCent]->Fill (psi2, b_event_weight);
        h_PbPbQ2Dist[iCent]->Fill (q2, b_event_weight);
        h_PbPbPsi2Dist[iCent]->Fill (psi2, b_event_weight);
      }
    }
    cout << "Done 1st Pb+Pb loop." << endl;



    // Normalize histograms
    if (h_PbPbFCalDist->Integral () > 0)
      h_PbPbFCalDist->Scale (1./h_PbPbFCalDist->Integral ());

    for (int iCent = 0; iCent < numFineCentBins; iCent++) {
      if (h_PbPbQ2Dist[iCent]->Integral () > 0)
        h_PbPbQ2Dist[iCent]->Scale (1./h_PbPbQ2Dist[iCent]->Integral ());
      if (h_PbPbPsi2Dist[iCent]->Integral () > 0)
        h_PbPbPsi2Dist[iCent]->Scale (1./h_PbPbPsi2Dist[iCent]->Integral ());
    }
    if (h_PbPbNchDist->Integral () > 0)
      h_PbPbNchDist->Scale (1./h_PbPbNchDist->Integral ());


    // Set weights to 1 in data
    if (strcmp (name, "data") == 0) {
      for (int ix = 1; ix <= h_PbPbFCal_weights->GetNbinsX (); ix++)
        h_PbPbFCal_weights->SetBinContent (ix, 1);
      for (int iCent = 0; iCent < numFineCentBins; iCent++) {
        for (int ix = 1; ix <= h_PbPbQ2_weights[iCent]->GetNbinsX (); ix++) {
          h_PbPbQ2_weights[iCent]->SetBinContent (ix, 1);
          h_PbPbPsi2_weights[iCent]->SetBinContent (ix, 1);
        }
      }

      for (int ix = 1; ix <= h_PbPbNch_weights->GetNbinsX (); ix++)
        h_PbPbNch_weights->SetBinContent (ix, 1);
    }
    // Calculate reweighting factors otherwise (requires obtaining raw distributions in data)
    else {
      // get reference (data) distributions
      TFile* ztrackFile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/eventWeightsFile.root", "read");
      referenceNchDist = (TH1D*)ztrackFile->Get ("h_PbPbNchDist_data");

      referenceFCalDist = (TH1D*) ztrackFile->Get ("h_PbPbFCalDist_data");
      for (int iCent = 0; iCent < numFineCentBins; iCent++) {
        referenceQ2Dist[iCent] = (TH1D*)ztrackFile->Get (Form ("h_PbPbQ2Dist_iCent%i_data", iCent));
        referencePsi2Dist[iCent] = (TH1D*)ztrackFile->Get (Form ("h_PbPbPsi2Dist_iCent%i_data", iCent));
      }

      // set FCal distribution weights, then calculate potential additional weights
      for (int ix = 1; ix <= h_PbPbFCal_weights->GetNbinsX (); ix++) {
        const double fcal_weight = (h_PbPbFCalDist->GetBinContent (ix) != 0 ? referenceFCalDist->GetBinContent (ix) / h_PbPbFCalDist->GetBinContent (ix) : 0);
        h_PbPbFCal_weights->SetBinContent (ix, fcal_weight);
      }
      for (int ix = 1; ix <= h_PbPbNch_weights->GetNbinsX (); ix++) {
        const double fcal_weight = (h_PbPbNchDist->GetBinContent (ix) != 0 ? referenceNchDist->GetBinContent (ix) / h_PbPbNchDist->GetBinContent (ix) : 0);
        h_PbPbNch_weights->SetBinContent (ix, fcal_weight);
      }

      // restore other histograms (so we don't fill it twice)
      for (int iCent = 0; iCent < numFineCentBins; iCent++) {
        for (int ix = 1; ix <= h_PbPbQ2Dist[iCent]->GetNbinsX (); ix++) {
          h_PbPbQ2Dist[iCent]->SetBinContent (ix, 0);
          h_PbPbQ2Dist[iCent]->SetBinError (ix, 0);
        }
        for (int ix = 1; ix <= h_PbPbPsi2Dist[iCent]->GetNbinsX (); ix++) {
          h_PbPbPsi2Dist[iCent]->SetBinContent (ix, 0);
          h_PbPbPsi2Dist[iCent]->SetBinError (ix, 0);
        }
      }
      

      for (int iEvt = 0; iEvt < nEvts; iEvt++) {
        if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
          cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
        PbPbTree->GetEntry (iEvt);

        short iCent = 0;
        while (iCent < numFineCentBins) {
          if (fcal_et < fineCentBins[iCent])
            break;
          else
            iCent++;
        }
        if (iCent < 1 || iCent > numFineCentBins-1)
          continue;

        //event_weight = b_event_weight * h_PbPbFCal_weights->GetBinContent (h_PbPbFCal_weights->FindBin (fcal_et));
        event_weight = b_event_weight;
        h_PbPbQ2Dist[iCent]->Fill (q2, event_weight);
        h_PbPbPsi2Dist[iCent]->Fill (psi2, event_weight);
      }
      cout << "Done 2nd Pb+Pb loop." << endl;

      // now set other weighting histograms
      for (int iCent = 0; iCent < numFineCentBins; iCent++) {
        if (h_PbPbQ2Dist[iCent]->Integral () > 0)
          h_PbPbQ2Dist[iCent]->Scale (1./h_PbPbQ2Dist[iCent]->Integral ());
        if (h_PbPbPsi2Dist[iCent]->Integral () > 0)
          h_PbPbPsi2Dist[iCent]->Scale (1./h_PbPbPsi2Dist[iCent]->Integral ());
      }
      
      for (int iCent = 0; iCent < numFineCentBins; iCent++) {
        for (int ix = 1; ix <= h_PbPbQ2_weights[iCent]->GetNbinsX (); ix++) {
          const double q2_weight = (h_PbPbQ2Dist[iCent]->GetBinContent (ix) != 0 ? referenceQ2Dist[iCent]->GetBinContent (ix) / h_PbPbQ2Dist[iCent]->GetBinContent (ix) : 0);
          const double q2_weight_err = (h_PbPbQ2Dist[iCent]->GetBinContent (ix) != 0 ? q2_weight * sqrt (pow (h_PbPbQ2Dist[iCent]->GetBinError (ix) / h_PbPbQ2Dist[iCent]->GetBinContent (ix), 2) + pow (referenceQ2Dist[iCent]->GetBinError (ix) / referenceQ2Dist[iCent]->GetBinContent (ix), 2)) : 0);
          h_PbPbQ2_weights[iCent]->SetBinContent (ix, q2_weight);
          h_PbPbQ2_weights[iCent]->SetBinError (ix, q2_weight_err);
        }
        for (int ix = 1; ix <= h_PbPbPsi2_weights[iCent]->GetNbinsX (); ix++) {
          const double psi2_weight = (h_PbPbPsi2Dist[iCent]->GetBinContent (ix) != 0 ? referencePsi2Dist[iCent]->GetBinContent (ix) / h_PbPbPsi2Dist[iCent]->GetBinContent (ix) : 0);
          const double psi2_weight_err = (h_PbPbPsi2Dist[iCent]->GetBinContent (ix) != 0 ? psi2_weight * sqrt (pow (h_PbPbPsi2Dist[iCent]->GetBinError (ix) / h_PbPbPsi2Dist[iCent]->GetBinContent (ix), 2) + pow (referencePsi2Dist[iCent]->GetBinError (ix) / referencePsi2Dist[iCent]->GetBinContent (ix), 2)) : 0);
          h_PbPbPsi2_weights[iCent]->SetBinContent (ix, psi2_weight);
          h_PbPbPsi2_weights[iCent]->SetBinError (ix, psi2_weight_err);
        }
      }
      
      ztrackFile->Close ();
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

    if (strcmp (name, "data") == 0) {
      for (int ix = 1; ix <= h_ppNch_weights->GetNbinsX (); ix++) {
        h_ppNch_weights->SetBinContent (ix, 1);
      }
    } else {
      TFile* ztrackFile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/eventWeightsFile.root", "read");
      referenceNchDist = (TH1D*)ztrackFile->Get ("h_ppNchDist_data");

      for (int ix = 1; ix <= h_ppNch_weights->GetNbinsX (); ix++) {
        const double nch_weight = (h_ppNchDist->GetBinContent (ix) != 0 ? referenceNchDist->GetBinContent (ix) / h_ppNchDist->GetBinContent (ix) : 0);
        h_ppNch_weights->SetBinContent (ix, nch_weight);
      }

      ztrackFile->Close ();
    }
  }


  TFile* eventWeightsFile = new TFile (outFileName, "recreate");

  SafeWrite (h_PbPbFCalDist);
  SafeWrite (h_PbPbFCal_weights);

  for (int iCent = 0; iCent < numFineCentBins; iCent++) {
    SafeWrite (h_PbPbQ2Dist[iCent]);
    SafeWrite (h_PbPbQ2_weights[iCent]);

    SafeWrite (h_PbPbPsi2Dist[iCent]);
    SafeWrite (h_PbPbPsi2_weights[iCent]);
  }
  

  SafeWrite (h_PbPbNchDist);
  SafeWrite (h_PbPbNch_weights);

  SafeWrite (h_ppNchDist);
  SafeWrite (h_ppNch_weights);

  eventWeightsFile->Close ();
}


void GenerateDataWeights () {
  GenerateWeights ("data", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/data18hi.root", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/eventWeightsFile.root");
}

void GeneratePbPbMCWeights () {
  GenerateWeights ("mc", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/data/mc_100/PbPbEventWeightsTree.root", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/data/mc_100/PbPbEventWeightsFile.root");
}

void GenerateppMCWeights () {
  GenerateWeights ("mc", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/data/mc_100/ppEventWeightsTree.root", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/data/mc_100/ppEventWeightsFile.root");
}

void GenerateHijingWeights () {
  GenerateWeights ("hijing", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/data/hijing_200/PbPbEventWeightsTree.root", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/data/hijing_200/PbPbEventWeightsFile.root");
}

//void Generate2015HijingWeights () {
//  GenerateWeights ("hijing", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/data/hijing_110/hijing2015_eventWeightsTree.root", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/data/hijing_110/hijing2015_eventWeightsFile.root");
//}

//void Generate2018HijingWeights () {
//  GenerateWeights ("hijing", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/data/hijing_110/hijing2018_eventWeightsTree.root", "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/data/hijing_110/hijing2018_eventWeightsFile.root");
//}

#endif
