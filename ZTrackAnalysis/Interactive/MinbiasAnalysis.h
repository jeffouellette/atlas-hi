#ifndef __MinbiasAnalysis_h__
#define __MinbiasAnalysis_h__

#include "Params.h"
#include "FullAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class MinbiasAnalysis : public FullAnalysis {

  private:
  TFile* zMixFile = nullptr;

  TTree* LoadEventMixingTree (const char* _treeName);

  public:
  MinbiasAnalysis (const char* _name = "minbias", const char* subDir = "Nominal", const bool _useHITight = false) : FullAnalysis () {
    name = _name;
    directory = Form ("MinbiasAnalysis/%s/", subDir);
    plotFill = true;
    plotSignal = false;
    useAltMarker = false;
    backgroundSubtracted = true;
    useHITight = _useHITight;
    LoadTrackingEfficiencies ();
    LoadTrackingPurities ();

    SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
    TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
    h_PbPbFCal_weights = (TH1D*) eventWeightsFile->Get ("h_PbPbFCal_weights_minbias");
    for (short iCent = 0; iCent < numFinerCentBins; iCent++) {
      h_PbPbQ2_weights[iCent] = (TH1D*) eventWeightsFile->Get (Form ("h_PbPbQ2_weights_iCent%i_minbias", iCent));
    }
    h_ppNch_weights = (TH1D*) eventWeightsFile->Get ("h_ppNch_weights_minbias");

    SetupDirectories (directory, "ZTrackAnalysis/");
  }

  void Execute (const char* inFileName = "outFile.root", const char* outFileName = "savedHists.root") override;

  void CombineHists () override;
  void ScaleHists () override;

  void GenerateWeights ();
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Loads the file for mixing Z's into minimum bias events
////////////////////////////////////////////////////////////////////////////////////////////////
TTree* MinbiasAnalysis :: LoadEventMixingTree (const char* _treeName) {
  if (! (zMixFile && zMixFile->IsOpen ())) {
    SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
    zMixFile = new TFile (Form ("%s/Nominal/zMixFile.root", rootPath.Data ()), "read");
    SetupDirectories (directory, "ZTrackAnalysis/");
  }
  return (TTree*) zMixFile->Get (_treeName);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (NEW VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  SetupDirectories (directory, "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*) inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*) inFile->Get ("ppZTrackTree");

  CreateHists ();

  int ntrk = 0;
  float event_weight = 1, fcal_weight = 1, q2_weight = 1, vz_weight = 1, nch_weight = 1;
  float fcal_et = 0, q2 = 0, vz = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];

  float z_pt = 0, z_phi = 0;

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    TTree* zMixTree = LoadEventMixingTree ("PbPbZTree");
    const int nMixEvts = zMixTree->GetEntries ();
    zMixTree->SetBranchAddress ("z_pt",   &z_pt);
    zMixTree->SetBranchAddress ("z_phi",  &z_phi);

    PbPbTree->SetBranchAddress ("fcal_et",  &fcal_et);
    PbPbTree->SetBranchAddress ("q2",       &q2);
    PbPbTree->SetBranchAddress ("vz",       &vz);
    PbPbTree->SetBranchAddress ("ntrk",     &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",   trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",  trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",  trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      zMixTree->GetEntry (iEvt % nMixEvts);

      const short iSpc = 0;
      const short iPhi = 0;

      short iPtZ = 0;
      while (iPtZ < nPtZBins) {
        if (z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
      }

      short iCent = 0;
      while (iCent < numCentBins) {
        if (fcal_et < centBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      short iFinerCent = 0;
      while (iFinerCent < numFinerCentBins) {
        if (fcal_et < finerCentBins[iFinerCent])
          break;
        else
          iFinerCent++;
      }
      if (iFinerCent < 1 || iFinerCent > numFinerCentBins-1)
        continue;

      fcal_weight = h_PbPbFCal_weights->GetBinContent (h_PbPbFCal_weights->FindBin (fcal_et));
      //q2_weight = h_PbPbQ2_weights[iFinerCent]->GetBinContent (h_PbPbQ2_weights[iFinerCent]->FindBin (q2));

      event_weight = fcal_weight * q2_weight * vz_weight;
      //event_weight = h_PbPb_event_reweights->GetBinContent (h_PbPb_event_reweights->FindBin (fcal_et, q2, vz));

      if (event_weight == 0)
        continue;

      h_fcal_et->Fill (fcal_et);
      //h_fcal_et_q2->Fill (fcal_et, q2);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);
      //h_fcal_et_q2_reweighted->Fill (fcal_et, q2, event_weight);

      h_q2[iFinerCent]->Fill (q2);
      h_q2_reweighted[iFinerCent]->Fill (q2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkEff = GetTrackingEfficiency (fcal_et, trk_pt[iTrk], trk_eta[iTrk], true);

        if (trkEff == 0)
          continue;

        TH1D* h = h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent];
        h->Fill (trk_pt[iTrk], event_weight / trkEff);

        h = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
        h->Fill (trk_pt[iTrk], event_weight / trkEff);

        const int iPtTrk = h->FindFixBin (trk_pt[iTrk]) - 1;
        if (iPtTrk >= 0 && iPtTrk < nPtTrkBins) {
          h = h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent];
          for (int iPhi = 1; iPhi <= h->GetNbinsX (); iPhi++)
            h->Fill (h->GetBinCenter (iPhi), event_weight / trkEff / h->GetNbinsX ());
        }

        h = h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
        h->Fill (trk_pt[iTrk] / z_pt, event_weight / trkEff);
      }

    } // end loop over Pb+Pb tree
    cout << "Done minbias Pb+Pb loop." << endl;
  }


  if (ppTree) {
    TTree* zMixTree = LoadEventMixingTree ("ppZTree");
    const int nMixEvts = zMixTree->GetEntries ();
    zMixTree->SetBranchAddress ("z_pt",   &z_pt);
    zMixTree->SetBranchAddress ("z_phi",  &z_phi);

    ppTree->SetBranchAddress ("vz",         &vz);
    ppTree->SetBranchAddress ("ntrk",       &ntrk);
    ppTree->SetBranchAddress ("trk_pt",     trk_pt);
    ppTree->SetBranchAddress ("trk_eta",    trk_eta);
    ppTree->SetBranchAddress ("trk_phi",    trk_phi);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      zMixTree->GetEntry (iEvt % nMixEvts);

      const short iSpc = 0;
      const short iPhi = 0;
      const short iCent = 0; // iCent = 0 for pp

      short iPtZ = 0;
      while (iPtZ < nPtZBins) {
        if (z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
      }

      //nch_weight = h_ppNch_weights->GetBinContent (h_ppNch_weights->FindBin (ntrk));

      event_weight = nch_weight;

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkEff = GetTrackingEfficiency (fcal_et, trk_pt[iTrk], trk_eta[iTrk], false);

        if (trkEff == 0)
          continue;

        TH1D* h = h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent];
        h->Fill (trk_pt[iTrk], event_weight / trkEff);

        h = h_z_trk_pt[iSpc][iPtZ][iPhi][iCent];
        h->Fill (trk_pt[iTrk], event_weight / trkEff);

        const int iPtTrk = h->FindFixBin (trk_pt[iTrk]) - 1;
        if (iPtTrk >= 0 && iPtTrk < nPtTrkBins) {
          h = h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent];
          for (int iPhi = 1; iPhi <= h->GetNbinsX (); iPhi++)
            h->Fill (h->GetBinCenter (iPhi), event_weight / trkEff / h->GetNbinsX ());
        }

        h = h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent];
        h->Fill (trk_pt[iTrk] / z_pt, event_weight / trkEff);
      }
    } // end loop over pp tree
    cout << "Done minbias pp loop." << endl;
  }

  //CombineHists ();
  //ScaleHists ();

  SaveHists (outFileName);

  inFile->Close ();
  SaferDelete (inFile);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: CombineHists () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {

        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          if (iSpc == 0 && iPhi == 0)
            continue;
          h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]->Add (h_z_trk_raw_pt[0][iPtZ][0][iCent]);
          h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Add (h_z_trk_pt[0][iPtZ][0][iCent]);
          h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Add (h_z_trk_xzh[0][iPtZ][0][iCent]);

        } // end loop over phi
        h_z_trk_zpt[iSpc][iPtZ][iCent]->Add (h_z_trk_raw_pt[iSpc][iPtZ][0][iCent]);
        h_z_trk_zxzh[iSpc][iPtZ][iCent]->Add (h_z_trk_xzh[iSpc][iPtZ][0][iCent]);

        if (iSpc == 0)
          continue;

        for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Add (h_z_trk_phi[0][iPtZ][iPtTrk][iCent]);
        }
        h_z_counts[iSpc][iPtZ][iCent]->Add (h_z_counts[0][iPtZ][iCent]);
      } // end loop over pT^Z
    } // end loop over species
  } // end loop over centralities
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Scale histograms for plotting, calculating signals, etc.
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: ScaleHists () {
  if (histsScaled || !histsLoaded)
    return;

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        TH1D* countsHist = h_z_counts[iSpc][iPtZ][iCent];
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          const double yieldNormFactor = countsHist->GetBinContent (1) * (pi);
          //RescaleWithError (h, yieldNormFactor, yieldNormFactorError);
          h_z_trk_raw_pt[iSpc][iPtZ][iPhi][iCent]->Scale (1./ countsHist->GetBinContent (1));
          if (yieldNormFactor > 0) {
            h_z_trk_pt[iSpc][iPtZ][iPhi][iCent]->Scale (1. / yieldNormFactor, "width");
            h_z_trk_xzh[iSpc][iPtZ][iPhi][iCent]->Scale (1. / yieldNormFactor, "width");
          }
        } // end loop over phi
        h_z_trk_zpt[iSpc][iPtZ][iCent]->Scale (1./ (countsHist->GetBinContent (1) * (pi)), "width");
        h_z_trk_zxzh[iSpc][iPtZ][iCent]->Scale (1./ (countsHist->GetBinContent (1) * (pi)), "width");

        const double normFactor = countsHist->GetBinContent (1);
        for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          //h_z_trk_phi[iPtTrk][iPtZ][iCent][iSpc] = h_z_trk_pt_phi[iPtZ][iCent][iSpc]->ProjectionX (Form ("h_z_trk_phi_iPtTrk%i_iPtZ%i_iCent%i_%s", iPtTrk, iPtZ, iCent, name.c_str ()), iPtTrk+1, iPtTrk+1);

          TH1D* h = h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent];
          h->Rebin (2);
          if (iPtTrk > 3)
            h->Rebin (2);
          if (iCent != 0)
            h->Rebin (2);
          if (normFactor > 0)
            h->Scale (1. / normFactor);
        }
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species

  histsScaled = true;
  return;
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (OLD VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: GenerateWeights () {
  const int gw_nFCalBins = 100;
  const double* gw_fcalBins = linspace (0, 5200, gw_nFCalBins);
  const int gw_nQ2Bins = 20;
  const double* gw_q2Bins = linspace (0, 0.3, gw_nQ2Bins);

  const int gw_nNchBins = 160;
  const double* gw_nchBins = linspace (-0.5, 160.5, gw_nNchBins);
  
  TH1D* _h_PbPbFCalDist = nullptr;
  TH1D* _h_PbPbFCal_weights = nullptr;
  TH1D* _h_PbPbQ2Dist[numFinerCentBins];
  TH1D* _h_PbPbQ2_weights[numFinerCentBins];

  TH1D* _h_ppNchDist = nullptr;
  TH1D* _h_ppNch_weights = nullptr;
  

  SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
  TFile* inFile = new TFile (Form ("%s/Nominal/eventWeightsTree.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "recreate");

  _h_PbPbFCalDist = new TH1D (Form ("h_PbPbFCalDist_%s", name.c_str ()), "", gw_nFCalBins, gw_fcalBins);
  _h_PbPbFCal_weights = new TH1D (Form ("h_PbPbFCal_weights_%s", name.c_str ()), "", gw_nFCalBins, gw_fcalBins);
  _h_PbPbFCalDist->Sumw2 ();
  _h_PbPbFCal_weights->Sumw2 ();
  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    _h_PbPbQ2Dist[iCent] = new TH1D (Form ("h_PbPbQ2Dist_iCent%i_%s", iCent, name.c_str ()), "", gw_nQ2Bins, gw_q2Bins);
    _h_PbPbQ2_weights[iCent] = new TH1D (Form ("h_PbPbQ2_weights_iCent%i_%s", iCent, name.c_str ()), "", gw_nQ2Bins, gw_q2Bins);
    _h_PbPbQ2Dist[iCent]->Sumw2 ();
    _h_PbPbQ2_weights[iCent]->Sumw2 ();
  }
  _h_ppNchDist = new TH1D (Form ("h_ppNchDist_%s", name.c_str ()), "", gw_nNchBins, gw_nchBins);
  _h_ppNchDist->Sumw2 ();
  _h_ppNch_weights = new TH1D (Form ("h_ppNch_weights_%s", name.c_str ()), "", gw_nNchBins, gw_nchBins);
  _h_ppNch_weights->Sumw2 ();

  //bool passes_toroid = true;
  float fcal_et = 0, q2 = 0, event_weight = 1;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("fcal_et",    &fcal_et);
    PbPbTree->SetBranchAddress ("q2",         &q2);
    PbPbTree->SetBranchStatus ("vz",          0);
    PbPbTree->SetBranchStatus ("psi2",        0);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);
      _h_PbPbFCalDist->Fill (fcal_et);
    }
    cout << "Done 1st Pb+Pb loop." << endl;

    if (_h_PbPbFCalDist->Integral () > 0)
      _h_PbPbFCalDist->Scale (1./_h_PbPbFCalDist->Integral ());

    SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
    TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
    TH1D* referenceFCalDist = (TH1D*)ztrackFile->Get ("h_PbPbFCalDist_data");

    TH1D* referenceQ2Dist[numFinerCentBins];
    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
      referenceQ2Dist[iCent] = (TH1D*)ztrackFile->Get (Form ("h_PbPbQ2Dist_iCent%i_data", iCent));
    }

    for (int ix = 1; ix <= _h_PbPbFCal_weights->GetNbinsX (); ix++) {
      const double fcal_weight = (_h_PbPbFCalDist->GetBinContent (ix) != 0 ? referenceFCalDist->GetBinContent (ix) / _h_PbPbFCalDist->GetBinContent (ix) : 0);
      _h_PbPbFCal_weights->SetBinContent (ix, fcal_weight);
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

      event_weight = _h_PbPbFCal_weights->GetBinContent (_h_PbPbFCal_weights->FindBin (fcal_et));

      _h_PbPbQ2Dist[iCent]->Fill (q2, event_weight);
    }
    cout << "Done 2nd Pb+Pb loop." << endl;

    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
      if (_h_PbPbQ2Dist[iCent]->Integral () > 0)
        _h_PbPbQ2Dist[iCent]->Scale (1./_h_PbPbQ2Dist[iCent]->Integral ());
    }

    for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
      for (int ix = 1; ix <= _h_PbPbQ2_weights[iCent]->GetNbinsX (); ix++) {
        const double q2_weight = (_h_PbPbQ2Dist[iCent]->GetBinContent (ix) != 0 ? referenceQ2Dist[iCent]->GetBinContent (ix) / _h_PbPbQ2Dist[iCent]->GetBinContent (ix) : 0);
        const double q2_weight_err = (_h_PbPbQ2Dist[iCent]->GetBinContent (ix) != 0 ? q2_weight * sqrt (pow (_h_PbPbQ2Dist[iCent]->GetBinError (ix) / _h_PbPbQ2Dist[iCent]->GetBinContent (ix), 2) + pow (referenceQ2Dist[iCent]->GetBinError (ix) / referenceQ2Dist[iCent]->GetBinContent (ix), 2)) : 0);
        _h_PbPbQ2_weights[iCent]->SetBinContent (ix, q2_weight);
        _h_PbPbQ2_weights[iCent]->SetBinError (ix, q2_weight_err);
      }
    }

    ztrackFile->Close ();
    SaferDelete (ztrackFile);

    SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
  }

  if (ppTree) {
    int ntrk = 0;
    //ppTree->SetBranchAddress ("ntrk", &ntrk);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      _h_ppNchDist->Fill (ntrk);//, event_weight);
    }
    cout << "Done pp loop." << endl;

    if (_h_ppNchDist->Integral () > 0)
      _h_ppNchDist->Scale (1./_h_ppNchDist->Integral ());

    if (name == "data") {
      for (int ix = 1; ix <= _h_ppNch_weights->GetNbinsX (); ix++) {
        _h_ppNch_weights->SetBinContent (ix, 1);
      }
    } else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceNchDist = (TH1D*)ztrackFile->Get ("h_ppNchDist_data");

      for (int ix = 1; ix <= _h_ppNch_weights->GetNbinsX (); ix++) {
        const double nch_weight = (_h_ppNchDist->GetBinContent (ix) != 0 ? referenceNchDist->GetBinContent (ix) / _h_ppNchDist->GetBinContent (ix) : 0);
        _h_ppNch_weights->SetBinContent (ix, nch_weight);
      }

      ztrackFile->Close ();

      if (name == "data")
        SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      else if (name == "mc")
        SetupDirectories ("MCAnalysis/", "ZTrackAnalysis/");
      else if (name == "minbias")
        SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
      else if (name == "truth")
        SetupDirectories ("TruthAnalysis/", "ZTrackAnalysis/");
    }
  }


  inFile->Close ();
  SaferDelete (inFile);

  eventWeightsFile->cd ();

  SafeWrite (_h_PbPbFCalDist);
  SafeWrite (_h_PbPbFCal_weights);
  for (int iCent = 0; iCent < numFinerCentBins; iCent++) {
    SafeWrite (_h_PbPbQ2Dist[iCent]);
    SafeWrite (_h_PbPbQ2_weights[iCent]);
  }
  SafeWrite (_h_ppNchDist);
  SafeWrite (_h_ppNch_weights);

  eventWeightsFile->Close ();
}


#endif
