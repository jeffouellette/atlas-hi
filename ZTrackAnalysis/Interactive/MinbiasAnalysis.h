#ifndef __MinbiasAnalysis_h__
#define __MinbiasAnalysis_h__

#include "Params.h"
#include "Analysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class MinbiasAnalysis : public Analysis {

  public:
  MinbiasAnalysis (const char* _name = "minbias", const char* subDir = "Nominal") : Analysis () {
    name = _name;
    directory = Form ("MinbiasAnalysis/%s/", subDir);
    plotFill = true;
    plotSignal = false;
    useAltMarker = false;
    backgroundSubtracted = true;
    LoadTrackingEfficiencies ();
    SetupDirectories (directory, "ZTrackAnalysis/");
  }

  void Execute () override;

  void CombineHists () override;
  void ScaleHists () override;

  void GenerateWeights ();

  void OldExecute ();
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (NEW VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: Execute () {
  SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
  h_PbPbFCal_weights = (TH1D*) eventWeightsFile->Get ("h_PbPbFCal_weights_minbias");
  h_PbPbQ2_weights = (TH1D*) eventWeightsFile->Get ("h_PbPbQ2_weights_minbias");
  h_PbPbVZ_weights = (TH1D*) eventWeightsFile->Get ("h_PbPbVZ_weights_minbias");
  h_ppVZ_weights = (TH1D*) eventWeightsFile->Get ("h_ppVZ_weights_minbias");
  //h_PbPb_event_reweights = (TH3D*)eventWeightsFile->Get ("h_PbPbEventReweights_minbias");
  //h_pp_event_reweights = (TH1D*)eventWeightsFile->Get ("h_ppEventReweights_minbias");

  SetupDirectories (directory, "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*) inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*) inFile->Get ("ppZTrackTree");

  CreateHists ();

  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, event_weight = 1, fcal_weight = 1, q2_weight = 1, vz_weight = 1;
  vector<float>* trk_yield = nullptr;

  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("fcal_et",    &fcal_et);
    PbPbTree->SetBranchAddress ("q2",         &q2);
    PbPbTree->SetBranchAddress ("vz",         &vz);
    PbPbTree->SetBranchAddress ("trk_yield",  &trk_yield);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      const short iSpc = 0;
      const short iPtZ = 0;
      const short iXZTrk = 0;
      short iCent = 0;
      while (iCent < numCentBins) {
        if (fcal_et < centBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      fcal_weight = h_PbPbFCal_weights->GetBinContent (h_PbPbFCal_weights->FindBin (fcal_et));
      q2_weight = h_PbPbQ2_weights->GetBinContent (h_PbPbQ2_weights->FindBin (q2));
      vz_weight = h_PbPbVZ_weights->GetBinContent (h_PbPbVZ_weights->FindBin (vz));

      event_weight = fcal_weight * q2_weight * vz_weight;
      //event_weight = h_PbPb_event_reweights->GetBinContent (h_PbPb_event_reweights->FindBin (fcal_et, q2, vz));

      h_fcal_et->Fill (fcal_et);
      //h_fcal_et_q2->Fill (fcal_et, q2);
      h_fcal_et_reweighted->Fill (fcal_et, fcal_weight);
      //h_fcal_et_q2_reweighted->Fill (fcal_et, q2, event_weight);

      h_q2->Fill (q2);
      h_q2_reweighted->Fill (q2, q2_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, vz_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);

      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        TH1D* h = h_z_trk_pt[iSpc][iPtZ][iXZTrk][idPhi][iCent];
        //for (int iPhiTrk = 0; iPhiTrk < 3; iPhiTrk++) {
          for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
            const float trkEff = GetTrackingEfficiency (fcal_et, h->GetBinCenter (iPtTrk), 0, true); // TODO temporary -- will need to update grid code when tracking efficiencies are available!
            if (trkEff == 0)
              continue;
            h->SetBinContent (iPtTrk+1, h->GetBinContent (iPtTrk+1) + trk_yield->at (iPtTrk+nPtTrkBins*idPhi) / trkEff);
            h->SetBinError (iPtTrk+1, sqrt ( pow (h->GetBinError (iPtTrk+1), 2) + pow (trk_yield->at (iPtTrk+nPtTrkBins*idPhi) / trkEff, 2))); 
          } // end loop over PtTrk bins
        //} // end loop over dPhiTrk (from grid)
      } // end loop over dPhi bins

    } // end loop over Pb+Pb tree
    cout << endl;
  }


  if (ppTree) {
    ppTree->SetBranchAddress ("vz",         &vz);
    ppTree->SetBranchAddress ("trk_yield",  &trk_yield);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      const short iSpc = 0;
      const short iPtZ = 0;
      const short iXZTrk = 0;
      const short iCent = 0; // iCent = 0 for pp

      vz_weight = h_PbPbVZ_weights->GetBinContent (h_ppVZ_weights->FindBin (vz));

      event_weight = vz_weight;

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, vz_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);

      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        TH1D* h = h_z_trk_pt[iSpc][iPtZ][iXZTrk][idPhi][iCent];
        for (int iPhiTrk = 0; iPhiTrk < 3; iPhiTrk++) {
          for (int iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
            const float trkEff = GetTrackingEfficiency (fcal_et, h->GetBinCenter (iPtTrk), 0, false); // TODO temporary -- will need to update grid code when tracking efficiencies are available!
            if (trkEff == 0)
              continue;
            h->SetBinContent (iPtTrk+1, h->GetBinContent (iPtTrk+1) + trk_yield->at (iPtTrk+7*iPhiTrk) / trkEff);
            h->SetBinError (iPtTrk+1, sqrt ( pow (h->GetBinError (iPtTrk+1), 2) + pow (trk_yield->at (iPtTrk+nPtTrkBins*iPhiTrk) / trkEff, 2))); 
          } // end loop over PtTrk bins
        } // end loop over dPhiTrk (from grid)
      } // end loop over dPhi bins
    } // end loop over pp tree
    cout << endl;
  }

  CombineHists ();
  ScaleHists ();

  SaveHists ();

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
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          if (iSpc == 0 && iPtZ == 0 && iXZTrk == 0)
            continue;
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]->Add (h_z_trk_pt[0][0][0][iPhi][iCent]);
          } // end loop over phi
        }
        if (iSpc == 0 && iPtZ == 0)
          continue;
        h_z_counts[iSpc][iPtZ][iCent]->Add (h_z_counts[0][0][iCent]);
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
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            TH1D* h = h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent];
            TH1D* countsHist = h_z_counts[iSpc][iPtZ][iCent];
            const double yieldNormFactor = countsHist->GetBinContent (1) * (pi/3.);
            //RescaleWithError (h, yieldNormFactor, yieldNormFactorError);
            if (yieldNormFactor > 0)
              h->Scale (1. / yieldNormFactor);
          } // end loop over phi
        } // end loop over xZTrk
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
  const int gw_nQ2Bins = 50;
  const double* gw_q2Bins = linspace (0, 1, gw_nQ2Bins);
  const int gw_nVertZBins = 50;
  const double* gw_vertZBins = linspace (-200, 200, gw_nVertZBins);
  
  TH1D* _h_PbPbFCalDist = nullptr;
  TH1D* _h_PbPbVZDist = nullptr;
  TH1D* _h_PbPbQ2Dist = nullptr;
  TH1D* _h_PbPbFCal_weights = nullptr;
  TH1D* _h_PbPbVZ_weights = nullptr;
  TH1D* _h_PbPbQ2_weights = nullptr;
  
  TH1D* _h_ppVZDist = nullptr;
  TH1D* _h_ppVZ_weights = nullptr;


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

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "recreate");

  _h_PbPbFCalDist = new TH1D (Form ("h_PbPbFCalDist_%s", name.c_str ()), "", gw_nFCalBins, gw_fcalBins);
  _h_PbPbFCal_weights = new TH1D (Form ("h_PbPbFCal_weights_%s", name.c_str ()), "", gw_nFCalBins, gw_fcalBins);
  _h_PbPbQ2Dist = new TH1D (Form ("h_PbPbQ2Dist_%s", name.c_str ()), "", gw_nQ2Bins, gw_q2Bins);
  _h_PbPbQ2_weights = new TH1D (Form ("h_PbPbQ2_weights_%s", name.c_str ()), "", gw_nQ2Bins, gw_q2Bins);
  _h_PbPbVZDist = new TH1D (Form ("h_PbPbVZDist_%s", name.c_str ()), "", gw_nVertZBins, gw_vertZBins);
  _h_PbPbVZ_weights = new TH1D (Form ("h_PbPbVZ_weights_%s", name.c_str ()), "", gw_nVertZBins, gw_vertZBins);
  _h_PbPbFCalDist->Sumw2 ();
  _h_PbPbFCal_weights->Sumw2 ();
  _h_PbPbQ2Dist->Sumw2 ();
  _h_PbPbQ2_weights->Sumw2 ();
  _h_PbPbVZDist->Sumw2 ();
  _h_PbPbVZ_weights->Sumw2 ();

  _h_ppVZDist = new TH1D (Form ("h_ppVZDist_%s", name.c_str ()), "", gw_nVertZBins, gw_vertZBins);
  _h_ppVZ_weights = new TH1D (Form ("h_ppVZ_weights_%s", name.c_str ()), "", gw_nVertZBins, gw_vertZBins);
  _h_ppVZDist->Sumw2 ();
  _h_ppVZ_weights->Sumw2 ();

  bool passes_toroid = true;
  int nvert = 0;
  float fcalA_et = 0, fcalC_et = 0, fcalA_et_Cos = 0, fcalC_et_Cos = 0, fcalA_et_Sin = 0, fcalC_et_Sin = 0;
  float fcal_et = 0, q2 = 0, vz = 0;
  vector<float>* vert_z = nullptr;
  vector<int>* vert_type = nullptr;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("fcal_et",    &fcal_et);
    PbPbTree->SetBranchAddress ("q2",         &q2);
    PbPbTree->SetBranchAddress ("vz",         &vz);
    PbPbTree->SetBranchStatus ("trk_yield",   0);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      _h_PbPbFCalDist->Fill (fcal_et);
      _h_PbPbQ2Dist->Fill (q2);
      _h_PbPbVZDist->Fill (vz);
    }
    cout << endl;

    if (_h_PbPbFCalDist->Integral () > 0)
      _h_PbPbFCalDist->Scale (1./_h_PbPbFCalDist->Integral ());
    if (_h_PbPbQ2Dist->Integral () > 0)
      _h_PbPbQ2Dist->Scale (1./_h_PbPbQ2Dist->Integral ());
    if (_h_PbPbVZDist->Integral () > 0)
      _h_PbPbVZDist->Scale (1./_h_PbPbVZDist->Integral ());

    if (name == "data") {
      for (int ix = 1; ix <= _h_PbPbFCal_weights->GetNbinsX (); ix++) {
        _h_PbPbFCal_weights->SetBinContent (ix, 1);
      }
      for (int ix = 1; ix <= _h_PbPbQ2_weights->GetNbinsX (); ix++) {
        _h_PbPbQ2_weights->SetBinContent (ix, 1);
      }
      for (int ix = 1; ix <= _h_PbPbVZ_weights->GetNbinsX (); ix++) {
        _h_PbPbVZ_weights->SetBinContent (ix, 1);
      }
    }

    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceFCalDist = (TH1D*)ztrackFile->Get ("h_PbPbFCalDist_data");
      TH1D* referenceVZDist = (TH1D*)ztrackFile->Get ("h_PbPbVZDist_data");
      TH1D* referenceQ2Dist = (TH1D*)ztrackFile->Get ("h_PbPbQ2Dist_data");

      for (int ix = 1; ix <= _h_PbPbFCal_weights->GetNbinsX (); ix++) {
        const double fcal_weight = (_h_PbPbFCalDist->GetBinContent (ix) != 0 ? referenceFCalDist->GetBinContent (ix) / _h_PbPbFCalDist->GetBinContent (ix) : 0);
        _h_PbPbFCal_weights->SetBinContent (ix, fcal_weight);
      }
      for (int ix = 1; ix <= _h_PbPbQ2_weights->GetNbinsX (); ix++) {
        const double q2_weight = (_h_PbPbQ2Dist->GetBinContent (ix) != 0 ? referenceQ2Dist->GetBinContent (ix) / _h_PbPbQ2Dist->GetBinContent (ix) : 0);
        _h_PbPbQ2_weights->SetBinContent (ix, q2_weight);
      }
      for (int ix = 1; ix <= _h_PbPbVZ_weights->GetNbinsX (); ix++) {
        const double vz_weight = (_h_PbPbVZDist->GetBinContent (ix) != 0 ? referenceVZDist->GetBinContent (ix) / _h_PbPbVZDist->GetBinContent (ix) : 0);
        _h_PbPbVZ_weights->SetBinContent (ix, vz_weight);
      }

      ztrackFile->Close ();
      SaferDelete (ztrackFile);
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
    ppTree->SetBranchAddress ("vz",         &vz);
    ppTree->SetBranchStatus ("trk_yield",   0);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      _h_ppVZDist->Fill (vz);
    }
    cout << endl;
    if (_h_ppVZDist->Integral () > 0)
      _h_ppVZDist->Scale (1./_h_ppVZDist->Integral ());

    if (name == "data") {
      for (int ix = 1; ix <= _h_ppVZ_weights->GetNbinsX (); ix++) {
        _h_ppVZ_weights->SetBinContent (ix, 1);
      }
    }
    else {
      SetupDirectories ("DataAnalysis/", "ZTrackAnalysis/");
      TFile* ztrackFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
      TH1D* referenceVZDist = (TH1D*)ztrackFile->Get ("h_ppVZDist_data");
      for (int ix = 1; ix <= _h_ppVZ_weights->GetNbinsX (); ix++) {
        const double vz_weight = (_h_ppVZDist->GetBinContent (ix) != 0. ? referenceVZDist->GetBinContent (ix) / _h_ppVZDist->GetBinContent (ix) : 0);
        _h_ppVZ_weights->SetBinContent (ix, vz_weight);
      }

      ztrackFile->Close ();
      SaferDelete (ztrackFile);
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

  inFile->Close ();
  SaferDelete (inFile);

  eventWeightsFile->cd ();

  SafeWrite (_h_PbPbFCalDist);
  SafeWrite (_h_PbPbVZDist);
  SafeWrite (_h_PbPbQ2Dist);
  SafeWrite (_h_PbPbFCal_weights);
  SafeWrite (_h_PbPbVZ_weights);
  SafeWrite (_h_PbPbQ2_weights);
  SafeWrite (_h_ppVZDist);
  SafeWrite (_h_ppVZ_weights);

  eventWeightsFile->Close ();
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately. (OLD VERSION)
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis :: OldExecute () {
  SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
  h_PbPbFCal_weights = (TH1D*) eventWeightsFile->Get ("h_PbPbFCal_weights_minbias");
  h_PbPbQ2_weights = (TH1D*) eventWeightsFile->Get ("h_PbPbQ2_weights_minbias");
  h_PbPbVZ_weights = (TH1D*) eventWeightsFile->Get ("h_PbPbVZ_weights_minbias");
  h_ppVZ_weights = (TH1D*) eventWeightsFile->Get ("h_ppVZ_weights_minbias");
  //h_PbPb_event_reweights = (TH3D*)eventWeightsFile->Get ("h_PbPbEventReweights_minbias");
  //h_pp_event_reweights = (TH1D*)eventWeightsFile->Get ("h_ppEventReweights_minbias");

  SetupDirectories (directory, "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_eta = 0, z_phi = 0, event_weight = 1, fcal_weight = 1, q2_weight = 1, vz_weight = 1;
  int ntrk = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;
  double** trkPtProj = Get2DArray <double> (numPhiBins, nPtTrkBins);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("fcal_et", &fcal_et);
    PbPbTree->SetBranchAddress ("q2",      &q2);
    PbPbTree->SetBranchAddress ("psi2",    &psi2);
    PbPbTree->SetBranchAddress ("vz",      &vz);
    PbPbTree->SetBranchAddress ("z_eta",   &z_eta);
    PbPbTree->SetBranchAddress ("z_phi",   &z_phi);
    PbPbTree->SetBranchAddress ("ntrk",    &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",  &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta", &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi", &trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      const short iSpc = 0;
      const short iPtZ = 0;
      short iCent = 0;
      while (iCent < numCentBins) {
        if (fcal_et < centBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      fcal_weight = h_PbPbFCal_weights->GetBinContent (h_PbPbFCal_weights->FindBin (fcal_et));
      q2_weight = h_PbPbQ2_weights->GetBinContent (h_PbPbQ2_weights->FindBin (q2));
      vz_weight = h_PbPbVZ_weights->GetBinContent (h_PbPbVZ_weights->FindBin (vz));

      event_weight = fcal_weight * q2_weight * vz_weight;
      //event_weight = h_PbPb_event_reweights->GetBinContent (h_PbPb_event_reweights->FindBin (fcal_et, q2, vz));

      h_fcal_et->Fill (fcal_et);
      //h_fcal_et_q2->Fill (fcal_et, q2);
      h_fcal_et_reweighted->Fill (fcal_et, fcal_weight);
      //h_fcal_et_q2_reweighted->Fill (fcal_et, q2, event_weight);

      h_q2->Fill (q2);
      h_q2_reweighted->Fill (q2, q2_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, vz_weight);

      float dphi = DeltaPhi (z_phi, psi2, false);
      if (dphi > pi/2)
        dphi = pi - dphi;

      h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);

      for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
          trkPtProj[iPhi][iPtTrk] = 0;
        }
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        const short iXZTrk = 0;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), true);
        if (trkEff == 0)
          continue;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight / trkEff);

        // Add to missing pT (requires dphi in +/-pi/2 to +/-pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        bool awaySide = false;
        if (dphi > pi/2) {
          dphi = pi-dphi;
          awaySide = true;
        }

        short iPtTrk = 0;
        while (iPtTrk < nPtTrkBins && trkpt > ptTrkBins[iPtTrk+1])
          iPtTrk++;
        // start at the 1st phi bin and integrate outwards until the track is no longer contained 
        // e.g. so 7pi/8->pi is a subset of pi/2->pi
        short iPhi = 0;
        while (iPhi < numPhiTrkBins && dphi > phiTrkBins[iPhi]) {
          if (awaySide)
            trkPtProj[iPhi][iPtTrk] += -trkpt * cos (dphi) / trkEff;
          else
            trkPtProj[iPhi][iPtTrk] += trkpt * cos (dphi) / trkEff;
          iPhi++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++)
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi])
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][idPhi][iCent]->Fill (trkpt, event_weight / trkEff);

        //// Study correlations (requires dphi in -pi/2 to 3pi/2)
        //dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        //if (dphi < -pi/2)
        //  dphi = dphi + 2*pi;

        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++)
          h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]->Fill (dphi, trkpt, event_weight / trkEff);
      } // end loop over tracks

      for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
          trkPtProj[iPhi][iPtTrk] = 0;
        }
      }
      
    } // end loop over Pb+Pb tree
    cout << endl;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("vz",      &vz);
    ppTree->SetBranchAddress ("z_eta",   &z_eta);
    ppTree->SetBranchAddress ("z_phi",   &z_phi);
    ppTree->SetBranchAddress ("ntrk",    &ntrk);
    ppTree->SetBranchAddress ("trk_pt",  &trk_pt);
    ppTree->SetBranchAddress ("trk_eta", &trk_eta);
    ppTree->SetBranchAddress ("trk_phi", &trk_phi);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      //if (fabs (vz) > 1.5)
      //  continue; // vertex cut

      const short iSpc = 0;
      const short iPtZ = 0;
      const short iCent = 0; // iCent = 0 for pp

      vz_weight = h_ppVZ_weights->GetBinContent (h_ppVZ_weights->FindBin (vz));

      event_weight = vz_weight;
      //event_weight = h_pp_event_reweights->GetBinContent (h_pp_event_reweights->FindBin (vz));

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, vz_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        const short iXZTrk = 0;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), false);
        if (trkEff == 0)
          continue;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight / trkEff);

        // Add to missing pT (requires dphi in -pi/2 to pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        bool awaySide = false;
        if (dphi > pi/2) {
          dphi = pi-dphi;
          awaySide = true;
        }

        short iPtTrk = 0;
        while (iPtTrk < nPtTrkBins && trkpt > ptTrkBins[iPtTrk+1])
          iPtTrk++;
        // start at the 1st phi bin and integrate outwards until the track is no longer contained 
        // e.g. so 7pi/8->pi is a subset of pi/2->pi
        short iPhi = 0;
        while (iPhi < numPhiTrkBins && dphi > phiTrkBins[iPhi]) {
          if (awaySide)
            trkPtProj[iPhi][iPtTrk] += -trkpt * cos (dphi) / trkEff;
          else
            trkPtProj[iPhi][iPtTrk] += trkpt * cos (dphi) / trkEff;
          iPhi++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++)
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi])
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][idPhi][iCent]->Fill (trkpt, event_weight / trkEff);

        //// Study correlations (requires dphi in -pi/2 to 3pi/2)
        //dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        //if (dphi < -pi/2)
        //  dphi = dphi + 2*pi;

        h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]->Fill (dphi, trkpt, event_weight / trkEff);
      } // end loop over tracks

      for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
          trkPtProj[iPhi][iPtTrk] = 0;
        }
      }
    } // end loop over pp tree
    cout << endl;
  }

  CombineHists ();
  ScaleHists ();

  SaveHists ();
  //LoadHists ();

  inFile->Close ();
  SaferDelete (inFile);

  Delete2DArray (trkPtProj, numPhiBins, nPtTrkBins);
}

#endif
