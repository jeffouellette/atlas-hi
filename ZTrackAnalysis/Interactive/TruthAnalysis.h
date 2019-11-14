#ifndef __TruthAnalysis_h__
#define __TruthAnalysis_h__

#include "Params.h"
#include "FullAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class TruthAnalysis : public FullAnalysis {

  public:
  TH1D*** h_z_jet_counts      = Get2DArray <TH1D*> (nPtZBins, 2); // iPtZ, hasHighPtTrk
  TH2D*   h_z_jet_pt          = nullptr;
  TH1D**  h_z_jet_xzj         = Get1DArray <TH1D*> (nPtZBins); // iPtZ
  TH2D*** h_z_jet_deta_dphi   = Get2DArray <TH2D*> (nPtZBins, 2); // iPtZ, hasHighPtTrk
  TH1D*** h_z_jet_dphi        = Get2DArray <TH1D*> (nPtZBins, 2); // iPtZ, hasHighPtTrk
  TH1D**  h_jet_jet_dphi      = Get1DArray <TH1D*> (nPtZBins); // iPtZ
  TH2D**  h_jet_trk_deta_dphi = Get1DArray <TH2D*> (nPtZBins); // iPtZ
  TH1D**  h_jet_trk_dphi      = Get1DArray <TH1D*> (nPtZBins); // iPtZ
  

  TruthAnalysis (const char* _name = "truth") : FullAnalysis () {
    name = _name;
    plotFill = false;
    useAltMarker = false;
    hasBkg = false;
    histsUnfolded = true;
    isMC = true;
  }

  void Execute (const char* inFileName, const char* outFileName) override;

  void CreateHists () override;
  void ScaleHists () override;
  void LoadHists (const char* histFileName = "savedHists.root", const bool _finishHists = true) override;
  void SaveHists (const char* histFileName = "savedHists.root") override;

  void PlotZJetPt ();
  void PlotxZJet ();
  void PlotZJetCorrelations ();
  void PlotJetTrkCorrelations ();
};


////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis :: CreateHists () {
  FullAnalysis :: CreateHists ();

  h_z_jet_pt = new TH2D ("h_z_jet_pt", "", 75, 0, 300, 75, 1, 300);
  h_z_jet_pt->Sumw2 ();

  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    h_z_jet_xzj[iPtZ] = new TH1D (Form ("h_z_jet_xzj_iPtZ%i", iPtZ), "", 100, 0, 1.5);
    h_z_jet_xzj[iPtZ]->Sumw2 ();

    for (int hasHighPtTrk = 0; hasHighPtTrk <= 1; hasHighPtTrk++) {
      h_z_jet_counts[iPtZ][hasHighPtTrk] = new TH1D (Form ("h_z_jet_counts_iPtZ%i_%s", iPtZ, hasHighPtTrk==1 ? "hasHighPtTrk":"inclusive"), "", 1, 0, 1);
      h_z_jet_counts[iPtZ][hasHighPtTrk]->Sumw2 ();

      h_z_jet_deta_dphi[iPtZ][hasHighPtTrk] = new TH2D (Form ("h_z_jet_deta_dphi_iPtZ%i_%s", iPtZ, hasHighPtTrk==1 ? "hasHighPtTrk":"inclusive"), "", 50, -5, 5, 18, -pi/2, 3*pi/2);
      h_z_jet_deta_dphi[iPtZ][hasHighPtTrk]->Sumw2 ();
    }

    h_jet_jet_dphi[iPtZ] = new TH1D (Form ("h_jet_jet_dphi_iPtZ%i", iPtZ), "", 18, -pi/2, 3*pi/2);
    h_jet_jet_dphi[iPtZ]->Sumw2 ();

    h_jet_trk_deta_dphi[iPtZ] = new TH2D (Form ("h_jet_trk_deta_dphi_iPtZ%i", iPtZ), "", 50, -5, 5, 80, -pi/2, 3*pi/2);
    h_jet_trk_deta_dphi[iPtZ]->Sumw2 ();
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Scale histograms for plotting, calculating signals, etc.
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis :: ScaleHists () {
  if (histsScaled || !histsLoaded)
    return;

  FullAnalysis :: ScaleHists ();

  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    h_z_jet_dphi[iPtZ][0] = h_z_jet_deta_dphi[iPtZ][0]->ProjectionY (Form ("h_z_jet_dphi_iPtZ%i_inclusive", iPtZ));
    h_z_jet_dphi[iPtZ][1] = h_z_jet_deta_dphi[iPtZ][1]->ProjectionY (Form ("h_z_jet_dphi_iPtZ%i_hasHighPtTrk", iPtZ));
    h_jet_trk_dphi[iPtZ] = h_jet_trk_deta_dphi[iPtZ]->ProjectionY (Form ("h_jet_trk_dphi_iPtZ%i", iPtZ));

    if (h_z_jet_counts[iPtZ][0]->GetBinContent (1) > 0) {
      h_z_jet_xzj[iPtZ]->Scale (1. / h_z_jet_counts[iPtZ][0]->GetBinContent (1), "width");
      h_z_jet_dphi[iPtZ][0]->Scale (1. / h_z_jet_counts[iPtZ][0]->GetBinContent (1));
      h_jet_trk_dphi[iPtZ]->Scale (1. / h_z_jet_counts[iPtZ][0]->GetBinContent (1));
    }

    h_jet_jet_dphi[iPtZ]->Scale (1. / h_z_jet_counts[iPtZ][1]->GetBinContent (1));

    if (h_z_jet_counts[iPtZ][1]->GetBinContent (1) > 0) {
      h_z_jet_dphi[iPtZ][1]->Scale (1. / h_z_jet_counts[iPtZ][1]->GetBinContent (1));
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis :: LoadHists (const char* histFileName, const bool _finishHists) {
  if (histsLoaded)
    return;

  FullAnalysis :: LoadHists (histFileName, false);

  if (!histFile) {
    SetupDirectories ("", "ZTrackAnalysis/");
    histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "read");
  }
  TDirectory* _gDirectory = gDirectory;
  if (!histFile->IsOpen ()) {
    cout << "Error in TruthAnalysis :: LoadHists: histFile not open after calling parent function, exiting." << endl;
    return;
  }

  h_z_jet_pt = (TH2D*)histFile->Get ("h_z_jet_pt");
  //h_z_jet_xzj = (TH2D*)histFile->Get ("h_z_jet_xzj");

  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    h_z_jet_counts[iPtZ][0] = (TH1D*)histFile->Get (Form ("h_z_jet_counts_iPtZ%i_inclusive", iPtZ));
    h_z_jet_counts[iPtZ][1] = (TH1D*)histFile->Get (Form ("h_z_jet_counts_iPtZ%i_hasHighPtTrk", iPtZ));

    h_z_jet_deta_dphi[iPtZ][0] = (TH2D*)histFile->Get (Form ("h_z_jet_deta_dphi_iPtZ%i_inclusive", iPtZ));
    h_z_jet_dphi[iPtZ][0] = (TH1D*)histFile->Get (Form ("h_z_jet_dphi_iPtZ%i_inclusive", iPtZ));

    h_z_jet_deta_dphi[iPtZ][1] = (TH2D*)histFile->Get (Form ("h_z_jet_deta_dphi_iPtZ%i_hasHighPtTrk", iPtZ));
    h_z_jet_dphi[iPtZ][1] = (TH1D*)histFile->Get (Form ("h_z_jet_dphi_iPtZ%i_hasHighPtTrk", iPtZ));

    h_jet_jet_dphi[iPtZ] = (TH1D*)histFile->Get (Form ("h_jet_jet_dphi_iPtZ%i", iPtZ));

    h_z_jet_xzj[iPtZ] = (TH1D*)histFile->Get (Form ("h_z_jet_xzj_iPtZ%i", iPtZ));

    h_jet_trk_deta_dphi[iPtZ] = (TH2D*)histFile->Get (Form ("h_jet_trk_deta_dphi_iPtZ%i", iPtZ));
    h_jet_trk_dphi[iPtZ] = (TH1D*)histFile->Get (Form ("h_jet_trk_dphi_iPtZ%i", iPtZ));
  }

  _gDirectory->cd ();

  histsLoaded = true;

  if (_finishHists) {
    FullAnalysis :: CombineHists ();
    TruthAnalysis :: ScaleHists ();
  }

  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis :: SaveHists (const char* histFileName) {
  FullAnalysis :: SaveHists (histFileName);

  if (!histFile) {
    SetupDirectories ("", "ZTrackAnalysis/");
    histFile = new TFile (Form ("%s/%s", rootPath.Data (), histFileName), "update");
    histFile->cd ();
  }

  SafeWrite (h_z_jet_pt);

  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    SafeWrite (h_z_jet_counts[iPtZ][0]);
    SafeWrite (h_z_jet_counts[iPtZ][1]);

    SafeWrite (h_z_jet_xzj[iPtZ]);

    SafeWrite (h_jet_jet_dphi[iPtZ]);

    SafeWrite (h_z_jet_deta_dphi[iPtZ][0]);
    SafeWrite (h_z_jet_dphi[iPtZ][0]);
    SafeWrite (h_z_jet_deta_dphi[iPtZ][1]);
    SafeWrite (h_z_jet_dphi[iPtZ][1]);

    SafeWrite (h_jet_trk_deta_dphi[iPtZ]);
    SafeWrite (h_jet_trk_dphi[iPtZ]);
  }

  histFile->Close ();
  histFile = nullptr;
  histsLoaded = false;
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1;//, fcal_weight = 1, q2_weight = 1, psi2_weight = 1, vz_weight = 1, nch_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0;
  float z_pt = 0, z_y = 0, z_phi = 0, z_eta = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0, njet = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;
  vector<float>* jet_pt = nullptr, *jet_eta = nullptr, *jet_phi = nullptr, *jet_e = nullptr;


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("isEE",         &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("vz",           &vz);
    PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
    PbPbTree->SetBranchAddress ("z_y",          &z_y);
    PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
    PbPbTree->SetBranchAddress ("z_m",          &z_m);
    PbPbTree->SetBranchAddress ("l1_pt",        &l1_pt);
    PbPbTree->SetBranchAddress ("l1_eta",       &l1_eta);
    PbPbTree->SetBranchAddress ("l1_phi",       &l1_phi);
    PbPbTree->SetBranchAddress ("l1_charge",    &l1_charge);
    PbPbTree->SetBranchAddress ("l2_pt",        &l2_pt);
    PbPbTree->SetBranchAddress ("l2_eta",       &l2_eta);
    PbPbTree->SetBranchAddress ("l2_phi",       &l2_phi);
    PbPbTree->SetBranchAddress ("l2_charge",    &l2_charge);
    PbPbTree->SetBranchAddress ("ntrk",         &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",       &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",      &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",      &trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

      const short iCent = GetCentBin (fcal_et);
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      //const short iFinerCent = GetFinerCentBin (fcal_et);
      //if (iFinerCent < 1 || iFinerCent > numFinerCentBins-1)
      //  continue;

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin

      //fcal_weight = h_PbPbFCal_weights[iPtZ]->GetBinContent (h_PbPbFCal_weights[iPtZ]->FindBin (fcal_et));
      //q2_weight = h_PbPbQ2_weights[iCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iCent][iPtZ]->FindBin (q2));
      //psi2_weight = h_PbPbPsi2_weights[iCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iCent][iPtZ]->FindBin (psi2));

      //event_weight = event_weight * fcal_weight * q2_weight * psi2_weight * vz_weight;
      if (event_weight == 0)
        continue;

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);

      h_q2[iCent]->Fill (q2);
      h_q2_reweighted[iCent]->Fill (q2, event_weight);
      h_psi2[iCent]->Fill (psi2);
      h_psi2_reweighted[iCent]->Fill (psi2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      TLorentzVector zvec;
      zvec.SetPxPyPzE (z_pt*cos(z_phi), z_pt*sin(z_phi), sqrt(z_pt*z_pt+z_m*z_m)*sinh(z_y), sqrt(z_pt*z_pt+z_m*z_m)*cosh(z_y));
      z_eta = zvec.Eta ();

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      h_z_y_phi[iCent][iSpc][iPtZ]->Fill (z_y, InTwoPi (z_phi), event_weight);
      h_z_eta[iCent][iSpc][iPtZ]->Fill (z_eta, event_weight);
      h_z_y[iCent][iSpc][iPtZ]->Fill (z_y, event_weight);
      int iReg = (fabs (z_y) > 1.00 ? 1 : 0); // barrel vs. endcaps
      h_z_m[iCent][iSpc][iReg]->Fill (z_m, event_weight);

      h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
      h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l1_eta, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l2_eta, event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);

      if (z_pt > 5) {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight);

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, event_weight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_trk_pt_dphi_unscaled[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, event_weight);
          }
        }
      } // end loop over tracks

    } // end loop over Pb+Pb tree
    cout << "Done truth-level Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("event_weight",  &event_weight);
    ppTree->SetBranchAddress ("isEE",          &isEE);
    ppTree->SetBranchAddress ("vz",            &vz);
    ppTree->SetBranchAddress ("z_pt",          &z_pt);
    ppTree->SetBranchAddress ("z_y",           &z_y);
    ppTree->SetBranchAddress ("z_phi",         &z_phi);
    ppTree->SetBranchAddress ("z_m",           &z_m);
    ppTree->SetBranchAddress ("l1_pt",         &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",        &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",        &l1_phi);
    ppTree->SetBranchAddress ("l1_charge",     &l1_charge);
    ppTree->SetBranchAddress ("l2_pt",         &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",        &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",        &l2_phi);
    ppTree->SetBranchAddress ("l2_charge",     &l2_charge);
    ppTree->SetBranchAddress ("ntrk",          &ntrk);
    ppTree->SetBranchAddress ("trk_pt",        &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",       &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",       &trk_phi);
    ppTree->SetBranchAddress ("truth_jet_n",   &njet);
    ppTree->SetBranchAddress ("truth_jet_pt",  &jet_pt);
    ppTree->SetBranchAddress ("truth_jet_eta", &jet_eta);
    ppTree->SetBranchAddress ("truth_jet_phi", &jet_phi);
    ppTree->SetBranchAddress ("truth_jet_e",   &jet_e);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin

      bool triggerEvent = false;
      //float triggerPhi = 0;
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        if (trk_pt->at (iTrk) > 20 && DeltaPhi (z_phi, trk_phi->at (iTrk)) < pi/2) {
          triggerEvent = true;
          //triggerPhi = trk_phi->at (iTrk);
        }
      }

      TLorentzVector ztlv;
      ztlv.SetPxPyPzE (z_pt * cos (z_phi), z_pt * sin (z_phi), sqrt (z_pt*z_pt + z_m*z_m) * sinh (z_y), sqrt (z_pt*z_pt + z_m*z_m) * cosh (z_y));
      z_eta = ztlv.Eta ();

      //nch_weight = h_ppNch_weights->GetBinContent (h_ppNch_weights->FindBin (ntrk));

      //event_weight = event_weight * vz_weight * nch_weight;
      if (event_weight == 0)
        continue;

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      TLorentzVector zvec;
      zvec.SetPxPyPzE (z_pt*cos(z_phi), z_pt*sin(z_phi), sqrt(z_pt*z_pt+z_m*z_m)*sinh(z_y), sqrt(z_pt*z_pt+z_m*z_m)*cosh(z_y));
      z_eta = zvec.Eta ();

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      h_z_y_phi[iCent][iSpc][iPtZ]->Fill (z_y, InTwoPi (z_phi), event_weight);
      h_z_eta[iCent][iSpc][iPtZ]->Fill (z_eta, event_weight);
      h_z_y[iCent][iSpc][iPtZ]->Fill (z_y, event_weight);
      int iReg = (fabs (z_y) > 1.00 ? 1 : 0); // barrel vs. endcaps
      h_z_m[iCent][iSpc][iReg]->Fill (z_m, event_weight);

      h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
      h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l1_eta, event_weight);
      h_lepton_eta[iCent][iSpc]->Fill (l2_eta, event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
      h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);

      if (z_pt > 5) {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

      h_z_jet_counts[iPtZ][0]->Fill (0.5, event_weight);
      h_z_jet_counts[iPtZ][0]->Fill (1.5);
      if (triggerEvent) {
        h_z_jet_counts[iPtZ][1]->Fill (0.5, event_weight);
        h_z_jet_counts[iPtZ][1]->Fill (1.5);
      }

      for (int iJet = 0; iJet < njet; iJet++) {
        if (DeltaR (l1_eta, jet_eta->at (iJet), l1_phi, jet_phi->at (iJet)) < 0.2)
          continue;
        if (DeltaR (l2_eta, jet_eta->at (iJet), l2_phi, jet_phi->at (iJet)) < 0.2)
          continue;

        h_z_jet_pt->Fill (z_pt, jet_pt->at (iJet), event_weight);
    
        float dphi = DeltaPhi (z_phi, jet_phi->at (iJet), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        if (dphi > 7*pi/8 && dphi < 9*pi/8)
          h_z_jet_xzj[iPtZ]->Fill (jet_pt->at (iJet) / z_pt, event_weight);

        if (jet_pt->at (iJet) > 20) {
          h_z_jet_deta_dphi[iPtZ][0]->Fill (z_eta - jet_eta->at (iJet), dphi, event_weight);
          if (triggerEvent) {
            h_z_jet_deta_dphi[iPtZ][1]->Fill (z_eta - jet_eta->at (iJet), dphi, event_weight);

            for (int iJet2 = 0; iJet2 < iJet; iJet2++) {
              if (DeltaR (l1_eta, jet_eta->at (iJet2), l1_phi, jet_phi->at (iJet2)) < 0.2)
                continue;
              if (DeltaR (l2_eta, jet_eta->at (iJet2), l2_phi, jet_phi->at (iJet2)) < 0.2)
                continue;
              if (jet_pt->at (iJet2) < 20)
                continue;

              dphi = DeltaPhi (jet_phi->at (iJet), jet_phi->at (iJet2), true);
              //dphi = DeltaPhi (jet_phi->at (iJet), triggerPhi, true);
              if (dphi < -pi/2) dphi = dphi+2*pi;
              h_jet_jet_dphi[iPtZ]->Fill (dphi, event_weight);
            }
          }
        }
      } // end loop over jets

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;

        for (int iJet = 0; iJet < njet; iJet++) {
          if (DeltaR (l1_eta, jet_eta->at (iJet), l1_phi, jet_phi->at (iJet)) < 0.2)
            continue;
          if (DeltaR (l2_eta, jet_eta->at (iJet), l2_phi, jet_phi->at (iJet)) < 0.2)
            continue;
          if (fabs (jet_eta->at (iJet) - trk_eta->at (iTrk)) < 2)
            continue;

          float dphi = DeltaPhi (jet_phi->at (iJet), trk_phi->at (iTrk), true);
          if (dphi < -pi/2)
            dphi = dphi + 2*pi;
          h_jet_trk_deta_dphi[iPtZ]->Fill (jet_eta->at (iJet) - trk_eta->at (iTrk), dphi, event_weight);
        } // end loop over jets

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight);

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, event_weight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_trk_pt_dphi_unscaled[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt / z_pt, event_weight);
          }
        }
      } // end loop over tracks

    } // end loop over pp tree
    cout << "Done truth-level pp loop." << endl;
  }

  //CombineHists ();
  //ScaleHists ();
  
  SaveHists (outFileName);

  inFile->Close ();
  if (inFile) { delete inFile; inFile = nullptr; }
}




////////////////////////////////////////////////////////////////////////////////////////////////
// Plots Z-Jet pT correlation.
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis::PlotZJetPt () {
  SetupDirectories ("", "ZTrackAnalysis/");

  TCanvas* c = new TCanvas ("ZJetPtCanvas", "", 800, 600);
  c->SetRightMargin (0.14);
  gPad->SetLogz ();

  h_z_jet_pt->GetXaxis ()->SetTitle ("Truth #it{p}_{T}^{ Z} [GeV]");
  h_z_jet_pt->GetYaxis ()->SetTitle ("Truth #it{p}_{T}^{ jet} [GeV]");
  h_z_jet_pt->GetZaxis ()->SetTitle ("Counts");

  h_z_jet_pt->Draw ("colz");
  c->SaveAs (Form ("%s/ZJetPtCorrelation.pdf", plotPath.Data ()));
}




void TruthAnalysis::PlotxZJet () {
  SetupDirectories ("", "ZTrackAnalysis/");

  TCanvas* c = new TCanvas ("xZJetCanvas", "", 800, 600);

  double max = 0, min = 1e30;
  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h = h_z_jet_xzj[iPtZ];
    if (h->GetMaximum () > max)  max = h->GetMaximum ();
    if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
  }
  max = 1.3*max;
  min = 0.7*min;

  for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h = h_z_jet_xzj[iPtZ];

    h->GetYaxis ()->SetRangeUser (min, max);

    h->GetXaxis ()->SetTitle ("#it{x}_{Z}^{jet}");
    h->GetYaxis ()->SetTitle ("1/N_{Z} dN/dx");

    h->SetMarkerColor (colors[iPtZ-1]);
    h->SetLineColor (colors[iPtZ-1]);

    h->Draw (iPtZ == 0 ? "hist" : "same hist");
    myText (0.65, 0.9-0.06*(iPtZ-1), colors[iPtZ-1], Form ("%g < #it{p}_{T}^{ Z} < %g [GeV]", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);
  }
  myText (0.3, 0.9, kBlack, "#it{p}_{T}^{ jet} > 5 GeV", 0.04);
  myText (0.3, 0.84, kBlack, "|#Delta#phi| > 7#pi/8", 0.04);
  c->SaveAs (Form ("%s/xZJetDists.pdf", plotPath.Data ()));

}




void TruthAnalysis::PlotZJetCorrelations () {
  TCanvas* c = new TCanvas ("ZJetCorrelationsCanvas", "", 800, 600);

  float max = 1.2 * fmax (h_z_jet_dphi[2][0]->GetMaximum (), h_z_jet_dphi[2][1]->GetMaximum ());
  float min = 0;

  TH1D* h = h_z_jet_dphi[2][0];

  h->GetYaxis ()->SetRangeUser (min, max);

  h->GetXaxis ()->SetTitle ("#Delta#phi");
  h->GetYaxis ()->SetTitle ("Y (#Delta#phi)");

  h->SetMarkerColor (kBlue+1);
  h->SetLineColor (kBlue+1);

  h->Draw ("e1");
  //myText (0.65, 0.9-0.06*(iPtZ-1), colors[iPtZ-1], Form ("%g < #it{p}_{T}^{ Z} < %g [GeV]", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);


  h = h_z_jet_dphi[2][1];

  h->GetYaxis ()->SetRangeUser (min, max);

  h->GetXaxis ()->SetTitle ("#Delta#phi");
  h->GetYaxis ()->SetTitle ("Y (#Delta#phi)");

  h->SetMarkerColor (kOrange+8);
  h->SetLineColor (kOrange+8);

  h->Draw ("same e1");


  h = h_jet_jet_dphi[2];

  h->GetYaxis ()->SetRangeUser (min, max);

  h->GetXaxis ()->SetTitle ("#Delta#phi");
  h->GetYaxis ()->SetTitle ("Y (#Delta#phi)");

  h->SetMarkerColor (kBlack);
  h->SetLineColor (kBlack);
  h->SetMarkerStyle (kOpenCircle);

  h->Draw ("same e1");

  myText (0.22, 0.9, kBlack, "#it{p}_{T}^{Z} > 25 GeV, #it{p}_{T}^{jet} > 20 GeV", 0.04);
  myText (0.22, 0.83, kBlue+1, "Inclusive events Z-jet #Delta#phi", 0.04);
  myText (0.22, 0.77, kOrange+8, "Triggered events Z-jet #Delta#phi", 0.04);
  myText (0.22, 0.71, kBlack, "Jet-jet #Delta#phi", 0.04);

  myText (0.62, 0.88, kBlack, "Pythia truth-level", 0.04);

  c->SaveAs (Form ("%s/ZJetCorrelations.pdf", plotPath.Data ()));
}




void TruthAnalysis::PlotJetTrkCorrelations () {
  SetupDirectories ("", "ZTrackAnalysis/");

  TCanvas* c = new TCanvas ("JetTrkCorrelationsCanvas", "", 800, 600);

  double max = 0, min = 1e30;
  for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h = h_jet_trk_dphi[iPtZ];
    if (h->GetMaximum () > max)  max = h->GetMaximum ();
    if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
  }
  max = (max > 0 ? 2*max : 1);
  min = (min > 0 ? 0.5*min : 0.1);

  gPad->SetLogy ();
  for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h = h_jet_trk_dphi[iPtZ];

    h->GetYaxis ()->SetRangeUser (min, max);

    h->GetXaxis ()->SetTitle ("Jet-Track #Delta#phi");
    h->GetYaxis ()->SetTitle ("1/N_{Z} dN_{assoc.}/d#Delta#phi");

    h->SetMarkerColor (colors[iPtZ-1]);
    h->SetLineColor (colors[iPtZ-1]);

    h->Draw (iPtZ == 0 ? "hist" : "same hist");
    myText (0.65, 0.9-0.06*(iPtZ-1), colors[iPtZ-1], Form ("%g < #it{p}_{T}^{ Z} < %g [GeV]", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);
  }
  myText (0.3, 0.9, kBlack, "#it{p}_{T}^{ jet} > 5 GeV", 0.04);
  c->SaveAs (Form ("%s/JetTrkCorrelations.pdf", plotPath.Data ()));
}

#endif
