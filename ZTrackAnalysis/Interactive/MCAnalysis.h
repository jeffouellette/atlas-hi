#ifndef __MCAnalysis_h__
#define __MCAnalysis_h__

#include "Params.h"
#include "FullAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class MCAnalysis : public FullAnalysis {

  public:

  bool takeNonTruthTracks = false;

  MCAnalysis (const char* _name = "mc") : FullAnalysis () {
    name = _name;
    plotFill = false;
    useAltMarker = false;
    isMC = true;
    eventWeightsFileName = "MCAnalysis/Nominal/eventWeightsFile.root";
  }

  void Execute (const char* inFileName, const char* outFileName) override;
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void MCAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  LoadEventWeights ();

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  if (!inFile || !inFile->IsOpen ()) {
    cout << "Error in Execute: File not open!" << endl;
    return;
  }
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1;//, fcal_weight = 1;//, q2_weight = 1, psi2_weight = 1;//, vz_weight = 1, nch_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0;
  float z_pt = 0, z_eta = 0, z_y = 0, z_phi = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  float l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];
  bool trk_truth_matched[10000];

  float** trk_counts = Get2DArray <float> (2, 6);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("isEE",       &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",    &fcal_et);
    PbPbTree->SetBranchAddress ("q2",         &q2);
    PbPbTree->SetBranchAddress ("psi2",       &psi2);
    PbPbTree->SetBranchAddress ("vz",         &vz);
    PbPbTree->SetBranchAddress ("z_pt",       &z_pt);
    PbPbTree->SetBranchAddress ("z_y",        &z_y);
    PbPbTree->SetBranchAddress ("z_phi",      &z_phi);
    PbPbTree->SetBranchAddress ("z_m",        &z_m);
    PbPbTree->SetBranchAddress ("l1_pt",      &l1_pt);
    PbPbTree->SetBranchAddress ("l1_eta",     &l1_eta);
    PbPbTree->SetBranchAddress ("l1_phi",     &l1_phi);
    PbPbTree->SetBranchAddress ("l1_charge",  &l1_charge);
    PbPbTree->SetBranchAddress ("l1_trk_pt",  &l1_trk_pt);
    PbPbTree->SetBranchAddress ("l1_trk_eta", &l1_trk_eta);
    PbPbTree->SetBranchAddress ("l1_trk_phi", &l1_trk_phi);
    PbPbTree->SetBranchAddress ("l2_pt",      &l2_pt);
    PbPbTree->SetBranchAddress ("l2_eta",     &l2_eta);
    PbPbTree->SetBranchAddress ("l2_phi",     &l2_phi);
    PbPbTree->SetBranchAddress ("l2_charge",  &l2_charge);
    PbPbTree->SetBranchAddress ("l2_trk_pt",  &l2_trk_pt);
    PbPbTree->SetBranchAddress ("l2_trk_eta", &l2_trk_eta);
    PbPbTree->SetBranchAddress ("l2_trk_phi", &l2_trk_phi);
    PbPbTree->SetBranchAddress ("ntrk",       &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",     trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",    trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",    trk_phi);
    PbPbTree->SetBranchAddress ("trk_truth_matched", trk_truth_matched);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      //if (fabs (vz) > 150)
      //  continue; // vertex cut

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

      const short iCent = GetCentBin (fcal_et);
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      const short iFineCent = GetFineCentBin (fcal_et);
      if (iFineCent < 1 || iFineCent > numFineCentBins-1)
        continue;

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;

      //{
      //  float dphi = DeltaPhi (z_phi, psi2, false);
      //  if (dphi > pi/2)
      //    dphi = pi - dphi;
      //  fcal_weight = h_PbPbFCal_weights[iSpc][iPtZ]->GetBinContent (h_PbPbFCal_weights[iSpc][iPtZ]->FindBin (fcal_et));
      //  //q2_weight = h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->FindBin (q2));
      //  //psi2_weight = h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->FindBin (dphi));

      //  event_weight *= fcal_weight;// * psi2_weight;
      //}

      if (event_weight == 0)
        continue;

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);

      h_q2[iFineCent]->Fill (q2);
      h_q2_reweighted[iFineCent]->Fill (q2, event_weight);
      h_psi2[iFineCent]->Fill (psi2);
      h_psi2_reweighted[iFineCent]->Fill (psi2, event_weight);
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

      h_lepton_trk_pt[iCent][iSpc]->Fill (l1_trk_pt, event_weight);
      h_lepton_trk_pt[iCent][iSpc]->Fill (l2_trk_pt, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (takeNonTruthTracks && trk_truth_matched[iTrk])
          continue;

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.03 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.03))
          continue;

        {
          float mindr = pi;
          float ptdiff = 0;
          float dr = DeltaR (trk_eta[iTrk], l1_trk_eta, trk_phi[iTrk], l1_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l1_trk_pt) / (trkpt + l1_trk_pt);
          }
          dr = DeltaR (trk_eta[iTrk], l2_trk_eta, trk_phi[iTrk], l2_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l2_trk_pt) / (trkpt + l2_trk_pt);
          }
          h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, ptdiff);
        }

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], true);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], true);
        if (trkEff == 0 || trkPur == 0)
          continue;
        const float trkWeight = trkPur / trkEff;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight * trkWeight);

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, event_weight * trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight * trkWeight);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (xhz, event_weight * trkWeight);
          }
        }

        if (3*pi/4 <= dphi) {
          if (trkpt < 60) {
            if (30 <= trkpt)      trk_counts[0][5] += trkWeight;
            else if (15 <= trkpt) trk_counts[0][4] += trkWeight;
            else if (8 <= trkpt)  trk_counts[0][3] += trkWeight;
            else if (4 <= trkpt)  trk_counts[0][2] += trkWeight;
            else if (2 <= trkpt)  trk_counts[0][1] += trkWeight;
            else if (1 <= trkpt)  trk_counts[0][0] += trkWeight;
          }
 
          if (1./60. <= xhz) {
            if (xhz <= 1./30.)      trk_counts[1][0] += trkWeight;
            else if (xhz <= 1./15.) trk_counts[1][1] += trkWeight;
            else if (xhz <= 1./8.)  trk_counts[1][2] += trkWeight;
            else if (xhz <= 1./4.)  trk_counts[1][3] += trkWeight;
            else if (xhz <= 1./2.)  trk_counts[1][4] += trkWeight;
            else if (xhz <= 1.)     trk_counts[1][5] += trkWeight;
          }
        }
      } // end loop over tracks

      // fill yield histograms and covariance matrices
      for (int i1 = 0; i1 < nPtTrkBins[iPtZ]; i1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->Fill (0.5*(ptTrkBins[iPtZ][i1]+ptTrkBins[iPtZ][i1+1]), event_weight*trk_counts[0][i1]);
        for (int i2 = 0 ; i2 < nPtTrkBins[iPtZ]; i2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->Fill (0.5*(ptTrkBins[iPtZ][i1]+ptTrkBins[iPtZ][i1+1]), 0.5*(ptTrkBins[iPtZ][i2]+ptTrkBins[iPtZ][i2+1]), event_weight*(trk_counts[0][i1])*(trk_counts[0][i2]));
      } // end loop over i1
      for (int i1 = 0; i1 < nXHZBins[iPtZ]; i1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Fill (0.5*(xHZBins[iPtZ][i1]+xHZBins[iPtZ][i1+1]), event_weight*trk_counts[1][i1]);
        for (int i2 = 0 ; i2 < nXHZBins[iPtZ]; i2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->Fill (0.5*(xHZBins[iPtZ][i1]+xHZBins[iPtZ][i1+1]), 0.5*(xHZBins[iPtZ][i2]+xHZBins[iPtZ][i2+1]), event_weight*(trk_counts[1][i1])*(trk_counts[1][i2]));
      } // end loop over i1

      // reset trk count measurements for next event
      for (int i = 0; i < 6; i++) {
        trk_counts[0][i] = 0;
        trk_counts[1][i] = 0;
      } // end loop over i

    } // end loop over Pb+Pb tree
    cout << "Done MC Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("event_weight", &event_weight);
    ppTree->SetBranchAddress ("isEE",       &isEE);
    ppTree->SetBranchAddress ("vz",         &vz);
    ppTree->SetBranchAddress ("z_pt",       &z_pt);
    ppTree->SetBranchAddress ("z_y",        &z_y);
    ppTree->SetBranchAddress ("z_phi",      &z_phi);
    ppTree->SetBranchAddress ("z_m",        &z_m);
    ppTree->SetBranchAddress ("l1_pt",      &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",     &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",     &l1_phi);
    ppTree->SetBranchAddress ("l1_charge",  &l1_charge);
    ppTree->SetBranchAddress ("l1_trk_pt",  &l1_trk_pt);
    ppTree->SetBranchAddress ("l1_trk_eta", &l1_trk_eta);
    ppTree->SetBranchAddress ("l1_trk_phi", &l1_trk_phi);
    ppTree->SetBranchAddress ("l2_pt",      &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",     &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",     &l2_phi);
    ppTree->SetBranchAddress ("l2_charge",  &l2_charge);
    ppTree->SetBranchAddress ("l2_trk_pt",  &l2_trk_pt);
    ppTree->SetBranchAddress ("l2_trk_eta", &l2_trk_eta);
    ppTree->SetBranchAddress ("l2_trk_phi", &l2_trk_phi);
    ppTree->SetBranchAddress ("ntrk",       &ntrk);
    ppTree->SetBranchAddress ("trk_pt",     trk_pt);
    ppTree->SetBranchAddress ("trk_eta",    trk_eta);
    ppTree->SetBranchAddress ("trk_phi",    trk_phi);
    ppTree->SetBranchAddress ("trk_truth_matched", trk_truth_matched);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      //if (fabs (vz) > 150)
      //  continue; // vertex cut

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1)
        continue;

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

      h_lepton_trk_pt[iCent][iSpc]->Fill (l1_trk_pt, event_weight);
      h_lepton_trk_pt[iCent][iSpc]->Fill (l2_trk_pt, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (takeNonTruthTracks && trk_truth_matched[iTrk])
          continue;

        if (trkpt < trk_min_pt)// || trk_max_pt < trkpt)
          continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.03 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.03))
          continue;

        {
          float mindr = pi;
          float ptdiff = 0;
          float dr = DeltaR (trk_eta[iTrk], l1_trk_eta, trk_phi[iTrk], l1_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l1_trk_pt) / (trkpt + l1_trk_pt);
          }
          dr = DeltaR (trk_eta[iTrk], l2_trk_eta, trk_phi[iTrk], l2_trk_phi);
          if (dr < mindr) {
            mindr = dr;
            ptdiff = 2. * fabs (trkpt - l2_trk_pt) / (trkpt + l2_trk_pt);
          }
          h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, ptdiff);
        }

        const double trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta[iTrk], false);
        const double trkPur = GetTrackingPurity (fcal_et, trkpt, trk_eta[iTrk], false);
        if (trkEff == 0 || trkPur == 0)
          continue;
        const float trkWeight = trkPur / trkEff;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight * trkWeight);

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        float dphi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins[iPtZ]; iPtTrk++) {
          if (ptTrkBins[iPtZ][iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtZ][iPtTrk+1])
            h_trk_dphi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, event_weight * trkWeight);
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt);
            h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight * trkWeight);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (xhz, event_weight * trkWeight);
          }
        }

        if (3*pi/4 <= dphi) {
          if (trkpt < 60) {
            if (30 <= trkpt)      trk_counts[0][5] += trkWeight;
            else if (15 <= trkpt) trk_counts[0][4] += trkWeight;
            else if (8 <= trkpt)  trk_counts[0][3] += trkWeight;
            else if (4 <= trkpt)  trk_counts[0][2] += trkWeight;
            else if (2 <= trkpt)  trk_counts[0][1] += trkWeight;
            else if (1 <= trkpt)  trk_counts[0][0] += trkWeight;
          }
 
          if (1./60. <= xhz) {
            if (xhz <= 1./30.)      trk_counts[1][0] += trkWeight;
            else if (xhz <= 1./15.) trk_counts[1][1] += trkWeight;
            else if (xhz <= 1./8.)  trk_counts[1][2] += trkWeight;
            else if (xhz <= 1./4.)  trk_counts[1][3] += trkWeight;
            else if (xhz <= 1./2.)  trk_counts[1][4] += trkWeight;
            else if (xhz <= 1.)     trk_counts[1][5] += trkWeight;
          }
        }
      } // end loop over tracks

      // fill yield histograms and covariance matrices
      for (int i1 = 0; i1 < nPtTrkBins[iPtZ]; i1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->Fill (0.5*(ptTrkBins[iPtZ][i1]+ptTrkBins[iPtZ][i1+1]), event_weight*trk_counts[0][i1]);
        for (int i2 = 0 ; i2 < nPtTrkBins[iPtZ]; i2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->Fill (0.5*(ptTrkBins[iPtZ][i1]+ptTrkBins[iPtZ][i1+1]), 0.5*(ptTrkBins[iPtZ][i2]+ptTrkBins[iPtZ][i2+1]), event_weight*(trk_counts[0][i1])*(trk_counts[0][i2]));
      } // end loop over i1
      for (int i1 = 0; i1 < nXHZBins[iPtZ]; i1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->Fill (0.5*(xHZBins[iPtZ][i1]+xHZBins[iPtZ][i1+1]), event_weight*trk_counts[1][i1]);
        for (int i2 = 0 ; i2 < nXHZBins[iPtZ]; i2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->Fill (0.5*(xHZBins[iPtZ][i1]+xHZBins[iPtZ][i1+1]), 0.5*(xHZBins[iPtZ][i2]+xHZBins[iPtZ][i2+1]), event_weight*(trk_counts[1][i1])*(trk_counts[1][i2]));
      } // end loop over i1

      // reset trk count measurements for next event
      for (int i = 0; i < 6; i++) {
        trk_counts[0][i] = 0;
        trk_counts[1][i] = 0;
      } // end loop over i

    } // end loop over pp tree
    cout << "Done MC pp loop." << endl;
  }

  Delete2DArray (trk_counts, 2, 6);

  SaveHists (outFileName);

  if (inFile) inFile->Close ();
  SaferDelete (inFile);
}


#endif
