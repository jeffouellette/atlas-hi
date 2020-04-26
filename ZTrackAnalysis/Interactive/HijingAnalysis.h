#ifndef __HijingAnalysis_h__
#define __HijingAnalysis_h__

#include "Params.h"
#include "FullAnalysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class HijingAnalysis : public FullAnalysis {

  public:

  bool takeNonTruthTracks = false;

  HijingAnalysis (const char* _name = "mc") : FullAnalysis () {
    name = _name;
    plotFill = false;
    useAltMarker = false;
    isMC = true;
    useHijingEffs = true;
    eventWeightsFileName = "MCAnalysis/Hijing/eventWeightsFile.root";
  }

  void Execute (const char* inFileName, const char* outFileName) override;
};




////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void HijingAnalysis :: Execute (const char* inFileName, const char* outFileName) {

  cout << "Arguments provided: " << endl;
  cout << "inFileName = " << inFileName << endl;
  cout << "outFileName = " << outFileName << endl;

  LoadEventWeights ();
  //eventPlaneCalibrator = EventPlaneCalibrator (Form ("%s/FCalCalibration/Nominal/data18hi.root", rootPath.Data ()));

  SetupDirectories ("", "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/%s", rootPath.Data (), inFileName), "read");
  if (!inFile || !inFile->IsOpen ()) {
    cout << "Error in Execute: File not open!" << endl;
    return;
  }
  cout << "Read input file from " << Form ("%s/%s", rootPath.Data (), inFileName) << endl;

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbMixedTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1, fcal_weight = 1, q2_weight = 1, psi2_weight = 1;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0;
  float z_pt = 0, z_eta = 0, z_y = 0, z_phi = 0, z_m = 0;
  float l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  float l1_trk_pt = 0, l1_trk_eta = 0, l1_trk_phi = 0, l2_trk_pt = 0, l2_trk_eta = 0, l2_trk_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  float trk_pt[10000], trk_eta[10000], trk_phi[10000];

  int***    trks_counts   = Get3DArray <int> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights1 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  float***  trks_weights2 = Get3DArray <float> (2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  int**     trks_counts_inPhi   = Get2DArray <int> (maxNPtchBins, 40);
  float**   trks_weights1_inPhi = Get2DArray <float> (maxNPtchBins, 40);
  float**   trks_weights2_inPhi = Get2DArray <float> (maxNPtchBins, 40);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("z_event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("isEE",       &isEE);
    PbPbTree->SetBranchAddress ("z_fcal_et",  &fcal_et);
    PbPbTree->SetBranchAddress ("z_q2",       &q2);
    PbPbTree->SetBranchAddress ("z_psi2",     &psi2);
    PbPbTree->SetBranchAddress ("z_vz",       &vz);
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
    PbPbTree->SetBranchAddress ("z_ntrk",     &ntrk);
    PbPbTree->SetBranchAddress ("z_trk_pt",   trk_pt);
    PbPbTree->SetBranchAddress ("z_trk_eta",  trk_eta);
    PbPbTree->SetBranchAddress ("z_trk_phi",  trk_phi);

    const int nEvts = PbPbTree->GetEntries () / mixingFraction;
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      if (fabs (vz) > 150) continue; // vertex cut

      //{
      //  CorrectQ2Vector (q2x_a, q2y_a, q2x_c, q2y_c);
      //  const float q2x = q2x_a + q2x_c;
      //  const float q2y = q2y_a + q2y_c;
      //  q2 = sqrt (q2x*q2x + q2y*q2y) / fcal_et;
      //  psi2 = 0.5 * atan2 (q2y, q2x);
      //}

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

      const short iCent = GetCentBin (fcal_et);
      if (iCent < 1 || iCent > numCentBins-1) continue;

      const short iFineCent = GetFineCentBin (fcal_et);
      if (iFineCent < 1 || iFineCent > numFineCentBins-1) continue;

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

      // do a reweighting procedure
      {
        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)  dphi = pi - dphi;
        if (useCentWgts)  fcal_weight = h_PbPbFCal_weights[iSpc][iPtZ]->GetBinContent (h_PbPbFCal_weights[iSpc][iPtZ]->FindBin (fcal_et));
        if (useQ2Wgts)    q2_weight   = h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbQ2_weights[iSpc][iFineCent][iPtZ]->FindBin (q2));
        if (usePsi2Wgts)  psi2_weight = h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->GetBinContent (h_PbPbPsi2_weights[iSpc][iFineCent][iPtZ]->FindBin (dphi));

        event_weight *= fcal_weight * q2_weight * psi2_weight;
      }

      if (event_weight == 0) continue;

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
        if (dphi > pi/2) dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

      h_lepton_trk_pt[iCent][iSpc]->Fill (l1_trk_pt, event_weight);
      h_lepton_trk_pt[iCent][iSpc]->Fill (l2_trk_pt, event_weight);

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (2.5, pow (event_weight, 2));

      if (iPtZ < 2) continue;

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);

      h_q2[iFineCent]->Fill (q2);
      h_q2_reweighted[iFineCent]->Fill (q2, event_weight);
      h_psi2[iFineCent]->Fill (psi2);
      h_psi2_reweighted[iFineCent]->Fill (psi2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt[iTrk];
        const float xhz = trkpt / z_pt;

        if (trkpt < trk_min_pt) continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.03 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.03)) continue;

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
        if (trkEff == 0 || trkPur == 0) continue;
        const float trkWeight = trkPur / trkEff;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight * trkWeight);

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (allPtchBins[0] <= trkpt) {
          iPtch = 0;
          while (iPtch < maxNPtchBins && allPtchBins[iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < maxNPtchBins) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (iPtch, iCent) && (-pi/2.)+(2.*pi/GetNdPhiBins (iPtch, iCent))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (allXhZBins[0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < maxNXhZBins && allXhZBins[iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < maxNPtchBins) {
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < maxNXhZBins) {
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if (3*pi/4 <= dPhi) {
          if (iPtch != -1 && iPtch < maxNPtchBins) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < maxNXhZBins) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks

      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
        for (int idPhi1 = 0; idPhi1 < GetNdPhiBins (iPtch, iCent); idPhi1++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi1+1, sqrt (pow (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinError (idPhi1+1), 2) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi1])));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (iPtch, iCent); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi1

      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < max (maxNPtchBins, maxNXhZBins); j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < maxNPtchBins; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

    } // end loop over Pb+Pb tree
    cout << "Done Hijing Pb+Pb loop." << endl;
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

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 100 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      if (fabs (vz) > 150) continue; // vertex cut

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      const short iPtZ = GetPtZBin (z_pt); // find z-pt bin
      if (iPtZ < 0 || iPtZ > nPtZBins-1) continue;

      if (event_weight == 0) continue;

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
        if (dphi > pi/2) dphi = pi - dphi;
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

        if (trkpt < trk_min_pt) continue;

        if (doLeptonRejVar && (DeltaR (l1_trk_eta, trk_eta[iTrk], l1_trk_phi, trk_phi[iTrk]) < 0.03 || DeltaR (l2_trk_eta, trk_eta[iTrk], l2_trk_phi, trk_phi[iTrk]) < 0.03)) continue;

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
        if (trkEff == 0 || trkPur == 0) continue;
        const float trkWeight = trkPur / trkEff;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight * trkWeight);

        // Study correlations (requires dPhi in -pi/2 to 3pi/2)
        float dPhi = DeltaPhi (z_phi, trk_phi[iTrk], true);
        if (dPhi < -pi/2) dPhi = dPhi + 2*pi;

        short iPtch = -1;
        if (allPtchBins[0] <= trkpt) {
          iPtch = 0;
          while (iPtch < maxNPtchBins && allPtchBins[iPtch+1] < trkpt) iPtch++;
        }

        if (iPtch != -1 && iPtch < maxNPtchBins) {
          short idPhi = 0;
          while (idPhi < GetNdPhiBins (iPtch, iCent) && (-pi/2.)+(2.*pi/GetNdPhiBins (iPtch, iCent))*(idPhi+1) < dPhi) idPhi++;

          trks_counts_inPhi[iPtch][idPhi]   += 1;
          trks_weights1_inPhi[iPtch][idPhi] += trkWeight;
          trks_weights2_inPhi[iPtch][idPhi] += pow (trkWeight, 2);
        }

        short iXhZ = -1;
        if (allXhZBins[0] <= xhz) {
          iXhZ = 0;
          while (iXhZ < maxNXhZBins && allXhZBins[iXhZ+1] < xhz) iXhZ++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dPhi = DeltaPhi (z_phi, trk_phi[iTrk], false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dPhi && dPhi <= phiHighBins[idPhi]) {
            if (iPtch != -1 && iPtch < maxNPtchBins) {
              trks_counts[0][iPtch][idPhi]    += 1;
              trks_weights1[0][iPtch][idPhi]  += trkWeight;
              trks_weights2[0][iPtch][idPhi]  += pow (trkWeight, 2);
            }
            if (iXhZ != -1 && iXhZ < maxNXhZBins) {
              trks_counts[1][iXhZ][idPhi]   += 1;
              trks_weights1[1][iXhZ][idPhi] += trkWeight;
              trks_weights2[1][iXhZ][idPhi] += pow (trkWeight, 2);
            }
          }
        } // end loop over idPhi
        if (3*pi/4 <= dPhi) {
          if (iPtch != -1 && iPtch < maxNPtchBins) {
            trks_counts[0][iPtch][numPhiBins]   += 1;
            trks_weights1[0][iPtch][numPhiBins] += trkWeight;
            trks_weights2[0][iPtch][numPhiBins] += pow (trkWeight, 2);
          }
          if (iXhZ != -1 && iXhZ < maxNXhZBins) {
            trks_counts[1][iXhZ][numPhiBins]    += 1;
            trks_weights1[1][iXhZ][numPhiBins]  += trkWeight;
            trks_weights2[1][iXhZ][numPhiBins]  += pow (trkWeight, 2);
          }
        }
      } // end loop over tracks

      // fill phi correlation histograms and covariance matrices
      for (int iPtch = 0; iPtch < maxNPtchBins; iPtch++) {
        for (int idPhi1 = 0; idPhi1 < GetNdPhiBins (iPtch, iCent); idPhi1++) {
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]));
          h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->SetBinError (idPhi1+1, sqrt (pow (h_trk_dphi[iSpc][iPtZ][iPtch][iCent]->GetBinError (idPhi1+1), 2) + pow (event_weight, 2) * (trks_weights2_inPhi[iPtch][idPhi1])));
          for (int idPhi2 = 0; idPhi2 < GetNdPhiBins (iPtch, iCent); idPhi2++)
            h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->SetBinContent (idPhi1+1, idPhi2+1, h2_trk_dphi_cov[iSpc][iPtZ][iPtch][iCent]->GetBinContent (idPhi1+1, idPhi2+1) + (event_weight) * (trks_weights1_inPhi[iPtch][idPhi1]) * (trks_weights1_inPhi[iPtch][idPhi2]));
        } // end loop over iPtch
      } // end loop over idPhi1

      // fill yield histograms binned in dPhi and covariance matrices
      for (int idPhi = 0; idPhi < numPhiBins; idPhi++) {
        for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
          h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi_raw[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + trks_counts[0][iPtch1][idPhi]);
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]));
          h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][idPhi])));
          for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
            h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][idPhi]) * (trks_weights1[0][iPtch2][idPhi]));
        } // end loop over iPtch1
        for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]));
          h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][idPhi])));
          for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
            h2_trk_xhz_dphi_cov[iSpc][iPtZ][idPhi][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_pt_dphi_cov[iSpc][iPtZ][idPhi][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][idPhi]) * (trks_weights1[1][iXhZ2][idPhi]));
        } // end loop over iXhZ1
      } // end loop over idPhi

      // fill yield histograms and covariance matrices (for dPhi integrated yield)
      for (int iPtch1 = 0; iPtch1 < maxNPtchBins; iPtch1++) {
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]));
        h_trk_pt_ptz[iSpc][iPtZ][iCent]->SetBinError   (iPtch1+1, sqrt (pow (h_trk_pt_ptz[iSpc][iPtZ][iCent]->GetBinError (iPtch1+1), 2) + pow (event_weight, 2) * (trks_weights2[0][iPtch1][numPhiBins])));
        for (int iPtch2 = 0; iPtch2 < maxNPtchBins; iPtch2++)
          h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iPtch1+1, iPtch2+1, h2_trk_pt_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iPtch1+1, iPtch2+1) + (event_weight) * (trks_weights1[0][iPtch1][numPhiBins]) * (trks_weights1[0][iPtch2][numPhiBins]));
      } // end loop over iPtch1
      for (int iXhZ1 = 0; iXhZ1 < maxNXhZBins; iXhZ1++) {
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]));
        h_trk_xhz_ptz[iSpc][iPtZ][iCent]->SetBinError   (iXhZ1+1, sqrt (pow (h_trk_xhz_ptz[iSpc][iPtZ][iCent]->GetBinError (iXhZ1+1), 2) + pow (event_weight, 2) * (trks_weights2[1][iXhZ1][numPhiBins])));
        for (int iXhZ2 = 0; iXhZ2 < maxNXhZBins; iXhZ2++)
          h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->SetBinContent (iXhZ1+1, iXhZ2+1, h2_trk_xhz_ptz_cov[iSpc][iPtZ][iCent]->GetBinContent (iXhZ1+1, iXhZ2+1) + (event_weight) * (trks_weights1[1][iXhZ1][numPhiBins]) * (trks_weights1[1][iXhZ2][numPhiBins]));
      } // end loop over iXhZ1

      // reset trk count measurements for next event
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < max (maxNPtchBins, maxNXhZBins); j++) {
          for (int k = 0; k <= numPhiBins; k++) {
            trks_counts[i][j][k] = 0;
            trks_weights1[i][j][k] = 0;
            trks_weights2[i][j][k] = 0;
          } // end loop over k
        } // end loop over j
      } // end loop over i
      for (int i = 0; i < maxNPtchBins; i++) {
        for (int j = 0; j < 40; j++) {
          trks_counts_inPhi[i][j] = 0;
          trks_weights1_inPhi[i][j] = 0;
          trks_weights2_inPhi[i][j] = 0;
        } // end loop over j
      } // end loop over i

    } // end loop over pp tree
    cout << "Done Hijing pp loop." << endl;
  }

  Delete3DArray (&trks_counts, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights1, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete3DArray (&trks_weights2, 2, max (maxNPtchBins, maxNXhZBins), numPhiBins+1);
  Delete2DArray (&trks_counts_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights1_inPhi, maxNPtchBins, 40);
  Delete2DArray (&trks_weights2_inPhi, maxNPtchBins, 40);

  SaveHists (outFileName);

  if (inFile) inFile->Close ();
  SaferDelete (&inFile);
}


#endif
