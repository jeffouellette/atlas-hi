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
  MCAnalysis (const char* _name = "mc", const char* subDir = "Nominal", const bool _useHITight = false) : FullAnalysis () {
    name = _name;
    directory = Form ("MCAnalysis/%s/", subDir);
    plotFill = true;
    useAltMarker = false;
    useHITight = _useHITight;
    LoadTrackingEfficiencies ();

    SetupDirectories ("MCAnalysis/", "ZTrackAnalysis/");
    TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
    h_PbPbFCal_weights = (TH1D*) eventWeightsFile->Get ("h_PbPbFCal_weights_mc");
    for (short iCent = 0; iCent < numFinerCentBins; iCent++) {
      h_PbPbQ2_weights[iCent] = (TH1D*) eventWeightsFile->Get (Form ("h_PbPbQ2_weights_iCent%i_mc", iCent));
    }
    h_ppNch_weights = (TH1D*) eventWeightsFile->Get ("h_ppNch_weights_mc");

    SetupDirectories (directory, "ZTrackAnalysis/");
  }

  void Execute () override;
};


////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void MCAnalysis :: Execute () {

  SetupDirectories (directory, "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1, fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_pt = 0, z_y = 0, z_phi = 0, z_m = 0, l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0, vz_weight = 1, q2_weight = 1, fcal_weight = 1, nch_weight = 1;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr, *l_trk_pt = nullptr, *l_trk_eta = nullptr, *l_trk_phi = nullptr;
  //double** trkPtProj = Get2DArray <double> (numPhiBins, nPtTrkBins);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("isEE",      &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",   &fcal_et);
    PbPbTree->SetBranchAddress ("q2",        &q2);
    PbPbTree->SetBranchAddress ("psi2",      &psi2);
    PbPbTree->SetBranchAddress ("vz",        &vz);
    PbPbTree->SetBranchAddress ("z_pt",      &z_pt);
    PbPbTree->SetBranchAddress ("z_y",       &z_y);
    PbPbTree->SetBranchAddress ("z_phi",     &z_phi);
    PbPbTree->SetBranchAddress ("z_m",       &z_m);
    PbPbTree->SetBranchAddress ("l1_pt",     &l1_pt);
    PbPbTree->SetBranchAddress ("l1_eta",    &l1_eta);
    PbPbTree->SetBranchAddress ("l1_phi",    &l1_phi);
    PbPbTree->SetBranchAddress ("l1_charge", &l1_charge);
    PbPbTree->SetBranchAddress ("l2_pt",     &l2_pt);
    PbPbTree->SetBranchAddress ("l2_eta",    &l2_eta);
    PbPbTree->SetBranchAddress ("l2_phi",    &l2_phi);
    PbPbTree->SetBranchAddress ("l2_charge", &l2_charge);
    PbPbTree->SetBranchAddress ("l_trk_pt",  &l_trk_pt);
    PbPbTree->SetBranchAddress ("l_trk_eta", &l_trk_eta);
    PbPbTree->SetBranchAddress ("l_trk_phi", &l_trk_phi);
    PbPbTree->SetBranchAddress ("ntrk",      &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",    &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",   &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",   &trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      //if (fabs (vz) > 1.5)
      //  continue; // vertex cut

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined

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

      short iPtZ = 0; // find z-pt bin
      while (iPtZ < nPtZBins) {
        if (z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
      }

      fcal_weight = h_PbPbFCal_weights->GetBinContent (h_PbPbFCal_weights->FindBin (fcal_et));
      //q2_weight = h_PbPbQ2_weights[iFinerCent]->GetBinContent (h_PbPbQ2_weights[iFinerCent]->FindBin (q2));

      event_weight *= fcal_weight * q2_weight * vz_weight;

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);

      h_q2[iFinerCent]->Fill (q2);
      h_q2_reweighted[iFinerCent]->Fill (q2, event_weight);
      h_PbPb_vz->Fill (vz);
      h_PbPb_vz_reweighted->Fill (vz, event_weight);

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      if (z_pt > zPtBins[1]) {
        int iReg = (fabs (z_y) > 1.00 ? 1 : 0); // barrel vs. endcaps
        h_z_m[iCent][iSpc][iReg]->Fill (z_m, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
        h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
        h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);
        h_z_y_phi[iCent][iSpc]->Fill (z_y, InTwoPi (z_phi), event_weight);

        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);

        for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
          h_lepton_trk_pt[iCent][iSpc]->Fill (l_trk_pt->at (iLTrk), event_weight);
        }
      }

      //for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
      //  for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
      //    trkPtProj[iPhi][iPtTrk] = 0;
      //  }
      //}

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        {
          float mindr = pi;
          //float phidiff = 0;
          float ptdiff = 0;
          for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
            const float dr = DeltaR (trk_eta->at (iTrk), l_trk_eta->at (iLTrk), trk_phi->at (iTrk), l_trk_phi->at (iLTrk));
            if (dr < mindr) {
              mindr = dr;
              ptdiff = 2. * fabs (trkpt - l_trk_pt->at (iLTrk)) / (trkpt + l_trk_pt->at (iLTrk));
              //phidiff = DeltaPhi (trk_phi->at (iTrk), l_trk_phi->at (iLTrk));
            }
          }
          h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, trkpt);
          //h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, phidiff);
        }

        const float zH = trkpt / z_pt;
        const short iZH = GetiZH (zH);
        if (iZH < 0 || iZH > nZHBins-1)
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), true);
        if (trkEff == 0)
          continue;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight / trkEff);

        //// Add to missing pT (requires dphi in +/-pi/2 to +/-pi)
        //float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        //bool awaySide = false;
        //if (dphi > pi/2) {
        //  dphi = pi-dphi;
        //  awaySide = true;
        //}

        //short iPtTrk = 0;
        //while (iPtTrk < nPtTrkBins && trkpt > ptTrkBins[iPtTrk+1])
        //  iPtTrk++;
        //// start at the 1st phi bin and integrate outwards until the track is no longer contained 
        //// e.g. so 7pi/8->pi is a subset of pi/2->pi
        //short iPhi = 0;
        //while (iPhi < numPhiTrkBins && dphi > phiTrkBins[iPhi]) {
        //  if (awaySide)
        //    trkPtProj[iPhi][iPtTrk] += -trkpt * cos (dphi) / trkEff;
        //  else
        //    trkPtProj[iPhi][iPtTrk] += trkpt * cos (dphi) / trkEff;
        //  iPhi++;
        //}

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight / trkEff);
            h_z_trk_xzh[iSpc][iPtZ][idPhi][iCent]->Fill (zH, event_weight / trkEff);
          }
        }

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          if (ptTrkBins[iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, event_weight / trkEff);
        }
        //h_z_trk_pt_phi[iPtZ][iCent][iSpc]->Fill (dphi, trkpt, event_weight / trkEff);
      } // end loop over tracks

      //for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
      //  for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
      //    h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
      //    trkPtProj[iPhi][iPtTrk] = 0;
      //  }
      //}
    } // end loop over Pb+Pb tree
    cout << "Done MC Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("event_weight", &event_weight);
    ppTree->SetBranchAddress ("isEE",      &isEE);
    ppTree->SetBranchAddress ("vz",        &vz);
    ppTree->SetBranchAddress ("z_pt",      &z_pt);
    ppTree->SetBranchAddress ("z_y",       &z_y);
    ppTree->SetBranchAddress ("z_phi",     &z_phi);
    ppTree->SetBranchAddress ("z_m",       &z_m);
    ppTree->SetBranchAddress ("l1_pt",     &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",    &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",    &l1_phi);
    ppTree->SetBranchAddress ("l1_charge", &l1_charge);
    ppTree->SetBranchAddress ("l2_pt",     &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",    &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",    &l2_phi);
    ppTree->SetBranchAddress ("l2_charge", &l2_charge);
    ppTree->SetBranchAddress ("l_trk_pt",  &l_trk_pt);
    ppTree->SetBranchAddress ("l_trk_eta", &l_trk_eta);
    ppTree->SetBranchAddress ("l_trk_phi", &l_trk_phi);
    ppTree->SetBranchAddress ("ntrk",      &ntrk);
    ppTree->SetBranchAddress ("trk_pt",    &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",   &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",   &trk_phi);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      ppTree->GetEntry (iEvt);

      //if (fabs (vz) > 1.5)
      //  continue; // vertex cut

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      short iPtZ = 0; // find z-pt bin
      while (iPtZ < nPtZBins) {
        if (z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
      }

      //nch_weight = h_ppNch_weights->GetBinContent (h_ppNch_weights->FindBin (ntrk));

      event_weight *= vz_weight * nch_weight;

      h_pp_vz->Fill (vz);
      h_pp_vz_reweighted->Fill (vz, event_weight);

      h_pp_nch->Fill (ntrk);
      h_pp_nch_reweighted->Fill (ntrk, event_weight);

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      if (z_pt > zPtBins[1]) {
        int iReg = (fabs (z_y) > 1.00 ? 1 : 0); // barrel vs. endcaps
        h_z_m[iCent][iSpc][iReg]->Fill (z_m, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
        h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l1_phi), event_weight);
        h_z_lepton_dphi[iCent][iSpc]->Fill (DeltaPhi (z_phi, l2_phi), event_weight);
        h_z_y_phi[iCent][iSpc]->Fill (z_y, InTwoPi (z_phi), event_weight);

        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);

        for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
          h_lepton_trk_pt[iCent][iSpc]->Fill (l_trk_pt->at (iLTrk), event_weight);
        }
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        {
          float mindr = pi;
          //float phidiff = 0;
          float ptdiff = 0;
          for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
            const float dr = DeltaR (trk_eta->at (iTrk), l_trk_eta->at (iLTrk), trk_phi->at (iTrk), l_trk_phi->at (iLTrk));
            if (dr < mindr) {
              mindr = dr;
              ptdiff = 2. * fabs (trkpt - l_trk_pt->at (iLTrk)) / (trkpt + l_trk_pt->at (iLTrk));
              //phidiff = DeltaPhi (trk_phi->at (iTrk), l_trk_phi->at (iLTrk));
            }
          }
          h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, trkpt);
          //h_lepton_trk_dr[iCent][iSpc]->Fill (mindr, phidiff);
        }

        const float zH = trkpt / z_pt;
        const short iZH = GetiZH (zH);
        if (iZH < 0 || iZH > nZHBins-1)
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), false);
        if (trkEff == 0)
          continue;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight / trkEff);

        //// Add to missing pT (requires dphi in -pi/2 to pi/2)
        //float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        //bool awaySide = false;
        //if (dphi > pi/2) {
        //  dphi = pi-dphi;
        //  awaySide = true;
        //}

        //short iPtTrk = 0;
        //while (iPtTrk < nPtTrkBins && trkpt > ptTrkBins[iPtTrk+1])
        //  iPtTrk++;
        //// start at the 1st phi bin and integrate outwards until the track is no longer contained 
        //// e.g. so 7pi/8->pi is a subset of pi/2->pi
        //short iPhi = 0;
        //while (iPhi < numPhiTrkBins && dphi > phiTrkBins[iPhi]) {
        //  if (awaySide)
        //    trkPtProj[iPhi][iPtTrk] += -trkpt * cos (dphi) / trkEff;
        //  else
        //    trkPtProj[iPhi][iPtTrk] += trkpt * cos (dphi) / trkEff;
        //  iPhi++;
        //}
        
        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        float dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight / trkEff);
            h_z_trk_xzh[iSpc][iPtZ][idPhi][iCent]->Fill (zH, event_weight / trkEff);
          }
        }

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          if (ptTrkBins[iPtTrk] <= trkpt && trkpt < ptTrkBins[iPtTrk+1])
            h_z_trk_phi[iSpc][iPtZ][iPtTrk][iCent]->Fill (dphi, event_weight / trkEff);
        }
        //h_z_trk_pt_phi[iPtZ][iCent][iSpc]->Fill (dphi, trkpt, event_weight / trkEff);
        
      } // end loop over tracks

      //for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
      //  for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
      //    h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
      //    trkPtProj[iPhi][iPtTrk] = 0;
      //  }
      //}
    } // end loop over pp tree
    cout << "Done MC pp loop." << endl;
  }

  CombineHists ();
  ScaleHists ();

  SaveHists ();
  //LoadHists ();

  inFile->Close ();
  if (inFile) { delete inFile; inFile = nullptr; }

  //Delete2DArray (trkPtProj, numPhiBins, nPtTrkBins);
}

#endif
