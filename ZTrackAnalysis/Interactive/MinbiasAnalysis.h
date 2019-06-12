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
    LoadTrackingEfficiencies ();
    SetupDirectories (directory, "ZTrackAnalysis/");
  }

  void Execute () override;

  void CombineHists () override;
};


////////////////////////////////////////////////////////////////////////////////////////////////
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis::CombineHists () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      h_z_phi[iCent][iSpc]->Add (h_z_phi[iCent][0]);
      h_z_pt[iCent][iSpc]->Add (h_z_pt[iCent][0]);
      h_z_m[iCent][iSpc]->Add (h_z_m[iCent][0]);
      h_trk_pt[iCent][iSpc]->Add (h_trk_pt[iCent][0]);
      //h_lepton_dr[iCent][iSpc]->Add (h_lepton_dr[iCent][0]);

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          if (iSpc == 0 && iPtZ == 0 && iXZTrk == 0)
            continue;
          h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]->Add (h_z_trk_pt_phi[0][0][iCent][0]);
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][iPhi][iCent]->Add (h_z_trk_pt[0][0][0][iPhi][iCent]);
          } // end loop over phi
        }
        for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          h_z_missing_pt[iSpc][iPtZ][iPhi][iCent]->Add (h_z_missing_pt[0][0][iPhi][iCent]);
        } // end loop over phi
        if (iSpc == 0 && iPtZ == 0)
          continue;
        h_z_counts[iSpc][iPtZ][iCent]->Add (h_z_counts[0][0][iCent]);
      } // end loop over pT^Z
    } // end loop over species
  } // end loop over centralities
  return;
}



////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately.
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis::Execute () {
  SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
  h_PbPb_event_reweights = (TH3D*)eventWeightsFile->Get ("h_PbPbEventReweights_minbias");
  h_pp_event_reweights = (TH1D*)eventWeightsFile->Get ("h_ppEventReweights_minbias");

  SetupDirectories (directory, "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_eta = 0, z_phi = 0, event_weight = 1;
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

      event_weight = h_PbPb_event_reweights->GetBinContent (h_PbPb_event_reweights->FindBin (fcal_et, q2, vz));

      h_fcal_et->Fill (fcal_et);
      h_fcal_et_q2->Fill (fcal_et, q2);
      h_fcal_et_reweighted->Fill (fcal_et, event_weight);
      h_fcal_et_q2_reweighted->Fill (fcal_et, q2, event_weight);

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

      event_weight = h_pp_event_reweights->GetBinContent (h_pp_event_reweights->FindBin (vz));

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

  Delete1DArray (trkPtProj, numPhiBins);
}

#endif
