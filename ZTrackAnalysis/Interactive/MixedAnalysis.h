#ifndef __MixedAnalysis_h__
#define __MixedAnalysis_h__
#include "ZTrackUtilities.h"

#include <ArrayTemplates.h>

#include <AtlasUtils.h>

#include <TEfficiency.h>
#include <TClass.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include <iostream>
#include <string>

using namespace atlashi;
using namespace std;

class MixedAnalysis : public PhysicsAnalysis {

  public:

  MixedAnalysis () { }

  MixedAnalysis (const char* _name, const char* subDir, const bool _useHITight = false) : PhysicsAnalysis (_name, subDir, _useHITight) { }

  public:
  virtual void Execute () override;

};





////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
// Designed to be overloaded. The default here is for analyzing data.
////////////////////////////////////////////////////////////////////////////////////////////////
void MixedAnalysis :: Execute () {
  SetupDirectories (directory.c_str (), "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1, fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_pt = 0, z_phi = 0;
  int ntrk = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;
  //double** trkPtProj = Get2DArray <double> (numPhiBins, nPtTrkBins);

  float rand_z_pt = 0, rand_z_phi = 0;

  TRandom3* rando = new TRandom3 ();
  rando->SetSeed ();

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("isEE",      &isEE);
    PbPbTree->SetBranchAddress ("fcal_et",   &fcal_et);
    PbPbTree->SetBranchAddress ("q2",        &q2);
    PbPbTree->SetBranchAddress ("psi2",      &psi2);
    PbPbTree->SetBranchAddress ("vz",        &vz);
    PbPbTree->SetBranchAddress ("z_pt",      &z_pt);
    PbPbTree->SetBranchAddress ("z_phi",     &z_phi);
    PbPbTree->SetBranchAddress ("ntrk",      &ntrk);
    PbPbTree->SetBranchAddress ("trk_pt",    &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",   &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",   &trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;

      PbPbTree->GetEntry ((int) (rando->Rndm () * nEvts));

      rand_z_pt = z_pt;
      rand_z_phi = z_phi;

      PbPbTree->GetEntry (iEvt);

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

      short iFineCent = 0;
      while (iFineCent < numFineCentBins) {
        if (fcal_et < fineCentBins[iFineCent])
          break;
        else
          iFineCent++;
      }
      if (iFineCent < 1 || iFineCent > numFineCentBins-1)
        continue;

      short iPtZ = 0; // find z-pt bin
      while (iPtZ < nPtZBins) {
        if (rand_z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
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

        const float zH = trkpt / rand_z_pt;
        const short iZH = GetiZH (zH);
        if (iZH < 0 || iZH > nZHBins-1)
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), true);
        if (trkEff == 0)
          continue;

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
        float dphi = DeltaPhi (rand_z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight / trkEff);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (zH, event_weight / trkEff);
          }
        }

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (rand_z_phi, trk_phi->at (iTrk), true);
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
    cout << "Done primary Pb+Pb loop." << endl;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (ppTree) {
    ppTree->SetBranchAddress ("isEE",      &isEE);
    ppTree->SetBranchAddress ("vz",        &vz);
    ppTree->SetBranchAddress ("z_pt",      &z_pt);
    ppTree->SetBranchAddress ("z_phi",     &z_phi);
    ppTree->SetBranchAddress ("ntrk",      &ntrk);
    ppTree->SetBranchAddress ("trk_pt",    &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",   &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",   &trk_phi);

    const int nEvts = ppTree->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;

      ppTree->GetEntry ((int) (rando->Rndm () * nEvts));

      rand_z_pt = z_pt;
      rand_z_phi = z_phi;

      ppTree->GetEntry (iEvt);

      const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
      const short iCent = 0; // iCent = 0 for pp

      short iPtZ = 0; // find z-pt bin
      while (iPtZ < nPtZBins) {
        if (rand_z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
      }

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      h_z_counts[iSpc][iPtZ][iCent]->Fill (1.5);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        const float zH = trkpt / rand_z_pt;
        const short iZH = GetiZH (zH);
        if (iZH < 0 || iZH > nZHBins-1)
          continue;

        const float trkEff = GetTrackingEfficiency (fcal_et, trkpt, trk_eta->at (iTrk), false);
        if (trkEff == 0)
          continue;

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
        float dphi = DeltaPhi (rand_z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++) {
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi]) {
            h_z_trk_raw_pt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight / trkEff);
            h_trk_xhz_dphi[iSpc][iPtZ][idPhi][iCent]->Fill (zH, event_weight / trkEff);
          }
        }

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (rand_z_phi, trk_phi->at (iTrk), true);
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
    cout << "Done primary pp loop." << endl;
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
