#ifndef __MCAnalysis_h__
#define __MCAnalysis_h__

#include "Params.h"
#include "Analysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class MCAnalysis : public Analysis {

  public:
  MCAnalysis () : Analysis () {
    name = "mc";
    directory = "MCAnalysis/";
    plotFill = true;
    useAltMarker = false;
    SetupDirectories (directory, "ZTrackAnalysis/");
  }

  void Execute ();
};


////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void MCAnalysis::Execute () {
  SetupDirectories (directory, "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");

  PbPbEventReweights = (TH3D*)eventWeightsFile->Get ("PbPbEventReweights_mc");
  ppEventReweights = (TH1D*)eventWeightsFile->Get ("ppEventReweights_mc");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1, fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_pt = 0, z_eta = 0, z_phi = 0, z_m = 0, l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;//, *l_trk_pt = nullptr, *l_trk_eta = nullptr, *l_trk_phi = nullptr;
  double** trkPtProj = Get2DArray <double> (numPhiBins, nPtTrkBins);


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
    PbPbTree->SetBranchAddress ("z_eta",     &z_eta);
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

      short iPtZ = 0; // find z-pt bin
      while (iPtZ < nPtZBins) {
        if (z_pt < zPtBins[iPtZ+1])
          break;
        else
          iPtZ++;
      }

      event_weight = PbPbEventReweights->GetBinContent (PbPbEventReweights->FindBin (fcal_et, q2, vz));

      FCalSpec->Fill (fcal_et, event_weight);
      FCalQ2Corr->Fill (fcal_et, q2, event_weight);

      ZPtSpecs[iCent][iSpc]->Fill (z_pt, event_weight);
      if (z_pt > zPtBins[1]) {
        ZMYields[iCent][iSpc]->Fill (z_m, event_weight);
        LeptonSpec[iCent][iSpc]->Fill (l1_pt, event_weight);
        LeptonSpec[iCent][iSpc]->Fill (l2_pt, event_weight);

        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        ZPhiYields[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

      for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
          trkPtProj[iPhi][iPtTrk] = 0;
        }
      }

      ZCounts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        const float xZTrk = trkpt / z_pt;
        const short iXZTrk = GetiXZTrk (xZTrk);
        if (iXZTrk < 0 || iXZTrk > nXZTrkBins-1)
          continue;

        TrackSpec[iCent][iSpc]->Fill (trkpt, event_weight);

        // Add to missing pT (requires dphi in +/-pi/2 to +/-pi)
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
            trkPtProj[iPhi][iPtTrk] += -trkpt * cos (dphi);
          else
            trkPtProj[iPhi][iPtTrk] += trkpt * cos (dphi);
          iPhi++;
        }

        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++)
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi])
            ZTracksPt[iSpc][iPtZ][iXZTrk][idPhi][iCent]->Fill (trkpt, event_weight);

        //// Study correlations (requires dphi in -pi/2 to 3pi/2)
        //dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        //if (dphi < -pi/2)
        //  dphi = dphi + 2*pi;

        ZTrackPtPhi[iPtZ][iXZTrk][iCent][iSpc]->Fill (dphi, trkpt, event_weight);
      } // end loop over tracks

      for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          ZMissingPt[iSpc][iPtZ][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
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
    ppTree->SetBranchAddress ("isEE",      &isEE);
    ppTree->SetBranchAddress ("vz",        &vz);
    ppTree->SetBranchAddress ("z_pt",      &z_pt);
    ppTree->SetBranchAddress ("z_eta",     &z_eta);
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

      event_weight = ppEventReweights->GetBinContent (ppEventReweights->FindBin (vz));

      ZPtSpecs[iCent][iSpc]->Fill (z_pt, event_weight);
      if (z_pt > zPtBins[1]) {
        ZMYields[iCent][iSpc]->Fill (z_m, event_weight);
        LeptonSpec[iCent][iSpc]->Fill (l1_pt, event_weight);
        LeptonSpec[iCent][iSpc]->Fill (l2_pt, event_weight);
      }

      ZCounts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        const float xZTrk = trkpt / z_pt;
        const short iXZTrk = GetiXZTrk (xZTrk);
        if (iXZTrk < 0 || iXZTrk > nXZTrkBins-1)
          continue;

        TrackSpec[iCent][iSpc]->Fill (trkpt, event_weight);

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
            trkPtProj[iPhi][iPtTrk] += -trkpt * cos (dphi);
          else
            trkPtProj[iPhi][iPtTrk] += trkpt * cos (dphi);
          iPhi++;
        }
        
        // Study track yield relative to Z-going direction (requires dphi in 0 to pi)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), false);
        for (short idPhi = 0; idPhi < numPhiBins; idPhi++)
          if (phiLowBins[idPhi] <= dphi && dphi <= phiHighBins[idPhi])
            ZTracksPt[iSpc][iPtZ][iXZTrk][idPhi][iCent]->Fill (trkpt, event_weight);

        //// Study correlations (requires dphi in -pi/2 to 3pi/2)
        //dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        //if (dphi < -pi/2)
        //  dphi = dphi + 2*pi;

        ZTrackPtPhi[iPtZ][iXZTrk][iCent][iSpc]->Fill (dphi, trkpt, event_weight);
        
      } // end loop over tracks

      for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          ZMissingPt[iSpc][iPtZ][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
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
  if (inFile) { delete inFile; inFile = nullptr; }

  Delete1DArray (trkPtProj, numPhiBins);
}

#endif
