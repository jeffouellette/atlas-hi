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
  MinbiasAnalysis () : Analysis () {
    name = "minbias";
    directory = "MinbiasAnalysis/";
    plotFill = true;
    plotSignal = false;
    useAltMarker = false;
    SetupDirectories (directory, "ZTrackAnalysis/");
  }

  void Execute ();

  void CombineHists () override;
};


////////////////////////////////////////////////////////////////////////////////////////////////
// Fill combined species histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis::CombineHists () {
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    for (short iSpc = 0; iSpc < 3; iSpc++) {
      ZPhiYields[iCent][iSpc]->Add (ZPhiYields[iCent][0]);
      ZPtSpecs[iCent][iSpc]->Add (ZPtSpecs[iCent][0]);
      ZMYields[iCent][iSpc]->Add (ZMYields[iCent][0]);
      TrackSpec[iCent][iSpc]->Add (TrackSpec[iCent][0]);
      //DRDists[iCent][iSpc]->Add (DRDists[iCent][0]);

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++) {
          if (iSpc == 0 && iPtZ == 0 && iXZTrk == 0)
            continue;
          ZTrackPtPhi[iPtZ][iXZTrk][iCent][iSpc]->Add (ZTrackPtPhi[0][0][iCent][0]);
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            ZTracksPt[iSpc][iPtZ][iXZTrk][iPhi][iCent]->Add (ZTracksPt[0][0][0][iPhi][iCent]);
          } // end loop over phi
        }
        for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          ZMissingPt[iSpc][iPtZ][iPhi][iCent]->Add (ZMissingPt[0][0][iPhi][iCent]);
        } // end loop over phi
        if (iSpc == 0 && iPtZ == 0)
          continue;
        ZCounts[iSpc][iPtZ][iCent]->Add (ZCounts[0][0][iCent]);
      } // end loop over pT^Z
    } // end loop over species
  } // end loop over centralities
  return;
}



////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately.
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis::Execute () {
  SetupDirectories (directory, "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
  PbPbEventReweights = (TH3D*)eventWeightsFile->Get ("PbPbEventReweights_minbias");
  ppEventReweights = (TH1D*)eventWeightsFile->Get ("ppEventReweights_minbias");

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

      event_weight = PbPbEventReweights->GetBinContent (PbPbEventReweights->FindBin (fcal_et, q2, vz));

      FCalSpec->Fill (fcal_et, event_weight);
      FCalQ2Corr->Fill (fcal_et, q2, event_weight);

      float dphi = DeltaPhi (z_phi, psi2, false);
      if (dphi > pi/2)
        dphi = pi - dphi;

      ZPhiYields[iCent][iSpc]->Fill (2*dphi, event_weight);

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

        const short iXZTrk = 0;

        TrackSpec[iCent][iSpc]->Fill (trkpt, event_weight);

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

        for (short iXZTrk = 0; iXZTrk < nXZTrkBins; iXZTrk++)
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

      event_weight = ppEventReweights->GetBinContent (ppEventReweights->FindBin (vz));

      ZCounts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);

      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        const short iXZTrk = 0;

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
  SaferDelete (inFile);

  Delete1DArray (trkPtProj, numPhiBins);
}

#endif
