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
    plotFill = false;
    useAltMarker = true;
    SetupDirectories (directory, "ZTrackAnalysis/");
  }
  void Execute ();

};



////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately.
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis::Execute () {
  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");
  PbPbEventReweights = (TH3D*)eventWeightsFile->Get ("PbPbEventReweights_minbias");

  CreateHists ();

  const bool isMC = false;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_eta = 0, z_phi = 0, event_weight = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;
  double** trkPtProj = Get2DArray <double> (numPhiBins, nPtTrkBins);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  if (PbPbTree) {
    PbPbTree->SetBranchAddress ("event_weight", &event_weight);
    PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
    PbPbTree->SetBranchAddress ("q2",           &q2);
    PbPbTree->SetBranchAddress ("psi2",         &psi2);
    PbPbTree->SetBranchAddress ("vz",           &vz);
    PbPbTree->SetBranchAddress ("z_eta",        &z_eta);
    PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
    PbPbTree->SetBranchAddress ("trk_pt",       &trk_pt);
    PbPbTree->SetBranchAddress ("trk_eta",      &trk_eta);
    PbPbTree->SetBranchAddress ("trk_phi",      &trk_phi);

    const int nEvts = PbPbTree->GetEntries ();
    //int nEvts = 10000;
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      PbPbTree->GetEntry (iEvt);

      //if (fabs (vz) > 1.5)
      //  continue; // vertex cut

      short iCent = 0;
      while (iCent < numCentBins) {
        if (fcal_et < centBins[iCent])
          break;
        else
          iCent++;
      }
      if (iCent < 1 || iCent > numCentBins-1)
        continue;

      FCalSpec->Fill (fcal_et, event_weight);
      FCalQ2Corr->Fill (fcal_et, q2, event_weight);
      event_weight *= PbPbEventReweights->GetBinContent (PbPbEventReweights->FindBin (fcal_et, psi2, vz));

      float dphi = DeltaPhi (z_phi, psi2, false);
      if (dphi > pi/2)
        dphi = pi - dphi;
      ZPhiYields[iCent]->Fill (2*dphi, event_weight);

      for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
          trkPtProj[iPhi][iPtTrk] = 0;
        }
      }

      for (short iSpc = 0; iSpc < 3; iSpc++)
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) 
          ZCounts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);

      for (int iTrk = 0; iTrk < trk_pt->size (); iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < 1)
          continue;

        for (short iSpc = 0; iSpc < 3; iSpc++)
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
        short idPhi = 0;
        while (idPhi < numPhiBins-1) {
          if (phiLowBins[idPhi] < dphi && dphi < phiHighBins[idPhi])
            break;
          else
            idPhi++;
        }

        if (idPhi >= 0 && idPhi < numPhiBins)
          for (short iSpc = 0; iSpc < 3; iSpc++)
            for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) 
              ZTracksPt[iSpc][iPtZ][idPhi][iCent]->Fill (trkpt, event_weight);

        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        for (short iSpc = 0; iSpc < 3; iSpc++)
          ZTrackPtPhi[iCent][iSpc]->Fill (dphi, trkpt, event_weight);
      } // end loop over tracks

      for (short iSpc = 0; iSpc < 3; iSpc++) {
        for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
            for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++)
              ZMissingPt[iSpc][iPtZ][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
            trkPtProj[iPhi][iPtTrk] = 0;
          }
        }
      }
    } // end loop over Pb+Pb tree
    cout << endl;
  }

  for (short iSpc = 0; iSpc < 3; iSpc++) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          TH1D* thisHist = ZTracksPt[iSpc][iPtZ][iPhi][iCent];
          TH1D* countsHist = ZCounts[iSpc][iPtZ][iCent];
          const double yieldNormFactor = countsHist->GetBinContent (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);
          //const double yieldNormFactorError = countsHist->GetBinError (1) * (phiHighBins[iPhi]-phiLowBins[iPhi]);

          //RescaleWithError (thisHist, yieldNormFactor, yieldNormFactorError);
          if (yieldNormFactor > 0)
            thisHist->Scale (1. / yieldNormFactor);
        } // end loop over phi
      } // end loop over pT^Z
    } // end loop over centralities
  } // end loop over species

  SaveHists ();
  //LoadHists ();

  inFile->Close ();
  SaferDelete (inFile);

  Delete1DArray (trkPtProj, numPhiBins);
}

#endif
