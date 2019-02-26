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
    SetupDirectories ("MinbiasAnalysis/", "ZTrackAnalysis/");
  }
  void Execute ();

};



////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over minbias trees and fills histograms appropriately.
////////////////////////////////////////////////////////////////////////////////////////////////
void MinbiasAnalysis::Execute () {
  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");

  TFile* fcalWeightsFile = new TFile ("FCalWeightsFile.root", "read");
  TH1D* fcalWeights = (TH1D*)fcalWeightsFile->Get ("FCalWeights");

  for (short iData = 0; iData < 2; iData++) {
    const char* data = iData == 0 ? "data" : "mc";
    PbPbEventInfoDist[iData] = new TH3D (Form ("PbPbEventInfoDist_%s_minbias", data), "", nFCalBins, fcalBins, nPsi2Bins, psi2Bins, nVertZBins, vertZBins);
    PbPbEventInfoDist[iData]->Sumw2 ();
    ppEventInfoDist[iData] = new TH1D (Form ("ppEventInfoDist_%s_minbias", data), "", nVertZBins, vertZBins);
    ppEventInfoDist[iData]->Sumw2 ();
    FCalSpec[iData] = new TH1D (Form ("FCalSpec_%s", data), "", 300, 0, 6000); 
    FCalSpec[iData]->Sumw2 ();
    FCalQ2Corr[iData] = new TH2D (Form ("FCalQ2Corr_%s_minbias", data), "", 300, 0, 6000, 150, 0, 300);
    FCalQ2Corr[iData]->Sumw2 ();
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      ZPhiYields[iData][iCent] = new TH1D (Form ("ZPhiYield_%s_iCent%i_minbias", data, iCent), "", 80, 0, pi);
      ZPhiYields[iData][iCent]->Sumw2 ();
      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        ZTrackPtPhi[iData][iCent][iSpc] = new TH2D (Form ("ZTrackPtPhi_%s_%s_iCent%i_minbias", data, spc, iCent), "", 80, -pi/2, 3*pi/2, nPtTrkBins, ptTrkBins);
        ZPtSpecs[iData][iCent][iSpc] = new TH1D (Form ("ZPtSpec_%s_%s_iCent%i_minbias", data, spc, iCent), "", 300, 0, 300);
        TrackSpec[iData][iCent][iSpc] = new TH1D (Form ("TrackSpec_%s_%s_iCent%i_minbias", data, spc, iCent), "", 100, 0, 100);
        
        ZTrackPtPhi[iData][iCent][iSpc]->Sumw2 ();
        ZPtSpecs[iData][iCent][iSpc]->Sumw2 ();
        TrackSpec[iData][iCent][iSpc]->Sumw2 ();
      }

      for (short iSpc = 0; iSpc < 3; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
            ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent] = new TH2D (Form ("ZMissingPt_%s_%s_iPtZ%i_iPhi%i_iCent%i_minbias", data, spc, iPtZ, iPhi, iCent), "", numZMissingPtBins, zMissingPtBins, nPtTrkBins, ptTrkBins);
            ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent]->Sumw2 ();
          }
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent] = new TH1D (Form ("ZTracksPt_%s_%s_iPtZ%i_iPhi%i_iCent%i_minbias", data, spc, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
            ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent]->Sumw2 ();
          }
          ZCounts[iSpc][iPtZ][iData][iCent] = new TH1D (Form ("ZCounts_%s_%s_iPtZ%i_iCent%i_minbias", data, spc, iPtZ, iCent), "", 1, 0, 1);
          ZCounts[iSpc][iPtZ][iData][iCent]->Sumw2 ();
        }
      }
    }
  }

  const bool isMC = false;
  float fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_eta = 0, z_phi = 0, event_weight = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;
  double** trkPtProj = Get2DArray <double> (numPhiBins, nPtTrkBins);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  PbPbTree->SetBranchAddress ("event_weight", &event_weight);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);
  //PbPbTree->SetBranchAddress ("vz",           &vz);
  PbPbTree->SetBranchAddress ("z_eta",        &z_eta);
  PbPbTree->SetBranchAddress ("z_phi",        &z_phi);
  PbPbTree->SetBranchAddress ("trk_pt",       &trk_pt);
  PbPbTree->SetBranchAddress ("trk_eta",      &trk_eta);
  PbPbTree->SetBranchAddress ("trk_phi",      &trk_phi);

  int nEvts = PbPbTree->GetEntries ();
  //int nEvts = 10000;
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    PbPbTree->GetEntry (iEvt);

    const short iData = 0; // 0 for not MC (data), 1 for MC

    short iCent = 0;
    while (iCent < numCentBins) {
      if (fcal_et < centBins[iCent])
        break;
      else
        iCent++;
    }
    if (iCent < 1 || iCent > numCentBins-1)
      continue;

    PbPbEventInfoDist[iData]->Fill (fcal_et, psi2, vz, event_weight);
    FCalSpec[iData]->Fill (fcal_et, event_weight);
    FCalQ2Corr[iData]->Fill (fcal_et, q2, event_weight);

    float dphi = DeltaPhi (z_phi, psi2, false);
    if (dphi > pi/2)
      dphi = pi - dphi;
    ZPhiYields[iData][iCent]->Fill (2*dphi, event_weight);

    for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
      for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
        trkPtProj[iPhi][iPtTrk] = 0;
      }
    }

    for (short iSpc = 0; iSpc < 3; iSpc++)
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) 
        ZCounts[iSpc][iPtZ][iData][iCent]->Fill (0.5, event_weight);

    for (int iTrk = 0; iTrk < trk_pt->size (); iTrk++) {
      const float trkpt = trk_pt->at (iTrk);

      if (trkpt < 1)
        continue;

      for (short iSpc = 0; iSpc < 3; iSpc++)
        TrackSpec[iData][iCent][iSpc]->Fill (trkpt, event_weight);

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
            ZTracksPt[iSpc][iPtZ][iData][idPhi][iCent]->Fill (trkpt, event_weight);

      // Study correlations (requires dphi in -pi/2 to 3pi/2)
      dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
      if (dphi < -pi/2)
        dphi = dphi + 2*pi;

      for (short iSpc = 0; iSpc < 3; iSpc++)
        ZTrackPtPhi[iData][iCent][iSpc]->Fill (dphi, trkpt, event_weight);
    } // end loop over tracks

    for (short iSpc = 0; iSpc < 3; iSpc++) {
      for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
        for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
          for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++)
            ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
          trkPtProj[iPhi][iPtTrk] = 0;
        }
      }
    }
  } // end loop over Pb+Pb tree
  cout << endl;

  SaveHists ();
  //LoadHists ("minbias");

  inFile->Close ();
  SaferDelete (inFile);

  fcalWeightsFile->Close ();
  SaferDelete (fcalWeightsFile);

  Delete1DArray (trkPtProj, numPhiBins);
}

#endif
