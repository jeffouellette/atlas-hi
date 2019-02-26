#ifndef __ZTrackAnalysis_h__
#define __ZTrackAnalysis_h__

#include "Params.h"
#include "Analysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class ZTrackAnalysis : public Analysis {

  public:
  ZTrackAnalysis () : Analysis () {
    SetupDirectories ("ZTrackAnalysis/", "ZTrackAnalysis/");
  }
  void Execute ();

};


////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void ZTrackAnalysis::Execute () {
  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  TFile* fcalWeightsFile = new TFile ("FCalWeightsFile.root", "read");
  TH1D* fcalWeights = (TH1D*)fcalWeightsFile->Get ("FCalWeights");

  for (short iData = 0; iData < 2; iData++) {
    const char* data = iData == 0 ? "data" : "mc";
    PbPbEventInfoDist[iData] = new TH3D (Form ("PbPbEventInfoDist_%s_ztrack", data), "", nFCalBins, fcalBins, nPsi2Bins, psi2Bins, nVertZBins, vertZBins);
    PbPbEventInfoDist[iData]->Sumw2 ();
    ppEventInfoDist[iData] = new TH1D (Form ("ppEventInfoDist_%s_ztrack", data), "", nVertZBins, vertZBins);
    ppEventInfoDist[iData]->Sumw2 ();
    FCalSpec[iData] = new TH1D (Form ("FCalSpec_%s_ztrack", data), "", 300, 0, 6000); 
    FCalSpec[iData]->Sumw2 ();
    FCalQ2Corr[iData] = new TH2D (Form ("FCalQ2Corr_%s_ztrack", data), "", 300, 0, 6000, 150, 0, 300);
    FCalQ2Corr[iData]->Sumw2 ();
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      ZPhiYields[iData][iCent] = new TH1D (Form ("ZPhiYield_%s_iCent%i_ztrack", data, iCent), "", 80, 0, pi);
      ZPhiYields[iData][iCent]->Sumw2 ();
      for (short iSpc = 0; iSpc < 2; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : "mumu");
        ZTrackPtPhi[iData][iCent][iSpc] = new TH2D (Form ("ZTrackPtPhi_%s_%s_iCent%i_ztrack", data, spc, iCent), "", 80, -pi/2, 3*pi/2, nPtTrkBins, ptTrkBins);
        ZPtSpecs[iData][iCent][iSpc] = new TH1D (Form ("ZPtSpec_%s_%s_iCent%i_ztrack", data, spc, iCent), "", 300, 0, 300);
        ZMYields[iData][iCent][iSpc] = new TH1D (Form ("ZMSpec_%s_%s_iCent%i_ztrack", data, spc, iCent), "", 40, 76, 106);
        TrackSpec[iData][iCent][iSpc] = new TH1D (Form ("TrackSpec_%s_%s_iCent%i_ztrack", data, spc, iCent), "", 100, 0, 100);
        //DRDists[iData][iCent][iSpc] = new TH2D (Form ("DRDist_%s_%s_iCent%i_ztrack", data, spc, iCent), "", 100, 0, 0.1, 80, 0, 0.8);
        
        ZTrackPtPhi[iData][iCent][iSpc]->Sumw2 ();
        ZPtSpecs[iData][iCent][iSpc]->Sumw2 ();
        ZMYields[iData][iCent][iSpc]->Sumw2 ();
        TrackSpec[iData][iCent][iSpc]->Sumw2 ();
        //DRDists[iData][iCent][iSpc]->Sumw2 ();
      }
      ElectronSpec[iData][iCent] = new TH1D (Form ("ElectronSpec_%s_iCent%i_ztrack", data, iCent), "", 250, 0, 250);
      MuonSpec[iData][iCent] = new TH1D (Form ("MuonSpec_%s_iCent%i_ztrack", data, iCent), "", 250, 0, 250);

      for (short iSpc = 0; iSpc < 2; iSpc++) {
        const char* spc = (iSpc == 0 ? "ee" : (iSpc == 1 ? "mumu" : "comb"));
        for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
          for (int iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
            ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent] = new TH2D (Form ("ZMissingPt_%s_%s_iPtZ%i_iPhi%i_iCent%i_ztrack", data, spc, iPtZ, iPhi, iCent), "", numZMissingPtBins, zMissingPtBins, nPtTrkBins, ptTrkBins);
            ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent]->Sumw2 ();
          }
          for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
            ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent] = new TH1D (Form ("ZTracksPt_%s_%s_iPtZ%i_iPhi%i_iCent%i_ztrack", data, spc, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
            ZTracksPt[iSpc][iPtZ][iData][iPhi][iCent]->Sumw2 ();
          }
          ZCounts[iSpc][iPtZ][iData][iCent] = new TH1D (Form ("ZCounts_%s_%s_iPtZ%i_iCent%i_ztrack", data, spc, iPtZ, iCent), "", 1, 0, 1);
          ZCounts[iSpc][iPtZ][iData][iCent]->Sumw2 ();
        }
      }
    }
  }

  bool isEE = false, isMC = false;
  float event_weight = 0, fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_pt = 0, z_eta = 0, z_phi = 0, z_m = 0, l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  int l1_charge = 0, l2_charge = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr;//, *l_trk_pt = nullptr, *l_trk_eta = nullptr, *l_trk_phi = nullptr;
  double** trkPtProj = Get2DArray <double> (numPhiBins, nPtTrkBins);


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over PbPb tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  PbPbTree->SetBranchAddress ("isEE",         &isEE);
  PbPbTree->SetBranchAddress ("isMC",         &isMC);
  PbPbTree->SetBranchAddress ("event_weight", &event_weight);
  PbPbTree->SetBranchAddress ("fcal_et",      &fcal_et);
  PbPbTree->SetBranchAddress ("q2",           &q2);
  PbPbTree->SetBranchAddress ("psi2",         &psi2);
  //PbPbTree->SetBranchAddress ("vz",         &vz);
  PbPbTree->SetBranchAddress ("z_pt",         &z_pt);
  PbPbTree->SetBranchAddress ("z_eta",        &z_eta);
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
  PbPbTree->SetBranchAddress ("trk_pt",       &trk_pt);
  PbPbTree->SetBranchAddress ("trk_eta",      &trk_eta);
  PbPbTree->SetBranchAddress ("trk_phi",      &trk_phi);
  //PbPbTree->SetBranchAddress ("l_trk_pt",     &l_trk_pt);
  //PbPbTree->SetBranchAddress ("l_trk_eta",    &l_trk_eta);
  //PbPbTree->SetBranchAddress ("l_trk_phi",    &l_trk_phi);

  int nEvts = PbPbTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    PbPbTree->GetEntry (iEvt);

    const short iData = isMC ? 1 : 0; // 0 for not MC (data), 1 for MC
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

    PbPbEventInfoDist[iData]->Fill (fcal_et, psi2, vz, event_weight);
    FCalSpec[iData]->Fill (fcal_et, event_weight);
    FCalQ2Corr[iData]->Fill (fcal_et, q2, event_weight);
    if (isMC)
      event_weight *= fcalWeights->GetBinContent (fcalWeights->FindBin (fcal_et));

    ZPtSpecs[iData][iCent][iSpc]->Fill (z_pt, event_weight);
    ZMYields[iData][iCent][iSpc]->Fill (z_m, event_weight);
    if (isEE) {
      ElectronSpec[iData][iCent]->Fill (l1_pt, event_weight);
      ElectronSpec[iData][iCent]->Fill (l2_pt, event_weight);
    }
    else {
      MuonSpec[iData][iCent]->Fill (l1_pt, event_weight);
      MuonSpec[iData][iCent]->Fill (l2_pt, event_weight);
    }

    float dphi = DeltaPhi (z_phi, psi2, false);
    if (dphi > pi/2)
      dphi = pi - dphi;
    ZPhiYields[iData][iCent]->Fill (2*dphi, event_weight);

    for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
      for (short iPhi = 0; iPhi < numPhiBins; iPhi++) {
        trkPtProj[iPhi][iPtTrk] = 0;
      }
    }

    ZCounts[iSpc][iPtZ][iData][iCent]->Fill (0.5, event_weight);
    for (int iTrk = 0; iTrk < trk_pt->size (); iTrk++) {
      const float trkpt = trk_pt->at (iTrk);

      if (trkpt < 1)
        continue;

      TrackSpec[iData][iCent][iSpc]->Fill (trkpt, event_weight);
      //float minDR = 2;
      //int minLTrk = -1;
      //for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
      //  float dR = DeltaR (l_trk_eta->at (iLTrk), trk_eta->at (iTrk), l_trk_phi->at (iLTrk), trk_phi->at (iTrk));
      //  if (dR < minDR) {
      //    minDR = dR;
      //    minLTrk = iLTrk;
      //  }
      //}
      //if (isEE && minLTrk != -1) 
      //  DRDists[iData][iCent][0]->Fill (minDR, fabs (l_trk_pt->at (minLTrk) - trk_pt->at (iTrk)) / l_trk_pt->at (minLTrk), event_weight);
      //else if (minLTrk != -1) 
      //  DRDists[iData][iCent][1]->Fill (minDR, fabs (l_trk_pt->at (minLTrk) - trk_pt->at (iTrk)) / l_trk_pt->at (minLTrk), event_weight);

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
        ZTracksPt[iSpc][iPtZ][iData][idPhi][iCent]->Fill (trkpt, event_weight);

      if (z_pt > 20) {
        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        ZTrackPtPhi[iData][iCent][iSpc]->Fill (dphi, trkpt, event_weight);
      }
    } // end loop over tracks

    for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
      for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
        trkPtProj[iPhi][iPtTrk] = 0;
      }
    }
  } // end loop over Pb+Pb tree
  cout << endl;


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over pp tree
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ppTree->SetBranchAddress ("isEE",         &isEE);
  ppTree->SetBranchAddress ("isMC",         &isMC);
  ppTree->SetBranchAddress ("event_weight", &event_weight);
  ppTree->SetBranchAddress ("psi2",         &psi2);
  //ppTree->SetBranchAddress ("vz",         &vz);
  ppTree->SetBranchAddress ("z_pt",         &z_pt);
  ppTree->SetBranchAddress ("z_eta",        &z_eta);
  ppTree->SetBranchAddress ("z_phi",        &z_phi);
  ppTree->SetBranchAddress ("z_m",          &z_m);
  ppTree->SetBranchAddress ("l1_pt",        &l1_pt);
  ppTree->SetBranchAddress ("l1_eta",       &l1_eta);
  ppTree->SetBranchAddress ("l1_phi",       &l1_phi);
  ppTree->SetBranchAddress ("l1_charge",    &l1_charge);
  ppTree->SetBranchAddress ("l2_pt",        &l2_pt);
  ppTree->SetBranchAddress ("l2_eta",       &l2_eta);
  ppTree->SetBranchAddress ("l2_phi",       &l2_phi);
  ppTree->SetBranchAddress ("l2_charge",    &l2_charge);
  ppTree->SetBranchAddress ("trk_pt",       &trk_pt);
  ppTree->SetBranchAddress ("trk_eta",      &trk_eta);
  ppTree->SetBranchAddress ("trk_phi",      &trk_phi);
  //ppTree->SetBranchAddress ("l_trk_pt",     &l_trk_pt);
  //ppTree->SetBranchAddress ("l_trk_eta",    &l_trk_eta);
  //ppTree->SetBranchAddress ("l_trk_phi",    &l_trk_phi);

  nEvts = ppTree->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    ppTree->GetEntry (iEvt);

    const short iData = isMC ? 1 : 0; // 0 for not MC (data), 1 for MC
    const short iSpc = isEE ? 0 : 1; // 0 for electrons, 1 for muons, 2 for combined
    const short iCent = 0; // iCent = 0 for pp

    short iPtZ = 0; // find z-pt bin
    while (iPtZ < nPtZBins) {
      if (z_pt < zPtBins[iPtZ+1])
        break;
      else
        iPtZ++;
    }

    ppEventInfoDist[iData]->Fill (vz, event_weight);

    ZPtSpecs[iData][iCent][iSpc]->Fill (z_pt, event_weight);
    ZMYields[iData][iCent][iSpc]->Fill (z_m, event_weight);
    if (isEE) {
      ElectronSpec[iData][iCent]->Fill (l1_pt, event_weight);
      ElectronSpec[iData][iCent]->Fill (l2_pt, event_weight);
    }
    else {
      MuonSpec[iData][iCent]->Fill (l1_pt, event_weight);
      MuonSpec[iData][iCent]->Fill (l2_pt, event_weight);
    }

    float dphi = DeltaPhi (z_phi, psi2, false);
    if (dphi > pi/2)
      dphi = pi - dphi;
    ZPhiYields[iData][iCent]->Fill (2*dphi, event_weight);

    ZCounts[iSpc][iPtZ][iData][iCent]->Fill (0.5, event_weight);
    for (int iTrk = 0; iTrk < trk_pt->size (); iTrk++) {
      const float trkpt = trk_pt->at (iTrk);

      if (trkpt < 1)
        continue;

      TrackSpec[iData][iCent][iSpc]->Fill (trkpt, event_weight);
      //float minDR = 2;
      //int minLTrk = -1;
      //if (!isMC) {
      //  for (int iLTrk = 0; iLTrk < l_trk_pt->size (); iLTrk++) {
      //    float dR = DeltaR (l_trk_eta->at (iLTrk), trk_eta->at (iTrk), l_trk_phi->at (iLTrk), trk_phi->at (iTrk));
      //    if (dR < minDR) {
      //      minDR = dR;
      //      minLTrk = iLTrk;
      //    }
      //  }
      //}
      //if (isEE && minLTrk != -1)
      //  DRDists[iData][iCent][0]->Fill (minDR, fabs (l_trk_pt->at (minLTrk) - trk_pt->at (iTrk)) / l_trk_pt->at (minLTrk), event_weight);
      //else if (minLTrk != -1)
      //  DRDists[iData][iCent][1]->Fill (minDR, fabs (l_trk_pt->at (minLTrk) - trk_pt->at (iTrk)) / l_trk_pt->at (minLTrk), event_weight);

      // Add to missing pT (requires dphi in -pi/2 to pi/2)
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
        ZTracksPt[iSpc][iPtZ][iData][idPhi][iCent]->Fill (trkpt, event_weight);

      if (z_pt > 20) {
        // Study correlations (requires dphi in -pi/2 to 3pi/2)
        dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        ZTrackPtPhi[iData][iCent][iSpc]->Fill (dphi, trkpt, event_weight);
      }
    } // end loop over tracks

    for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
      for (short iPtTrk = 0; iPtTrk < nPtTrkBins; iPtTrk++) {
        ZMissingPt[iSpc][iPtZ][iData][iPhi][iCent]->Fill (trkPtProj[iPhi][iPtTrk], 0.5*(ptTrkBins[iPtTrk]+ptTrkBins[iPtTrk+1]), event_weight);
        trkPtProj[iPhi][iPtTrk] = 0;
      }
    }
  } // end loop over pp tree


  for (short iData = 0; iData < 2; iData++) {
    const char* data = iData == 0 ? "data" : "mc";
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      ZTrackPtPhi[iData][iCent][2] = new TH2D (Form ("ZTrackPtPhi_%s_comb_iCent%i_ztrack", data, iCent), "", 80, -pi/2, 3*pi/2, nPtTrkBins, ptTrkBins);
      ZPtSpecs[iData][iCent][2] = new TH1D (Form ("ZPtSpec_%s_comb_iCent%i_ztrack", data, iCent), "", 300, 0, 300);
      ZMYields[iData][iCent][2] = new TH1D (Form ("ZMSpec_%s_comb_iCent%i_ztrack", data, iCent), "", 40, 76, 106);
      TrackSpec[iData][iCent][2] = new TH1D (Form ("TrackSpec_%s_comb_iCent%i_ztrack", data, iCent), "", 100, 0, 100);
      //DRDists[iData][iCent][2] = new TH2D (Form ("DRDist_%s_comb_iCent%i_ztrack", data, iCent), "", 100, 0, 0.1, 80, 0, 0.8);

      ZTrackPtPhi[iData][iCent][2]->Sumw2 ();
      ZPtSpecs[iData][iCent][2]->Sumw2 ();
      ZMYields[iData][iCent][2]->Sumw2 ();
      TrackSpec[iData][iCent][2]->Sumw2 ();
      //DRDists[iData][iCent][2]->Sumw2 ();

      ZTrackPtPhi[iData][iCent][2]->Add (ZTrackPtPhi[iData][iCent][0]);
      ZTrackPtPhi[iData][iCent][2]->Add (ZTrackPtPhi[iData][iCent][1]);
      ZPtSpecs[iData][iCent][2]->Add (ZPtSpecs[iData][iCent][0]);
      ZPtSpecs[iData][iCent][2]->Add (ZPtSpecs[iData][iCent][1]);
      ZMYields[iData][iCent][2]->Add (ZMYields[iData][iCent][0]);
      ZMYields[iData][iCent][2]->Add (ZMYields[iData][iCent][1]);
      TrackSpec[iData][iCent][2]->Add (TrackSpec[iData][iCent][0]);
      TrackSpec[iData][iCent][2]->Add (TrackSpec[iData][iCent][1]);
      //DRDists[iData][iCent][2]->Add (DRDists[iData][iCent][0]);
      //DRDists[iData][iCent][2]->Add (DRDists[iData][iCent][1]);

      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        for (short iPhi = 0; iPhi < numPhiTrkBins; iPhi++) {
          ZMissingPt[2][iPtZ][iData][iPhi][iCent] = new TH2D (Form ("ZMissingPt_%s_comb_iPtZ%i_iPhi%i_iCent%i_ztrack", data, iPtZ, iPhi, iCent), "", numZMissingPtBins, zMissingPtBins, nPtTrkBins, ptTrkBins);
          ZMissingPt[2][iPtZ][iData][iPhi][iCent]->Sumw2 ();
          ZMissingPt[2][iPtZ][iData][iPhi][iCent]->Add (ZMissingPt[0][iPtZ][iData][iPhi][iCent]);
          ZMissingPt[2][iPtZ][iData][iPhi][iCent]->Add (ZMissingPt[1][iPtZ][iData][iPhi][iCent]);
        }
        for (int iPhi = 0; iPhi < numPhiBins; iPhi++) {
          ZTracksPt[2][iPtZ][iData][iPhi][iCent] = new TH1D (Form ("ZTracksPt_%s_comb_iPtZ%i_iPhi%i_iCent%i_ztrack", data, iPtZ, iPhi, iCent), "", nPtTrkBins, ptTrkBins);
          ZTracksPt[2][iPtZ][iData][iPhi][iCent]->Sumw2 ();
          ZTracksPt[2][iPtZ][iData][iPhi][iCent]->Add (ZTracksPt[0][iPtZ][iData][iPhi][iCent]);
          ZTracksPt[2][iPtZ][iData][iPhi][iCent]->Add (ZTracksPt[1][iPtZ][iData][iPhi][iCent]);
        }
        ZCounts[2][iPtZ][iData][iCent] = new TH1D (Form ("ZCounts_%s_comb_iPtZ%i_iCent%i_ztrack", data, iPtZ, iCent), "", 1, 0, 1);
        ZCounts[2][iPtZ][iData][iCent]->Sumw2 ();
        ZCounts[2][iPtZ][iData][iCent]->Add (ZCounts[0][iPtZ][iData][iCent]);
        ZCounts[2][iPtZ][iData][iCent]->Add (ZCounts[1][iPtZ][iData][iCent]);
      }
    }
  }
  cout << endl;

  SaveHists ();
  //LoadHists ("ztrack");

  inFile->Close ();
  if (inFile) { delete inFile; inFile = nullptr; }

  fcalWeightsFile->Close ();
  if (fcalWeightsFile) { delete fcalWeightsFile; fcalWeightsFile = nullptr; }

  Delete1DArray (trkPtProj, numPhiBins);
}

#endif
