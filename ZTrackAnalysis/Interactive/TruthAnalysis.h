#ifndef __TruthAnalysis_h__
#define __TruthAnalysis_h__

#include "Params.h"
#include "Analysis.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>

#include <iostream>

using namespace std;
using namespace atlashi;

class TruthAnalysis : public Analysis {

  public:
  TH1D**  ZJetCounts      = Get1DArray <TH1D*> (nPtZBins); // iPtZ 
  TH2D*   ZJetPt          = nullptr;
  TH2D*   ZJetXZJet       = nullptr;
  TH1D**  xZJetDists      = Get1DArray <TH1D*> (nPtZBins); // iPtZ
  TH2D**  ZJetdEtadPhi    = Get1DArray <TH2D*> (nPtZBins); // iPtZ
  TH1D**  ZJetdPhi        = Get1DArray <TH1D*> (nPtZBins); // iPtZ
  TH2D**  JetTrkdEtadPhi  = Get1DArray <TH2D*> (nPtZBins); // iPtZ
  TH1D**  JetTrkdPhi      = Get1DArray <TH1D*> (nPtZBins); // iPtZ
  

  TruthAnalysis () : Analysis () {
    name = "truth";
    directory = "TruthAnalysis/";
    plotFill = false;
    useAltMarker = false;
    SetupDirectories (directory, "ZTrackAnalysis/");
  }

  void Execute ();

  void CreateHists () override;
  void ScaleHists () override;
  void LoadHists () override;
  void SaveHists () override;  

  void PlotZJetPt ();
  void PlotxZJet ();
  void PlotZJetCorrelations ();
  void PlotJetTrkCorrelations ();
};


////////////////////////////////////////////////////////////////////////////////////////////////
// Create new histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis::CreateHists () {
  Analysis::CreateHists ();

  ZJetPt = new TH2D ("ZJetPt", "", 75, 0, 300, 75, 1, 300);
  ZJetPt->Sumw2 ();

  ZJetXZJet = new TH2D ("ZJetXZJet", "", 75, 0, 300, 200, 0, 4);

  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    ZJetCounts[iPtZ] = new TH1D (Form ("ZJetCounts_iPtZ%i", iPtZ), "", 1, 0, 1);
    ZJetCounts[iPtZ]->Sumw2 ();
    xZJetDists[iPtZ] = new TH1D (Form ("xZJetDists_iPtZ%i", iPtZ), "", 100, 0, 1.5);
    xZJetDists[iPtZ]->Sumw2 ();
    ZJetdEtadPhi[iPtZ] = new TH2D (Form ("ZJetdEtadPhi_iPtZ%i", iPtZ), "", 50, -5, 5, 80, -pi/2, 3*pi/2);
    ZJetdEtadPhi[iPtZ]->Sumw2 ();
    JetTrkdEtadPhi[iPtZ] = new TH2D (Form ("JetTrkdEtadPhi_iPtZ%i", iPtZ), "", 50, -5, 5, 80, -pi/2, 3*pi/2);
    JetTrkdEtadPhi[iPtZ]->Sumw2 ();
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Scale histograms for plotting, calculating signals, etc.
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis::ScaleHists () {
  if (histsScaled || !histsLoaded)
    return;

  Analysis::ScaleHists ();

  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    ZJetdPhi[iPtZ] = ZJetdEtadPhi[iPtZ]->ProjectionY (Form ("ZJetdPhi_iPtZ%i", iPtZ));
    JetTrkdPhi[iPtZ] = JetTrkdEtadPhi[iPtZ]->ProjectionY (Form ("JetTrkdPhi_iPtZ%i", iPtZ));

    if (ZJetCounts[iPtZ]->GetBinContent (1) > 0) {
      xZJetDists[iPtZ]->Scale (1. / ZJetCounts[iPtZ]->GetBinContent (1), "width");
      ZJetdPhi[iPtZ]->Scale (1. / ZJetCounts[iPtZ]->GetBinContent (1), "width");
      JetTrkdPhi[iPtZ]->Scale (1. / ZJetCounts[iPtZ]->GetBinContent (1), "width");
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Load pre-filled histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis::LoadHists () {
  if (histsLoaded)
    return;

  Analysis::LoadHists ();
  if (!histFile) {
    SetupDirectories (directory, "ZTrackAnalysis/");
    histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "read");
  }

  ZJetPt = (TH2D*)histFile->Get ("ZJetPt");
  ZJetXZJet = (TH2D*)histFile->Get ("ZJetXZJet");

  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    ZJetCounts[iPtZ] = (TH1D*)histFile->Get (Form ("ZJetCounts_iPtZ%i", iPtZ));
    xZJetDists[iPtZ] = (TH1D*)histFile->Get (Form ("xZJetDists_iPtZ%i", iPtZ));

    ZJetdEtadPhi[iPtZ] = (TH2D*)histFile->Get (Form ("ZJetdEtadPhi_iPtZ%i", iPtZ));
    ZJetdPhi[iPtZ] = (TH1D*)histFile->Get (Form ("ZJetdPhi_iPtZ%i", iPtZ));

    JetTrkdEtadPhi[iPtZ] = (TH2D*)histFile->Get (Form ("JetTrkdEtadPhi_iPtZ%i", iPtZ));
    JetTrkdPhi[iPtZ] = (TH1D*)histFile->Get (Form ("JetTrkdPhi_iPtZ%i", iPtZ));
  }

  histsLoaded = true;
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Save histograms
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis::SaveHists () {
  Analysis::SaveHists ();

  if (!histFile) {
    SetupDirectories (directory, "ZTrackAnalysis/");
    histFile = new TFile (Form ("%s/savedHists.root", rootPath.Data ()), "recreate");
    histFile->cd ();
  }

  SafeWrite (ZJetPt);
  SafeWrite (ZJetXZJet);

  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    SafeWrite (ZJetCounts[iPtZ]);
    SafeWrite (xZJetDists[iPtZ]);

    SafeWrite (ZJetdEtadPhi[iPtZ]);
    SafeWrite (ZJetdPhi[iPtZ]);

    SafeWrite (JetTrkdEtadPhi[iPtZ]);
    SafeWrite (JetTrkdPhi[iPtZ]);
  }

  histFile->Close ();
  histFile = nullptr;
  histsLoaded = false;
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Main macro. Loops over Pb+Pb and pp trees and fills histograms appropriately, then saves them.
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis::Execute () {
  SetupDirectories (directory, "ZTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TTree* PbPbTree = (TTree*)inFile->Get ("PbPbZTrackTree");
  TTree* ppTree = (TTree*)inFile->Get ("ppZTrackTree");

  TFile* eventWeightsFile = new TFile (Form ("%s/eventWeightsFile.root", rootPath.Data ()), "read");

  h_PbPbEventReweights = (TH3D*)eventWeightsFile->Get ("h_PbPbEventReweights_truth");
  h_ppEventReweights = (TH1D*)eventWeightsFile->Get ("h_ppEventReweights_truth");

  CreateHists ();

  bool isEE = false;
  float event_weight = 1, fcal_et = 0, q2 = 0, psi2 = 0, vz = 0, z_pt = 0, z_eta = 0, z_phi = 0, z_m = 0, l1_pt = 0, l1_eta = 0, l1_phi = 0, l2_pt = 0, l2_eta = 0, l2_phi = 0;
  int l1_charge = 0, l2_charge = 0, ntrk = 0, njet = 0;
  vector<float>* trk_pt = nullptr, *trk_eta = nullptr, *trk_phi = nullptr, *jet_pt = nullptr, *jet_eta = nullptr, *jet_phi = nullptr, *jet_e = nullptr;
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

      event_weight = h_PbPbEventReweights->GetBinContent (h_PbPbEventReweights->FindBin (fcal_et, psi2, vz));

      h_fcal_et->Fill (fcal_et, event_weight);
      h_fcal_et_q2->Fill (fcal_et, q2, event_weight);

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      if (z_pt > zPtBins[1]) {
        h_z_m[iCent][iSpc]->Fill (z_m, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);

        float dphi = DeltaPhi (z_phi, psi2, false);
        if (dphi > pi/2)
          dphi = pi - dphi;
        h_z_phi[iCent][iSpc]->Fill (2*dphi, event_weight);
      }

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

        const float xZTrk = trkpt / z_pt;
        const short iXZTrk = GetiXZTrk (xZTrk);
        if (iXZTrk < 0 || iXZTrk > nXZTrkBins-1)
          continue;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight);

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
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][idPhi][iCent]->Fill (trkpt, event_weight);

        //// Study correlations (requires dphi in -pi/2 to 3pi/2)
        //dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        //if (dphi < -pi/2)
        //  dphi = dphi + 2*pi;

        h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]->Fill (dphi, trkpt, event_weight);
        
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
    ppTree->SetBranchAddress ("isEE",          &isEE);
    ppTree->SetBranchAddress ("vz",            &vz);
    ppTree->SetBranchAddress ("z_pt",          &z_pt);
    ppTree->SetBranchAddress ("z_eta",         &z_eta);
    ppTree->SetBranchAddress ("z_phi",         &z_phi);
    ppTree->SetBranchAddress ("z_m",           &z_m);
    ppTree->SetBranchAddress ("l1_pt",         &l1_pt);
    ppTree->SetBranchAddress ("l1_eta",        &l1_eta);
    ppTree->SetBranchAddress ("l1_phi",        &l1_phi);
    ppTree->SetBranchAddress ("l1_charge",     &l1_charge);
    ppTree->SetBranchAddress ("l2_pt",         &l2_pt);
    ppTree->SetBranchAddress ("l2_eta",        &l2_eta);
    ppTree->SetBranchAddress ("l2_phi",        &l2_phi);
    ppTree->SetBranchAddress ("l2_charge",     &l2_charge);
    ppTree->SetBranchAddress ("ntrk",          &ntrk);
    ppTree->SetBranchAddress ("trk_pt",        &trk_pt);
    ppTree->SetBranchAddress ("trk_eta",       &trk_eta);
    ppTree->SetBranchAddress ("trk_phi",       &trk_phi);
    ppTree->SetBranchAddress ("truth_jet_n",   &njet);
    ppTree->SetBranchAddress ("truth_jet_pt",  &jet_pt);
    ppTree->SetBranchAddress ("truth_jet_eta", &jet_eta);
    ppTree->SetBranchAddress ("truth_jet_phi", &jet_phi);
    ppTree->SetBranchAddress ("truth_jet_e",   &jet_e);

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

      event_weight = h_ppEventReweights->GetBinContent (h_ppEventReweights->FindBin (vz));

      h_z_pt[iCent][iSpc]->Fill (z_pt, event_weight);
      if (z_pt > zPtBins[1]) {
        h_z_m[iCent][iSpc]->Fill (z_m, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l1_pt, event_weight);
        h_lepton_pt[iCent][iSpc]->Fill (l2_pt, event_weight);
      }

      ZJetCounts[iPtZ]->Fill (0.5, event_weight);
      for (int iJet = 0; iJet < njet; iJet++) {
        if (DeltaR (l1_eta, jet_eta->at (iJet), l1_phi, jet_phi->at (iJet)) < 0.2)
          continue;
        if (DeltaR (l2_eta, jet_eta->at (iJet), l2_phi, jet_phi->at (iJet)) < 0.2)
          continue;

        ZJetPt->Fill (z_pt, jet_pt->at (iJet), event_weight);
    
        float dphi = DeltaPhi (z_phi, jet_phi->at (iJet), true);
        if (dphi < -pi/2)
          dphi = dphi + 2*pi;

        if (dphi > 7*pi/8 && dphi < 9*pi/8)
          xZJetDists[iPtZ]->Fill (jet_pt->at (iJet) / z_pt, event_weight);

        ZJetdEtadPhi[iPtZ]->Fill (z_eta - jet_eta->at (iJet), dphi, event_weight);
      } // end loop over jets

      h_z_counts[iSpc][iPtZ][iCent]->Fill (0.5, event_weight);
      for (int iTrk = 0; iTrk < ntrk; iTrk++) {
        const float trkpt = trk_pt->at (iTrk);

        if (trkpt < trk_min_pt)
          continue;

        for (int iJet = 0; iJet < njet; iJet++) {
          if (DeltaR (l1_eta, jet_eta->at (iJet), l1_phi, jet_phi->at (iJet)) < 0.2)
            continue;
          if (DeltaR (l2_eta, jet_eta->at (iJet), l2_phi, jet_phi->at (iJet)) < 0.2)
            continue;
          if (fabs (jet_eta->at (iJet) - trk_eta->at (iTrk)) < 2)
            continue;

          float dphi = DeltaPhi (jet_phi->at (iJet), trk_phi->at (iTrk), true);
          if (dphi < -pi/2)
            dphi = dphi + 2*pi;
          JetTrkdEtadPhi[iPtZ]->Fill (jet_eta->at (iJet) - trk_eta->at (iTrk), dphi, event_weight);
        } // end loop over jets

        const float xZTrk = trkpt / z_pt;
        const short iXZTrk = GetiXZTrk (xZTrk);
        if (iXZTrk < 0 || iXZTrk > nXZTrkBins-1)
          continue;

        h_trk_pt[iCent][iSpc]->Fill (trkpt, event_weight);

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
            h_z_trk_pt[iSpc][iPtZ][iXZTrk][idPhi][iCent]->Fill (trkpt, event_weight);

        //// Study correlations (requires dphi in -pi/2 to 3pi/2)
        //dphi = DeltaPhi (z_phi, trk_phi->at (iTrk), true);
        //if (dphi < -pi/2)
        //  dphi = dphi + 2*pi;

        h_z_trk_pt_phi[iPtZ][iXZTrk][iCent][iSpc]->Fill (dphi, trkpt, event_weight);
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

  inFile->Close ();
  if (inFile) { delete inFile; inFile = nullptr; }

  Delete1DArray (trkPtProj, numPhiBins);
}


////////////////////////////////////////////////////////////////////////////////////////////////
// Plots Z-Jet pT correlation.
////////////////////////////////////////////////////////////////////////////////////////////////
void TruthAnalysis::PlotZJetPt () {
  SetupDirectories (directory, "ZTrackAnalysis/");

  TCanvas* c = new TCanvas ("ZJetPtCanvas", "", 800, 600);
  c->SetRightMargin (0.14);
  gPad->SetLogz ();

  ZJetPt->GetXaxis ()->SetTitle ("Truth #it{p}_{T}^{ Z} [GeV]");
  ZJetPt->GetYaxis ()->SetTitle ("Truth #it{p}_{T}^{ jet} [GeV]");
  ZJetPt->GetZaxis ()->SetTitle ("Counts");

  ZJetPt->Draw ("colz");
  c->SaveAs (Form ("%s/ZJetPtCorrelation.pdf", plotPath.Data ()));
}


void TruthAnalysis::PlotxZJet () {
  SetupDirectories (directory, "ZTrackAnalysis/");

  TCanvas* c = new TCanvas ("xZJetCanvas", "", 800, 600);

  double max = 0, min = 1e30;
  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h = xZJetDists[iPtZ];
    if (h->GetMaximum () > max)  max = h->GetMaximum ();
    if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
  }
  max = 1.3*max;
  min = 0.7*min;

  for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h = xZJetDists[iPtZ];

    h->GetYaxis ()->SetRangeUser (min, max);

    h->GetXaxis ()->SetTitle ("#it{x}_{Z}^{jet}");
    h->GetYaxis ()->SetTitle ("1/N_{Z} dN/dx");

    h->SetMarkerColor (colors[iPtZ-1]);
    h->SetLineColor (colors[iPtZ-1]);

    h->Draw (iPtZ == 0 ? "hist" : "same hist");
    myText (0.65, 0.9-0.06*(iPtZ-1), colors[iPtZ-1], Form ("%g < #it{p}_{T}^{ Z} < %g [GeV]", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);
  }
  myText (0.3, 0.9, kBlack, "#it{p}_{T}^{ jet} > 5 GeV", 0.04);
  myText (0.3, 0.84, kBlack, "|#Delta#phi| > 7#pi/8", 0.04);
  c->SaveAs (Form ("%s/xZJetDists.pdf", plotPath.Data ()));

}


void TruthAnalysis::PlotZJetCorrelations () {
  TCanvas* c = new TCanvas ("ZJetCorrelationsCanvas", "", 800, 600);

  double max = 0, min = 1e30;
  for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h = ZJetdPhi[iPtZ];
    if (h->GetMaximum () > max)  max = h->GetMaximum ();
    if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
  }
  max = 1.3*max;
  min = 0.7*min;

  for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h = ZJetdPhi[iPtZ];

    h->GetYaxis ()->SetRangeUser (min, max);

    h->GetXaxis ()->SetTitle ("Z-Jet #Delta#phi");
    h->GetYaxis ()->SetTitle ("1/N_{Z} dN/d#Delta#phi");

    h->SetMarkerColor (colors[iPtZ-1]);
    h->SetLineColor (colors[iPtZ-1]);

    h->Draw (iPtZ == 0 ? "hist" : "same hist");
    myText (0.65, 0.9-0.06*(iPtZ-1), colors[iPtZ-1], Form ("%g < #it{p}_{T}^{ Z} < %g [GeV]", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);
  }
  myText (0.3, 0.9, kBlack, "#it{p}_{T}^{ jet} > 5 GeV", 0.04);
  c->SaveAs (Form ("%s/ZJetCorrelations.pdf", plotPath.Data ()));
}


void TruthAnalysis::PlotJetTrkCorrelations () {
  SetupDirectories (directory, "ZTrackAnalysis/");

  TCanvas* c = new TCanvas ("JetTrkCorrelationsCanvas", "", 800, 600);

  double max = 0, min = 1e30;
  for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h = JetTrkdPhi[iPtZ];
    if (h->GetMaximum () > max)  max = h->GetMaximum ();
    if (h->GetMinimum (0) < min) min = h->GetMinimum (0);
  }
  max = (max > 0 ? 2*max : 1);
  min = (min > 0 ? 0.5*min : 0.1);

  gPad->SetLogy ();
  for (short iPtZ = 1; iPtZ < nPtZBins; iPtZ++) {
    TH1D* h = JetTrkdPhi[iPtZ];

    h->GetYaxis ()->SetRangeUser (min, max);

    h->GetXaxis ()->SetTitle ("Jet-Track #Delta#phi");
    h->GetYaxis ()->SetTitle ("1/N_{Z} dN_{assoc.}/d#Delta#phi");

    h->SetMarkerColor (colors[iPtZ-1]);
    h->SetLineColor (colors[iPtZ-1]);

    h->Draw (iPtZ == 0 ? "hist" : "same hist");
    myText (0.65, 0.9-0.06*(iPtZ-1), colors[iPtZ-1], Form ("%g < #it{p}_{T}^{ Z} < %g [GeV]", zPtBins[iPtZ], zPtBins[iPtZ+1]), 0.04);
  }
  myText (0.3, 0.9, kBlack, "#it{p}_{T}^{ jet} > 5 GeV", 0.04);
  c->SaveAs (Form ("%s/JetTrkCorrelations.pdf", plotPath.Data ()));
}

#endif
