#ifndef __CalcPerEventYields_C__
#define __CalcPerEventYields_C__

#include <ArrayTemplates.h>
#include <GlobalParams.h>
#include <Utilities.h>

#include "../Params.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

void CalcPerEventYields () {

  SetupDirectories ("", "ZTrackAnalysis/");

  int part_n = 0;
  float part_pt[10000];
  float part_eta[10000];
  float part_phi[10000];
  float part_e[10000];
  int part_pdgid[10000];
  int part_barcode[10000];
  int part_charge[10000];

  float z_pt   = 0;
  float z_eta  = 0;
  float z_phi  = 0;
  float z_m    = 0;
  bool  isEE   = false;

  TH1D**** h_z_trk_pth = Get3DArray <TH1D*> (2, 3, nPtZBins); // pp/PbPb, species, ptZ
  TH1D**** h_z_trk_xhz = Get3DArray <TH1D*> (2, 3, nPtZBins); // pp/PbPb, species, ptZ
  int***   z_counts    = Get3DArray <int>   (2, 3, nPtZBins); // pp/PbPb, species, ptZ

  double* pthBins = logspace (0.5, 60, 25);
  double* xhzBins = logspace (0.01, 1, 25);

  for (int iCollSys : {0, 1}) {
    for (int iPtZ = nPtZBins-3; iPtZ < nPtZBins; iPtZ++) {
      for (int iSpc = 0; iSpc < 3; iSpc++) {
        h_z_trk_pth[iCollSys][iSpc][iPtZ] = new TH1D (Form ("h_z_trk_pth_%s_iSpc%i_iPtZ%i", iCollSys == 0 ? "vacuum":"medium", iSpc, iPtZ), ";#it{p}_{T}^{ ch} [GeV];Counts", 25, pthBins);
        h_z_trk_xhz[iCollSys][iSpc][iPtZ] = new TH1D (Form ("h_z_trk_xhz_%s_iSpc%i_iPtZ%i", iCollSys == 0 ? "vacuum":"medium", iSpc, iPtZ), ";#it{x}_{hZ};Counts", 25, xhzBins);
      }
    }
  }

  TFile* mediumFile = new TFile (Form ("%s/Jewel/mediumFile.root", rootPath.Data ()), "read");
  TTree* mediumTree = (TTree*) mediumFile->Get ("tree");

  mediumTree->SetBranchAddress ("part_n",        &part_n);
  mediumTree->SetBranchAddress ("part_pt",       part_pt);
  mediumTree->SetBranchAddress ("part_eta",      part_eta);
  mediumTree->SetBranchAddress ("part_phi",      part_phi);
  mediumTree->SetBranchAddress ("part_e",        part_e);
  mediumTree->SetBranchAddress ("part_pdgid",    part_pdgid);
  mediumTree->SetBranchAddress ("part_barcode",  part_barcode);
  mediumTree->SetBranchAddress ("part_charge",   part_charge);

  mediumTree->SetBranchAddress ("z_pt",   &z_pt);
  mediumTree->SetBranchAddress ("z_eta",  &z_eta);
  mediumTree->SetBranchAddress ("z_phi",  &z_phi);
  mediumTree->SetBranchAddress ("z_m",    &z_m);
  mediumTree->SetBranchAddress ("isEE",   &isEE);

  const int nMediumEvts = mediumTree->GetEntries ();

  for (int iEvt = 0; iEvt < nMediumEvts; iEvt++) {
    mediumTree->GetEntry (iEvt);

    if (z_pt < 15)
      continue;

    const short iPtZ = GetPtZBin (z_pt);
    const short iSpc = (isEE ? 0 : 1);

    for (int iPart = 0; iPart < part_n; iPart++) {
      if (DeltaPhi (z_phi, part_phi[iPart]) < 3*pi/4)
        continue;

      h_z_trk_pth[1][iSpc][iPtZ]->Fill (part_pt[iPart]);
      h_z_trk_xhz[1][iSpc][iPtZ]->Fill (part_pt[iPart] / z_pt);
    } // end ch. hadron loop

    z_counts[1][iSpc][iPtZ]++;
  } // end event loop

  mediumFile->Close ();

  TFile* vacuumFile = new TFile (Form ("%s/Jewel/vacuumFile.root", rootPath.Data ()), "read");
  TTree* vacuumTree = (TTree*) vacuumFile->Get ("tree");
  
  vacuumTree->SetBranchAddress ("part_n",        &part_n);
  vacuumTree->SetBranchAddress ("part_pt",       part_pt);
  vacuumTree->SetBranchAddress ("part_eta",      part_eta);
  vacuumTree->SetBranchAddress ("part_phi",      part_phi);
  vacuumTree->SetBranchAddress ("part_e",        part_e);
  vacuumTree->SetBranchAddress ("part_pdgid",    part_pdgid);
  vacuumTree->SetBranchAddress ("part_barcode",  part_barcode);
  vacuumTree->SetBranchAddress ("part_charge",   part_charge);

  vacuumTree->SetBranchAddress ("z_pt",   &z_pt);
  vacuumTree->SetBranchAddress ("z_eta",  &z_eta);
  vacuumTree->SetBranchAddress ("z_phi",  &z_phi);
  vacuumTree->SetBranchAddress ("z_m",    &z_m);
  vacuumTree->SetBranchAddress ("isEE",   &isEE);

  const int nVacuumEvts = vacuumTree->GetEntries ();

  for (int iEvt = 0; iEvt < nVacuumEvts; iEvt++) {
    vacuumTree->GetEntry (iEvt);

    if (z_pt < 15)
      continue;

    const short iPtZ = GetPtZBin (z_pt);
    const short iSpc = (isEE ? 0 : 1);

    for (int iPart = 0; iPart < part_n; iPart++) {
      if (DeltaPhi (z_phi, part_phi[iPart]) < 3*pi/4)
        continue;

      h_z_trk_pth[0][iSpc][iPtZ]->Fill (part_pt[iPart]);
      h_z_trk_xhz[0][iSpc][iPtZ]->Fill (part_pt[iPart] / z_pt);
    } // end ch. hadron loop

    z_counts[0][iSpc][iPtZ]++;
  } // end event loop


  // postprocessing
  for (int iCollSys : {0, 1}) {
    for (int iPtZ = nPtZBins-3; iPtZ < nPtZBins; iPtZ++) {
      h_z_trk_pth[iCollSys][2][iPtZ]->Add (h_z_trk_pth[iCollSys][0][iPtZ]);
      h_z_trk_pth[iCollSys][2][iPtZ]->Add (h_z_trk_pth[iCollSys][1][iPtZ]);
      h_z_trk_xhz[iCollSys][2][iPtZ]->Add (h_z_trk_xhz[iCollSys][0][iPtZ]);
      h_z_trk_xhz[iCollSys][2][iPtZ]->Add (h_z_trk_xhz[iCollSys][1][iPtZ]);
      z_counts[iCollSys][2][iPtZ] = z_counts[iCollSys][0][iPtZ] + z_counts[iCollSys][1][iPtZ];
    }
  }

  TFile* outFile = new TFile (Form ("%s/Jewel/hists.root", rootPath.Data ()), "recreate");
  for (int iCollSys : {0, 1}) {
    for (int iPtZ = nPtZBins-3; iPtZ < nPtZBins; iPtZ++) {
      for (int iSpc = 0; iSpc < 3; iSpc++) {

        if (z_counts[iCollSys][iSpc][iPtZ] > 0) {
          h_z_trk_pth[iCollSys][iSpc][iPtZ]->Scale (1./ (z_counts[iCollSys][iSpc][iPtZ] * (phiHighBins[numPhiBins-1] - phiLowBins[1])), "width");
          h_z_trk_xhz[iCollSys][iSpc][iPtZ]->Scale (1./ (z_counts[iCollSys][iSpc][iPtZ] * (phiHighBins[numPhiBins-1] - phiLowBins[1])), "width");
        }

        h_z_trk_pth[iCollSys][iSpc][iPtZ]->Write ();
        h_z_trk_xhz[iCollSys][iSpc][iPtZ]->Write ();
      }
    }
  }


  outFile->Close ();
  

} // end main

#endif
