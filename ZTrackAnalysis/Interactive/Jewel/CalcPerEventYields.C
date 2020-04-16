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
  TH1D**** h_z_trk_pth_var = Get3DArray <TH1D*> (2, 3, nPtZBins); // pp/PbPb, species, ptZ
  TH1D**** h_z_trk_xhz = Get3DArray <TH1D*> (2, 3, nPtZBins); // pp/PbPb, species, ptZ
  TH1D**** h_z_trk_xhz_var = Get3DArray <TH1D*> (2, 3, nPtZBins); // pp/PbPb, species, ptZ
  int***   z_counts    = Get3DArray <int>   (2, 3, nPtZBins); // pp/PbPb, species, ptZ

  const int nBinsPtch = 18;
  const int nBinsXhZ = 15;
  double** pthBins = Get1DArray <double*> (nPtZBins);
  double** xhzBins = Get1DArray <double*> (nPtZBins);

  TFile* outFile = new TFile (Form ("%s/Jewel/hists.root", rootPath.Data ()), "recreate");
  for (int iPtZ = nPtZBins-3; iPtZ < nPtZBins; iPtZ++) {
    const double maxPtch = pTchBins[iPtZ][nPtchBins[iPtZ]];
    const double minPtch = pTchBins[iPtZ][0];
    const double maxXhZ = xhZBins[iPtZ][nXhZBins[iPtZ]];
    const double minXhZ = xhZBins[iPtZ][0];
    pthBins[iPtZ] = logspace (minPtch, maxPtch, nBinsPtch);
    xhzBins[iPtZ] = logspace (minXhZ, maxXhZ, nBinsXhZ);
    for (int iCollSys : {0, 1}) {
      for (int iSpc = 0; iSpc < 3; iSpc++) {
        h_z_trk_pth[iCollSys][iSpc][iPtZ] = new TH1D (Form ("h_z_trk_pth_%s_iSpc%i_iPtZ%i", iCollSys == 0 ? "vacuum":"medium", iSpc, iPtZ), ";#it{p}_{T}^{ ch} [GeV];Counts", nBinsPtch, pthBins[iPtZ]);
        h_z_trk_pth_var[iCollSys][iSpc][iPtZ] = new TH1D (Form ("h_z_trk_pth_var_%s_iSpc%i_iPtZ%i", iCollSys == 0 ? "vacuum":"medium", iSpc, iPtZ), ";#it{p}_{T}^{ ch} [GeV];Counts", nBinsPtch, pthBins[iPtZ]);
        h_z_trk_xhz[iCollSys][iSpc][iPtZ] = new TH1D (Form ("h_z_trk_xhz_%s_iSpc%i_iPtZ%i", iCollSys == 0 ? "vacuum":"medium", iSpc, iPtZ), ";#it{x}_{hZ};Counts", nBinsXhZ, xhzBins[iPtZ]);
        h_z_trk_xhz_var[iCollSys][iSpc][iPtZ] = new TH1D (Form ("h_z_trk_xhz_var_%s_iSpc%i_iPtZ%i", iCollSys == 0 ? "vacuum":"medium", iSpc, iPtZ), ";#it{x}_{hZ};Counts", nBinsXhZ, xhzBins[iPtZ]);
      }
    }
  }

  TFile* mediumFile = new TFile (Form ("%s/Jewel/z_jet_medium_nevt_1760000.root", rootPath.Data ()), "read");
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

  int* yields_pth = Get1DArray <int> (nBinsPtch);
  int* yields_xhz = Get1DArray <int> (nBinsXhZ);

  for (int iEvt = 0; iEvt < nMediumEvts; iEvt++) {
    mediumTree->GetEntry (iEvt);

    for (int i = 0; i < nBinsPtch; i++)
      yields_pth[i] = 0;
    for (int i = 0; i < nBinsXhZ; i++)
      yields_xhz[i] = 0;

    if (z_pt < 15)
      continue;

    const short iPtZ = GetPtZBin (z_pt);
    const short iSpc = (isEE ? 0 : 1);

    for (int iPart = 0; iPart < part_n; iPart++) {
      if (DeltaPhi (z_phi, part_phi[iPart]) < 3*pi/4)
        continue;
      const double pth = part_pt[iPart];
      const double xhz = pth / z_pt;

      if (pth >= pthBins[iPtZ][0]) {
        int iPth = 0;
        while (iPth < nBinsPtch && pth >= pthBins[iPtZ][iPth+1]) iPth++;
        if (0 <= iPth && iPth < nBinsPtch)
          yields_pth[iPth]++;
      }
    
      if (xhz >= xhzBins[iPtZ][0]) {
        int iXhZ = 0;
        while (iXhZ < nBinsXhZ && xhz >= xhzBins[iPtZ][iXhZ+1]) iXhZ++; 
        if (0 <= iXhZ && iXhZ < nBinsXhZ)
          yields_xhz[iXhZ]++;
      }
    } // end ch. hadron loop

    
    for (int i = 0; i < nBinsPtch; i++) {
      h_z_trk_pth[1][iSpc][iPtZ]->SetBinContent (i+1, h_z_trk_pth[1][iSpc][iPtZ]->GetBinContent (i+1) + yields_pth[i]);
      h_z_trk_pth_var[1][iSpc][iPtZ]->SetBinContent (i+1, h_z_trk_pth_var[1][iSpc][iPtZ]->GetBinContent (i+1) + pow (yields_pth[i], 2));
    }
    for (int i = 0; i < nBinsXhZ; i++) {
      h_z_trk_xhz[1][iSpc][iPtZ]->SetBinContent (i+1, h_z_trk_xhz[1][iSpc][iPtZ]->GetBinContent (i+1) + yields_xhz[i]);
      h_z_trk_xhz_var[1][iSpc][iPtZ]->SetBinContent (i+1, h_z_trk_xhz_var[1][iSpc][iPtZ]->GetBinContent (i+1) + pow (yields_xhz[i], 2));
    }

    z_counts[1][iSpc][iPtZ]++;
  } // end event loop

  mediumFile->Close ();



  TFile* vacuumFile = new TFile (Form ("%s/Jewel/z_jet_vacuum_nevt_1760000.root", rootPath.Data ()), "read");
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

    for (int i = 0; i < nBinsPtch; i++)
      yields_pth[i] = 0;
    for (int i = 0; i < nBinsXhZ; i++)
      yields_xhz[i] = 0;

    if (z_pt < 15)
      continue;

    const short iPtZ = GetPtZBin (z_pt);
    const short iSpc = (isEE ? 0 : 1);

    for (int iPart = 0; iPart < part_n; iPart++) {
      if (DeltaPhi (z_phi, part_phi[iPart]) < 3*pi/4)
        continue;
      const double pth = part_pt[iPart];
      const double xhz = pth / z_pt;

      if (pth >= pthBins[iPtZ][0]) {
        int iPth = 0;
        while (iPth < nBinsPtch && pth >= pthBins[iPtZ][iPth+1]) iPth++;
        if (0 <= iPth && iPth < nBinsPtch)
          yields_pth[iPth]++;
      }
    
      if (xhz >= xhzBins[iPtZ][0]) {
        int iXhZ = 0;
        while (iXhZ < nBinsXhZ && xhz >= xhzBins[iPtZ][iXhZ+1]) iXhZ++; 
        if (0 <= iXhZ && iXhZ < nBinsXhZ)
          yields_xhz[iXhZ]++;
      }
    } // end ch. hadron loop

    
    for (int i = 0; i < nBinsPtch; i++) {
      h_z_trk_pth[0][iSpc][iPtZ]->SetBinContent (i+1, h_z_trk_pth[0][iSpc][iPtZ]->GetBinContent (i+1) + yields_pth[i]);
      h_z_trk_pth_var[0][iSpc][iPtZ]->SetBinContent (i+1, h_z_trk_pth_var[0][iSpc][iPtZ]->GetBinContent (i+1) + pow (yields_pth[i], 2));
    }
    for (int i = 0; i < nBinsXhZ; i++) {
      h_z_trk_xhz[0][iSpc][iPtZ]->SetBinContent (i+1, h_z_trk_xhz[0][iSpc][iPtZ]->GetBinContent (i+1) + yields_xhz[i]);
      h_z_trk_xhz_var[0][iSpc][iPtZ]->SetBinContent (i+1, h_z_trk_xhz_var[0][iSpc][iPtZ]->GetBinContent (i+1) + pow (yields_xhz[i], 2));
    }

    z_counts[0][iSpc][iPtZ]++;
  } // end event loop


  for (int iPtZ = nPtZBins-3; iPtZ < nPtZBins; iPtZ++) {
    delete[] pthBins[iPtZ];
    pthBins[iPtZ] = nullptr;
    delete[] xhzBins[iPtZ];
    xhzBins[iPtZ] = nullptr;
  }
  Delete1DArray (pthBins, nPtZBins);
  Delete1DArray (xhzBins, nPtZBins);


  // postprocessing
  for (int iCollSys : {0, 1}) {
    for (int iPtZ = nPtZBins-3; iPtZ < nPtZBins; iPtZ++) {
      h_z_trk_pth[iCollSys][2][iPtZ]->Add (h_z_trk_pth[iCollSys][0][iPtZ]);
      h_z_trk_pth[iCollSys][2][iPtZ]->Add (h_z_trk_pth[iCollSys][1][iPtZ]);
      h_z_trk_xhz[iCollSys][2][iPtZ]->Add (h_z_trk_xhz[iCollSys][0][iPtZ]);
      h_z_trk_xhz[iCollSys][2][iPtZ]->Add (h_z_trk_xhz[iCollSys][1][iPtZ]);
      h_z_trk_pth_var[iCollSys][2][iPtZ]->Add (h_z_trk_pth_var[iCollSys][0][iPtZ]);
      h_z_trk_pth_var[iCollSys][2][iPtZ]->Add (h_z_trk_pth_var[iCollSys][1][iPtZ]);
      h_z_trk_xhz_var[iCollSys][2][iPtZ]->Add (h_z_trk_xhz_var[iCollSys][0][iPtZ]);
      h_z_trk_xhz_var[iCollSys][2][iPtZ]->Add (h_z_trk_xhz_var[iCollSys][1][iPtZ]);
      z_counts[iCollSys][2][iPtZ] = z_counts[iCollSys][0][iPtZ] + z_counts[iCollSys][1][iPtZ];
    }
  }

  outFile->cd ();
  for (int iCollSys : {0, 1}) {
    for (int iPtZ = nPtZBins-3; iPtZ < nPtZBins; iPtZ++) {
      for (int iSpc = 0; iSpc < 3; iSpc++) {

        const int nZ = z_counts[iCollSys][iSpc][iPtZ];

        if (nZ > 0) {
          h_z_trk_pth[iCollSys][iSpc][iPtZ]->Scale (1./ nZ);// * (phiHighBins[numPhiBins-1] - phiLowBins[1])), "width");
          h_z_trk_xhz[iCollSys][iSpc][iPtZ]->Scale (1./ nZ);// * (phiHighBins[numPhiBins-1] - phiLowBins[1])), "width");
          for (int iX = 1; iX <= nBinsPtch; iX++)
            h_z_trk_pth_var[iCollSys][iSpc][iPtZ]->SetBinContent (iX, h_z_trk_pth_var[iCollSys][iSpc][iPtZ]->GetBinContent (iX) - nZ * pow (h_z_trk_pth[iCollSys][iSpc][iPtZ]->GetBinContent (iX), 2));
          for (int iX = 1; iX <= nBinsXhZ; iX++)
            h_z_trk_xhz_var[iCollSys][iSpc][iPtZ]->SetBinContent (iX, h_z_trk_xhz_var[iCollSys][iSpc][iPtZ]->GetBinContent (iX) - nZ * pow (h_z_trk_xhz[iCollSys][iSpc][iPtZ]->GetBinContent (iX), 2));
          h_z_trk_pth_var[iCollSys][iSpc][iPtZ]->Scale (1./ (nZ-1));
          h_z_trk_xhz_var[iCollSys][iSpc][iPtZ]->Scale (1./ (nZ-1));

          for (int iX = 1; iX <= nBinsPtch; iX++)
            h_z_trk_pth[iCollSys][iSpc][iPtZ]->SetBinError (iX, sqrt (h_z_trk_pth_var[iCollSys][iSpc][iPtZ]->GetBinContent (iX) / nZ));
          for (int iX = 1; iX <= nBinsXhZ; iX++)
            h_z_trk_xhz[iCollSys][iSpc][iPtZ]->SetBinError (iX, sqrt (h_z_trk_xhz_var[iCollSys][iSpc][iPtZ]->GetBinContent (iX) / nZ));

          h_z_trk_pth[iCollSys][iSpc][iPtZ]->Scale (1./ (phiHighBins[numPhiBins-1] - phiLowBins[1]), "width");
          h_z_trk_xhz[iCollSys][iSpc][iPtZ]->Scale (1./ (phiHighBins[numPhiBins-1] - phiLowBins[1]), "width");
        }

        h_z_trk_pth[iCollSys][iSpc][iPtZ]->Write ();
        h_z_trk_xhz[iCollSys][iSpc][iPtZ]->Write ();
      }
    }
  }

  outFile->Close ();
  

} // end main

#endif
