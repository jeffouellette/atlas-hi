#ifndef __ReweightMuons_C__
#define __ReweightMuons_C__

#include "Params.h"

#include <ArrayTemplates.h>

using namespace atlashi;

void ReweightMuons () {

  TH1D*** h_z_y_ee = Get2DArray <TH1D*> (numCentBins, nPtZBins);
  TH1D*** h_z_y_mumu = Get2DArray <TH1D*> (numCentBins, nPtZBins);

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      h_z_y_ee[iCent][iPtZ] = new TH1D (Form ("h_z_y_ee_iCent%i_iPtZ%i", iCent, iPtZ), "", 20, -2.5, 2.5);
      h_z_y_mumu[iCent][iPtZ] = new TH1D (Form ("h_z_y_mumu_iCent%i_iPtZ%i", iCent, iPtZ), "", 20, -2.5, 2.5);
      h_z_y_ee[iCent][iPtZ]->Sumw2 ();
      h_z_y_mumu[iCent][iPtZ]->Sumw2 ();
    }
  }

  bool isEE = 0;
  float fcal_et = 0, z_y = 0, z_pt = 0;

  TFile* inFile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/data18hi.root", "read");
  TTree* inTree = (TTree*) inFile->Get ("PbPbZTrackTree");

  inTree->SetBranchAddress ("isEE", &isEE);
  inTree->SetBranchAddress ("fcal_et", &fcal_et);
  inTree->SetBranchAddress ("z_y", &z_y);
  inTree->SetBranchAddress ("z_pt", &z_pt);

  for (int i = 0; i < inTree->GetEntries (); i++) {
    inTree->GetEntry (i);

    const short iCent = GetCentBin (fcal_et);
    const short iPtZ = GetPtZBin (z_pt);

    if (iCent < 0 || iCent > numCentBins-1 || iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;

    if (isEE)
      h_z_y_ee[iCent][iPtZ]->Fill (z_y);
    else
      h_z_y_mumu[iCent][iPtZ]->Fill (z_y);
  }


  inTree = (TTree*) inFile->Get ("ppZTrackTree");

  inTree->SetBranchAddress ("isEE", &isEE);
  inTree->SetBranchAddress ("z_y", &z_y);
  inTree->SetBranchAddress ("z_pt", &z_pt);

  for (int i = 0; i < inTree->GetEntries (); i++) {
    inTree->GetEntry (i);

    const short iCent = 0;
    const short iPtZ = GetPtZBin (z_pt);

    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;

    if (isEE)
      h_z_y_ee[iCent][iPtZ]->Fill (z_y);
    else
      h_z_y_mumu[iCent][iPtZ]->Fill (z_y);
  }

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      h_z_y_mumu[iCent][iPtZ]->Divide (h_z_y_ee[iCent][iPtZ]);
      h_z_y_ee[iCent][iPtZ]->Divide (h_z_y_ee[iCent][iPtZ]);
    }
  }

  TFile* outFile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/zWeights.root", "recreate");

  for (int iCent = 0; iCent < numCentBins; iCent++) {
    for (int iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
      h_z_y_mumu[iCent][iPtZ]->Write ("", TObject :: kOverwrite);
      h_z_y_ee[iCent][iPtZ]->Write ("", TObject :: kOverwrite);
    }
  }
  outFile->Close ();


}

#endif
