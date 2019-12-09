#ifndef __MakeMixingCatTrees_C__
#define __MakeMixingCatTrees_C__

#include "Params.h"

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>


const vector<int> runs = {
  365502,
  365512,
  365573,
  365602,
  365627,
  365678,
  365681,
  365709,
  365752,
  365834,
  365914,
  365932,
  366011,
  366029,
  366092,
  366142,
  366268,
  366337,
  366383,
  366413,
  366476,
  366526,
  366528,
  366627,
  366691,
  366754,
  366805,
  366860,
  366878,
  366919,
  366931,
  366994,
  367023,
  367099,
  367134,
  367165,
  367170,
  367233,
  367273,
  367318,
  367321,
  367363,
  367364,
  367365,
  367384
};


const float fileCentBins[13] = {
    66.402, // 80%
   296.17,  // 60%
   885.172, // 40%
  1378.92,  // 30%
  1690.47,  // 25%
  2055.77,  // 20%
  2484.75,  // 15%
  2995.94,  // 10%
  3229.67,  //  8%
  3485.57,  //  6%
  3767.,    //  4%
  4083.38,  //  2%
  5000      //  0%,   entry in array is numFileCentBins-1
};
const int numFileCentBins = sizeof (fileCentBins) / sizeof (fileCentBins[0]);

short GetFileCentBin (const float fcal_et) {
  short i = 0;
  while (i < numFileCentBins) {
    if (fcal_et < fileCentBins[i])
      break;
    i++;
  }
  return i;
}


const short numFileQ2Bins = 10;
const double* fileQ2Bins = linspace (0, 0.2, numFileQ2Bins);

short GetFileQ2Bin (const float q2) {
  short i = 0;
  while (i < numFileQ2Bins) {
    if (q2 < fileQ2Bins[i+1])
      break;
    i++;
  }
  return i;
}


void MakeMixingCatTrees (const char* path, TChain* c, const float fcal_et_shift = 0) {

  gSystem->Exec (Form ("mkdir -p %s/MixingTrees", path));

  TFile**** farr = Get3DArray <TFile*> (numFileCentBins, numFileQ2Bins, numPsi2Bins);
  TTree**** tarr = Get3DArray <TTree*> (numFileCentBins, numFileQ2Bins, numPsi2Bins);

  for (int iCent = 1; iCent < numFileCentBins; iCent++) {
    for (int iQ2 = 0; iQ2 < numFileQ2Bins; iQ2++) {
      for (int iPsi2 = 0; iPsi2 < numPsi2Bins; iPsi2++) {

        farr[iCent][iQ2][iPsi2] = new TFile (Form ("%s/MixingTrees/tree_iCent%i_iQ2%i_iPsi2%i.root", path, iCent, iQ2, iPsi2), "update");
        tarr[iCent][iQ2][iPsi2] = (TTree*) farr[iCent][iQ2][iPsi2]->Get ("PbPbZTrackTree");
        if (!tarr[iCent][iQ2][iPsi2])
          tarr[iCent][iQ2][iPsi2] = (TTree*) c->CloneTree (0);
        tarr[iCent][iQ2][iPsi2]->SetDirectory (farr[iCent][iQ2][iPsi2]);
      }
    }
  }

  float fcal_et = 0, q2 = 0, psi2 = 0;
  c->SetBranchAddress ("fcal_et", &fcal_et);
  c->SetBranchAddress ("q2",      &q2);
  c->SetBranchAddress ("psi2",    &psi2);

  const int nEvts = c->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    c->GetEntry (iEvt);

    //fcal_et = fcal_et + fcal_et_shift;

    const short iCent = GetFileCentBin (fcal_et);
    if (iCent < 1 || iCent > numFileCentBins-1)
      continue;

    const short iQ2 = GetFileQ2Bin (q2);
    if (iQ2 < 0 || iQ2 > numFileQ2Bins-1)
      continue;

    const short iPsi2 = GetPsi2Bin (psi2);
    if (iPsi2 < 0 || iPsi2 > numPsi2Bins-1)
      continue;

    tarr[iCent][iQ2][iPsi2]->Fill ();
  }
  cout << endl;

  
  for (int iCent = 1; iCent < numFileCentBins; iCent++) {
    for (int iQ2 = 0; iQ2 < numFileQ2Bins; iQ2++) {
      for (int iPsi2 = 0; iPsi2 < numPsi2Bins; iPsi2++) {
        farr[iCent][iQ2][iPsi2]->cd ();
        tarr[iCent][iQ2][iPsi2]->Write ("", TObject :: kOverwrite);
        farr[iCent][iQ2][iPsi2]->Close ();
      }
    }
  }

  //Delete1DArray <TFile*> (farr, numFileCentBins);
  //Delete1DArray <TTree*> (tarr, numFileCentBins);

}


void MakeZTaggedMixingCatTrees () {
  TChain* c = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  const char* path = "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal";

  TFile**** farr = Get3DArray <TFile*> (numFileCentBins, 1, 1);
  TTree**** tarr = Get3DArray <TTree*> (numFileCentBins, 1, 1);
  //TFile**** farr = Get3DArray <TFile*> (numFileCentBins, numFileQ2Bins, numPsi2Bins);
  //TTree**** tarr = Get3DArray <TTree*> (numFileCentBins, numFileQ2Bins, numPsi2Bins);

  for (const auto& group : runGroups) {
    c->Reset ();
    for (int run : *(group.second)) {
      c->Add (Form ("%s/%i.root", path, run));
    }
    gSystem->Exec (Form ("mkdir -p %s/%s_sorted", path, group.first.c_str ()));

  //gSystem->Exec (Form ("mkdir -p %s/MixingTrees", path));

    for (int iCent = 1; iCent < numFileCentBins; iCent++) {
      for (int iQ2 = 0; iQ2 < 1; iQ2++) {
        for (int iPsi2 = 0; iPsi2 < 1; iPsi2++) {
      //for (int iQ2 = 0; iQ2 < numFileQ2Bins; iQ2++) {
      //  for (int iPsi2 = 0; iPsi2 < numPsi2Bins; iPsi2++) {
          farr[iCent][iQ2][iPsi2] = new TFile (Form ("%s/%s_sorted/tree_iCent%i_iQ2%i_iPsi2%i.root", path, group.first.c_str (), iCent, iQ2, iPsi2), "update");
          //farr[iCent][iQ2][iPsi2] = new TFile (Form ("%s/MixingTrees/tree_iCent%i_iQ2%i_iPsi2%i.root", path, iCent, iQ2, iPsi2), "update");
          tarr[iCent][iQ2][iPsi2] = (TTree*) farr[iCent][iQ2][iPsi2]->Get ("PbPbZTrackTree");
          if (!tarr[iCent][iQ2][iPsi2])
            tarr[iCent][iQ2][iPsi2] = (TTree*) c->CloneTree (0);
          tarr[iCent][iQ2][iPsi2]->SetDirectory (farr[iCent][iQ2][iPsi2]);
        }
      }
    }

    float fcal_et = 0, q2 = 0, psi2 = 0;
    c->SetBranchAddress ("fcal_et", &fcal_et);
    c->SetBranchAddress ("q2",      &q2);
    c->SetBranchAddress ("psi2",    &psi2);

    const int nEvts = c->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      c->GetEntry (iEvt);

      //fcal_et = fcal_et + fcal_et_shift;

      const short iCent = GetFileCentBin (fcal_et);
      if (iCent < 1 || iCent > numFileCentBins-1)
        continue;

      const short iQ2 = 0;
      const short iPsi2 = 0;

      //const short iQ2 = GetFileQ2Bin (q2);
      //if (iQ2 < 0 || iQ2 > numFileQ2Bins-1)
      //  continue;

      //const short iPsi2 = GetPsi2Bin (psi2);
      //if (iPsi2 < 0 || iPsi2 > numPsi2Bins-1)
      //  continue;

      tarr[iCent][iQ2][iPsi2]->Fill ();
    }
    cout << endl;
    c->Reset ();
    
    for (int iCent = 1; iCent < numFileCentBins; iCent++) {
      for (int iQ2 = 0; iQ2 < 1; iQ2++) {
        for (int iPsi2 = 0; iPsi2 < 1; iPsi2++) {
      //for (int iQ2 = 0; iQ2 < numFileQ2Bins; iQ2++) {
      //  for (int iPsi2 = 0; iPsi2 < numPsi2Bins; iPsi2++) {
          farr[iCent][iQ2][iPsi2]->cd ();
          tarr[iCent][iQ2][iPsi2]->Write ("", TObject :: kOverwrite);
          farr[iCent][iQ2][iPsi2]->Close ();
        }
      }
    }

  }
}




void MakeMCMixingCatTrees () {
  TChain* c = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  const char* path = "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal";

  c->Add (Form ("%s/PbPb18_*Zee.root", path));
  c->Add (Form ("%s/PbPb18_*Zmumu.root", path));
  MakeMixingCatTrees (path, c, -24.78);

  //gSystem->Exec (Form ("hadd -f %s/MixingTrees/tree_pp.root %s/pp_Zee.root %s/pp_Zmumu.root", path, path, path));
}




void MakeMinbiasMixingCatTrees () {
  TChain* c = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  const char* path = "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MinbiasAnalysis/Nominal";

  TFile**** farr = Get3DArray <TFile*> (numFileCentBins, 1, 1);
  TTree**** tarr = Get3DArray <TTree*> (numFileCentBins, 1, 1);
  //TFile**** farr = Get3DArray <TFile*> (numFileCentBins, numFileQ2Bins, numPsi2Bins);
  //TTree**** tarr = Get3DArray <TTree*> (numFileCentBins, numFileQ2Bins, numPsi2Bins);

  for (const auto& group : runGroups) {
    c->Reset ();
    for (int run : *(group.second)) {
      c->Add (Form ("%s/%i_cc.root", path, run));
      c->Add (Form ("%s/%i_pc.root", path, run));
    }

    gSystem->Exec (Form ("mkdir -p %s/%s_sorted", path, group.first.c_str ()));

    for (int iCent = 1; iCent < numFileCentBins; iCent++) {
      for (int iQ2 = 0; iQ2 < 1; iQ2++) {
        for (int iPsi2 = 0; iPsi2 < 1; iPsi2++) {
      //for (int iQ2 = 0; iQ2 < numFileQ2Bins; iQ2++) {
      //  for (int iPsi2 = 0; iPsi2 < numPsi2Bins; iPsi2++) {
          farr[iCent][iQ2][iPsi2] = new TFile (Form ("%s/%s_sorted/tree_iCent%i_iQ2%i_iPsi2%i.root", path, group.first.c_str (), iCent, iQ2, iPsi2), "update");
          tarr[iCent][iQ2][iPsi2] = (TTree*) farr[iCent][iQ2][iPsi2]->Get ("PbPbZTrackTree");
          if (!tarr[iCent][iQ2][iPsi2])
            tarr[iCent][iQ2][iPsi2] = (TTree*) c->CloneTree (0);
          tarr[iCent][iQ2][iPsi2]->SetDirectory (farr[iCent][iQ2][iPsi2]);
        }
      }
    }

    float fcal_et = 0, q2 = 0, psi2 = 0;
    c->SetBranchAddress ("fcal_et", &fcal_et);
    c->SetBranchAddress ("q2",      &q2);
    c->SetBranchAddress ("psi2",    &psi2);

    const int nEvts = c->GetEntries ();
    for (int iEvt = 0; iEvt < nEvts; iEvt++) {
      if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
        cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
      c->GetEntry (iEvt);

      //fcal_et = fcal_et + fcal_et_shift;

      const short iCent = GetFileCentBin (fcal_et);
      if (iCent < 1 || iCent > numFileCentBins-1)
        continue;

      const short iQ2 = 0;
      const short iPsi2 = 0;

      //const short iQ2 = GetFileQ2Bin (q2);
      //if (iQ2 < 0 || iQ2 > numFileQ2Bins-1)
      //  continue;

      //const short iPsi2 = GetPsi2Bin (psi2);
      //if (iPsi2 < 0 || iPsi2 > numPsi2Bins-1)
      //  continue;

      tarr[iCent][iQ2][iPsi2]->Fill ();
    }
    cout << endl;
    c->Reset ();
  
    for (int iCent = 1; iCent < numFileCentBins; iCent++) {
      for (int iQ2 = 0; iQ2 < 1; iQ2++) {
        for (int iPsi2 = 0; iPsi2 < 1; iPsi2++) {
      //for (int iQ2 = 0; iQ2 < numFileQ2Bins; iQ2++) {
      //  for (int iPsi2 = 0; iPsi2 < numPsi2Bins; iPsi2++) {
          farr[iCent][iQ2][iPsi2]->cd ();
          tarr[iCent][iQ2][iPsi2]->Write ("", TObject :: kOverwrite);
          farr[iCent][iQ2][iPsi2]->Close ();
        }
      }
    }

  }
}




void MakeAllCatTrees () {
  //MakeMCMixingCatTrees ();
  MakeMinbiasMixingCatTrees ();
}


#endif
