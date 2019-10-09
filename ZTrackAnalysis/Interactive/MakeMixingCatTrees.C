#ifndef __MakeMixingCatTrees_C__
#define __MakeMixingCatTrees_C__

#include "Params.h"

void MakeMixingCatTrees (const char* path, TChain* c) {

  gSystem->Exec (Form ("mkdir -p %s/MixingTrees", path));

  TTree** tarr = Get1DArray <TTree*> (numSuperFineCentBins);
  TFile** farr = Get1DArray <TFile*> (numSuperFineCentBins);

  for (int iCent = 0; iCent < numSuperFineCentBins; iCent++) {
    farr[iCent] = new TFile (Form ("%s/MixingTrees/tree_iCent%i.root", path, iCent), "recreate");
    tarr[iCent] = (TTree*) c->CloneTree (0);
    //tarr[iCent]->SetName (Form ("mixingTree_iCent%i", iCent));
    tarr[iCent]->SetDirectory (farr[iCent]);
  }

  float fcal_et = 0;
  c->SetBranchAddress ("fcal_et", &fcal_et);

  const int nEvts = c->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {
    if (nEvts > 0 && iEvt % (nEvts / 100) == 0)
      cout << iEvt / (nEvts / 100) << "\% done...\r" << flush;
    c->GetEntry (iEvt);

    const short iCent = GetSuperFineCentBin (fcal_et);
    if (iCent < 1 || iCent > numSuperFineCentBins-1)
      continue;

    tarr[iCent]->Fill ();
  }

  
  for (int iCent = 1; iCent < numSuperFineCentBins; iCent++) {
    farr[iCent]->cd ();
    tarr[iCent]->Write ("", TObject :: kOverwrite);
    farr[iCent]->Close ();
  }

  //Delete1DArray <TFile*> (farr, numSuperFineCentBins);
  //Delete1DArray <TTree*> (tarr, numSuperFineCentBins);

}

void MakeMCMixingCatTrees () {
  TChain* c = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  const char* path = "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MCAnalysis/Nominal";
  c->Add (Form ("%s/PbPb18_*Zee.root", path));
  c->Add (Form ("%s/PbPb18_*Zmumu.root", path));
  MakeMixingCatTrees (path, c);

  //gSystem->Exec (Form ("hadd -f %s/MixingTrees/tree_iCent0.root pp_Zee.root pp_Zmumu.root", path));
}

void MakeMinbiasMixingCatTrees () {
  TChain* c = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  const char* path = "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MinbiasAnalysis/Nominal";
  c->Add (Form ("%s/36*_pc.root", path));
  c->Add (Form ("%s/36*_cc.root", path));
  MakeMixingCatTrees (path, c);

  //gSystem->Exec (Form ("hadd -f %s/MixingTrees/tree_iCent0.root 340644.root 340683.root 340697.root 340718.root 340814.root 340849.root 340850.root 340910.root 340925.root 340973.root 341027.root 341123.root 341184.root", path));
}

void MakeAllCatTrees () {
  //MakeMCMixingCatTrees ();
  MakeMinbiasMixingCatTrees ();
}
#endif
