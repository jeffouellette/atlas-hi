#ifndef __MakeMixingCatTrees_C__
#define __MakeMixingCatTrees_C__

#include "Params.h"

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

void MakeMixingCatTrees (const char* path, TChain* c, const float fcal_et_shift = 0) {

  gSystem->Exec (Form ("mkdir -p %s/MixingTrees", path));

  TTree** tarr = Get1DArray <TTree*> (numSuperFineCentBins);
  TFile** farr = Get1DArray <TFile*> (numSuperFineCentBins);

  for (int iCent = 0; iCent < numSuperFineCentBins; iCent++) {
    farr[iCent] = new TFile (Form ("%s/MixingTrees/tree_iCent%i.root", path, iCent), "update");
    tarr[iCent] = (TTree*) farr[iCent]->Get ("PbPbZTrackTree");
    if (!tarr[iCent])
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

    //fcal_et = fcal_et + fcal_et_shift;

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
  MakeMixingCatTrees (path, c, -24.78);

  gSystem->Exec (Form ("hadd -f %s/MixingTrees/tree_iCent0.root %s/pp_Zee.root %s/pp_Zmumu.root", path, path, path));
}

void MakeMinbiasMixingCatTrees () {
  TChain* c = new TChain ("PbPbZTrackTree", "PbPbZTrackTree");
  const char* path = "/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MinbiasAnalysis/Nominal";

  vector<string> files = {"GroupA.root", "GroupB.root", "GroupC.root", "GroupD.root", "GroupE.root", "GroupF.root", "GroupG.root", "GroupH.root", "GroupI.root", "GroupJ.root", "GroupK.root"};
  for (string f : files) {
    c->Add (Form ("%s/%s", path, f.c_str ()));
    MakeMixingCatTrees (path, c);
    c->Reset ();
  }
  //c->Add (Form ("%s/36*_pc.root", path));
  //c->Add (Form ("%s/36*_cc.root", path));
  //MakeMixingCatTrees (path, c);

  gSystem->Exec (Form ("hadd -f %s/MixingTrees/tree_iCent0.root %s/340644.root %s/340683.root %s/340697.root %s/340718.root %s/340814.root %s/340849.root %s/340850.root %s/340910.root %s/340925.root %s/340973.root %s/341027.root %s/341123.root %s/341184.root", path, path, path, path, path, path, path, path, path, path, path, path, path, path));
}

void MakeAllCatTrees () {
  //MakeMCMixingCatTrees ();
  MakeMinbiasMixingCatTrees ();
}
#endif
