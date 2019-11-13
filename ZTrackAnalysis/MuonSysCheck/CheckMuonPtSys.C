#include "../Interactive/Params.h"
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>

void CheckMuonPtSys () {
  TFile* inFile = new TFile ("mc16_5TeV.root", "read");

  TTree* inTree = (TTree*) inFile->Get ("bush");

  int mn = 0;
  float mpt[40];
  float meta[40];
  float mphi[40];
  vector <vector <double>>* mptsys = nullptr;
  vector <vector <double>>* metasys = nullptr;
  vector <vector <double>>* mphisys = nullptr;

  inTree->SetBranchAddress ("muon_n", &mn);
  inTree->SetBranchAddress ("muon_pt", mpt);
  inTree->SetBranchAddress ("muon_eta", meta);
  inTree->SetBranchAddress ("muon_phi", mphi);
  inTree->SetBranchAddress ("muon_pt_sys", &mptsys);
  inTree->SetBranchAddress ("muon_eta_sys", &metasys);
  inTree->SetBranchAddress ("muon_phi_sys", &mphisys);

  TH1D* h_dmpt[10];
  for (int isys = 0; isys < 10; isys++) {
    h_dmpt[isys] = new TH1D (Form ("h_dmpt_sys%i", isys), ";#Delta#it{p}_{T}^{#mu} [GeV];Counts", 40, -2, 2);
    h_dmpt[isys]->Sumw2 ();
  }

  for (int i = 0; i < inTree->GetEntries (); i++) {
    inTree->GetEntry (i);

    for (int im = 0; im < mn; im++) {
      for (int isys = 0; isys < 10; isys++) {
        h_dmpt[isys]->Fill (mpt[im] - (*mptsys)[im][isys]);
      }
    }
  }

  TCanvas* c = new TCanvas ("c", "", 600, 600);

  for (int isys = 0; isys < 10; isys++) {
    h_dmpt[isys]->SetLineColor (colors[isys]);
    h_dmpt[isys]->SetMarkerColor (colors[isys]);
    h_dmpt[isys]->Draw (isys == 0 ? "hist" : "same hist");
  }
  
}
