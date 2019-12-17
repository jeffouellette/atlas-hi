#include "Params.h"

#include <TH2.h>
#include <TFile.h>
#include <TTree.h>

#include <string>

void CalculateMixingRatios () {

  const string names[11] = {"GroupA.root", "GroupB.root", "GroupC.root", "GroupD.root", "GroupE.root", "GroupF.root", "GroupG.root", "GroupH.root", "GroupI.root", "GroupJ.root", "GroupK.root"};

  TH1D* h_fcal_mb[11] = {};
  TH1D* h_fcal_data[11] = {};
  TH1D* h_fcal_ratios[11] = {};

  for (int i = 0; i < 11; i++) {
    h_fcal_mb[i] = new TH1D (Form ("h_fcal_mb_%i", i), "#Sigma#it{E}_{T}^{FCal} [GeV]", numSuperFineCentBins-1, superFineCentBins);
    h_fcal_mb[i]->Sumw2 ();
    h_fcal_data[i] = new TH1D (Form ("h_fcal_data_%i", i), "#Sigma#it{E}_{T}^{FCal} [GeV]", numSuperFineCentBins-1, superFineCentBins);
    h_fcal_data[i]->Sumw2 ();
  }


  TDirectory* _cd = (TDirectory*) gDirectory;
  for (int i = 0; i < 11; i++) {

    TFile* mbFile = new TFile (Form ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MinbiasAnalysis/Nominal/%s", names[i].c_str ()), "read");
    TTree* mbTree = (TTree*) mbFile->Get ("PbPbZTrackTree");

    _cd->cd ();

    float fcal_et = 0;
    mbTree->SetBranchAddress ("fcal_et", &fcal_et);
    for (int j = 0; j < mbTree->GetEntries (); j++) {
      mbTree->GetEntry (j);
      h_fcal_mb[i]->Fill (fcal_et);
    }

    mbFile->Close ();
  }

  for (int i = 0; i < 11; i++) {

    TFile* dataFile = new TFile (Form ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/%s", names[i].c_str ()), "read");
    TTree* dataTree = (TTree*) dataFile->Get ("PbPbZTrackTree");

    _cd->cd ();

    float fcal_et = 0;
    dataTree->SetBranchAddress ("fcal_et", &fcal_et);
    for (int j = 0; j < dataTree->GetEntries (); j++) {
      dataTree->GetEntry (j);
      h_fcal_data[i]->Fill (fcal_et);
    }

    dataFile->Close ();

    h_fcal_mb[i]->SetLineColor (kGray);
    h_fcal_data[i]->SetLineColor (colors[i+1]);

    h_fcal_ratios[i] = (TH1D*) h_fcal_mb[i]->Clone (Form ("h_fcal_ratios_%i", i));
    h_fcal_ratios[i]->Divide (h_fcal_data[i]);
  }

  for (int i = 0; i < 11; i++) {
    h_fcal_mb[i]->SetLineColor (kGray);
    h_fcal_data[i]->SetLineColor (colors[i+1]);

    h_fcal_ratios[i] = (TH1D*) h_fcal_mb[i]->Clone (Form ("h_fcal_ratios_%i", i));
    h_fcal_ratios[i]->Divide (h_fcal_data[i]);
  }

  TFile* outFile = new TFile ("MixingRatios.root", "recreate");
  for (int i = 0; i < 11; i++) {
    h_fcal_mb[i]->Write ();
    h_fcal_data[i]->Write ();
    h_fcal_ratios[i]->Write ();
  }
  outFile->Close ();

  for (int i = 0; i < 11; i++) {
    h_fcal_mb[i]->Draw (i == 0 ? "hist" : "same hist");
    h_fcal_data[i]->Draw ("same hist");
    h_fcal_ratios[i]->Draw ("same hist");
  }
}
