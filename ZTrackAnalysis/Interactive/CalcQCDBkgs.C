#ifndef __CalcQCDBkgs_C__
#define __CalcQCDBkgs_C__

#include "Params.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

void CalcQCDBkgs () {

  TFile* dataFile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/BkgEstimator/Nominal/data18hi.root", "read");
  TTree* dataTree = (TTree*) dataFile->Get ("PbPbZTree");

  TFile* mcFile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/BkgEstimator/Nominal/mc.root", "read");
  TTree* mcTree = (TTree*) mcFile->Get ("PbPbZTree");

  bool isQCDBkg = false;
  bool isSameSign = false;
  bool isEE = false;
  float z_m = 0;
  float event_weight = 0;

  dataTree->SetBranchAddress ("isQCDBkg",     &isQCDBkg);
  dataTree->SetBranchAddress ("isSameSign",   &isSameSign);
  dataTree->SetBranchAddress ("isEE",         &isEE);
  dataTree->SetBranchAddress ("z_m",          &z_m);
  dataTree->SetBranchAddress ("event_weight", &event_weight);

  mcTree->SetBranchAddress ("isQCDBkg",   &isQCDBkg);
  mcTree->SetBranchAddress ("isSameSign", &isSameSign);
  mcTree->SetBranchAddress ("isEE",       &isEE);
  mcTree->SetBranchAddress ("z_m",        &z_m);
  mcTree->SetBranchAddress ("event_weight", &event_weight);

  float n_os_mumu = 0;
  float n_ss_mumu = 0;

  TH1D* h_zm_ee_os_data = new TH1D ("h_zm_ee_os_data", ";#it{m}_{#it{ee}} [GeV];Events", 75, 60, 135);
  TH1D* h_zm_ee_ss_data = new TH1D ("h_zm_ee_ss_data", ";#it{m}_{#it{ee}} [GeV];Events", 75, 60, 135);
  TH1D* h_zm_ee_qcd_data = new TH1D ("h_zm_ee_qcd_data", ";#it{m}_{#it{ee}} [GeV];Events", 75, 60, 135);
  TH1D* h_zm_ee_os_mc = new TH1D ("h_zm_ee_os_mc", ";#it{m}_{#it{ee}} [GeV];Events", 75, 60, 135);
  TH1D* h_zm_ee_ss_mc = new TH1D ("h_zm_ee_ss_mc", ";#it{m}_{#it{ee}} [GeV];Events", 75, 60, 135);
  TH1D* h_zm_ee_qcd_mc = new TH1D ("h_zm_ee_qcd_mc", ";#it{m}_{#it{ee}} [GeV];Events", 75, 60, 135);

  for (int iEvt = 0; iEvt < dataTree->GetEntries (); iEvt++) {
    dataTree->GetEntry (iEvt);

    if (!isEE) {
      if (76 < z_m && z_m < 106) {
        if (isSameSign) n_ss_mumu += event_weight;
        else            n_os_mumu += event_weight;
      }
    }

    else {
      if (isQCDBkg)         h_zm_ee_qcd_data->Fill (z_m, event_weight);
      else if (isSameSign)  h_zm_ee_ss_data->Fill (z_m, event_weight);
      else if (!isSameSign) h_zm_ee_os_data->Fill (z_m, event_weight);
    }
  }

  const float muonQCDBkgRate = n_ss_mumu / (n_ss_mumu + n_os_mumu);
  cout << "QCD bkg in Z->mumu = " << muonQCDBkgRate * 100 << "%" << endl;


  for (int iEvt = 0; iEvt < mcTree->GetEntries (); iEvt++) {
    mcTree->GetEntry (iEvt);

    if (isEE) {
      if (isQCDBkg)         h_zm_ee_qcd_mc->Fill (z_m, event_weight);
      else if (isSameSign)  h_zm_ee_ss_mc->Fill (z_m, event_weight);
      else if (!isSameSign) h_zm_ee_os_mc->Fill (z_m, event_weight);
    }
  }


  const float mcNormFactor = h_zm_ee_os_data->Integral () / h_zm_ee_os_mc->Integral ();
  h_zm_ee_qcd_mc->Scale (mcNormFactor);
  h_zm_ee_ss_mc->Scale (mcNormFactor);
  h_zm_ee_os_mc->Scale (mcNormFactor);


  TH1D* h_zm_ee_ss_data_minus_mc = (TH1D*) h_zm_ee_ss_data->Clone ("h_zm_ee_ss_data_minus_mc");
  h_zm_ee_ss_data_minus_mc->Add (h_zm_ee_ss_mc, -1);
  
  const float qcdNormFactor = h_zm_ee_ss_data_minus_mc->Integral (h_zm_ee_ss_data_minus_mc->FindBin (60), h_zm_ee_ss_data_minus_mc->FindBin (70)) / h_zm_ee_qcd_data->Integral (h_zm_ee_qcd_data->FindBin (60), h_zm_ee_qcd_data->FindBin (70));  
  h_zm_ee_qcd_data->Scale (qcdNormFactor);
  
  const float electronQCDBkgRate = h_zm_ee_qcd_data->Integral (h_zm_ee_qcd_data->FindBin (76), h_zm_ee_qcd_data->FindBin (106) - 1) / h_zm_ee_os_data->Integral (h_zm_ee_os_data->FindBin (76), h_zm_ee_os_data->FindBin (106) - 1);
  cout << "QCD bkg in Z->ee   = " << electronQCDBkgRate * 100 << "%" << endl;
  
}

#endif
