#ifndef __CalcQCDBkgs_C__
#define __CalcQCDBkgs_C__

#include "Params.h"

#include "AtlasUtils.h"

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>

#include <fstream>

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

  mcTree->SetBranchAddress ("isQCDBkg",     &isQCDBkg);
  mcTree->SetBranchAddress ("isSameSign",   &isSameSign);
  mcTree->SetBranchAddress ("isEE",         &isEE);
  mcTree->SetBranchAddress ("z_m",          &z_m);
  mcTree->SetBranchAddress ("event_weight", &event_weight);

  float n_os_mumu = 0;
  float n_ss_mumu = 0;

  TH1D* h_zm_mumu_os_data = new TH1D ("h_zm_mumu_os_data", ";m_{#it{#mu#mu}} [GeV];Counts / GeV", 90, 60, 150);
  h_zm_mumu_os_data->Sumw2 ();
  TH1D* h_zm_mumu_ss_data = new TH1D ("h_zm_mumu_ss_data", ";m_{#it{#mu#mu}} [GeV];Counts / GeV", 90, 60, 150);
  h_zm_mumu_ss_data->Sumw2 ();
  TH1D* h_zm_mumu_os_mc = new TH1D ("h_zm_mumu_os_mc", ";m_{#it{#mu#mu}} [GeV];Counts / GeV", 90, 60, 150);
  h_zm_mumu_os_mc->Sumw2 ();
  TH1D* h_zm_mumu_ss_mc = new TH1D ("h_zm_mumu_ss_mc", ";m_{#it{#mu#mu}} [GeV];Counts / GeV", 90, 60, 150);
  h_zm_mumu_ss_mc->Sumw2 ();

  TH1D* h_zm_ee_os_data = new TH1D ("h_zm_ee_os_data", ";m_{#it{ee}} [GeV];Counts / GeV", 90, 60, 150);
  h_zm_ee_os_data->Sumw2 ();
  TH1D* h_zm_ee_ss_data = new TH1D ("h_zm_ee_ss_data", ";m_{#it{ee}} [GeV];Counts / GeV", 90, 60, 150);
  h_zm_ee_ss_data->Sumw2 ();
  TH1D* h_zm_ee_qcd_data = new TH1D ("h_zm_ee_qcd_data", ";m_{#it{ee}} [GeV];Counts / GeV", 90, 60, 150);
  h_zm_ee_qcd_data->Sumw2 ();
  TH1D* h_zm_ee_os_mc = new TH1D ("h_zm_ee_os_mc", ";m_{#it{ee}} [GeV];Counts / GeV", 90, 60, 150);
  h_zm_ee_os_mc->Sumw2 ();
  TH1D* h_zm_ee_ss_mc = new TH1D ("h_zm_ee_ss_mc", ";m_{#it{ee}} [GeV];Counts / GeV", 90, 60, 150);
  h_zm_ee_ss_mc->Sumw2 ();
  //TH1D* h_zm_ee_qcd_mc = new TH1D ("h_zm_ee_qcd_mc", ";#it{m}_{#it{ee}} [GeV];Counts / GeV", 90, 60, 150);
  //h_zm_ee_qcd_mc->Sumw2 ();

  for (int iEvt = 0; iEvt < dataTree->GetEntries (); iEvt++) {
    dataTree->GetEntry (iEvt);

    if (!isEE) {
      if (isSameSign) h_zm_mumu_ss_data->Fill (z_m);
      else            h_zm_mumu_os_data->Fill (z_m);
      if (76 < z_m && z_m < 106) {
        if (isSameSign) n_ss_mumu += 1;
        else            n_os_mumu += 1;
      }
    }

    else {
      if (isQCDBkg)         h_zm_ee_qcd_data->Fill (z_m);
      else if (isSameSign)  h_zm_ee_ss_data->Fill (z_m);
      else if (!isSameSign) h_zm_ee_os_data->Fill (z_m);
    }
  }

  for (int iEvt = 0; iEvt < mcTree->GetEntries (); iEvt++) {
    mcTree->GetEntry (iEvt);

    if (isEE) {
      if (isSameSign) h_zm_ee_ss_mc->Fill (z_m);
      else            h_zm_ee_os_mc->Fill (z_m);
    }
    else {
      if (isSameSign) h_zm_mumu_ss_mc->Fill (z_m);
      else            h_zm_mumu_os_mc->Fill (z_m);
    }
  }




  ofstream output;
  output.open ("QCDbkgs.out");

  h_zm_mumu_ss_data->Rebin (2);
  h_zm_mumu_ss_mc->Rebin (2);
  const float mumuMCNormFactor = h_zm_mumu_os_data->Integral () / h_zm_mumu_os_mc->Integral ();
  h_zm_mumu_ss_mc->Scale (mumuMCNormFactor);
  h_zm_mumu_os_mc->Scale (mumuMCNormFactor);

  const float muonQCDBkgRate = n_ss_mumu / (n_ss_mumu + n_os_mumu);
  cout << "QCD bkg in Z->mumu = " << muonQCDBkgRate * 100 << "%" << endl;
  output << "QCD bkg in Z->mumu = " << muonQCDBkgRate * 100 << "%" << endl;


  h_zm_ee_ss_data->Rebin (2);
  h_zm_ee_qcd_data->Rebin (2);
  h_zm_ee_ss_mc->Rebin (2);
  const float eeMCNormFactor = h_zm_ee_os_data->Integral () / h_zm_ee_os_mc->Integral ();
  h_zm_ee_ss_mc->Scale (eeMCNormFactor);
  h_zm_ee_os_mc->Scale (eeMCNormFactor);

  TH1D* h_zm_ee_ss_data_minus_mc = (TH1D*) h_zm_ee_ss_data->Clone ("h_zm_ee_ss_data_minus_mc");
  h_zm_ee_ss_data_minus_mc->Add (h_zm_ee_ss_mc, -1);
  
  const float qcdNormFactor = h_zm_ee_ss_data_minus_mc->Integral (h_zm_ee_ss_data_minus_mc->FindBin (60), h_zm_ee_ss_data_minus_mc->FindBin (70)) / h_zm_ee_qcd_data->Integral (h_zm_ee_qcd_data->FindBin (60), h_zm_ee_qcd_data->FindBin (70));  
  h_zm_ee_qcd_data->Scale (qcdNormFactor);
  
  const float electronQCDBkgRate = h_zm_ee_qcd_data->Integral (h_zm_ee_qcd_data->FindBin (76), h_zm_ee_qcd_data->FindBin (106) - 1) / h_zm_ee_os_data->Integral (h_zm_ee_os_data->FindBin (76), h_zm_ee_os_data->FindBin (106) - 1);
  cout << "QCD bkg in Z->ee   = " << electronQCDBkgRate * 100 << "%" << endl;
  output << "QCD bkg in Z->ee   = " << electronQCDBkgRate * 100 << "%" << endl;

  TH1D* h_zm_ee_os_mc_plus_qcd = (TH1D*) h_zm_ee_os_mc->Clone ("h_zm_ee_os_mc_plus_qcd");
  for (int iX = 1; iX <= h_zm_ee_os_mc_plus_qcd->GetNbinsX (); iX++) {
    h_zm_ee_os_mc_plus_qcd->SetBinContent (iX, h_zm_ee_os_mc_plus_qcd->GetBinContent (iX) + 0.5 * h_zm_ee_qcd_data->GetBinContent (h_zm_ee_qcd_data->FindBin (h_zm_ee_os_mc_plus_qcd->GetBinCenter (iX))));
  }

  output.close ();




  h_zm_mumu_os_mc->Scale (1., "width");
  h_zm_mumu_ss_mc->Scale (1., "width");
  h_zm_mumu_os_data->Scale (1., "width");
  h_zm_mumu_ss_data->Scale (1., "width");

  h_zm_ee_os_mc->Scale (1., "width");
  h_zm_ee_ss_mc->Scale (1., "width");
  h_zm_ee_os_data->Scale (1., "width");
  h_zm_ee_ss_data->Scale (1., "width");
  h_zm_ee_qcd_data->Scale (1., "width");
  h_zm_ee_ss_data_minus_mc->Scale (1., "width");
  h_zm_ee_os_mc_plus_qcd->Scale (1., "width");




  TCanvas* c_zm_mumu = new TCanvas ("c_zm_mumu", "", 600, 600);
  gPad->SetLogy ();
  h_zm_mumu_os_mc->SetFillColor (kYellow-9);
  h_zm_mumu_ss_mc->SetLineColor (kAzure-6);
  h_zm_mumu_ss_mc->SetLineStyle (2);
  h_zm_mumu_os_data->SetMarkerColor (kBlack);
  h_zm_mumu_os_data->SetLineColor (kBlack);
  h_zm_mumu_os_data->SetMarkerStyle (kFullCircle);
  h_zm_mumu_os_data->SetMarkerSize (0.8);
  h_zm_mumu_ss_data->SetMarkerColor (kRed+1);
  h_zm_mumu_ss_data->SetLineColor (kRed+1);
  h_zm_mumu_ss_data->SetMarkerStyle (kFullCircle);
  h_zm_mumu_ss_data->SetMarkerSize (0.8);

  h_zm_mumu_os_mc->GetYaxis ()->SetRangeUser (0.01, 3e4);
  h_zm_mumu_os_mc->GetXaxis ()->SetLabelSize (0.9*h_zm_mumu_os_mc->GetXaxis ()->GetLabelSize ());
  h_zm_mumu_os_mc->GetYaxis ()->SetLabelSize (0.9*h_zm_mumu_os_mc->GetYaxis ()->GetLabelSize ());

  h_zm_mumu_os_mc->Draw ("hist");
  h_zm_mumu_ss_mc->Draw ("hist same");
  h_zm_mumu_os_data->Draw ("e1 same");
  h_zm_mumu_ss_data->Draw ("e1 same");

  myText (0.20, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.042);
  myText (0.20, 0.84, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.036);
  myText (0.20, 0.79, kBlack, "2018 data, 1.4 nb^{-1}", 0.036);

  myMarkerText  (0.67, 0.88, kBlack, kFullCircle, "Opposite-sign data", 0.8, 0.032);
  myOnlyBoxText (0.67, 0.84, 1, kYellow-9, kBlack, 1, "Opposite-sign MC", 0.032, 1001, 1);
  myMarkerText  (0.67, 0.80, kRed+1, kFullCircle, "Same-sign data", 0.8, 0.032);
  myLineText    (0.67, 0.76, kAzure-6, 2, "Same-sign MC", 1, 0.032);

  c_zm_mumu->SaveAs ("Plots/ZmumuQCDBkg.pdf");
  




  TCanvas* c_zm_ee = new TCanvas ("c_zm_ee", "", 600, 600);
  gPad->SetLogy ();
  h_zm_ee_os_mc->SetFillColor (kYellow-9);
  h_zm_ee_ss_mc->SetLineColor (kAzure-6);
  h_zm_ee_ss_mc->SetLineStyle (2);
  h_zm_ee_os_data->SetMarkerColor (kBlack);
  h_zm_ee_os_data->SetLineColor (kBlack);
  h_zm_ee_os_data->SetMarkerStyle (kFullCircle);
  h_zm_ee_os_data->SetMarkerSize (0.8);
  h_zm_ee_ss_data->SetMarkerColor (kRed+1);
  h_zm_ee_ss_data->SetLineColor (kRed+1);
  h_zm_ee_ss_data->SetMarkerStyle (kFullCircle);
  h_zm_ee_ss_data->SetMarkerSize (0.8);
  h_zm_ee_qcd_data->SetLineColor (kGreen+3);
  h_zm_ee_ss_data_minus_mc->SetMarkerStyle (kOpenCircle);
  h_zm_ee_ss_data_minus_mc->SetMarkerSize (0.8);
  h_zm_ee_ss_data_minus_mc->SetMarkerColor (kMagenta);
  h_zm_ee_ss_data_minus_mc->SetLineColor (kMagenta);
  h_zm_ee_os_mc_plus_qcd->SetLineColor (kRed+1);
  h_zm_ee_os_mc_plus_qcd->SetLineStyle (2);

  h_zm_ee_os_mc->GetYaxis ()->SetRangeUser (0.08, 2e4);
  h_zm_ee_os_mc->GetXaxis ()->SetLabelSize (0.9*h_zm_ee_os_mc->GetXaxis ()->GetLabelSize ());
  h_zm_ee_os_mc->GetYaxis ()->SetLabelSize (0.9*h_zm_ee_os_mc->GetYaxis ()->GetLabelSize ());

  h_zm_ee_os_mc->Draw ("hist");
  h_zm_ee_ss_mc->Draw ("hist same");
  h_zm_ee_os_data->Draw ("e1 same");
  h_zm_ee_ss_data->Draw ("e1 same");
  h_zm_ee_qcd_data->Draw ("hist same");
  h_zm_ee_ss_data_minus_mc->Draw ("e1 same");
  h_zm_ee_os_mc_plus_qcd->Draw ("hist same");

  myText (0.20, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.042);
  myText (0.20, 0.84, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.036);
  myText (0.20, 0.79, kBlack, "2018 data, 1.7 nb^{-1}", 0.036);

  myMarkerText  (0.67, 0.88, kBlack, kFullCircle, "Opposite-sign data", 0.8, 0.032);
  myOnlyBoxText (0.67, 0.84, 1, kYellow-9, kBlack, 1, "Opposite-sign MC", 0.032, 1001, 1);
  myMarkerText  (0.67, 0.80, kRed+1, kFullCircle, "Same-sign data", 0.8, 0.032);
  myLineText    (0.67, 0.76, kAzure-6, 2, "Same-sign MC", 1, 0.032);
  myMarkerText  (0.67, 0.72, kMagenta, kOpenCircle, "SS data #minus MC", 0.8, 0.032);
  myLineText    (0.67, 0.68, kGreen+3, 1, "QCD template", 1, 0.032);
  myLineText    (0.67, 0.64, kRed+1, 2, "MC Sig. + QCD Bkg.", 1, 0.032);

  c_zm_ee->SaveAs ("Plots/ZeeQCDBkg.pdf");

  
}

#endif
