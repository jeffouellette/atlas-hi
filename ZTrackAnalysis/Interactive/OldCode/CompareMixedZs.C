#ifndef __CompareMixedZs_C__
#define __CompareMixedZs_C__

#include "Params.h"

void CompareMixedZs () {

  TFile* f_nominal = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/DataAnalysis/Nominal/zMixFile.root", "read");
  TFile* f_electronlhmed = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/DataAnalysis/Variations/ElectronLHMediumWPVariation/zMixFile.root", "read");
  TFile* f_muontight = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/DataAnalysis/Variations/MuonTightWPVariation/zMixFile.root", "read");

  TTree* t_nominal = (TTree*) f_nominal->Get ("ppZTree");
  TTree* t_electronlhmed = (TTree*) f_electronlhmed->Get ("ppZTree");
  TTree* t_muontight = (TTree*) f_muontight->Get ("ppZTree");

  TH1D* h_electronlhloose_z_pt = new TH1D ("h_electronlhloose_z_pt", "#it{p}_{T}^{Z} [GeV];#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{p}_{T}}", 39, 5, 200);
  h_electronlhloose_z_pt->Sumw2 ();
  TH1D* h_muonmedium_z_pt = new TH1D ("h_muonmedium_z_pt", "#it{p}_{T}^{Z} [GeV];#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{p}_{T}}", 39, 5, 200);
  h_muonmedium_z_pt->Sumw2 ();
  TH1D* h_electronlhmed_z_pt = new TH1D ("h_electronlhmed_z_pt", "#it{p}_{T}^{Z} [GeV];#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{p}_{T}}", 39, 5, 200);
  h_electronlhmed_z_pt->Sumw2 ();
  TH1D* h_muontight_z_pt = new TH1D ("h_muontight_z_pt", ";#it{p}_{T}^{Z} [GeV];#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{p}_{T}}", 39, 5, 200);
  h_muontight_z_pt->Sumw2 ();
  
  TH1D* h_electronlhloose_z_y = new TH1D ("h_electronlhloose_z_y", ";#it{y}_{Z};#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{y}}", 25, -2.5, 2.5);
  h_electronlhloose_z_y->Sumw2 ();
  TH1D* h_muonmedium_z_y = new TH1D ("h_muonmedium_z_y", ";#it{y}_{Z};#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{y}}", 25, -2.5, 2.5);
  h_muonmedium_z_y->Sumw2 ();
  TH1D* h_electronlhmed_z_y = new TH1D ("h_electronlhmed_z_y", ";#it{y}_{Z};#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{y}}", 25, -2.5, 2.5);
  h_electronlhmed_z_y->Sumw2 ();
  TH1D* h_muontight_z_y = new TH1D ("h_muontight_z_y", ";#it{y}_{Z};#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{y}}", 25, -2.5, 2.5);
  h_muontight_z_y->Sumw2 ();

  TH1D* h_electronlhloose_z_phi = new TH1D ("h_electronlhloose_z_phi", ";#phi_{Z};#frac{1}{N_{Z}} #frac{dN_{Z}}{d#phi}", 24, -pi, pi);
  h_electronlhloose_z_phi->Sumw2 ();
  TH1D* h_muonmedium_z_phi = new TH1D ("h_muonmedium_z_phi", ";#phi_{Z};#frac{1}{N_{Z}} #frac{dN_{Z}}{d#phi}", 24, -pi, pi);
  h_muonmedium_z_phi->Sumw2 ();
  TH1D* h_electronlhmed_z_phi = new TH1D ("h_electronlhmed_z_phi", ";#phi_{Z};#frac{1}{N_{Z}} #frac{dN_{Z}}{d#phi}", 24, -pi, pi);
  h_electronlhmed_z_phi->Sumw2 ();
  TH1D* h_muontight_z_phi = new TH1D ("h_muontight_z_phi", ";#phi_{Z};#frac{1}{N_{Z}} #frac{dN_{Z}}{d#phi}", 24, -pi, pi);
  h_muontight_z_phi->Sumw2 ();

  TH1D* h_electronlhloose_z_m = new TH1D ("h_electronlhloose_z_m", ";#it{m}_{Z} [GeV];#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{m}}", 60, 76, 106);
  h_electronlhloose_z_m->Sumw2 ();
  TH1D* h_muonmedium_z_m = new TH1D ("h_muonmedium_z_m", ";#it{m}_{Z} [GeV];#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{m}}", 60, 76, 106);
  h_muonmedium_z_m->Sumw2 ();
  TH1D* h_electronlhmed_z_m = new TH1D ("h_electronlhmed_z_m", ";#it{m}_{Z} [GeV];#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{m}}", 60, 76, 106);
  h_electronlhmed_z_m->Sumw2 ();
  TH1D* h_muontight_z_m = new TH1D ("h_muontight_z_m", ";#it{m}_{Z} [GeV];#frac{1}{N_{Z}} #frac{dN_{Z}}{d#it{m}}", 60, 76, 106);
  h_muontight_z_m->Sumw2 ();

  t_nominal->Draw ("z_pt>>h_electronlhloose_z_pt", "isEE");
  t_nominal->Draw ("z_pt>>h_muonmedium_z_pt", "!isEE");
  t_electronlhmed->Draw ("z_pt>>h_electronlhmed_z_pt", "isEE");
  t_muontight->Draw ("z_pt>>h_muontight_z_pt", "!isEE");
  t_nominal->Draw ("z_y>>h_electronlhloose_z_y", "isEE");
  t_nominal->Draw ("z_y>>h_muonmedium_z_y", "!isEE");
  t_electronlhmed->Draw ("z_y>>h_electronlhmed_z_y", "isEE");
  t_muontight->Draw ("z_y>>h_muontight_z_y", "!isEE");
  t_nominal->Draw ("z_phi>>h_electronlhloose_z_phi", "isEE");
  t_nominal->Draw ("z_phi>>h_muonmedium_z_phi", "!isEE");
  t_electronlhmed->Draw ("z_phi>>h_electronlhmed_z_phi", "isEE");
  t_muontight->Draw ("z_phi>>h_muontight_z_phi", "!isEE");
  t_nominal->Draw ("z_m>>h_electronlhloose_z_m", "isEE");
  t_nominal->Draw ("z_m>>h_muonmedium_z_m", "!isEE");
  t_electronlhmed->Draw ("z_m>>h_electronlhmed_z_m", "isEE");
  t_muontight->Draw ("z_m>>h_muontight_z_m", "!isEE");

  h_electronlhloose_z_pt->Scale (1./h_electronlhloose_z_pt->Integral (), "width");
  h_muonmedium_z_pt->Scale (1./h_muonmedium_z_pt->Integral (), "width");
  h_electronlhmed_z_pt->Scale (1./h_electronlhmed_z_pt->Integral (), "width");
  h_muontight_z_pt->Scale (1./h_muontight_z_pt->Integral (), "width");

  h_electronlhloose_z_y->Scale (1./h_electronlhloose_z_y->Integral (), "width");
  h_muonmedium_z_y->Scale (1./h_muonmedium_z_y->Integral (), "width");
  h_electronlhmed_z_y->Scale (1./h_electronlhmed_z_y->Integral (), "width");
  h_muontight_z_y->Scale (1./h_muontight_z_y->Integral (), "width");

  h_electronlhloose_z_phi->Scale (1./h_electronlhloose_z_phi->Integral (), "width");
  h_muonmedium_z_phi->Scale (1./h_muonmedium_z_phi->Integral (), "width");
  h_electronlhmed_z_phi->Scale (1./h_electronlhmed_z_phi->Integral (), "width");
  h_muontight_z_phi->Scale (1./h_muontight_z_phi->Integral (), "width");

  h_electronlhloose_z_m->Scale (1./h_electronlhloose_z_m->Integral (), "width");
  h_muonmedium_z_m->Scale (1./h_muonmedium_z_m->Integral (), "width");
  h_electronlhmed_z_m->Scale (1./h_electronlhmed_z_m->Integral (), "width");
  h_muontight_z_m->Scale (1./h_muontight_z_m->Integral (), "width");


  TCanvas* c_z_pt = new TCanvas ("c_z_pt", "", 800, 800);
  TPad* p_up_z_pt = new TPad ("p_up_z_pt", "", 0.0, 0.4, 1.0, 1.0);
  TPad* p_down_z_pt = new TPad ("p_down_z_pt", "", 0.0, 0.0, 1.0, 0.4);
  p_up_z_pt->SetBottomMargin (0);
  p_down_z_pt->SetTopMargin (0);
  p_down_z_pt->SetBottomMargin (0.25);
  p_up_z_pt->Draw ();
  p_down_z_pt->Draw ();

  h_electronlhloose_z_pt->SetLineColor (kRed+1);
  h_electronlhloose_z_pt->SetMarkerColor (kRed+1);
  h_electronlhloose_z_pt->SetMarkerStyle (kFullCircle);
  h_muonmedium_z_pt->SetLineColor (kAzure+2);
  h_muonmedium_z_pt->SetMarkerColor (kAzure+2);
  h_muonmedium_z_pt->SetMarkerStyle (kFullCircle);
  h_electronlhmed_z_pt->SetLineColor (kRed+1);
  h_electronlhmed_z_pt->SetMarkerColor (kRed+1);
  h_electronlhmed_z_pt->SetMarkerStyle (kOpenCircle);
  h_muontight_z_pt->SetLineColor (kAzure+2);
  h_muontight_z_pt->SetMarkerColor (kAzure+2);
  h_muontight_z_pt->SetMarkerStyle (kOpenCircle);

  p_up_z_pt->cd ();
  h_electronlhloose_z_pt->Draw ("e1");
  h_muonmedium_z_pt->Draw ("e1 same");
  h_electronlhmed_z_pt->DrawCopy ("e1 same");
  h_muontight_z_pt->DrawCopy ("e1 same");
  myMarkerTextNoLine (0.65, 0.88, kAzure+2, kFullCircle, "Medium Muons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.80, kRed+1, kFullCircle, "LHLoose Electrons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.72, kAzure+2, kOpenCircle, "Tight Muons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.64, kRed+1, kOpenCircle, "LHMedium Electrons", 1.4, 0.05);

  p_down_z_pt->cd ();
  h_electronlhmed_z_pt->Divide (h_electronlhloose_z_pt);
  h_muontight_z_pt->Divide (h_muonmedium_z_pt);
  h_electronlhmed_z_pt->Draw ("e1");
  h_muontight_z_pt->Draw ("e1 same");


  TCanvas* c_z_y = new TCanvas ("c_z_y", "", 800, 800);
  TPad* p_up_z_y = new TPad ("p_up_z_y", "", 0.0, 0.4, 1.0, 1.0);
  TPad* p_down_z_y = new TPad ("p_down_z_y", "", 0.0, 0.0, 1.0, 0.4);
  p_up_z_y->SetBottomMargin (0);
  p_down_z_y->SetTopMargin (0);
  p_down_z_y->SetBottomMargin (0.25);
  p_up_z_y->Draw ();
  p_down_z_y->Draw ();
  h_electronlhloose_z_y->SetLineColor (kRed+1);
  h_electronlhloose_z_y->SetMarkerColor (kRed+1);
  h_electronlhloose_z_y->SetMarkerStyle (kFullCircle);
  h_muonmedium_z_y->SetLineColor (kAzure+2);
  h_muonmedium_z_y->SetMarkerColor (kAzure+2);
  h_muonmedium_z_y->SetMarkerStyle (kFullCircle);
  h_electronlhmed_z_y->SetLineColor (kRed+1);
  h_electronlhmed_z_y->SetMarkerColor (kRed+1);
  h_electronlhmed_z_y->SetMarkerStyle (kOpenCircle);
  h_muontight_z_y->SetLineColor (kAzure+2);
  h_muontight_z_y->SetMarkerColor (kAzure+2);
  h_muontight_z_y->SetMarkerStyle (kOpenCircle);

  p_up_z_y->cd ();
  h_electronlhloose_z_y->Draw ("e1");
  h_muonmedium_z_y->Draw ("e1 same");
  h_electronlhmed_z_y->DrawCopy ("e1 same");
  h_muontight_z_y->DrawCopy ("e1 same");
  myMarkerTextNoLine (0.65, 0.88, kAzure+2, kFullCircle, "Medium Muons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.80, kRed+1, kFullCircle, "LHLoose Electrons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.72, kAzure+2, kOpenCircle, "Tight Muons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.64, kRed+1, kOpenCircle, "LHMedium Electrons", 1.4, 0.05);

  p_down_z_y->cd ();
  h_electronlhmed_z_y->Divide (h_electronlhloose_z_y);
  h_muontight_z_y->Divide (h_muonmedium_z_y);
  h_electronlhmed_z_y->Draw ("e1");
  h_muontight_z_y->Draw ("e1 same");


  TCanvas* c_z_phi = new TCanvas ("c_z_phi", "", 800, 800);
  TPad* p_up_z_phi = new TPad ("p_up_z_phi", "", 0.0, 0.4, 1.0, 1.0);
  TPad* p_down_z_phi = new TPad ("p_down_z_phi", "", 0.0, 0.0, 1.0, 0.4);
  p_up_z_phi->SetBottomMargin (0);
  p_down_z_phi->SetTopMargin (0);
  p_down_z_phi->SetBottomMargin (0.25);
  p_up_z_phi->Draw ();
  p_down_z_phi->Draw ();
  h_electronlhloose_z_phi->SetLineColor (kRed+1);
  h_electronlhloose_z_phi->SetMarkerColor (kRed+1);
  h_electronlhloose_z_phi->SetMarkerStyle (kFullCircle);
  h_muonmedium_z_phi->SetLineColor (kAzure+2);
  h_muonmedium_z_phi->SetMarkerColor (kAzure+2);
  h_muonmedium_z_phi->SetMarkerStyle (kFullCircle);
  h_electronlhmed_z_phi->SetLineColor (kRed+1);
  h_electronlhmed_z_phi->SetMarkerColor (kRed+1);
  h_electronlhmed_z_phi->SetMarkerStyle (kOpenCircle);
  h_muontight_z_phi->SetLineColor (kAzure+2);
  h_muontight_z_phi->SetMarkerColor (kAzure+2);
  h_muontight_z_phi->SetMarkerStyle (kOpenCircle);

  p_up_z_phi->cd ();
  h_electronlhloose_z_phi->Draw ("e1");
  h_muonmedium_z_phi->Draw ("e1 same");
  h_electronlhmed_z_phi->DrawCopy ("e1 same");
  h_muontight_z_phi->DrawCopy ("e1 same");
  myMarkerTextNoLine (0.65, 0.88, kAzure+2, kFullCircle, "Medium Muons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.80, kRed+1, kFullCircle, "LHLoose Electrons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.72, kAzure+2, kOpenCircle, "Tight Muons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.64, kRed+1, kOpenCircle, "LHMedium Electrons", 1.4, 0.05);

  p_down_z_phi->cd ();
  h_electronlhmed_z_phi->Divide (h_electronlhloose_z_phi);
  h_muontight_z_phi->Divide (h_muonmedium_z_phi);
  h_electronlhmed_z_phi->Draw ("e1");
  h_muontight_z_phi->Draw ("e1 same");


  TCanvas* c_z_m = new TCanvas ("c_z_m", "", 800, 800);
  TPad* p_up_z_m = new TPad ("p_up_z_m", "", 0.0, 0.4, 1.0, 1.0);
  TPad* p_down_z_m = new TPad ("p_down_z_m", "", 0.0, 0.0, 1.0, 0.4);
  p_up_z_m->SetBottomMargin (0);
  p_down_z_m->SetTopMargin (0);
  p_down_z_m->SetBottomMargin (0.25);
  p_up_z_m->Draw ();
  p_down_z_m->Draw ();
  h_electronlhloose_z_m->SetLineColor (kRed+1);
  h_electronlhloose_z_m->SetMarkerColor (kRed+1);
  h_electronlhloose_z_m->SetMarkerStyle (kFullCircle);
  h_muonmedium_z_m->SetLineColor (kAzure+2);
  h_muonmedium_z_m->SetMarkerColor (kAzure+2);
  h_muonmedium_z_m->SetMarkerStyle (kFullCircle);
  h_electronlhmed_z_m->SetLineColor (kRed+1);
  h_electronlhmed_z_m->SetMarkerColor (kRed+1);
  h_electronlhmed_z_m->SetMarkerStyle (kOpenCircle);
  h_muontight_z_m->SetLineColor (kAzure+2);
  h_muontight_z_m->SetMarkerColor (kAzure+2);
  h_muontight_z_m->SetMarkerStyle (kOpenCircle);

  p_up_z_m->cd ();
  h_electronlhloose_z_m->Draw ("e1");
  h_muonmedium_z_m->Draw ("e1 same");
  h_electronlhmed_z_m->DrawCopy ("e1 same");
  h_muontight_z_m->DrawCopy ("e1 same");
  myMarkerTextNoLine (0.65, 0.88, kAzure+2, kFullCircle, "Medium Muons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.80, kRed+1, kFullCircle, "LHLoose Electrons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.72, kAzure+2, kOpenCircle, "Tight Muons", 1.4, 0.05);
  myMarkerTextNoLine (0.65, 0.64, kRed+1, kOpenCircle, "LHMedium Electrons", 1.4, 0.05);

  p_down_z_m->cd ();
  h_electronlhmed_z_m->Divide (h_electronlhloose_z_m);
  h_muontight_z_m->Divide (h_muonmedium_z_m);
  h_electronlhmed_z_m->Draw ("e1");
  h_muontight_z_m->Draw ("e1 same");
}

#endif

