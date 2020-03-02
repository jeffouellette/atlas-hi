#ifndef __PlotLeptonEffs_C__
#define __PlotLeptonEffs_C__

#include "Params.h"
//#include "ZTrackUtilities.h"

#include <ArrayTemplates.h>

#include <TEfficiency.h>

#include <iostream>

typedef TGraphAsymmErrors TGAE;

TEfficiency* t_muonTrigEff_pt_pp;
TEfficiency* t_muonTrigEff_pt_PbPb18;
TEfficiency* t_muonTrigEff_pt_PbPb15;
TEfficiency* t_muonTrigEff_fcal_PbPb18;
TEfficiency* t_muonTrigEff_fcal_PbPb15;
TEfficiency* t_muonTrigEff_eta_phi_pp;
TEfficiency* t_muonTrigEff_eta_phi_PbPb18;
TEfficiency* t_muonTrigEff_eta_phi_PbPb15;
TEfficiency* t_muonTrigEff_eta_pp;
TEfficiency* t_muonTrigEff_eta_PbPb18;
TEfficiency* t_muonTrigEff_eta_PbPb15;

//TEfficiency* t_muonMed_IDEff_pt_pp[3];
//TEfficiency* t_muonMed_IDEff_pt_PbPb18[3];
//TEfficiency* t_muonMed_IDEff_pt_PbPb15[3];
//TEfficiency* t_muonMed_IDEff_eta_pp[3];
//TEfficiency* t_muonMed_IDEff_eta_PbPb18[3];
//TEfficiency* t_muonMed_IDEff_eta_PbPb15[3];
//TEfficiency* t_muonMed_IDEff_fcal_PbPb18[3];
//TEfficiency* t_muonMed_IDEff_fcal_PbPb15[3];
//
//TEfficiency* t_muonID_MSEff_pt_pp[3];
//TEfficiency* t_muonID_MSEff_pt_PbPb18[3];
//TEfficiency* t_muonID_MSEff_pt_PbPb15[3];
//TEfficiency* t_muonID_MSEff_eta_pp[3];
//TEfficiency* t_muonID_MSEff_eta_PbPb18[3];
//TEfficiency* t_muonID_MSEff_eta_PbPb15[3];
//TEfficiency* t_muonID_MSEff_fcal_PbPb18[3];
//TEfficiency* t_muonID_MSEff_fcal_PbPb15[3];

TEfficiency* t_muonIDEff_pt_pp;
TEfficiency* t_muonIDEff_pt_PbPb18;
TEfficiency* t_muonIDEff_pt_PbPb15;
TEfficiency* t_muonIDEff_fcal_PbPb18;
TEfficiency* t_muonIDEff_fcal_PbPb15;
TEfficiency* t_muonIDEff_eta_phi_pp;
TEfficiency* t_muonIDEff_eta_phi_PbPb18;
TEfficiency* t_muonIDEff_eta_phi_PbPb15;
TEfficiency* t_muonIDEff_eta_pp;
TEfficiency* t_muonIDEff_eta_PbPb18;
TEfficiency* t_muonIDEff_eta_PbPb15;

TEfficiency* t_electronTrigEff_pt_pp;
TEfficiency* t_electronTrigEff_pt_PbPb18;
TEfficiency* t_electronTrigEff_pt_PbPb15;
TEfficiency* t_electronTrigEff_eta_pp;
TEfficiency* t_electronTrigEff_eta_PbPb18;
TEfficiency* t_electronTrigEff_eta_PbPb15;
TEfficiency* t_electronTrigEff_pt_eta_pp;
TEfficiency* t_electronTrigEff_pt_eta_PbPb18;
TEfficiency* t_electronTrigEff_pt_eta_PbPb15;
TEfficiency* t_electronTrigEff_fcal_PbPb18;
TEfficiency* t_electronTrigEff_fcal_PbPb15;

TEfficiency* t_electronIDEff_pt_pp;
TEfficiency* t_electronIDEff_pt_PbPb18;
TEfficiency* t_electronIDEff_pt_PbPb15;
TEfficiency* t_electronIDEff_eta_pp;
TEfficiency* t_electronIDEff_eta_PbPb18;
TEfficiency* t_electronIDEff_eta_PbPb15;
TEfficiency* t_electronIDEff_pt_eta_pp;
TEfficiency* t_electronIDEff_pt_eta_PbPb18;
TEfficiency* t_electronIDEff_pt_eta_PbPb15;
TEfficiency* t_electronIDEff_fcal_PbPb18;
TEfficiency* t_electronIDEff_fcal_PbPb15;


bool isPbPb = true;

void PlotLeptonEffs () {

  SetupDirectories ("TagAndProbe/", "ZTrackAnalysis/");
  TFile* inFile = new TFile (Form ("%s/Nominal/outFile.root", rootPath.Data ()), "read");

  TH1D* h_muonTrigEffNum_pt[3];
  TH1D* h_muonTrigEffDen_pt[3];
  TH1D* h_muonTrigEffNum_eta[3];
  TH1D* h_muonTrigEffDen_eta[3];
  TH2D* h2_muonTrigEffNum_eta_phi[3];
  TH2D* h2_muonTrigEffDen_eta_phi[3];
  TH1D* h_muonTrigEffNum_fcal[2];
  TH1D* h_muonTrigEffDen_fcal[2];
  TH1D* h_muonIDEffNum_pt[3];
  TH1D* h_muonIDEffDen_pt[3];
  TH1D* h_muonIDEffNum_eta[3];
  TH1D* h_muonIDEffDen_eta[3];
  TH2D* h2_muonIDEffNum_eta_phi[3];
  TH2D* h2_muonIDEffDen_eta_phi[3];
  TH1D* h_muonIDEffNum_fcal[2];
  TH1D* h_muonIDEffDen_fcal[2];

  h_muonTrigEffNum_pt[0] = (TH1D*) inFile->Get ("h_muonTrigEffNum_pt_pp");
  h_muonTrigEffDen_pt[0] = (TH1D*) inFile->Get ("h_muonTrigEffDen_pt_pp");
  h_muonTrigEffNum_pt[1] = (TH1D*) inFile->Get ("h_muonTrigEffNum_pt_PbPb18");
  h_muonTrigEffDen_pt[1] = (TH1D*) inFile->Get ("h_muonTrigEffDen_pt_PbPb18");
  //h_muonTrigEffNum_pt[2] = (TH1D*) inFile->Get ("h_muonTrigEffNum_pt_PbPb15");
  //h_muonTrigEffDen_pt[2] = (TH1D*) inFile->Get ("h_muonTrigEffDen_pt_PbPb15");
  h_muonIDEffNum_pt[0]   = (TH1D*) inFile->Get ("h_muonIDEffNum_pt_pp");
  h_muonIDEffDen_pt[0]   = (TH1D*) inFile->Get ("h_muonIDEffDen_pt_pp");
  h_muonIDEffNum_pt[1]   = (TH1D*) inFile->Get ("h_muonIDEffNum_pt_PbPb18");
  h_muonIDEffDen_pt[1]   = (TH1D*) inFile->Get ("h_muonIDEffDen_pt_PbPb18");
  //h_muonIDEffNum_pt[2]   = (TH1D*) inFile->Get ("h_muonIDEffNum_pt_PbPb15");
  //h_muonIDEffDen_pt[2]   = (TH1D*) inFile->Get ("h_muonIDEffDen_pt_PbPb15");

  t_muonTrigEff_pt_pp = new TEfficiency (*(h_muonTrigEffNum_pt[0]), *(h_muonTrigEffDen_pt[0]));
  t_muonTrigEff_pt_PbPb18 = new TEfficiency (*(h_muonTrigEffNum_pt[1]), *(h_muonTrigEffDen_pt[1]));
  //t_muonTrigEff_pt_PbPb15 = new TEfficiency (*(h_muonTrigEffNum_pt[2]), *(h_muonTrigEffDen_pt[2]));
  t_muonIDEff_pt_pp = new TEfficiency (*(h_muonIDEffNum_pt[0]), *(h_muonIDEffDen_pt[0]));
  t_muonIDEff_pt_PbPb18 = new TEfficiency (*(h_muonIDEffNum_pt[1]), *(h_muonIDEffDen_pt[1]));
  //t_muonIDEff_pt_PbPb15 = new TEfficiency (*(h_muonIDEffNum_pt[2]), *(h_muonIDEffDen_pt[2]));

  h2_muonTrigEffNum_eta_phi[0] = (TH2D*) inFile->Get ("h2_muonTrigEffNum_eta_phi_pp");
  h2_muonTrigEffDen_eta_phi[0] = (TH2D*) inFile->Get ("h2_muonTrigEffDen_eta_phi_pp");
  h2_muonTrigEffNum_eta_phi[1] = (TH2D*) inFile->Get ("h2_muonTrigEffNum_eta_phi_PbPb18");
  h2_muonTrigEffDen_eta_phi[1] = (TH2D*) inFile->Get ("h2_muonTrigEffDen_eta_phi_PbPb18");
  //h2_muonTrigEffNum_eta_phi[2] = (TH2D*) inFile->Get ("h2_muonTrigEffNum_eta_phi_PbPb15");
  //h2_muonTrigEffDen_eta_phi[2] = (TH2D*) inFile->Get ("h2_muonTrigEffDen_eta_phi_PbPb15");
  h2_muonIDEffNum_eta_phi[0]   = (TH2D*) inFile->Get ("h2_muonIDEffNum_eta_phi_pp");
  h2_muonIDEffDen_eta_phi[0]   = (TH2D*) inFile->Get ("h2_muonIDEffDen_eta_phi_pp");
  h2_muonIDEffNum_eta_phi[1]   = (TH2D*) inFile->Get ("h2_muonIDEffNum_eta_phi_PbPb18");
  h2_muonIDEffDen_eta_phi[1]   = (TH2D*) inFile->Get ("h2_muonIDEffDen_eta_phi_PbPb18");
  //h2_muonIDEffNum_eta_phi[2]   = (TH2D*) inFile->Get ("h2_muonIDEffNum_eta_phi_PbPb15");
  //h2_muonIDEffDen_eta_phi[2]   = (TH2D*) inFile->Get ("h2_muonIDEffDen_eta_phi_PbPb15");

  h_muonTrigEffNum_eta[0] = (TH1D*) h2_muonTrigEffNum_eta_phi[0]->ProjectionX ("h_muonTrigEffNum_eta_pp");
  h_muonTrigEffDen_eta[0] = (TH1D*) h2_muonTrigEffDen_eta_phi[0]->ProjectionX ("h_muonTrigEffDen_eta_pp");
  h_muonTrigEffNum_eta[1] = (TH1D*) h2_muonTrigEffNum_eta_phi[1]->ProjectionX ("h_muonTrigEffNum_eta_PbPb18");
  h_muonTrigEffDen_eta[1] = (TH1D*) h2_muonTrigEffDen_eta_phi[1]->ProjectionX ("h_muonTrigEffDen_eta_PbPb18");
  //h_muonTrigEffNum_eta[2] = (TH1D*) h2_muonTrigEffNum_eta_phi[2]->ProjectionX ("h_muonTrigEffNum_eta_PbPb15");
  //h_muonTrigEffDen_eta[2] = (TH1D*) h2_muonTrigEffDen_eta_phi[2]->ProjectionX ("h_muonTrigEffDen_eta_PbPb15");
  h_muonIDEffNum_eta[0]   = (TH1D*) h2_muonIDEffNum_eta_phi[0]->ProjectionX ("h_muonIDEffNum_eta_pp");
  h_muonIDEffDen_eta[0]   = (TH1D*) h2_muonIDEffDen_eta_phi[0]->ProjectionX ("h_muonIDEffDen_eta_pp");
  h_muonIDEffNum_eta[1]   = (TH1D*) h2_muonIDEffNum_eta_phi[1]->ProjectionX ("h_muonIDEffNum_eta_PbPb18");
  h_muonIDEffDen_eta[1]   = (TH1D*) h2_muonIDEffDen_eta_phi[1]->ProjectionX ("h_muonIDEffDen_eta_PbPb18");
  //h_muonIDEffNum_eta[2]   = (TH1D*) h2_muonIDEffNum_eta_phi[2]->ProjectionX ("h_muonIDEffNum_eta_PbPb15");
  //h_muonIDEffDen_eta[2]   = (TH1D*) h2_muonIDEffDen_eta_phi[2]->ProjectionX ("h_muonIDEffDen_eta_PbPb15");

  t_muonTrigEff_eta_phi_pp = new TEfficiency (*(h2_muonTrigEffNum_eta_phi[0]), *(h2_muonTrigEffDen_eta_phi[0]));
  t_muonTrigEff_eta_phi_PbPb18 = new TEfficiency (*(h2_muonTrigEffNum_eta_phi[1]), *(h2_muonTrigEffDen_eta_phi[1]));
  //t_muonTrigEff_eta_phi_PbPb15 = new TEfficiency (*(h2_muonTrigEffNum_eta_phi[2]), *(h2_muonTrigEffDen_eta_phi[2]));
  t_muonTrigEff_eta_pp = new TEfficiency (*(h_muonTrigEffNum_eta[0]), *(h_muonTrigEffDen_eta[0]));
  t_muonTrigEff_eta_PbPb18 = new TEfficiency (*(h_muonTrigEffNum_eta[1]), *(h_muonTrigEffDen_eta[1]));
  //t_muonTrigEff_eta_PbPb15 = new TEfficiency (*(h_muonTrigEffNum_eta[2]), *(h_muonTrigEffDen_eta[2]));
  t_muonIDEff_eta_phi_pp = new TEfficiency (*(h2_muonIDEffNum_eta_phi[0]), *(h2_muonIDEffDen_eta_phi[0]));
  t_muonIDEff_eta_phi_PbPb18 = new TEfficiency (*(h2_muonIDEffNum_eta_phi[1]), *(h2_muonIDEffDen_eta_phi[1]));
  //t_muonIDEff_eta_phi_PbPb15 = new TEfficiency (*(h2_muonIDEffNum_eta_phi[2]), *(h2_muonIDEffDen_eta_phi[2]));
  t_muonIDEff_eta_pp = new TEfficiency (*(h_muonIDEffNum_eta[0]), *(h_muonIDEffDen_eta[0]));
  t_muonIDEff_eta_PbPb18 = new TEfficiency (*(h_muonIDEffNum_eta[1]), *(h_muonIDEffDen_eta[1]));
  //t_muonIDEff_eta_PbPb15 = new TEfficiency (*(h_muonIDEffNum_eta[2]), *(h_muonIDEffDen_eta[2]));

  h_muonTrigEffNum_fcal[0] = (TH1D*) inFile->Get ("h_muonTrigEffNum_fcal_PbPb18");
  h_muonTrigEffDen_fcal[0] = (TH1D*) inFile->Get ("h_muonTrigEffDen_fcal_PbPb18");
  //h_muonTrigEffNum_fcal[1] = (TH1D*) inFile->Get ("h_muonTrigEffNum_fcal_PbPb15");
  //h_muonTrigEffDen_fcal[1] = (TH1D*) inFile->Get ("h_muonTrigEffDen_fcal_PbPb15");
  t_muonTrigEff_fcal_PbPb18 = new TEfficiency (*(h_muonTrigEffNum_fcal[0]), *(h_muonTrigEffDen_fcal[0]));
  //t_muonTrigEff_fcal_PbPb15 = new TEfficiency (*(h_muonTrigEffNum_fcal[1]), *(h_muonTrigEffDen_fcal[1]));
  h_muonIDEffNum_fcal[0] = (TH1D*) inFile->Get ("h_muonIDEffNum_fcal_PbPb18");
  h_muonIDEffDen_fcal[0] = (TH1D*) inFile->Get ("h_muonIDEffDen_fcal_PbPb18");
  //h_muonIDEffNum_fcal[1] = (TH1D*) inFile->Get ("h_muonIDEffNum_fcal_PbPb15");
  //h_muonIDEffDen_fcal[1] = (TH1D*) inFile->Get ("h_muonIDEffDen_fcal_PbPb15");
  t_muonIDEff_fcal_PbPb18 = new TEfficiency (*(h_muonIDEffNum_fcal[0]), *(h_muonIDEffDen_fcal[0]));
  //t_muonIDEff_fcal_PbPb15 = new TEfficiency (*(h_muonIDEffNum_fcal[1]), *(h_muonIDEffDen_fcal[1]));

  //TH1D* h_muonMed_IDEffNum_pt[3][3];
  //TH1D* h_muonMed_IDEffDen_pt[3][3];
  //TH1D* h_muonMed_IDEffNum_eta[3][3];
  //TH1D* h_muonMed_IDEffDen_eta[3][3];
  //TH1D* h_muonMed_IDEffNum_fcal[2][3];
  //TH1D* h_muonMed_IDEffDen_fcal[2][3];

  //TH1D* h_muonID_MSEffNum_pt[3][3];
  //TH1D* h_muonID_MSEffDen_pt[3][3];
  //TH1D* h_muonID_MSEffNum_eta[3][3];
  //TH1D* h_muonID_MSEffDen_eta[3][3];
  //TH1D* h_muonID_MSEffNum_fcal[2][3];
  //TH1D* h_muonID_MSEffDen_fcal[2][3];

  //for (int iSign : {0, 1}) {
  //  h_muonMed_IDEffNum_pt[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_pt_pp_sign%i", iSign));
  //  h_muonMed_IDEffDen_pt[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_pt_pp_sign%i", iSign));
  //  h_muonMed_IDEffNum_pt[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_pt_PbPb18_sign%i", iSign));
  //  h_muonMed_IDEffDen_pt[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_pt_PbPb18_sign%i", iSign));
  //  h_muonMed_IDEffNum_pt[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_pt_PbPb15_sign%i", iSign));
  //  h_muonMed_IDEffDen_pt[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_pt_PbPb15_sign%i", iSign));

  //  t_muonMed_IDEff_pt_pp[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_pt[0][iSign]), *(h_muonMed_IDEffDen_pt[0][iSign]));
  //  t_muonMed_IDEff_pt_PbPb18[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_pt[1][iSign]), *(h_muonMed_IDEffDen_pt[1][iSign]));
  //  t_muonMed_IDEff_pt_PbPb15[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_pt[2][iSign]), *(h_muonMed_IDEffDen_pt[2][iSign]));

  //  h_muonMed_IDEffNum_eta[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_eta_pp_sign%i", iSign));
  //  h_muonMed_IDEffDen_eta[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_eta_pp_sign%i", iSign));
  //  h_muonMed_IDEffNum_eta[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_eta_PbPb18_sign%i", iSign));
  //  h_muonMed_IDEffDen_eta[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_eta_PbPb18_sign%i", iSign));
  //  h_muonMed_IDEffNum_eta[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_eta_PbPb15_sign%i", iSign));
  //  h_muonMed_IDEffDen_eta[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_eta_PbPb15_sign%i", iSign));

  //  t_muonMed_IDEff_eta_pp[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_eta[0][iSign]), *(h_muonMed_IDEffDen_eta[0][iSign]));
  //  t_muonMed_IDEff_eta_PbPb18[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_eta[1][iSign]), *(h_muonMed_IDEffDen_eta[1][iSign]));
  //  t_muonMed_IDEff_eta_PbPb15[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_eta[2][iSign]), *(h_muonMed_IDEffDen_eta[2][iSign]));

  //  h_muonMed_IDEffNum_fcal[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_fcal_PbPb18_sign%i", iSign));
  //  h_muonMed_IDEffDen_fcal[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_fcal_PbPb18_sign%i", iSign));
  //  h_muonMed_IDEffNum_fcal[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_fcal_PbPb15_sign%i", iSign));
  //  h_muonMed_IDEffDen_fcal[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_fcal_PbPb15_sign%i", iSign));
  //  t_muonMed_IDEff_fcal_PbPb18[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_fcal[0][iSign]), *(h_muonMed_IDEffDen_fcal[0][iSign]));
  //  t_muonMed_IDEff_fcal_PbPb15[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_fcal[1][iSign]), *(h_muonMed_IDEffDen_fcal[1][iSign]));


  //  h_muonID_MSEffNum_pt[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_pt_pp_sign%i", iSign));
  //  h_muonID_MSEffDen_pt[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_pt_pp_sign%i", iSign));
  //  h_muonID_MSEffNum_pt[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_pt_PbPb18_sign%i", iSign));
  //  h_muonID_MSEffDen_pt[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_pt_PbPb18_sign%i", iSign));
  //  h_muonID_MSEffNum_pt[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_pt_PbPb15_sign%i", iSign));
  //  h_muonID_MSEffDen_pt[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_pt_PbPb15_sign%i", iSign));

  //  t_muonID_MSEff_pt_pp[iSign] = new TEfficiency (*(h_muonID_MSEffNum_pt[0][iSign]), *(h_muonID_MSEffDen_pt[0][iSign]));
  //  t_muonID_MSEff_pt_PbPb18[iSign] = new TEfficiency (*(h_muonID_MSEffNum_pt[1][iSign]), *(h_muonID_MSEffDen_pt[1][iSign]));
  //  t_muonID_MSEff_pt_PbPb15[iSign] = new TEfficiency (*(h_muonID_MSEffNum_pt[2][iSign]), *(h_muonID_MSEffDen_pt[2][iSign]));

  //  h_muonID_MSEffNum_eta[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_eta_pp_sign%i", iSign));
  //  h_muonID_MSEffDen_eta[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_eta_pp_sign%i", iSign));
  //  h_muonID_MSEffNum_eta[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_eta_PbPb18_sign%i", iSign));
  //  h_muonID_MSEffDen_eta[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_eta_PbPb18_sign%i", iSign));
  //  h_muonID_MSEffNum_eta[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_eta_PbPb15_sign%i", iSign));
  //  h_muonID_MSEffDen_eta[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_eta_PbPb15_sign%i", iSign));

  //  t_muonID_MSEff_eta_pp[iSign] = new TEfficiency (*(h_muonID_MSEffNum_eta[0][iSign]), *(h_muonID_MSEffDen_eta[0][iSign]));
  //  t_muonID_MSEff_eta_PbPb18[iSign] = new TEfficiency (*(h_muonID_MSEffNum_eta[1][iSign]), *(h_muonID_MSEffDen_eta[1][iSign]));
  //  t_muonID_MSEff_eta_PbPb15[iSign] = new TEfficiency (*(h_muonID_MSEffNum_eta[2][iSign]), *(h_muonID_MSEffDen_eta[2][iSign]));

  //  h_muonID_MSEffNum_fcal[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_fcal_PbPb18_sign%i", iSign));
  //  h_muonID_MSEffDen_fcal[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_fcal_PbPb18_sign%i", iSign));
  //  h_muonID_MSEffNum_fcal[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_fcal_PbPb15_sign%i", iSign));
  //  h_muonID_MSEffDen_fcal[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_fcal_PbPb15_sign%i", iSign));
  //  t_muonID_MSEff_fcal_PbPb18[iSign] = new TEfficiency (*(h_muonID_MSEffNum_fcal[0][iSign]), *(h_muonID_MSEffDen_fcal[0][iSign]));
  //  t_muonID_MSEff_fcal_PbPb15[iSign] = new TEfficiency (*(h_muonID_MSEffNum_fcal[1][iSign]), *(h_muonID_MSEffDen_fcal[1][iSign]));
  //}

  //{
  //  h_muonMed_IDEffNum_pt[0][2] = (TH1D*) h_muonMed_IDEffNum_pt[0][1]->Clone ("h_muonMed_IDEffNum_pt_pp_sign2");
  //  h_muonMed_IDEffDen_pt[0][2] = (TH1D*) h_muonMed_IDEffDen_pt[0][1]->Clone ("h_muonMed_IDEffDen_pt_pp_sign2");
  //  h_muonMed_IDEffNum_pt[1][2] = (TH1D*) h_muonMed_IDEffNum_pt[1][1]->Clone ("h_muonMed_IDEffNum_pt_PbPb18_sign2");
  //  h_muonMed_IDEffDen_pt[1][2] = (TH1D*) h_muonMed_IDEffDen_pt[1][1]->Clone ("h_muonMed_IDEffDen_pt_PbPb18_sign2");
  //  h_muonMed_IDEffNum_pt[2][2] = (TH1D*) h_muonMed_IDEffNum_pt[2][1]->Clone ("h_muonMed_IDEffNum_pt_PbPb15_sign2");
  //  h_muonMed_IDEffDen_pt[2][2] = (TH1D*) h_muonMed_IDEffDen_pt[2][1]->Clone ("h_muonMed_IDEffDen_pt_PbPb15_sign2");
  //  h_muonMed_IDEffNum_pt[0][2]->Add (h_muonMed_IDEffNum_pt[0][0], -1);
  //  h_muonMed_IDEffDen_pt[0][2]->Add (h_muonMed_IDEffDen_pt[0][0], -1);
  //  h_muonMed_IDEffNum_pt[1][2]->Add (h_muonMed_IDEffNum_pt[1][0], -1);
  //  h_muonMed_IDEffDen_pt[1][2]->Add (h_muonMed_IDEffDen_pt[1][0], -1);
  //  h_muonMed_IDEffNum_pt[2][2]->Add (h_muonMed_IDEffNum_pt[2][0], -1);
  //  h_muonMed_IDEffDen_pt[2][2]->Add (h_muonMed_IDEffDen_pt[2][0], -1);

  //  t_muonMed_IDEff_pt_pp[2] = new TEfficiency (*(h_muonMed_IDEffNum_pt[0][2]), *(h_muonMed_IDEffDen_pt[0][2]));
  //  t_muonMed_IDEff_pt_PbPb18[2] = new TEfficiency (*(h_muonMed_IDEffNum_pt[1][2]), *(h_muonMed_IDEffDen_pt[1][2]));
  //  t_muonMed_IDEff_pt_PbPb15[2] = new TEfficiency (*(h_muonMed_IDEffNum_pt[2][2]), *(h_muonMed_IDEffDen_pt[2][2]));

  //  h_muonMed_IDEffNum_eta[0][2] = (TH1D*) h_muonMed_IDEffNum_eta[0][1]->Clone ("h_muonMed_IDEffNum_eta_pp_sign2");
  //  h_muonMed_IDEffDen_eta[0][2] = (TH1D*) h_muonMed_IDEffDen_eta[0][1]->Clone ("h_muonMed_IDEffDen_eta_pp_sign2");
  //  h_muonMed_IDEffNum_eta[1][2] = (TH1D*) h_muonMed_IDEffNum_eta[1][1]->Clone ("h_muonMed_IDEffNum_eta_PbPb18_sign2");
  //  h_muonMed_IDEffDen_eta[1][2] = (TH1D*) h_muonMed_IDEffDen_eta[1][1]->Clone ("h_muonMed_IDEffDen_eta_PbPb18_sign2");
  //  h_muonMed_IDEffNum_eta[2][2] = (TH1D*) h_muonMed_IDEffNum_eta[2][1]->Clone ("h_muonMed_IDEffNum_eta_PbPb15_sign2");
  //  h_muonMed_IDEffDen_eta[2][2] = (TH1D*) h_muonMed_IDEffDen_eta[2][1]->Clone ("h_muonMed_IDEffDen_eta_PbPb15_sign2");
  //  h_muonMed_IDEffNum_eta[0][2]->Add (h_muonMed_IDEffNum_eta[0][0], -1);
  //  h_muonMed_IDEffDen_eta[0][2]->Add (h_muonMed_IDEffDen_eta[0][0], -1);
  //  h_muonMed_IDEffNum_eta[1][2]->Add (h_muonMed_IDEffNum_eta[1][0], -1);
  //  h_muonMed_IDEffDen_eta[1][2]->Add (h_muonMed_IDEffDen_eta[1][0], -1);
  //  h_muonMed_IDEffNum_eta[2][2]->Add (h_muonMed_IDEffNum_eta[2][0], -1);
  //  h_muonMed_IDEffDen_eta[2][2]->Add (h_muonMed_IDEffDen_eta[2][0], -1);

  //  t_muonMed_IDEff_eta_pp[2] = new TEfficiency (*(h_muonMed_IDEffNum_eta[0][2]), *(h_muonMed_IDEffDen_eta[0][2]));
  //  t_muonMed_IDEff_eta_PbPb18[2] = new TEfficiency (*(h_muonMed_IDEffNum_eta[1][2]), *(h_muonMed_IDEffDen_eta[1][2]));
  //  t_muonMed_IDEff_eta_PbPb15[2] = new TEfficiency (*(h_muonMed_IDEffNum_eta[2][2]), *(h_muonMed_IDEffDen_eta[2][2]));

  //  h_muonMed_IDEffNum_fcal[0][2] = (TH1D*) h_muonMed_IDEffNum_fcal[0][1]->Clone ("h_muonMed_IDEffNum_fcal_PbPb18_sign2");
  //  h_muonMed_IDEffDen_fcal[0][2] = (TH1D*) h_muonMed_IDEffDen_fcal[0][1]->Clone ("h_muonMed_IDEffDen_fcal_PbPb18_sign2");
  //  h_muonMed_IDEffNum_fcal[1][2] = (TH1D*) h_muonMed_IDEffNum_fcal[1][1]->Clone ("h_muonMed_IDEffNum_fcal_PbPb15_sign2");
  //  h_muonMed_IDEffDen_fcal[1][2] = (TH1D*) h_muonMed_IDEffDen_fcal[1][1]->Clone ("h_muonMed_IDEffDen_fcal_PbPb15_sign2");
  //  h_muonMed_IDEffNum_fcal[0][2]->Add (h_muonMed_IDEffNum_fcal[0][0], -1);
  //  h_muonMed_IDEffDen_fcal[0][2]->Add (h_muonMed_IDEffDen_fcal[0][0], -1);
  //  h_muonMed_IDEffNum_fcal[1][2]->Add (h_muonMed_IDEffNum_fcal[1][0], -1);
  //  h_muonMed_IDEffDen_fcal[1][2]->Add (h_muonMed_IDEffDen_fcal[1][0], -1);

  //  t_muonMed_IDEff_fcal_PbPb18[2] = new TEfficiency (*(h_muonMed_IDEffNum_fcal[0][2]), *(h_muonMed_IDEffDen_fcal[0][2]));
  //  t_muonMed_IDEff_fcal_PbPb15[2] = new TEfficiency (*(h_muonMed_IDEffNum_fcal[1][2]), *(h_muonMed_IDEffDen_fcal[1][2]));


  //  h_muonID_MSEffNum_pt[0][2] = (TH1D*) h_muonID_MSEffNum_pt[0][1]->Clone ("h_muonID_MSEffNum_pt_pp_sign2");
  //  h_muonID_MSEffDen_pt[0][2] = (TH1D*) h_muonID_MSEffDen_pt[0][1]->Clone ("h_muonID_MSEffDen_pt_pp_sign2");
  //  h_muonID_MSEffNum_pt[1][2] = (TH1D*) h_muonID_MSEffNum_pt[1][1]->Clone ("h_muonID_MSEffNum_pt_PbPb18_sign2");
  //  h_muonID_MSEffDen_pt[1][2] = (TH1D*) h_muonID_MSEffDen_pt[1][1]->Clone ("h_muonID_MSEffDen_pt_PbPb18_sign2");
  //  h_muonID_MSEffNum_pt[2][2] = (TH1D*) h_muonID_MSEffNum_pt[2][1]->Clone ("h_muonID_MSEffNum_pt_PbPb15_sign2");
  //  h_muonID_MSEffDen_pt[2][2] = (TH1D*) h_muonID_MSEffDen_pt[2][1]->Clone ("h_muonID_MSEffDen_pt_PbPb15_sign2");
  //  h_muonID_MSEffNum_pt[0][2]->Add (h_muonID_MSEffNum_pt[0][0], -1);
  //  h_muonID_MSEffDen_pt[0][2]->Add (h_muonID_MSEffDen_pt[0][0], -1);
  //  h_muonID_MSEffNum_pt[1][2]->Add (h_muonID_MSEffNum_pt[1][0], -1);
  //  h_muonID_MSEffDen_pt[1][2]->Add (h_muonID_MSEffDen_pt[1][0], -1);
  //  h_muonID_MSEffNum_pt[2][2]->Add (h_muonID_MSEffNum_pt[2][0], -1);
  //  h_muonID_MSEffDen_pt[2][2]->Add (h_muonID_MSEffDen_pt[2][0], -1);

  //  t_muonID_MSEff_pt_pp[2] = new TEfficiency (*(h_muonID_MSEffNum_pt[0][2]), *(h_muonID_MSEffDen_pt[0][2]));
  //  t_muonID_MSEff_pt_PbPb18[2] = new TEfficiency (*(h_muonID_MSEffNum_pt[1][2]), *(h_muonID_MSEffDen_pt[1][2]));
  //  t_muonID_MSEff_pt_PbPb15[2] = new TEfficiency (*(h_muonID_MSEffNum_pt[2][2]), *(h_muonID_MSEffDen_pt[2][2]));

  //  h_muonID_MSEffNum_eta[0][2] = (TH1D*) h_muonID_MSEffNum_eta[0][1]->Clone ("h_muonID_MSEffNum_eta_pp_sign2");
  //  h_muonID_MSEffDen_eta[0][2] = (TH1D*) h_muonID_MSEffDen_eta[0][1]->Clone ("h_muonID_MSEffDen_eta_pp_sign2");
  //  h_muonID_MSEffNum_eta[1][2] = (TH1D*) h_muonID_MSEffNum_eta[1][1]->Clone ("h_muonID_MSEffNum_eta_PbPb18_sign2");
  //  h_muonID_MSEffDen_eta[1][2] = (TH1D*) h_muonID_MSEffDen_eta[1][1]->Clone ("h_muonID_MSEffDen_eta_PbPb18_sign2");
  //  h_muonID_MSEffNum_eta[2][2] = (TH1D*) h_muonID_MSEffNum_eta[2][1]->Clone ("h_muonID_MSEffNum_eta_PbPb15_sign2");
  //  h_muonID_MSEffDen_eta[2][2] = (TH1D*) h_muonID_MSEffDen_eta[2][1]->Clone ("h_muonID_MSEffDen_eta_PbPb15_sign2");
  //  h_muonID_MSEffNum_eta[0][2]->Add (h_muonID_MSEffNum_eta[0][0], -1);
  //  h_muonID_MSEffDen_eta[0][2]->Add (h_muonID_MSEffDen_eta[0][0], -1);
  //  h_muonID_MSEffNum_eta[1][2]->Add (h_muonID_MSEffNum_eta[1][0], -1);
  //  h_muonID_MSEffDen_eta[1][2]->Add (h_muonID_MSEffDen_eta[1][0], -1);
  //  h_muonID_MSEffNum_eta[2][2]->Add (h_muonID_MSEffNum_eta[2][0], -1);
  //  h_muonID_MSEffDen_eta[2][2]->Add (h_muonID_MSEffDen_eta[2][0], -1);

  //  t_muonID_MSEff_eta_pp[2] = new TEfficiency (*(h_muonID_MSEffNum_eta[0][2]), *(h_muonID_MSEffDen_eta[0][2]));
  //  t_muonID_MSEff_eta_PbPb18[2] = new TEfficiency (*(h_muonID_MSEffNum_eta[1][2]), *(h_muonID_MSEffDen_eta[1][2]));
  //  t_muonID_MSEff_eta_PbPb15[2] = new TEfficiency (*(h_muonID_MSEffNum_eta[2][2]), *(h_muonID_MSEffDen_eta[2][2]));

  //  h_muonID_MSEffNum_fcal[0][2] = (TH1D*) h_muonID_MSEffNum_fcal[0][1]->Clone ("h_muonID_MSEffNum_fcal_PbPb18_sign2");
  //  h_muonID_MSEffDen_fcal[0][2] = (TH1D*) h_muonID_MSEffDen_fcal[0][1]->Clone ("h_muonID_MSEffDen_fcal_PbPb18_sign2");
  //  h_muonID_MSEffNum_fcal[1][2] = (TH1D*) h_muonID_MSEffNum_fcal[1][1]->Clone ("h_muonID_MSEffNum_fcal_PbPb15_sign2");
  //  h_muonID_MSEffDen_fcal[1][2] = (TH1D*) h_muonID_MSEffDen_fcal[1][1]->Clone ("h_muonID_MSEffDen_fcal_PbPb15_sign2");
  //  h_muonID_MSEffNum_fcal[0][2]->Add (h_muonID_MSEffNum_fcal[0][0], -1);
  //  h_muonID_MSEffDen_fcal[0][2]->Add (h_muonID_MSEffDen_fcal[0][0], -1);
  //  h_muonID_MSEffNum_fcal[1][2]->Add (h_muonID_MSEffNum_fcal[1][0], -1);
  //  h_muonID_MSEffDen_fcal[1][2]->Add (h_muonID_MSEffDen_fcal[1][0], -1);

  //  t_muonID_MSEff_fcal_PbPb18[2] = new TEfficiency (*(h_muonID_MSEffNum_fcal[0][2]), *(h_muonID_MSEffDen_fcal[0][2]));
  //  t_muonID_MSEff_fcal_PbPb15[2] = new TEfficiency (*(h_muonID_MSEffNum_fcal[1][2]), *(h_muonID_MSEffDen_fcal[1][2]));
  //}



  TCanvas* c_1 = new TCanvas ("c_1", "", 1500, 750);
  TPad* p_1 = new TPad ("p_1", "", 0.0, 0.0, 0.5, 1.0); // 0.5, 0.4, 0.5, 1.0
  TPad* p_2 = new TPad ("p_2", "", 0.5, 0.0, 1.0, 1.0); // 0.5, 0.4, 1.0, 1.0
  //TPad* p_3 = new TPad ("p_3", "", 0.0, 0.0, 0.5, 0.4);
  //TPad* p_4 = new TPad ("p_4", "", 0.5, 0.0, 1.0, 0.4);
  //p_1->SetBottomMargin (0);
  //p_2->SetBottomMargin (0);
  //p_3->SetTopMargin (0);
  //p_3->SetBottomMargin (0.25);
  //p_4->SetTopMargin (0);
  //p_4->SetBottomMargin (0.25);
  p_1->Draw ();
  p_2->Draw ();
  //p_3->Draw ();
  //p_4->Draw ();

  TLine* line = new TLine (20, 1, 85, 1);
  line->SetLineColor (kPink-8);
  line->SetLineStyle (2);

  //TGAE* g_muonTrigEffPtpp = TEff2TGAE (t_muonTrigEff_pt_pp);
  TGAE* g_muonTrigEffPtPbPb18 = TEff2TGAE (t_muonTrigEff_pt_PbPb18);

  //g_muonTrigEffPtpp->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  g_muonTrigEffPtPbPb18->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  //g_muonTrigEffPtpp->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEffPtPbPb18->GetYaxis ()->SetTitle ("Muon trigger efficiency");

  //g_muonTrigEffPtpp->GetYaxis ()->SetRangeUser (0.0, 1.25);
  g_muonTrigEffPtPbPb18->GetYaxis ()->SetRangeUser (0.0, 1.25);
  //g_muonTrigEffPtpp->GetXaxis ()->SetRangeUser (20, 85);
  g_muonTrigEffPtPbPb18->GetXaxis ()->SetRangeUser (20, 85);

  //g_muonTrigEffPtpp->SetMarkerStyle (kOpenSquare);
  //g_muonTrigEffPtpp->SetMarkerColor (kAzure+2);
  //g_muonTrigEffPtpp->SetLineColor (kAzure+2);
  g_muonTrigEffPtPbPb18->SetMarkerStyle (kOpenCircle);
  g_muonTrigEffPtPbPb18->SetMarkerColor (kRed+1);
  g_muonTrigEffPtPbPb18->SetLineColor (kRed+1);

  p_1->cd ();
  //g_muonTrigEffPtpp->Draw ("AP");
  g_muonTrigEffPtPbPb18->Draw ("AP");

  line->DrawLine (20, 1, 85, 1);

  myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.22, 0.82, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.66, 0.885, kBlack, "HLT_mu14", 0.04);

  //TGAE* g_muonTrigEffEtapp = TEff2TGAE (t_muonTrigEff_eta_pp);
  TGAE* g_muonTrigEffEtaPbPb18 = TEff2TGAE (t_muonTrigEff_eta_PbPb18);

  //g_muonTrigEffEtapp->GetXaxis ()->SetTitle ("#eta");
  g_muonTrigEffEtaPbPb18->GetXaxis ()->SetTitle ("#eta");
  //g_muonTrigEffEtapp->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEffEtaPbPb18->GetYaxis ()->SetTitle ("Muon trigger efficiency");

  //g_muonTrigEffEtapp->GetYaxis ()->SetRangeUser (0.0, 1.25);
  g_muonTrigEffEtaPbPb18->GetYaxis ()->SetRangeUser (0.0, 1.25);
  //g_muonTrigEffEtapp->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  g_muonTrigEffEtaPbPb18->GetXaxis ()->SetRangeUser (-2.5, 2.5);

  //g_muonTrigEffEtapp->SetMarkerStyle (kOpenSquare);
  //g_muonTrigEffEtapp->SetMarkerColor (kAzure+2);
  //g_muonTrigEffEtapp->SetLineColor (kAzure+2);
  g_muonTrigEffEtaPbPb18->SetMarkerStyle (kOpenCircle);
  g_muonTrigEffEtaPbPb18->SetMarkerColor (kRed+1);
  g_muonTrigEffEtaPbPb18->SetLineColor (kRed+1);

  p_2->cd ();
  //g_muonTrigEffEtapp->Draw ("AP");
  g_muonTrigEffEtaPbPb18->Draw ("AP");

  line->DrawLine (-2.5, 1, 2.5, 1);

  myText (0.61, 0.89, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.042);
  //myMarkerTextNoLine (0.22, 0.893, kAzure+2, kOpenSquare, "2017 #it{pp}, 258 pb^{-1}", 1.5, 0.042);
  myMarkerTextNoLine (0.22, 0.893, kRed+1, kOpenCircle, "2018 Pb+Pb, 1.4 nb^{-1}", 1.5, 0.042);

  c_1->SaveAs ("../Plots/LeptonPerformance/MuonTrigEffs/pp_PbPb18_pt_eta.pdf");



  //TGAE* g_muonIDEffPtpp = TEff2TGAE (t_muonIDEff_pt_pp);
  TGAE* g_muonIDEffPtPbPb18 = TEff2TGAE (t_muonIDEff_pt_PbPb18);

  //g_muonIDEffPtpp->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  g_muonIDEffPtPbPb18->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  //g_muonIDEffPtpp->GetYaxis ()->SetTitle ("Muon ID efficiency");
  g_muonIDEffPtPbPb18->GetYaxis ()->SetTitle ("Muon ID efficiency");

  //g_muonIDEffPtpp->GetYaxis ()->SetRangeUser (0.0, 1.25);
  g_muonIDEffPtPbPb18->GetYaxis ()->SetRangeUser (0.0, 1.25);
  //g_muonIDEffPtpp->GetXaxis ()->SetRangeUser (20, 85);
  g_muonIDEffPtPbPb18->GetXaxis ()->SetRangeUser (20, 85);

  //g_muonIDEffPtpp->SetMarkerStyle (kOpenSquare);
  //g_muonIDEffPtpp->SetMarkerColor (kAzure+2);
  //g_muonIDEffPtpp->SetLineColor (kAzure+2);
  g_muonIDEffPtPbPb18->SetMarkerStyle (kOpenCircle);
  g_muonIDEffPtPbPb18->SetMarkerColor (kRed+1);
  g_muonIDEffPtPbPb18->SetLineColor (kRed+1);

  p_1->cd ();
  p_1->Clear ();
  //g_muonIDEffPtpp->Draw ("AP");
  g_muonIDEffPtPbPb18->Draw ("AP");

  line->DrawLine (20, 1, 85, 1);

  myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.22, 0.82, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.22, 0.20, kBlack, "Medium muons", 0.04);

  //TGAE* g_muonIDEffEtapp = TEff2TGAE (t_muonIDEff_eta_pp);
  TGAE* g_muonIDEffEtaPbPb18 = TEff2TGAE (t_muonIDEff_eta_PbPb18);

  //g_muonIDEffEtapp->GetXaxis ()->SetTitle ("#eta");
  g_muonIDEffEtaPbPb18->GetXaxis ()->SetTitle ("#eta");
  //g_muonIDEffEtapp->GetYaxis ()->SetTitle ("Muon ID efficiency");
  g_muonIDEffEtaPbPb18->GetYaxis ()->SetTitle ("Muon ID efficiency");

  //g_muonIDEffEtapp->GetYaxis ()->SetRangeUser (0.0, 1.25);
  g_muonIDEffEtaPbPb18->GetYaxis ()->SetRangeUser (0.0, 1.25);
  //g_muonIDEffEtapp->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  g_muonIDEffEtaPbPb18->GetXaxis ()->SetRangeUser (-2.5, 2.5);

  //g_muonIDEffEtapp->SetMarkerStyle (kOpenSquare);
  //g_muonIDEffEtapp->SetMarkerColor (kAzure+2);
  //g_muonIDEffEtapp->SetLineColor (kAzure+2);
  g_muonIDEffEtaPbPb18->SetMarkerStyle (kOpenCircle);
  g_muonIDEffEtaPbPb18->SetMarkerColor (kRed+1);
  g_muonIDEffEtaPbPb18->SetLineColor (kRed+1);

  p_2->cd ();
  p_2->Clear ();
  //g_muonIDEffEtapp->Draw ("AP");
  g_muonIDEffEtaPbPb18->Draw ("AP");

  line->DrawLine (-2.5, 1, 2.5, 1);

  myText (0.61, 0.89, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.042);
  //myMarkerTextNoLine (0.22, 0.893, kAzure+2, kOpenSquare, "2017 #it{pp}, 258 pb^{-1}", 1.5, 0.042);
  myMarkerTextNoLine (0.22, 0.893, kRed+1, kOpenCircle, "2018 Pb+Pb, 1.4 nb^{-1}", 1.5, 0.042);

  c_1->SaveAs ("../Plots/LeptonPerformance/MuonIDEffs/pp_PbPb18_pt_eta.pdf");



  TCanvas* c_2 = new TCanvas ("c_2", "", 800, 600);
  FormatTH2Canvas (c_2, true);

  TH1* h2_muonTrigEff_eta_phi_pp = t_muonTrigEff_eta_phi_pp->GetCopyPassedHisto ();
  h2_muonTrigEff_eta_phi_pp->Divide (t_muonTrigEff_eta_phi_pp->GetTotalHistogram ());
  h2_muonTrigEff_eta_phi_pp->GetXaxis ()->SetTitle ("#eta_{#mu}");
  h2_muonTrigEff_eta_phi_pp->GetYaxis ()->SetTitle ("#phi_{#mu}");
  h2_muonTrigEff_eta_phi_pp->GetYaxis ()->SetTitleOffset (1.1);
  h2_muonTrigEff_eta_phi_pp->GetZaxis ()->SetRangeUser (0, 1);
  h2_muonTrigEff_eta_phi_pp->GetZaxis ()->SetTitle ("Muon trigger efficiency");
  h2_muonTrigEff_eta_phi_pp->Draw ("colz");
  myText (0.23, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.23, 0.84, kBlack, "#it{pp}, 5.02 TeV", 0.045);
  myText (0.23, 0.78, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.045);
  myText (0.23, 0.72, kBlack, "HLT_mu14", 0.045);
  c_2->SaveAs ("../Plots/LeptonPerformance/MuonTrigEffs/pp_eta_phi.pdf");

  TH1* h2_muonTrigEff_eta_phi_PbPb18 = t_muonTrigEff_eta_phi_PbPb18->GetCopyPassedHisto ();
  h2_muonTrigEff_eta_phi_PbPb18->Divide (t_muonTrigEff_eta_phi_PbPb18->GetTotalHistogram ());
  h2_muonTrigEff_eta_phi_PbPb18->GetXaxis ()->SetTitle ("#eta_{#mu}");
  h2_muonTrigEff_eta_phi_PbPb18->GetYaxis ()->SetTitle ("#phi_{#mu}");
  h2_muonTrigEff_eta_phi_PbPb18->GetYaxis ()->SetTitleOffset (1.1);
  h2_muonTrigEff_eta_phi_PbPb18->GetZaxis ()->SetRangeUser (0, 1);
  h2_muonTrigEff_eta_phi_PbPb18->GetZaxis ()->SetTitle ("Muon trigger efficiency");
  h2_muonTrigEff_eta_phi_PbPb18->Draw ("colz");
  myText (0.23, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.23, 0.84, kBlack, "Pb+Pb, 5.02 TeV", 0.045);
  myText (0.23, 0.78, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.045);
  myText (0.23, 0.72, kBlack, "HLT_mu14", 0.045);
  c_2->SaveAs ("../Plots/LeptonPerformance/MuonTrigEffs/PbPb18_eta_phi.pdf");

  TH1* h2_muonIDEff_eta_phi_pp = t_muonIDEff_eta_phi_pp->GetCopyPassedHisto ();
  h2_muonIDEff_eta_phi_pp->Divide (t_muonIDEff_eta_phi_pp->GetTotalHistogram ());
  h2_muonIDEff_eta_phi_pp->GetXaxis ()->SetTitle ("#eta_{#mu}");
  h2_muonIDEff_eta_phi_pp->GetYaxis ()->SetTitle ("#phi_{#mu}");
  h2_muonIDEff_eta_phi_pp->GetYaxis ()->SetTitleOffset (1.1);
  h2_muonIDEff_eta_phi_pp->GetZaxis ()->SetRangeUser (0, 1);
  h2_muonIDEff_eta_phi_pp->GetZaxis ()->SetTitle ("Muon ID efficiency");
  h2_muonIDEff_eta_phi_pp->Draw ("colz");
  myText (0.23, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.23, 0.84, kBlack, "#it{pp}, 5.02 TeV", 0.045);
  myText (0.23, 0.78, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.045);
  myText (0.23, 0.72, kBlack, "Medium muons", 0.045);
  c_2->SaveAs ("../Plots/LeptonPerformance/MuonIDEffs/pp_eta_phi.pdf");

  TH1* h2_muonIDEff_eta_phi_PbPb18 = t_muonIDEff_eta_phi_PbPb18->GetCopyPassedHisto ();
  h2_muonIDEff_eta_phi_PbPb18->Divide (t_muonIDEff_eta_phi_PbPb18->GetTotalHistogram ());
  h2_muonIDEff_eta_phi_PbPb18->GetXaxis ()->SetTitle ("#eta_{#mu}");
  h2_muonIDEff_eta_phi_PbPb18->GetYaxis ()->SetTitle ("#phi_{#mu}");
  h2_muonIDEff_eta_phi_PbPb18->GetYaxis ()->SetTitleOffset (1.1);
  h2_muonIDEff_eta_phi_PbPb18->GetZaxis ()->SetRangeUser (0, 1);
  h2_muonIDEff_eta_phi_PbPb18->GetZaxis ()->SetTitle ("Muon ID efficiency");
  h2_muonIDEff_eta_phi_PbPb18->Draw ("colz");
  myText (0.23, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.23, 0.84, kBlack, "Pb+Pb, 5.02 TeV", 0.045);
  myText (0.23, 0.78, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.045);
  myText (0.23, 0.72, kBlack, "Medium muons", 0.045);
  c_2->SaveAs ("../Plots/LeptonPerformance/MuonIDEffs/PbPb18_eta_phi.pdf");


  TCanvas* c_3 = new TCanvas ("c_3", "", 800, 600);
  FormatTH2Canvas (c_3, false);
  c_3->cd ();
  TH1D* h_muonTrigEff_fcal_PbPb18 = (TH1D*) h_muonTrigEffNum_fcal[0]->Clone ("h_muonTrigEff_fcal_PbPb18");
  BinomialDivide (h_muonTrigEff_fcal_PbPb18, h_muonTrigEffDen_fcal[0]);
  //TGAE* g_muonTrigEff_fcal_PbPb18 = TEff2TGAE (t_muonTrigEff_fcal_PbPb18);
  TGAE* g_muonTrigEff_fcal_PbPb18 = make_graph (h_muonTrigEff_fcal_PbPb18);
  g_muonTrigEff_fcal_PbPb18->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [TeV]");
  g_muonTrigEff_fcal_PbPb18->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEff_fcal_PbPb18->GetXaxis ()->SetRangeUser (-0.5, 5);
  g_muonTrigEff_fcal_PbPb18->GetYaxis ()->SetRangeUser (0.0, 1.1);
  g_muonTrigEff_fcal_PbPb18->SetMarkerStyle (kOpenCircle);
  g_muonTrigEff_fcal_PbPb18->Draw ("AP");

  TF1* f_muonTrigEff_fcal_PbPb18 = new TF1 ("f_muonTrigEff_fcal_PbPb18", "[0]+[1]*x", 0, 5);
  g_muonTrigEff_fcal_PbPb18->Fit (f_muonTrigEff_fcal_PbPb18, "RNQ0");
  f_muonTrigEff_fcal_PbPb18->SetLineStyle (2);
  f_muonTrigEff_fcal_PbPb18->SetLineColor (kGreen+2);
  f_muonTrigEff_fcal_PbPb18->SetLineWidth (2);
  f_muonTrigEff_fcal_PbPb18->Draw ("same");

  line->DrawLine (-0.5, 1, 5, 1);

  myText (0.2, 0.35, kBlack, Form ("Slope = %s TeV^{-1}", FormatMeasurement (f_muonTrigEff_fcal_PbPb18->GetParameter (1), f_muonTrigEff_fcal_PbPb18->GetParError (1), 2)), 0.04);
  myText (0.2, 0.30, kBlack, Form ("y-int. = %s", FormatMeasurement (f_muonTrigEff_fcal_PbPb18->GetParameter (0), f_muonTrigEff_fcal_PbPb18->GetParError (0), 2)), 0.04);
  myText (0.2, 0.25, kBlack, Form ("#chi^{2} / dof = %.2f / %i", f_muonTrigEff_fcal_PbPb18->GetChisquare (), f_muonTrigEff_fcal_PbPb18->GetNDF ()), 0.04);

  myText (0.62, 0.4, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.62, 0.35, kBlack, "2018 Pb+Pb, 1.4 nb^{-1}", 0.04);
  myText (0.62, 0.3, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.62, 0.25, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.04);
  myText (0.62, 0.20, kBlack, "Medium muons", 0.04);
  myMarkerTextNoLine (0.65, 0.900, kBlack, kOpenCircle, "HLT_mu14", 1.25, 0.04);
  c_3->SaveAs ("../Plots/LeptonPerformance/MuonTrigEffs/PbPb18_cent.pdf");
  SaferDelete (f_muonTrigEff_fcal_PbPb18);


  TH1D* h_electronTrigEffNum_pt[3];
  TH1D* h_electronTrigEffDen_pt[3];
  TH1D* h_electronTrigEffNum_eta[3];
  TH1D* h_electronTrigEffDen_eta[3];
  TH2D* h2_electronTrigEffNum_pt_eta[3];
  TH2D* h2_electronTrigEffDen_pt_eta[3];

  TH1D* h_electronTrigEffNum_fcal[2];
  TH1D* h_electronTrigEffDen_fcal[2];

  TH1D* h_electronIDEffNum_pt[3];
  TH1D* h_electronIDEffDen_pt[3];
  TH1D* h_electronIDEffNum_eta[3];
  TH1D* h_electronIDEffDen_eta[3];
  TH2D* h2_electronIDEffNum_pt_eta[3];
  TH2D* h2_electronIDEffDen_pt_eta[3];

  TH1D* h_electronIDEffNum_fcal[2];
  TH1D* h_electronIDEffDen_fcal[2];

  h_electronTrigEffNum_pt[0] = (TH1D*) inFile->Get ("h_electronTrigEffNum_pt_pp");
  h_electronTrigEffDen_pt[0] = (TH1D*) inFile->Get ("h_electronTrigEffDen_pt_pp");
  h_electronTrigEffNum_pt[1] = (TH1D*) inFile->Get ("h_electronTrigEffNum_pt_PbPb18");
  h_electronTrigEffDen_pt[1] = (TH1D*) inFile->Get ("h_electronTrigEffDen_pt_PbPb18");
  //h_electronTrigEffNum_pt[2] = (TH1D*) inFile->Get ("h_electronTrigEffNum_pt_PbPb15");
  //h_electronTrigEffDen_pt[2] = (TH1D*) inFile->Get ("h_electronTrigEffDen_pt_PbPb15");
  h_electronTrigEffNum_eta[0] = (TH1D*) inFile->Get ("h_electronTrigEffNum_eta_pp");
  h_electronTrigEffDen_eta[0] = (TH1D*) inFile->Get ("h_electronTrigEffDen_eta_pp");
  h_electronTrigEffNum_eta[1] = (TH1D*) inFile->Get ("h_electronTrigEffNum_eta_PbPb18");
  h_electronTrigEffDen_eta[1] = (TH1D*) inFile->Get ("h_electronTrigEffDen_eta_PbPb18");
  //h_electronTrigEffNum_eta[2] = (TH1D*) inFile->Get ("h_electronTrigEffNum_eta_PbPb15");
  //h_electronTrigEffDen_eta[2] = (TH1D*) inFile->Get ("h_electronTrigEffDen_eta_PbPb15");
  h2_electronTrigEffNum_pt_eta[0] = (TH2D*) inFile->Get ("h2_electronTrigEffNum_pt_eta_pp");
  h2_electronTrigEffDen_pt_eta[0] = (TH2D*) inFile->Get ("h2_electronTrigEffDen_pt_eta_pp");
  h2_electronTrigEffNum_pt_eta[1] = (TH2D*) inFile->Get ("h2_electronTrigEffNum_pt_eta_PbPb18");
  h2_electronTrigEffDen_pt_eta[1] = (TH2D*) inFile->Get ("h2_electronTrigEffDen_pt_eta_PbPb18");
  //h2_electronTrigEffNum_pt_eta[2] = (TH2D*) inFile->Get ("h2_electronTrigEffNum_pt_eta_PbPb15");
  //h2_electronTrigEffDen_pt_eta[2] = (TH2D*) inFile->Get ("h2_electronTrigEffDen_pt_eta_PbPb15");

  h_electronTrigEffNum_fcal[0] = (TH1D*) inFile->Get ("h_electronTrigEffNum_fcal_PbPb18");
  h_electronTrigEffDen_fcal[0] = (TH1D*) inFile->Get ("h_electronTrigEffDen_fcal_PbPb18");
  //h_electronTrigEffNum_fcal[1] = (TH1D*) inFile->Get ("h_electronTrigEffNum_fcal_PbPb15");
  //h_electronTrigEffDen_fcal[1] = (TH1D*) inFile->Get ("h_electronTrigEffDen_fcal_PbPb15");
  h_electronIDEffNum_pt[0] = (TH1D*) inFile->Get ("h_electronIDEffNum_pt_pp");
  h_electronIDEffDen_pt[0] = (TH1D*) inFile->Get ("h_electronIDEffDen_pt_pp");
  h_electronIDEffNum_pt[1] = (TH1D*) inFile->Get ("h_electronIDEffNum_pt_PbPb18");
  h_electronIDEffDen_pt[1] = (TH1D*) inFile->Get ("h_electronIDEffDen_pt_PbPb18");
  //h_electronIDEffNum_pt[2] = (TH1D*) inFile->Get ("h_electronIDEffNum_pt_PbPb15");
  //h_electronIDEffDen_pt[2] = (TH1D*) inFile->Get ("h_electronIDEffDen_pt_PbPb15");
  h_electronIDEffNum_eta[0] = (TH1D*) inFile->Get ("h_electronIDEffNum_eta_pp");
  h_electronIDEffDen_eta[0] = (TH1D*) inFile->Get ("h_electronIDEffDen_eta_pp");
  h_electronIDEffNum_eta[1] = (TH1D*) inFile->Get ("h_electronIDEffNum_eta_PbPb18");
  h_electronIDEffDen_eta[1] = (TH1D*) inFile->Get ("h_electronIDEffDen_eta_PbPb18");
  //h_electronIDEffNum_eta[2] = (TH1D*) inFile->Get ("h_electronIDEffNum_eta_PbPb15");
  //h_electronIDEffDen_eta[2] = (TH1D*) inFile->Get ("h_electronIDEffDen_eta_PbPb15");
  h2_electronIDEffNum_pt_eta[0] = (TH2D*) inFile->Get ("h2_electronIDEffNum_pt_eta_pp");
  h2_electronIDEffDen_pt_eta[0] = (TH2D*) inFile->Get ("h2_electronIDEffDen_pt_eta_pp");
  h2_electronIDEffNum_pt_eta[1] = (TH2D*) inFile->Get ("h2_electronIDEffNum_pt_eta_PbPb18");
  h2_electronIDEffDen_pt_eta[1] = (TH2D*) inFile->Get ("h2_electronIDEffDen_pt_eta_PbPb18");
  //h2_electronIDEffNum_pt_eta[2] = (TH2D*) inFile->Get ("h2_electronIDEffNum_pt_eta_PbPb15");
  //h2_electronIDEffDen_pt_eta[2] = (TH2D*) inFile->Get ("h2_electronIDEffDen_pt_eta_PbPb15");

  h_electronIDEffNum_fcal[0] = (TH1D*) inFile->Get ("h_electronIDEffNum_fcal_PbPb18");
  h_electronIDEffDen_fcal[0] = (TH1D*) inFile->Get ("h_electronIDEffDen_fcal_PbPb18");
  //h_electronIDEffNum_fcal[1] = (TH1D*) inFile->Get ("h_electronIDEffNum_fcal_PbPb15");
  //h_electronIDEffDen_fcal[1] = (TH1D*) inFile->Get ("h_electronIDEffDen_fcal_PbPb15");

  t_electronTrigEff_pt_pp = new TEfficiency (*(h_electronTrigEffNum_pt[0]), *(h_electronTrigEffDen_pt[0]));
  t_electronTrigEff_pt_PbPb18 = new TEfficiency (*(h_electronTrigEffNum_pt[1]), *(h_electronTrigEffDen_pt[1]));
  //t_electronTrigEff_pt_PbPb15 = new TEfficiency (*(h_electronTrigEffNum_pt[2]), *(h_electronTrigEffDen_pt[2]));
  t_electronTrigEff_eta_pp = new TEfficiency (*(h_electronTrigEffNum_eta[0]), *(h_electronTrigEffDen_eta[0]));
  t_electronTrigEff_eta_PbPb18 = new TEfficiency (*(h_electronTrigEffNum_eta[1]), *(h_electronTrigEffDen_eta[1]));
  //t_electronTrigEff_eta_PbPb15 = new TEfficiency (*(h_electronTrigEffNum_eta[2]), *(h_electronTrigEffDen_eta[2]));
  t_electronTrigEff_pt_eta_pp = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[0]), *(h2_electronTrigEffDen_pt_eta[0]));
  t_electronTrigEff_pt_eta_PbPb18 = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[1]), *(h2_electronTrigEffDen_pt_eta[1]));
  //t_electronTrigEff_pt_eta_PbPb15 = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[2]), *(h2_electronTrigEffDen_pt_eta[2]));

  t_electronTrigEff_fcal_PbPb18 = new TEfficiency (*(h_electronTrigEffNum_fcal[0]), *(h_electronTrigEffDen_fcal[0]));
  //t_electronTrigEff_fcal_PbPb15 = new TEfficiency (*(h_electronTrigEffNum_fcal[1]), *(h_electronTrigEffDen_fcal[1]));

  t_electronIDEff_pt_pp = new TEfficiency (*(h_electronIDEffNum_pt[0]), *(h_electronIDEffDen_pt[0]));
  t_electronIDEff_pt_PbPb18 = new TEfficiency (*(h_electronIDEffNum_pt[1]), *(h_electronIDEffDen_pt[1]));
  //t_electronIDEff_pt_PbPb15 = new TEfficiency (*(h_electronIDEffNum_pt[2]), *(h_electronIDEffDen_pt[2]));
  t_electronIDEff_eta_pp = new TEfficiency (*(h_electronIDEffNum_eta[0]), *(h_electronIDEffDen_eta[0]));
  t_electronIDEff_eta_PbPb18 = new TEfficiency (*(h_electronIDEffNum_eta[1]), *(h_electronIDEffDen_eta[1]));
  //t_electronIDEff_eta_PbPb15 = new TEfficiency (*(h_electronIDEffNum_eta[2]), *(h_electronIDEffDen_eta[2]));
  t_electronIDEff_pt_eta_pp = new TEfficiency (*(h2_electronIDEffNum_pt_eta[0]), *(h2_electronIDEffDen_pt_eta[0]));
  t_electronIDEff_pt_eta_PbPb18 = new TEfficiency (*(h2_electronIDEffNum_pt_eta[1]), *(h2_electronIDEffDen_pt_eta[1]));
  //t_electronIDEff_pt_eta_PbPb15 = new TEfficiency (*(h2_electronIDEffNum_pt_eta[2]), *(h2_electronIDEffDen_pt_eta[2]));

  t_electronIDEff_fcal_PbPb18 = new TEfficiency (*(h_electronIDEffNum_fcal[0]), *(h_electronIDEffDen_fcal[0]));
  //t_electronIDEff_fcal_PbPb15 = new TEfficiency (*(h_electronIDEffNum_fcal[1]), *(h_electronIDEffDen_fcal[1]));


  TGAE* g_electronTrigEffPtpp = TEff2TGAE (t_electronTrigEff_pt_pp);
  TGAE* g_electronTrigEffPtPbPb18 = TEff2TGAE (t_electronTrigEff_pt_PbPb18);
  //TGAE* g_electronTrigEffPtPbPb15 = TEff2TGAE (t_electronTrigEff_pt_PbPb15);

  g_electronTrigEffPtpp->GetXaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  g_electronTrigEffPtPbPb18->GetXaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  g_electronTrigEffPtpp->GetYaxis ()->SetTitle ("Electron trigger efficiency");
  g_electronTrigEffPtPbPb18->GetYaxis ()->SetTitle ("Electron trigger efficiency");

  g_electronTrigEffPtpp->GetYaxis ()->SetRangeUser (0.0, 1.25);
  g_electronTrigEffPtPbPb18->GetYaxis ()->SetRangeUser (0.0, 1.25);
  g_electronTrigEffPtpp->GetXaxis ()->SetRangeUser (20, 100);
  g_electronTrigEffPtPbPb18->GetXaxis ()->SetRangeUser (20, 100);

  g_electronTrigEffPtpp->SetMarkerStyle (kOpenSquare);
  g_electronTrigEffPtpp->SetMarkerColor (kAzure+2);
  g_electronTrigEffPtpp->SetLineColor (kAzure+2);
  g_electronTrigEffPtPbPb18->SetMarkerStyle (kOpenCircle);
  g_electronTrigEffPtPbPb18->SetMarkerColor (kRed+1);
  g_electronTrigEffPtPbPb18->SetLineColor (kRed+1);

  p_1->cd ();
  p_1->Clear ();
  g_electronTrigEffPtpp->Draw ("AP");
  g_electronTrigEffPtPbPb18->Draw ("P");

  line->DrawLine (20, 1, 100, 1);

  myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.22, 0.82, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  //myText (0.50, 0.885, kBlack, "HLT_e15_lhloose(_ion)_L1EM12", 0.04);

  TGAE* g_electronTrigEffEtapp = TEff2TGAE (t_electronTrigEff_eta_pp);
  TGAE* g_electronTrigEffEtaPbPb18 = TEff2TGAE (t_electronTrigEff_eta_PbPb18);

  g_electronTrigEffEtapp->GetXaxis ()->SetTitle ("#eta");
  g_electronTrigEffEtaPbPb18->GetXaxis ()->SetTitle ("#eta");
  g_electronTrigEffEtapp->GetYaxis ()->SetTitle ("Electron trigger efficiency");
  g_electronTrigEffEtaPbPb18->GetYaxis ()->SetTitle ("Electron trigger efficiency");

  g_electronTrigEffEtapp->GetYaxis ()->SetRangeUser (0.0, 1.25);
  g_electronTrigEffEtaPbPb18->GetYaxis ()->SetRangeUser (0.0, 1.25);
  g_electronTrigEffEtapp->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  g_electronTrigEffEtaPbPb18->GetXaxis ()->SetRangeUser (-2.5, 2.5);

  g_electronTrigEffEtapp->SetMarkerStyle (kOpenSquare);
  g_electronTrigEffEtapp->SetMarkerColor (kAzure+2);
  g_electronTrigEffEtapp->SetLineColor (kAzure+2);
  g_electronTrigEffEtaPbPb18->SetMarkerStyle (kOpenCircle);
  g_electronTrigEffEtaPbPb18->SetMarkerColor (kRed+1);
  g_electronTrigEffEtaPbPb18->SetLineColor (kRed+1);

  p_2->cd ();
  p_2->Clear ();
  g_electronTrigEffEtapp->Draw ("AP");
  g_electronTrigEffEtaPbPb18->Draw ("P");

  line->DrawLine (-2.5, 1, 2.5, 1);

  myText (0.61, 0.89, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.042);
  myMarkerTextNoLine (0.22, 0.893, kAzure+2, kOpenSquare, "2017 #it{pp}, 258 pb^{-1}", 1.5, 0.042);
  myMarkerTextNoLine (0.22, 0.835, kRed+1, kOpenCircle, "2018 Pb+Pb, 1.7 nb^{-1}", 1.5, 0.042);

  c_1->SaveAs ("../Plots/LeptonPerformance/ElectronTrigEffs/pp_PbPb18_pt_eta.pdf");


  TGAE* g_electronIDEffPtpp = TEff2TGAE (t_electronIDEff_pt_pp);
  TGAE* g_electronIDEffPtPbPb18 = TEff2TGAE (t_electronIDEff_pt_PbPb18);

  g_electronIDEffPtpp->GetXaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  g_electronIDEffPtPbPb18->GetXaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  g_electronIDEffPtpp->GetYaxis ()->SetTitle ("Electron ID efficiency");
  g_electronIDEffPtPbPb18->GetYaxis ()->SetTitle ("Electron ID efficiency");

  g_electronIDEffPtpp->GetYaxis ()->SetRangeUser (0., 1.25);
  g_electronIDEffPtPbPb18->GetYaxis ()->SetRangeUser (0., 1.25);
  g_electronIDEffPtpp->GetXaxis ()->SetRangeUser (20, 100);
  g_electronIDEffPtPbPb18->GetXaxis ()->SetRangeUser (20, 100);

  g_electronIDEffPtpp->SetMarkerStyle (kOpenSquare);
  g_electronIDEffPtpp->SetMarkerColor (kAzure+2);
  g_electronIDEffPtpp->SetLineColor (kAzure+2);
  g_electronIDEffPtPbPb18->SetMarkerStyle (kOpenCircle);
  g_electronIDEffPtPbPb18->SetMarkerColor (kRed+1);
  g_electronIDEffPtPbPb18->SetLineColor (kRed+1);

  p_1->cd ();
  p_1->Clear ();
  g_electronIDEffPtpp->Draw ("AP");
  g_electronIDEffPtPbPb18->Draw ("P");

  line->DrawLine (20, 1, 100, 1);

  myText (0.22, 0.88, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.22, 0.82, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.22, 0.20, kBlack, "LHLoose(_HI) electrons", 0.04);

  TGAE* g_electronIDEffEtapp = TEff2TGAE (t_electronIDEff_eta_pp);
  TGAE* g_electronIDEffEtaPbPb18 = TEff2TGAE (t_electronIDEff_eta_PbPb18);

  g_electronIDEffEtapp->GetXaxis ()->SetTitle ("#eta");
  g_electronIDEffEtaPbPb18->GetXaxis ()->SetTitle ("#eta");
  g_electronIDEffEtapp->GetYaxis ()->SetTitle ("Electron ID efficiency");
  g_electronIDEffEtaPbPb18->GetYaxis ()->SetTitle ("Electron ID efficiency");

  g_electronIDEffEtapp->GetYaxis ()->SetRangeUser (0., 1.25);
  g_electronIDEffEtaPbPb18->GetYaxis ()->SetRangeUser (0., 1.25);
  g_electronIDEffEtapp->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  g_electronIDEffEtaPbPb18->GetXaxis ()->SetRangeUser (-2.5, 2.5);

  g_electronIDEffEtapp->SetMarkerStyle (kOpenSquare);
  g_electronIDEffEtapp->SetMarkerColor (kAzure+2);
  g_electronIDEffEtapp->SetLineColor (kAzure+2);
  g_electronIDEffEtaPbPb18->SetMarkerStyle (kOpenCircle);
  g_electronIDEffEtaPbPb18->SetMarkerColor (kRed+1);
  g_electronIDEffEtaPbPb18->SetLineColor (kRed+1);

  p_2->cd ();
  p_2->Clear ();
  g_electronIDEffEtapp->Draw ("AP");
  g_electronIDEffEtaPbPb18->Draw ("P");

  line->DrawLine (-2.5, 1, 2.5, 1);

  myText (0.61, 0.89, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.042);
  myMarkerTextNoLine (0.22, 0.893, kAzure+2, kOpenSquare, "2017 #it{pp}, 258 pb^{-1}", 1.5, 0.042);
  myMarkerTextNoLine (0.22, 0.835, kRed+1, kOpenCircle, "2018 Pb+Pb, 1.7 nb^{-1}", 1.5, 0.042);

  c_1->SaveAs ("../Plots/LeptonPerformance/ElectronIDEffs/pp_PbPb18_pt_eta.pdf");


  FormatTH2Canvas (c_2, true);
  c_2->cd ();
  TH1* h2_electronTrigEff_pt_eta_pp = t_electronTrigEff_pt_eta_pp->GetCopyPassedHisto ();
  h2_electronTrigEff_pt_eta_pp->Divide (t_electronTrigEff_pt_eta_pp->GetTotalHistogram ());
  h2_electronTrigEff_pt_eta_pp->GetXaxis ()->SetTitle ("#eta_{e}");
  h2_electronTrigEff_pt_eta_pp->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronTrigEff_pt_eta_pp->GetYaxis ()->SetTitleOffset (1.1);
  h2_electronTrigEff_pt_eta_pp->GetZaxis ()->SetRangeUser (0, 1);
  h2_electronTrigEff_pt_eta_pp->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronTrigEff_pt_eta_pp->GetZaxis ()->SetTitle ("Electron trigger efficiency");
  h2_electronTrigEff_pt_eta_pp->Draw ("colz");
  myText (0.23, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.23, 0.84, kBlack, "2017 #it{pp}, 5.02 TeV", 0.045);
  //myText (0.23, 0.78, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.045);
  myText (0.23, 0.78, kBlack, "HLT_e15_lhloose_L1EM12", 0.045);
  c_2->SaveAs ("../Plots/LeptonPerformance/ElectronTrigEffs/pp_pt_eta.pdf");

  TH1* h2_electronTrigEff_pt_eta_PbPb18 = t_electronTrigEff_pt_eta_PbPb18->GetCopyPassedHisto ();
  h2_electronTrigEff_pt_eta_PbPb18->Divide (t_electronTrigEff_pt_eta_PbPb18->GetTotalHistogram ());
  h2_electronTrigEff_pt_eta_PbPb18->GetXaxis ()->SetTitle ("#eta_{e}");
  h2_electronTrigEff_pt_eta_PbPb18->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronTrigEff_pt_eta_PbPb18->GetYaxis ()->SetTitleOffset (1.1);
  h2_electronTrigEff_pt_eta_PbPb18->GetZaxis ()->SetRangeUser (0, 1);
  h2_electronTrigEff_pt_eta_PbPb18->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronTrigEff_pt_eta_PbPb18->GetZaxis ()->SetTitle ("Electron trigger efficiency");
  h2_electronTrigEff_pt_eta_PbPb18->Draw ("colz");
  myText (0.23, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.23, 0.84, kBlack, "2018 Pb+Pb, 5.02 TeV", 0.045);
  //myText (0.23, 0.78, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.045);
  myText (0.23, 0.78, kBlack, "HLT_e15_lhloose_ion_L1EM12", 0.045);
  c_2->SaveAs ("../Plots/LeptonPerformance/ElectronTrigEffs/PbPb18_pt_eta.pdf");

  TH1* h2_electronIDEff_pt_eta_pp = t_electronIDEff_pt_eta_pp->GetCopyPassedHisto ();
  h2_electronIDEff_pt_eta_pp->Divide (t_electronIDEff_pt_eta_pp->GetTotalHistogram ());
  h2_electronIDEff_pt_eta_pp->GetXaxis ()->SetTitle ("#eta_{e}");
  h2_electronIDEff_pt_eta_pp->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronIDEff_pt_eta_pp->GetYaxis ()->SetTitleOffset (1.1);
  h2_electronIDEff_pt_eta_pp->GetZaxis ()->SetRangeUser (0, 1);
  h2_electronIDEff_pt_eta_pp->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronIDEff_pt_eta_pp->GetZaxis ()->SetTitle ("Electron ID efficiency");
  h2_electronIDEff_pt_eta_pp->Draw ("colz");
  myText (0.23, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.23, 0.84, kBlack, "2017 #it{pp}, 5.02 TeV", 0.045);
  //myText (0.23, 0.78, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.045);
  myText (0.23, 0.78, kBlack, "LHLoose electrons", 0.045);
  c_2->SaveAs ("../Plots/LeptonPerformance/ElectronIDEffs/pp_pt_eta.pdf");

  TH1* h2_electronIDEff_pt_eta_PbPb18 = t_electronIDEff_pt_eta_PbPb18->GetCopyPassedHisto ();
  h2_electronIDEff_pt_eta_PbPb18->Divide (t_electronIDEff_pt_eta_PbPb18->GetTotalHistogram ());
  h2_electronIDEff_pt_eta_PbPb18->GetXaxis ()->SetTitle ("#eta_{e}");
  h2_electronIDEff_pt_eta_PbPb18->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronIDEff_pt_eta_PbPb18->GetYaxis ()->SetTitleOffset (1.1);
  h2_electronIDEff_pt_eta_PbPb18->GetZaxis ()->SetRangeUser (0, 1);
  h2_electronIDEff_pt_eta_PbPb18->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronIDEff_pt_eta_PbPb18->GetZaxis ()->SetTitle ("Electron ID efficiency");
  h2_electronIDEff_pt_eta_PbPb18->Draw ("colz");
  myText (0.23, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.23, 0.84, kBlack, "2018 Pb+Pb, 5.02 TeV", 0.045);
  //myText (0.23, 0.78, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.045);
  myText (0.23, 0.78, kBlack, "LHLoose_HI electrons", 0.045);
  c_2->SaveAs ("../Plots/LeptonPerformance/ElectronIDEffs/PbPb18_pt_eta.pdf");


  FormatTH2Canvas (c_3, false);
  c_3->cd ();
  TH1D* h_electronTrigEff_fcal_PbPb18 = (TH1D*) h_electronTrigEffNum_fcal[0]->Clone ("h_electronTrigEff_fcal_PbPb18");
  BinomialDivide (h_electronTrigEff_fcal_PbPb18, h_electronTrigEffDen_fcal[0]);
  //TGAE* g_electronTrigEff_fcal_PbPb18 = TEff2TGAE (t_electronTrigEff_fcal_PbPb18);
  TGAE* g_electronTrigEff_fcal_PbPb18 = make_graph (h_electronTrigEff_fcal_PbPb18);
  g_electronTrigEff_fcal_PbPb18->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [TeV]");
  g_electronTrigEff_fcal_PbPb18->GetYaxis ()->SetTitle ("Electron trigger efficiency");
  g_electronTrigEff_fcal_PbPb18->GetXaxis ()->SetRangeUser (-0.5, 5);
  g_electronTrigEff_fcal_PbPb18->GetYaxis ()->SetRangeUser (0.0, 1.10);
  g_electronTrigEff_fcal_PbPb18->SetMarkerStyle (kOpenCircle);
  g_electronTrigEff_fcal_PbPb18->Draw ("AP");

  TF1* f_electronTrigEff_fcal_PbPb18 = new TF1 ("f_electronTrigEff_fcal_PbPb18", "[0]+[1]*x", 0, 5);
  f_electronTrigEff_fcal_PbPb18->SetParameter (0, 1);
  f_electronTrigEff_fcal_PbPb18->SetParameter (1, 0);
  g_electronTrigEff_fcal_PbPb18->Fit (f_electronTrigEff_fcal_PbPb18, "RNQ0");
  f_electronTrigEff_fcal_PbPb18->SetLineStyle (2);
  f_electronTrigEff_fcal_PbPb18->SetLineColor (kGreen+2);
  f_electronTrigEff_fcal_PbPb18->SetLineWidth (2);
  f_electronTrigEff_fcal_PbPb18->Draw ("same");

  line->DrawLine (-0.5, 1, 5, 1);

  myText (0.2, 0.35, kBlack, Form ("Slope = %s TeV^{-1}", FormatMeasurement (f_electronTrigEff_fcal_PbPb18->GetParameter (1), f_electronTrigEff_fcal_PbPb18->GetParError (1), 2)), 0.04);
  myText (0.2, 0.30, kBlack, Form ("y-int. = %s", FormatMeasurement (f_electronTrigEff_fcal_PbPb18->GetParameter (0), f_electronTrigEff_fcal_PbPb18->GetParError (0), 2)), 0.04);
  myText (0.2, 0.25, kBlack, Form ("#chi^{2} / dof = %.2f / %i", f_electronTrigEff_fcal_PbPb18->GetChisquare (), f_electronTrigEff_fcal_PbPb18->GetNDF ()), 0.04);

  myText (0.62, 0.4, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.62, 0.35, kBlack, "2018 Pb+Pb, 1.7 nb^{-1}", 0.04);
  myText (0.62, 0.3, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.62, 0.25, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.04);
  myText (0.62, 0.20, kBlack, "LHLoose_HI electrons", 0.04);
  myMarkerTextNoLine (0.55, 0.90, kBlack, kOpenCircle, "HLT_e15_lhloose_ion_L1EM12", 1.25, 0.04);
  c_3->SaveAs ("../Plots/LeptonPerformance/ElectronTrigEffs/PbPb18_cent.pdf");
  SaferDelete (f_electronTrigEff_fcal_PbPb18);



  TH2D* h2_zmumuTrigEffNum_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffNum_pt_y_pp");
  TH2D* h2_zmumuTrigEffDen_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffDen_pt_y_pp");
  TH2D* h2_zmumuTrigEffNum_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffNum_pt_y_PbPb18");
  TH2D* h2_zmumuTrigEffDen_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffDen_pt_y_PbPb18");
  //TH2D* h2_zmumuTrigEffNum_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffNum_pt_y_PbPb15");
  //TH2D* h2_zmumuTrigEffDen_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffDen_pt_y_PbPb15");

  TH2D* h2_zmumuTrigEff_pt_y_pp = (TH2D*) h2_zmumuTrigEffNum_pt_y_pp->Clone ("h2_zmumuTrigEffNum_pt_y_pp");
  h2_zmumuTrigEff_pt_y_pp->Divide (h2_zmumuTrigEffDen_pt_y_pp);
  TH2D* h2_zmumuTrigEff_pt_y_PbPb18 = (TH2D*) h2_zmumuTrigEffNum_pt_y_PbPb18->Clone ("h2_zmumuTrigEffNum_pt_y_PbPb18");
  h2_zmumuTrigEff_pt_y_PbPb18->Divide (h2_zmumuTrigEffDen_pt_y_PbPb18);
  //TH2D* h2_zmumuTrigEff_pt_y_PbPb15 = (TH2D*) h2_zmumuTrigEffNum_pt_y_PbPb15->Clone ("h2_zmumuTrigEffNum_pt_y_PbPb15");
  //h2_zmumuTrigEff_pt_y_PbPb15->Divide (h2_zmumuTrigEffDen_pt_y_PbPb15);

  TH2D* h2_zeeTrigEffNum_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zeeTrigEffNum_pt_y_pp");
  TH2D* h2_zeeTrigEffDen_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zeeTrigEffDen_pt_y_pp");
  TH2D* h2_zeeTrigEffNum_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zeeTrigEffNum_pt_y_PbPb18");
  TH2D* h2_zeeTrigEffDen_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zeeTrigEffDen_pt_y_PbPb18");
  //TH2D* h2_zeeTrigEffNum_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zeeTrigEffNum_pt_y_PbPb15");
  //TH2D* h2_zeeTrigEffDen_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zeeTrigEffDen_pt_y_PbPb15");

  TH2D* h2_zeeTrigEff_pt_y_pp = (TH2D*) h2_zeeTrigEffNum_pt_y_pp->Clone ("h2_zeeTrigEffNum_pt_y_pp");
  h2_zeeTrigEff_pt_y_pp->Divide (h2_zeeTrigEffDen_pt_y_pp);
  TH2D* h2_zeeTrigEff_pt_y_PbPb18 = (TH2D*) h2_zeeTrigEffNum_pt_y_PbPb18->Clone ("h2_zeeTrigEffNum_pt_y_PbPb18");
  h2_zeeTrigEff_pt_y_PbPb18->Divide (h2_zeeTrigEffDen_pt_y_PbPb18);
  //TH2D* h2_zeeTrigEff_pt_y_PbPb15 = (TH2D*) h2_zeeTrigEffNum_pt_y_PbPb15->Clone ("h2_zeeTrigEffNum_pt_y_PbPb15");
  //h2_zeeTrigEff_pt_y_PbPb15->Divide (h2_zeeTrigEffDen_pt_y_PbPb15);


  TCanvas* c_4 = new TCanvas ("c_4", "", 1000, 800);
  FormatTH2Canvas (c_4, false);
  c_4->SetLogy ();

  h2_zmumuTrigEff_pt_y_pp->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zmumuTrigEff_pt_y_pp->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zmumuTrigEff_pt_y_pp->GetZaxis ()->SetTitle ("Z#rightarrow#mu#mu trigger efficiency");
  h2_zmumuTrigEff_pt_y_pp->GetXaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_pp->GetYaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_pp->GetZaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_pp->GetZaxis ()->SetRangeUser (0.75, 1);
  h2_zmumuTrigEff_pt_y_pp->Draw ("lego2");
  myText (0.08, 0.95, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.08, 0.90, kBlack, "2017 #it{pp}, #sqrt{s} = 5.02 TeV, 258 pb^{-1}", 0.04);
  myText (0.08, 0.85, kBlack, "HLT_mu14", 0.04);
  c_4->SaveAs ("../Plots/LeptonPerformance/ZmumuTrigEffs/pp_pt_y.pdf");

  TCanvas* c_5 = new TCanvas ("c_5", "", 1000, 800);
  FormatTH2Canvas (c_5, false);
  c_5->SetLogy ();

  h2_zmumuTrigEff_pt_y_PbPb18->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zmumuTrigEff_pt_y_PbPb18->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zmumuTrigEff_pt_y_PbPb18->GetZaxis ()->SetTitle ("Z#rightarrow#mu#mu trigger efficiency");
  h2_zmumuTrigEff_pt_y_PbPb18->GetXaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_PbPb18->GetYaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_PbPb18->GetZaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_PbPb18->GetZaxis ()->SetRangeUser (0.75, 1);
  h2_zmumuTrigEff_pt_y_PbPb18->Draw ("lego2");
  myText (0.08, 0.95, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.08, 0.90, kBlack, "2018 Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4 nb^{-1}", 0.04);
  myText (0.08, 0.85, kBlack, "HLT_mu14", 0.04);
  c_5->SaveAs ("../Plots/LeptonPerformance/ZmumuTrigEffs/PbPb18_pt_y.pdf");

  TCanvas* c_6 = new TCanvas ("c_6", "", 1000, 800);
  FormatTH2Canvas (c_6, false);
  c_6->SetLogy ();

  h2_zeeTrigEff_pt_y_pp->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zeeTrigEff_pt_y_pp->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zeeTrigEff_pt_y_pp->GetZaxis ()->SetTitle ("Z#rightarrow ee trigger efficiency");
  h2_zeeTrigEff_pt_y_pp->GetXaxis ()->SetLabelSize (0.04);
  h2_zeeTrigEff_pt_y_pp->GetYaxis ()->SetLabelSize (0.04);
  h2_zeeTrigEff_pt_y_pp->GetZaxis ()->SetLabelSize (0.04);
  h2_zeeTrigEff_pt_y_pp->GetZaxis ()->SetRangeUser (0.75, 1);
  h2_zeeTrigEff_pt_y_pp->Draw ("lego2");
  myText (0.08, 0.95, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.08, 0.90, kBlack, "2017 #it{pp}, #sqrt{s} = 5.02 TeV, 258 pb^{-1}", 0.04);
  myText (0.08, 0.85, kBlack, "HLT_e15_lhloose_L1EM12", 0.04);
  c_6->SaveAs ("../Plots/LeptonPerformance/ZeeTrigEffs/pp_pt_y.pdf");

  TCanvas* c_7 = new TCanvas ("c_7", "", 1000, 800);
  FormatTH2Canvas (c_7, false);
  c_7->SetLogy ();

  h2_zeeTrigEff_pt_y_PbPb18->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zeeTrigEff_pt_y_PbPb18->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zeeTrigEff_pt_y_PbPb18->GetZaxis ()->SetTitle ("Z#rightarrow ee trigger efficiency");
  h2_zeeTrigEff_pt_y_PbPb18->GetXaxis ()->SetLabelSize (0.04);
  h2_zeeTrigEff_pt_y_PbPb18->GetYaxis ()->SetLabelSize (0.04);
  h2_zeeTrigEff_pt_y_PbPb18->GetZaxis ()->SetLabelSize (0.04);
  h2_zeeTrigEff_pt_y_PbPb18->GetZaxis ()->SetRangeUser (0.75, 1);
  h2_zeeTrigEff_pt_y_PbPb18->Draw ("lego2");
  myText (0.08, 0.95, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.08, 0.90, kBlack, "2018 Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4 nb^{-1}", 0.04);
  myText (0.08, 0.85, kBlack, "HLT_e15_lhloose_ion_L1EM12", 0.04);
  c_7->SaveAs ("../Plots/LeptonPerformance/ZeeTrigEffs/PbPb18_pt_y.pdf");



  TH2D* h2_zmumuIDEffNum_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zmumuIDEffNum_pt_y_pp");
  TH2D* h2_zmumuIDEffDen_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zmumuIDEffDen_pt_y_pp");
  TH2D* h2_zmumuIDEffNum_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zmumuIDEffNum_pt_y_PbPb18");
  TH2D* h2_zmumuIDEffDen_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zmumuIDEffDen_pt_y_PbPb18");
  //TH2D* h2_zmumuIDEffNum_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zmumuIDEffNum_pt_y_PbPb15");
  //TH2D* h2_zmumuIDEffDen_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zmumuIDEffDen_pt_y_PbPb15");

  TH2D* h2_zmumuIDEff_pt_y_pp = (TH2D*) h2_zmumuIDEffNum_pt_y_pp->Clone ("h2_zmumuIDEffNum_pt_y_pp");
  h2_zmumuIDEff_pt_y_pp->Divide (h2_zmumuIDEffDen_pt_y_pp);
  TH2D* h2_zmumuIDEff_pt_y_PbPb18 = (TH2D*) h2_zmumuIDEffNum_pt_y_PbPb18->Clone ("h2_zmumuIDEffNum_pt_y_PbPb18");
  h2_zmumuIDEff_pt_y_PbPb18->Divide (h2_zmumuIDEffDen_pt_y_PbPb18);
  //TH2D* h2_zmumuIDEff_pt_y_PbPb15 = (TH2D*) h2_zmumuIDEffNum_pt_y_PbPb15->Clone ("h2_zmumuIDEffNum_pt_y_PbPb15");
  //h2_zmumuIDEff_pt_y_PbPb15->Divide (h2_zmumuIDEffDen_pt_y_PbPb15);

  TH2D* h2_zeeIDEffNum_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zeeIDEffNum_pt_y_pp");
  TH2D* h2_zeeIDEffDen_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zeeIDEffDen_pt_y_pp");
  TH2D* h2_zeeIDEffNum_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zeeIDEffNum_pt_y_PbPb18");
  TH2D* h2_zeeIDEffDen_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zeeIDEffDen_pt_y_PbPb18");
  //TH2D* h2_zeeIDEffNum_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zeeIDEffNum_pt_y_PbPb15");
  //TH2D* h2_zeeIDEffDen_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zeeIDEffDen_pt_y_PbPb15");

  TH2D* h2_zeeIDEff_pt_y_pp = (TH2D*) h2_zeeIDEffNum_pt_y_pp->Clone ("h2_zeeIDEffNum_pt_y_pp");
  h2_zeeIDEff_pt_y_pp->Divide (h2_zeeIDEffDen_pt_y_pp);
  TH2D* h2_zeeIDEff_pt_y_PbPb18 = (TH2D*) h2_zeeIDEffNum_pt_y_PbPb18->Clone ("h2_zeeIDEffNum_pt_y_PbPb18");
  h2_zeeIDEff_pt_y_PbPb18->Divide (h2_zeeIDEffDen_pt_y_PbPb18);
  //TH2D* h2_zeeIDEff_pt_y_PbPb15 = (TH2D*) h2_zeeIDEffNum_pt_y_PbPb15->Clone ("h2_zeeIDEffNum_pt_y_PbPb15");
  //h2_zeeIDEff_pt_y_PbPb15->Divide (h2_zeeIDEffDen_pt_y_PbPb15);


  FormatTH2Canvas (c_3, false);
  c_3->cd ();
  TH1D* h_muonIDEff_fcal_PbPb18 = (TH1D*) h_muonIDEffNum_fcal[0]->Clone ("h_muonIDEff_fcal_PbPb18");
  BinomialDivide (h_muonIDEff_fcal_PbPb18, h_muonIDEffDen_fcal[0]);
  //TGAE* g_muonIDEff_fcal_PbPb18 = TEff2TGAE (t_muonIDEff_fcal_PbPb18);
  TGAE* g_muonIDEff_fcal_PbPb18 = make_graph (h_muonIDEff_fcal_PbPb18);
  g_muonIDEff_fcal_PbPb18->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [TeV]");
  g_muonIDEff_fcal_PbPb18->GetYaxis ()->SetTitle ("Muon ID efficiency");
  g_muonIDEff_fcal_PbPb18->GetXaxis ()->SetRangeUser (-0.5, 5);
  g_muonIDEff_fcal_PbPb18->GetYaxis ()->SetRangeUser (0.0, 1.10);
  g_muonIDEff_fcal_PbPb18->SetMarkerStyle (kOpenCircle);
  g_muonIDEff_fcal_PbPb18->Draw ("AP");

  TF1* f_muonIDEff_fcal_PbPb18 = new TF1 ("f_muonIDEff_fcal_PbPb18", "[0]+[1]*x", 0, 5);
  f_muonIDEff_fcal_PbPb18->SetParameter (0, 1);
  f_muonIDEff_fcal_PbPb18->SetParameter (1, 0);
  g_muonIDEff_fcal_PbPb18->Fit (f_muonIDEff_fcal_PbPb18, "RNQ0");
  f_muonIDEff_fcal_PbPb18->SetLineStyle (2);
  f_muonIDEff_fcal_PbPb18->SetLineColor (kGreen+2);
  f_muonIDEff_fcal_PbPb18->SetLineWidth (2);
  f_muonIDEff_fcal_PbPb18->Draw ("same");

  line->DrawLine (-0.5, 1, 5, 1);

  myText (0.2, 0.35, kBlack, Form ("Slope = %s TeV^{-1}", FormatMeasurement (f_muonIDEff_fcal_PbPb18->GetParameter (1), f_muonIDEff_fcal_PbPb18->GetParError (1), 2)), 0.04);
  myText (0.2, 0.30, kBlack, Form ("y-int. = %s", FormatMeasurement (f_muonIDEff_fcal_PbPb18->GetParameter (0), f_muonIDEff_fcal_PbPb18->GetParError (0), 2)), 0.04);
  myText (0.2, 0.25, kBlack, Form ("#chi^{2} / dof = %.2f / %i", f_muonIDEff_fcal_PbPb18->GetChisquare (), f_muonIDEff_fcal_PbPb18->GetNDF ()), 0.04);

  myText (0.62, 0.4, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.62, 0.35, kBlack, "2018 Pb+Pb, 1.4 nb^{-1}", 0.04);
  myText (0.62, 0.3, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.62, 0.25, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.04);
  myText (0.62, 0.20, kBlack, "Medium muons", 0.04);
  c_3->SaveAs ("../Plots/LeptonPerformance/MuonIDEffs/PbPb18_cent.pdf");
  SaferDelete (f_muonIDEff_fcal_PbPb18);


  FormatTH2Canvas (c_3, false);
  c_3->cd ();
  TH1D* h_electronIDEff_fcal_PbPb18 = (TH1D*) h_electronIDEffNum_fcal[0]->Clone ("h_electronIDEff_fcal_PbPb18");
  BinomialDivide (h_electronIDEff_fcal_PbPb18, h_electronIDEffDen_fcal[0]);
  //TGAE* g_electronIDEff_fcal_PbPb18 = TEff2TGAE (t_electronIDEff_fcal_PbPb18);
  TGAE* g_electronIDEff_fcal_PbPb18 = make_graph (h_electronIDEff_fcal_PbPb18);
  g_electronIDEff_fcal_PbPb18->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [TeV]");
  g_electronIDEff_fcal_PbPb18->GetYaxis ()->SetTitle ("Electron ID efficiency");
  g_electronIDEff_fcal_PbPb18->GetXaxis ()->SetRangeUser (-0.5, 5);
  g_electronIDEff_fcal_PbPb18->GetYaxis ()->SetRangeUser (0.0, 1.10);
  g_electronIDEff_fcal_PbPb18->SetMarkerStyle (kOpenCircle);
  g_electronIDEff_fcal_PbPb18->Draw ("AP");

  line->DrawLine (-0.5, 1, 5, 1);

  //TF1* f_electronIDEff_fcal_PbPb18 = new TF1 ("f_electronIDEff_fcal_PbPb18", "[0]+[1]*x+[2]*x^2+[3]*x^3", 0, 5);
  //f_electronIDEff_fcal_PbPb18->SetParameter (0, 1);
  //f_electronIDEff_fcal_PbPb18->SetParameter (1, 0);
  //f_electronIDEff_fcal_PbPb18->SetParameter (2, 0);
  //f_electronIDEff_fcal_PbPb18->SetParameter (3, 0);
  //g_electronIDEff_fcal_PbPb18->Fit (f_electronIDEff_fcal_PbPb18, "RNQ0");
  //f_electronIDEff_fcal_PbPb18->SetLineStyle (2);
  //f_electronIDEff_fcal_PbPb18->SetLineColor (kGreen+2);
  //f_electronIDEff_fcal_PbPb18->SetLineWidth (2);
  //f_electronIDEff_fcal_PbPb18->Draw ("same");

  //myText (0.2, 0.35, kBlack, Form ("Slope = %s TeV^{-1}", FormatMeasurement (f_electronIDEff_fcal_PbPb18->GetParameter (1), f_electronIDEff_fcal_PbPb18->GetParError (1), 2)), 0.04);
  //myText (0.2, 0.30, kBlack, Form ("y-int. = %s", FormatMeasurement (f_electronIDEff_fcal_PbPb18->GetParameter (0), f_electronIDEff_fcal_PbPb18->GetParError (0), 2)), 0.04);
  //myText (0.2, 0.25, kBlack, Form ("#chi^{2} / dof = %.2f / %i", f_electronIDEff_fcal_PbPb18->GetChisquare (), f_electronIDEff_fcal_PbPb18->GetNDF ()), 0.04);

  myText (0.62, 0.4, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.62, 0.35, kBlack, "2018 Pb+Pb, 1.7 nb^{-1}", 0.04);
  myText (0.62, 0.3, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.62, 0.25, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.04);
  myText (0.62, 0.20, kBlack, "LHLoose_HI electrons", 0.04);
  c_3->SaveAs ("../Plots/LeptonPerformance/ElectronIDEffs/PbPb18_cent.pdf");
  //SaferDelete (f_electronIDEff_fcal_PbPb18);


  TCanvas* c_8 = new TCanvas ("c_8", "", 1000, 800);
  FormatTH2Canvas (c_8, false);
  c_8->SetLogy ();

  h2_zmumuIDEff_pt_y_pp->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zmumuIDEff_pt_y_pp->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zmumuIDEff_pt_y_pp->GetZaxis ()->SetTitle ("Z#rightarrow#mu#mu ID efficiency");
  h2_zmumuIDEff_pt_y_pp->GetXaxis ()->SetLabelSize (0.04);
  h2_zmumuIDEff_pt_y_pp->GetYaxis ()->SetLabelSize (0.04);
  h2_zmumuIDEff_pt_y_pp->GetZaxis ()->SetLabelSize (0.04);
  h2_zmumuIDEff_pt_y_pp->GetZaxis ()->SetRangeUser (0., 1);
  h2_zmumuIDEff_pt_y_pp->Draw ("lego2");
  myText (0.08, 0.95, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.08, 0.90, kBlack, "2017 #it{pp}, #sqrt{s} = 5.02 TeV, 258 pb^{-1}", 0.04);
  myText (0.08, 0.85, kBlack, "Medium muons", 0.04);
  c_8->SaveAs ("../Plots/LeptonPerformance/ZmumuIDEffs/pp_pt_y.pdf");

  TCanvas* c_9 = new TCanvas ("c_9", "", 1000, 800);
  FormatTH2Canvas (c_9, false);
  c_9->SetLogy ();

  h2_zmumuIDEff_pt_y_PbPb18->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zmumuIDEff_pt_y_PbPb18->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zmumuIDEff_pt_y_PbPb18->GetZaxis ()->SetTitle ("Z#rightarrow#mu#mu ID efficiency");
  h2_zmumuIDEff_pt_y_PbPb18->GetXaxis ()->SetLabelSize (0.04);
  h2_zmumuIDEff_pt_y_PbPb18->GetYaxis ()->SetLabelSize (0.04);
  h2_zmumuIDEff_pt_y_PbPb18->GetZaxis ()->SetLabelSize (0.04);
  h2_zmumuIDEff_pt_y_PbPb18->GetZaxis ()->SetRangeUser (0., 1);
  h2_zmumuIDEff_pt_y_PbPb18->Draw ("lego2");
  myText (0.08, 0.95, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.08, 0.90, kBlack, "2018 Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4 nb^{-1}", 0.04);
  myText (0.08, 0.85, kBlack, "Medium muons", 0.04);
  c_9->SaveAs ("../Plots/LeptonPerformance/ZmumuIDEffs/PbPb18_pt_y.pdf");

  TCanvas* c_10 = new TCanvas ("c_10", "", 1000, 800);
  FormatTH2Canvas (c_10, false);
  c_10->SetLogy ();

  h2_zeeIDEff_pt_y_pp->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zeeIDEff_pt_y_pp->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zeeIDEff_pt_y_pp->GetZaxis ()->SetTitle ("Z#rightarrow ee ID efficiency");
  h2_zeeIDEff_pt_y_pp->GetXaxis ()->SetLabelSize (0.04);
  h2_zeeIDEff_pt_y_pp->GetYaxis ()->SetLabelSize (0.04);
  h2_zeeIDEff_pt_y_pp->GetZaxis ()->SetLabelSize (0.04);
  h2_zeeIDEff_pt_y_pp->GetZaxis ()->SetRangeUser (0., 1);
  h2_zeeIDEff_pt_y_pp->Draw ("lego2");
  myText (0.08, 0.95, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.08, 0.90, kBlack, "2017 #it{pp}, #sqrt{s} = 5.02 TeV, 258 pb^{-1}", 0.04);
  myText (0.08, 0.85, kBlack, "LHLoose electrons", 0.04);
  c_10->SaveAs ("../Plots/LeptonPerformance/ZeeIDEffs/pp_pt_y.pdf");

  TCanvas* c_11 = new TCanvas ("c_11", "", 1000, 800);
  FormatTH2Canvas (c_11, false);
  c_11->SetLogy ();

  h2_zeeIDEff_pt_y_PbPb18->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zeeIDEff_pt_y_PbPb18->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zeeIDEff_pt_y_PbPb18->GetZaxis ()->SetTitle ("Z#rightarrow ee ID efficiency");
  h2_zeeIDEff_pt_y_PbPb18->GetXaxis ()->SetLabelSize (0.04);
  h2_zeeIDEff_pt_y_PbPb18->GetYaxis ()->SetLabelSize (0.04);
  h2_zeeIDEff_pt_y_PbPb18->GetZaxis ()->SetLabelSize (0.04);
  h2_zeeIDEff_pt_y_PbPb18->GetZaxis ()->SetRangeUser (0., 1);
  h2_zeeIDEff_pt_y_PbPb18->Draw ("lego2");
  myText (0.08, 0.95, kBlack, "#bf{#it{ATLAS}} Internal", 0.045);
  myText (0.08, 0.90, kBlack, "2018 Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV, 1.4 nb^{-1}", 0.04);
  myText (0.08, 0.85, kBlack, "LHLoose_HI electrons", 0.04);
  c_11->SaveAs ("../Plots/LeptonPerformance/ZeeIDEffs/PbPb18_pt_y.pdf");

}

#endif
