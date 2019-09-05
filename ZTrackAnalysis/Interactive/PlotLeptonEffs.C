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

TEfficiency* t_muonMed_IDEff_pt_pp[3];
TEfficiency* t_muonMed_IDEff_pt_PbPb18[3];
TEfficiency* t_muonMed_IDEff_pt_PbPb15[3];
TEfficiency* t_muonMed_IDEff_eta_pp[3];
TEfficiency* t_muonMed_IDEff_eta_PbPb18[3];
TEfficiency* t_muonMed_IDEff_eta_PbPb15[3];
TEfficiency* t_muonMed_IDEff_fcal_PbPb18[3];
TEfficiency* t_muonMed_IDEff_fcal_PbPb15[3];

TEfficiency* t_muonID_MSEff_pt_pp[3];
TEfficiency* t_muonID_MSEff_pt_PbPb18[3];
TEfficiency* t_muonID_MSEff_pt_PbPb15[3];
TEfficiency* t_muonID_MSEff_eta_pp[3];
TEfficiency* t_muonID_MSEff_eta_PbPb18[3];
TEfficiency* t_muonID_MSEff_eta_PbPb15[3];
TEfficiency* t_muonID_MSEff_fcal_PbPb18[3];
TEfficiency* t_muonID_MSEff_fcal_PbPb15[3];

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

  h_muonTrigEffNum_pt[0] = (TH1D*) inFile->Get ("h_muonTrigEffNum_pt_pp");
  h_muonTrigEffDen_pt[0] = (TH1D*) inFile->Get ("h_muonTrigEffDen_pt_pp");
  h_muonTrigEffNum_pt[1] = (TH1D*) inFile->Get ("h_muonTrigEffNum_pt_PbPb18");
  h_muonTrigEffDen_pt[1] = (TH1D*) inFile->Get ("h_muonTrigEffDen_pt_PbPb18");
  h_muonTrigEffNum_pt[2] = (TH1D*) inFile->Get ("h_muonTrigEffNum_pt_PbPb15");
  h_muonTrigEffDen_pt[2] = (TH1D*) inFile->Get ("h_muonTrigEffDen_pt_PbPb15");

  t_muonTrigEff_pt_pp = new TEfficiency (*(h_muonTrigEffNum_pt[0]), *(h_muonTrigEffDen_pt[0]));
  t_muonTrigEff_pt_PbPb18 = new TEfficiency (*(h_muonTrigEffNum_pt[1]), *(h_muonTrigEffDen_pt[1]));
  t_muonTrigEff_pt_PbPb15 = new TEfficiency (*(h_muonTrigEffNum_pt[2]), *(h_muonTrigEffDen_pt[2]));

  h2_muonTrigEffNum_eta_phi[0] = (TH2D*) inFile->Get ("h2_muonTrigEffNum_eta_phi_pp");
  h2_muonTrigEffDen_eta_phi[0] = (TH2D*) inFile->Get ("h2_muonTrigEffDen_eta_phi_pp");
  h2_muonTrigEffNum_eta_phi[1] = (TH2D*) inFile->Get ("h2_muonTrigEffNum_eta_phi_PbPb18");
  h2_muonTrigEffDen_eta_phi[1] = (TH2D*) inFile->Get ("h2_muonTrigEffDen_eta_phi_PbPb18");
  h2_muonTrigEffNum_eta_phi[2] = (TH2D*) inFile->Get ("h2_muonTrigEffNum_eta_phi_PbPb15");
  h2_muonTrigEffDen_eta_phi[2] = (TH2D*) inFile->Get ("h2_muonTrigEffDen_eta_phi_PbPb15");

  h_muonTrigEffNum_eta[0] = (TH1D*) h2_muonTrigEffNum_eta_phi[0]->ProjectionX ("h_muonTrigEffNum_eta_pp");
  h_muonTrigEffDen_eta[0] = (TH1D*) h2_muonTrigEffDen_eta_phi[0]->ProjectionX ("h_muonTrigEffDen_eta_pp");
  h_muonTrigEffNum_eta[1] = (TH1D*) h2_muonTrigEffNum_eta_phi[1]->ProjectionX ("h_muonTrigEffNum_eta_PbPb18");
  h_muonTrigEffDen_eta[1] = (TH1D*) h2_muonTrigEffDen_eta_phi[1]->ProjectionX ("h_muonTrigEffDen_eta_PbPb18");
  h_muonTrigEffNum_eta[2] = (TH1D*) h2_muonTrigEffNum_eta_phi[2]->ProjectionX ("h_muonTrigEffNum_eta_PbPb15");
  h_muonTrigEffDen_eta[2] = (TH1D*) h2_muonTrigEffDen_eta_phi[2]->ProjectionX ("h_muonTrigEffDen_eta_PbPb15");

  t_muonTrigEff_eta_phi_pp = new TEfficiency (*(h2_muonTrigEffNum_eta_phi[0]), *(h2_muonTrigEffDen_eta_phi[0]));
  t_muonTrigEff_eta_phi_PbPb18 = new TEfficiency (*(h2_muonTrigEffNum_eta_phi[1]), *(h2_muonTrigEffDen_eta_phi[1]));
  t_muonTrigEff_eta_phi_PbPb15 = new TEfficiency (*(h2_muonTrigEffNum_eta_phi[2]), *(h2_muonTrigEffDen_eta_phi[2]));
  t_muonTrigEff_eta_pp = new TEfficiency (*(h_muonTrigEffNum_eta[0]), *(h_muonTrigEffDen_eta[0]));
  t_muonTrigEff_eta_PbPb18 = new TEfficiency (*(h_muonTrigEffNum_eta[1]), *(h_muonTrigEffDen_eta[1]));
  t_muonTrigEff_eta_PbPb15 = new TEfficiency (*(h_muonTrigEffNum_eta[2]), *(h_muonTrigEffDen_eta[2]));

  h_muonTrigEffNum_fcal[0] = (TH1D*) inFile->Get ("h_muonTrigEffNum_fcal_PbPb18");
  h_muonTrigEffDen_fcal[0] = (TH1D*) inFile->Get ("h_muonTrigEffDen_fcal_PbPb18");
  h_muonTrigEffNum_fcal[1] = (TH1D*) inFile->Get ("h_muonTrigEffNum_fcal_PbPb15");
  h_muonTrigEffDen_fcal[1] = (TH1D*) inFile->Get ("h_muonTrigEffDen_fcal_PbPb15");
  t_muonTrigEff_fcal_PbPb18 = new TEfficiency (*(h_muonTrigEffNum_fcal[0]), *(h_muonTrigEffDen_fcal[0]));
  t_muonTrigEff_fcal_PbPb15 = new TEfficiency (*(h_muonTrigEffNum_fcal[1]), *(h_muonTrigEffDen_fcal[1]));

  TH1D* h_muonMed_IDEffNum_pt[3][3];
  TH1D* h_muonMed_IDEffDen_pt[3][3];
  TH1D* h_muonMed_IDEffNum_eta[3][3];
  TH1D* h_muonMed_IDEffDen_eta[3][3];
  TH1D* h_muonMed_IDEffNum_fcal[2][3];
  TH1D* h_muonMed_IDEffDen_fcal[2][3];

  TH1D* h_muonID_MSEffNum_pt[3][3];
  TH1D* h_muonID_MSEffDen_pt[3][3];
  TH1D* h_muonID_MSEffNum_eta[3][3];
  TH1D* h_muonID_MSEffDen_eta[3][3];
  TH1D* h_muonID_MSEffNum_fcal[2][3];
  TH1D* h_muonID_MSEffDen_fcal[2][3];

  for (int iSign : {0, 1}) {
    h_muonMed_IDEffNum_pt[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_pt_pp_sign%i", iSign));
    h_muonMed_IDEffDen_pt[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_pt_pp_sign%i", iSign));
    h_muonMed_IDEffNum_pt[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_pt_PbPb18_sign%i", iSign));
    h_muonMed_IDEffDen_pt[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_pt_PbPb18_sign%i", iSign));
    h_muonMed_IDEffNum_pt[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_pt_PbPb15_sign%i", iSign));
    h_muonMed_IDEffDen_pt[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_pt_PbPb15_sign%i", iSign));

    t_muonMed_IDEff_pt_pp[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_pt[0][iSign]), *(h_muonMed_IDEffDen_pt[0][iSign]));
    t_muonMed_IDEff_pt_PbPb18[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_pt[1][iSign]), *(h_muonMed_IDEffDen_pt[1][iSign]));
    t_muonMed_IDEff_pt_PbPb15[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_pt[2][iSign]), *(h_muonMed_IDEffDen_pt[2][iSign]));

    h_muonMed_IDEffNum_eta[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_eta_pp_sign%i", iSign));
    h_muonMed_IDEffDen_eta[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_eta_pp_sign%i", iSign));
    h_muonMed_IDEffNum_eta[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_eta_PbPb18_sign%i", iSign));
    h_muonMed_IDEffDen_eta[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_eta_PbPb18_sign%i", iSign));
    h_muonMed_IDEffNum_eta[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_eta_PbPb15_sign%i", iSign));
    h_muonMed_IDEffDen_eta[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_eta_PbPb15_sign%i", iSign));

    t_muonMed_IDEff_eta_pp[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_eta[0][iSign]), *(h_muonMed_IDEffDen_eta[0][iSign]));
    t_muonMed_IDEff_eta_PbPb18[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_eta[1][iSign]), *(h_muonMed_IDEffDen_eta[1][iSign]));
    t_muonMed_IDEff_eta_PbPb15[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_eta[2][iSign]), *(h_muonMed_IDEffDen_eta[2][iSign]));

    h_muonMed_IDEffNum_fcal[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_fcal_PbPb18_sign%i", iSign));
    h_muonMed_IDEffDen_fcal[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_fcal_PbPb18_sign%i", iSign));
    h_muonMed_IDEffNum_fcal[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffNum_fcal_PbPb15_sign%i", iSign));
    h_muonMed_IDEffDen_fcal[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonMed_IDEffDen_fcal_PbPb15_sign%i", iSign));
    t_muonMed_IDEff_fcal_PbPb18[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_fcal[0][iSign]), *(h_muonMed_IDEffDen_fcal[0][iSign]));
    t_muonMed_IDEff_fcal_PbPb15[iSign] = new TEfficiency (*(h_muonMed_IDEffNum_fcal[1][iSign]), *(h_muonMed_IDEffDen_fcal[1][iSign]));


    h_muonID_MSEffNum_pt[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_pt_pp_sign%i", iSign));
    h_muonID_MSEffDen_pt[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_pt_pp_sign%i", iSign));
    h_muonID_MSEffNum_pt[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_pt_PbPb18_sign%i", iSign));
    h_muonID_MSEffDen_pt[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_pt_PbPb18_sign%i", iSign));
    h_muonID_MSEffNum_pt[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_pt_PbPb15_sign%i", iSign));
    h_muonID_MSEffDen_pt[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_pt_PbPb15_sign%i", iSign));

    t_muonID_MSEff_pt_pp[iSign] = new TEfficiency (*(h_muonID_MSEffNum_pt[0][iSign]), *(h_muonID_MSEffDen_pt[0][iSign]));
    t_muonID_MSEff_pt_PbPb18[iSign] = new TEfficiency (*(h_muonID_MSEffNum_pt[1][iSign]), *(h_muonID_MSEffDen_pt[1][iSign]));
    t_muonID_MSEff_pt_PbPb15[iSign] = new TEfficiency (*(h_muonID_MSEffNum_pt[2][iSign]), *(h_muonID_MSEffDen_pt[2][iSign]));

    h_muonID_MSEffNum_eta[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_eta_pp_sign%i", iSign));
    h_muonID_MSEffDen_eta[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_eta_pp_sign%i", iSign));
    h_muonID_MSEffNum_eta[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_eta_PbPb18_sign%i", iSign));
    h_muonID_MSEffDen_eta[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_eta_PbPb18_sign%i", iSign));
    h_muonID_MSEffNum_eta[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_eta_PbPb15_sign%i", iSign));
    h_muonID_MSEffDen_eta[2][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_eta_PbPb15_sign%i", iSign));

    t_muonID_MSEff_eta_pp[iSign] = new TEfficiency (*(h_muonID_MSEffNum_eta[0][iSign]), *(h_muonID_MSEffDen_eta[0][iSign]));
    t_muonID_MSEff_eta_PbPb18[iSign] = new TEfficiency (*(h_muonID_MSEffNum_eta[1][iSign]), *(h_muonID_MSEffDen_eta[1][iSign]));
    t_muonID_MSEff_eta_PbPb15[iSign] = new TEfficiency (*(h_muonID_MSEffNum_eta[2][iSign]), *(h_muonID_MSEffDen_eta[2][iSign]));

    h_muonID_MSEffNum_fcal[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_fcal_PbPb18_sign%i", iSign));
    h_muonID_MSEffDen_fcal[0][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_fcal_PbPb18_sign%i", iSign));
    h_muonID_MSEffNum_fcal[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffNum_fcal_PbPb15_sign%i", iSign));
    h_muonID_MSEffDen_fcal[1][iSign] = (TH1D*) inFile->Get (Form ("h_muonID_MSEffDen_fcal_PbPb15_sign%i", iSign));
    t_muonID_MSEff_fcal_PbPb18[iSign] = new TEfficiency (*(h_muonID_MSEffNum_fcal[0][iSign]), *(h_muonID_MSEffDen_fcal[0][iSign]));
    t_muonID_MSEff_fcal_PbPb15[iSign] = new TEfficiency (*(h_muonID_MSEffNum_fcal[1][iSign]), *(h_muonID_MSEffDen_fcal[1][iSign]));
  }

  {
    h_muonMed_IDEffNum_pt[0][2] = (TH1D*) h_muonMed_IDEffNum_pt[0][1]->Clone ("h_muonMed_IDEffNum_pt_pp_sign2");
    h_muonMed_IDEffDen_pt[0][2] = (TH1D*) h_muonMed_IDEffDen_pt[0][1]->Clone ("h_muonMed_IDEffDen_pt_pp_sign2");
    h_muonMed_IDEffNum_pt[1][2] = (TH1D*) h_muonMed_IDEffNum_pt[1][1]->Clone ("h_muonMed_IDEffNum_pt_PbPb18_sign2");
    h_muonMed_IDEffDen_pt[1][2] = (TH1D*) h_muonMed_IDEffDen_pt[1][1]->Clone ("h_muonMed_IDEffDen_pt_PbPb18_sign2");
    h_muonMed_IDEffNum_pt[2][2] = (TH1D*) h_muonMed_IDEffNum_pt[2][1]->Clone ("h_muonMed_IDEffNum_pt_PbPb15_sign2");
    h_muonMed_IDEffDen_pt[2][2] = (TH1D*) h_muonMed_IDEffDen_pt[2][1]->Clone ("h_muonMed_IDEffDen_pt_PbPb15_sign2");
    h_muonMed_IDEffNum_pt[0][2]->Add (h_muonMed_IDEffNum_pt[0][0], -1);
    h_muonMed_IDEffDen_pt[0][2]->Add (h_muonMed_IDEffDen_pt[0][0], -1);
    h_muonMed_IDEffNum_pt[1][2]->Add (h_muonMed_IDEffNum_pt[1][0], -1);
    h_muonMed_IDEffDen_pt[1][2]->Add (h_muonMed_IDEffDen_pt[1][0], -1);
    h_muonMed_IDEffNum_pt[2][2]->Add (h_muonMed_IDEffNum_pt[2][0], -1);
    h_muonMed_IDEffDen_pt[2][2]->Add (h_muonMed_IDEffDen_pt[2][0], -1);

    t_muonMed_IDEff_pt_pp[2] = new TEfficiency (*(h_muonMed_IDEffNum_pt[0][2]), *(h_muonMed_IDEffDen_pt[0][2]));
    t_muonMed_IDEff_pt_PbPb18[2] = new TEfficiency (*(h_muonMed_IDEffNum_pt[1][2]), *(h_muonMed_IDEffDen_pt[1][2]));
    t_muonMed_IDEff_pt_PbPb15[2] = new TEfficiency (*(h_muonMed_IDEffNum_pt[2][2]), *(h_muonMed_IDEffDen_pt[2][2]));

    h_muonMed_IDEffNum_eta[0][2] = (TH1D*) h_muonMed_IDEffNum_eta[0][1]->Clone ("h_muonMed_IDEffNum_eta_pp_sign2");
    h_muonMed_IDEffDen_eta[0][2] = (TH1D*) h_muonMed_IDEffDen_eta[0][1]->Clone ("h_muonMed_IDEffDen_eta_pp_sign2");
    h_muonMed_IDEffNum_eta[1][2] = (TH1D*) h_muonMed_IDEffNum_eta[1][1]->Clone ("h_muonMed_IDEffNum_eta_PbPb18_sign2");
    h_muonMed_IDEffDen_eta[1][2] = (TH1D*) h_muonMed_IDEffDen_eta[1][1]->Clone ("h_muonMed_IDEffDen_eta_PbPb18_sign2");
    h_muonMed_IDEffNum_eta[2][2] = (TH1D*) h_muonMed_IDEffNum_eta[2][1]->Clone ("h_muonMed_IDEffNum_eta_PbPb15_sign2");
    h_muonMed_IDEffDen_eta[2][2] = (TH1D*) h_muonMed_IDEffDen_eta[2][1]->Clone ("h_muonMed_IDEffDen_eta_PbPb15_sign2");
    h_muonMed_IDEffNum_eta[0][2]->Add (h_muonMed_IDEffNum_eta[0][0], -1);
    h_muonMed_IDEffDen_eta[0][2]->Add (h_muonMed_IDEffDen_eta[0][0], -1);
    h_muonMed_IDEffNum_eta[1][2]->Add (h_muonMed_IDEffNum_eta[1][0], -1);
    h_muonMed_IDEffDen_eta[1][2]->Add (h_muonMed_IDEffDen_eta[1][0], -1);
    h_muonMed_IDEffNum_eta[2][2]->Add (h_muonMed_IDEffNum_eta[2][0], -1);
    h_muonMed_IDEffDen_eta[2][2]->Add (h_muonMed_IDEffDen_eta[2][0], -1);

    t_muonMed_IDEff_eta_pp[2] = new TEfficiency (*(h_muonMed_IDEffNum_eta[0][2]), *(h_muonMed_IDEffDen_eta[0][2]));
    t_muonMed_IDEff_eta_PbPb18[2] = new TEfficiency (*(h_muonMed_IDEffNum_eta[1][2]), *(h_muonMed_IDEffDen_eta[1][2]));
    t_muonMed_IDEff_eta_PbPb15[2] = new TEfficiency (*(h_muonMed_IDEffNum_eta[2][2]), *(h_muonMed_IDEffDen_eta[2][2]));

    h_muonMed_IDEffNum_fcal[0][2] = (TH1D*) h_muonMed_IDEffNum_fcal[0][1]->Clone ("h_muonMed_IDEffNum_fcal_PbPb18_sign2");
    h_muonMed_IDEffDen_fcal[0][2] = (TH1D*) h_muonMed_IDEffDen_fcal[0][1]->Clone ("h_muonMed_IDEffDen_fcal_PbPb18_sign2");
    h_muonMed_IDEffNum_fcal[1][2] = (TH1D*) h_muonMed_IDEffNum_fcal[1][1]->Clone ("h_muonMed_IDEffNum_fcal_PbPb15_sign2");
    h_muonMed_IDEffDen_fcal[1][2] = (TH1D*) h_muonMed_IDEffDen_fcal[1][1]->Clone ("h_muonMed_IDEffDen_fcal_PbPb15_sign2");
    h_muonMed_IDEffNum_fcal[0][2]->Add (h_muonMed_IDEffNum_fcal[0][0], -1);
    h_muonMed_IDEffDen_fcal[0][2]->Add (h_muonMed_IDEffDen_fcal[0][0], -1);
    h_muonMed_IDEffNum_fcal[1][2]->Add (h_muonMed_IDEffNum_fcal[1][0], -1);
    h_muonMed_IDEffDen_fcal[1][2]->Add (h_muonMed_IDEffDen_fcal[1][0], -1);

    t_muonMed_IDEff_fcal_PbPb18[2] = new TEfficiency (*(h_muonMed_IDEffNum_fcal[0][2]), *(h_muonMed_IDEffDen_fcal[0][2]));
    t_muonMed_IDEff_fcal_PbPb15[2] = new TEfficiency (*(h_muonMed_IDEffNum_fcal[1][2]), *(h_muonMed_IDEffDen_fcal[1][2]));


    h_muonID_MSEffNum_pt[0][2] = (TH1D*) h_muonID_MSEffNum_pt[0][1]->Clone ("h_muonID_MSEffNum_pt_pp_sign2");
    h_muonID_MSEffDen_pt[0][2] = (TH1D*) h_muonID_MSEffDen_pt[0][1]->Clone ("h_muonID_MSEffDen_pt_pp_sign2");
    h_muonID_MSEffNum_pt[1][2] = (TH1D*) h_muonID_MSEffNum_pt[1][1]->Clone ("h_muonID_MSEffNum_pt_PbPb18_sign2");
    h_muonID_MSEffDen_pt[1][2] = (TH1D*) h_muonID_MSEffDen_pt[1][1]->Clone ("h_muonID_MSEffDen_pt_PbPb18_sign2");
    h_muonID_MSEffNum_pt[2][2] = (TH1D*) h_muonID_MSEffNum_pt[2][1]->Clone ("h_muonID_MSEffNum_pt_PbPb15_sign2");
    h_muonID_MSEffDen_pt[2][2] = (TH1D*) h_muonID_MSEffDen_pt[2][1]->Clone ("h_muonID_MSEffDen_pt_PbPb15_sign2");
    h_muonID_MSEffNum_pt[0][2]->Add (h_muonID_MSEffNum_pt[0][0], -1);
    h_muonID_MSEffDen_pt[0][2]->Add (h_muonID_MSEffDen_pt[0][0], -1);
    h_muonID_MSEffNum_pt[1][2]->Add (h_muonID_MSEffNum_pt[1][0], -1);
    h_muonID_MSEffDen_pt[1][2]->Add (h_muonID_MSEffDen_pt[1][0], -1);
    h_muonID_MSEffNum_pt[2][2]->Add (h_muonID_MSEffNum_pt[2][0], -1);
    h_muonID_MSEffDen_pt[2][2]->Add (h_muonID_MSEffDen_pt[2][0], -1);

    t_muonID_MSEff_pt_pp[2] = new TEfficiency (*(h_muonID_MSEffNum_pt[0][2]), *(h_muonID_MSEffDen_pt[0][2]));
    t_muonID_MSEff_pt_PbPb18[2] = new TEfficiency (*(h_muonID_MSEffNum_pt[1][2]), *(h_muonID_MSEffDen_pt[1][2]));
    t_muonID_MSEff_pt_PbPb15[2] = new TEfficiency (*(h_muonID_MSEffNum_pt[2][2]), *(h_muonID_MSEffDen_pt[2][2]));

    h_muonID_MSEffNum_eta[0][2] = (TH1D*) h_muonID_MSEffNum_eta[0][1]->Clone ("h_muonID_MSEffNum_eta_pp_sign2");
    h_muonID_MSEffDen_eta[0][2] = (TH1D*) h_muonID_MSEffDen_eta[0][1]->Clone ("h_muonID_MSEffDen_eta_pp_sign2");
    h_muonID_MSEffNum_eta[1][2] = (TH1D*) h_muonID_MSEffNum_eta[1][1]->Clone ("h_muonID_MSEffNum_eta_PbPb18_sign2");
    h_muonID_MSEffDen_eta[1][2] = (TH1D*) h_muonID_MSEffDen_eta[1][1]->Clone ("h_muonID_MSEffDen_eta_PbPb18_sign2");
    h_muonID_MSEffNum_eta[2][2] = (TH1D*) h_muonID_MSEffNum_eta[2][1]->Clone ("h_muonID_MSEffNum_eta_PbPb15_sign2");
    h_muonID_MSEffDen_eta[2][2] = (TH1D*) h_muonID_MSEffDen_eta[2][1]->Clone ("h_muonID_MSEffDen_eta_PbPb15_sign2");
    h_muonID_MSEffNum_eta[0][2]->Add (h_muonID_MSEffNum_eta[0][0], -1);
    h_muonID_MSEffDen_eta[0][2]->Add (h_muonID_MSEffDen_eta[0][0], -1);
    h_muonID_MSEffNum_eta[1][2]->Add (h_muonID_MSEffNum_eta[1][0], -1);
    h_muonID_MSEffDen_eta[1][2]->Add (h_muonID_MSEffDen_eta[1][0], -1);
    h_muonID_MSEffNum_eta[2][2]->Add (h_muonID_MSEffNum_eta[2][0], -1);
    h_muonID_MSEffDen_eta[2][2]->Add (h_muonID_MSEffDen_eta[2][0], -1);

    t_muonID_MSEff_eta_pp[2] = new TEfficiency (*(h_muonID_MSEffNum_eta[0][2]), *(h_muonID_MSEffDen_eta[0][2]));
    t_muonID_MSEff_eta_PbPb18[2] = new TEfficiency (*(h_muonID_MSEffNum_eta[1][2]), *(h_muonID_MSEffDen_eta[1][2]));
    t_muonID_MSEff_eta_PbPb15[2] = new TEfficiency (*(h_muonID_MSEffNum_eta[2][2]), *(h_muonID_MSEffDen_eta[2][2]));

    h_muonID_MSEffNum_fcal[0][2] = (TH1D*) h_muonID_MSEffNum_fcal[0][1]->Clone ("h_muonID_MSEffNum_fcal_PbPb18_sign2");
    h_muonID_MSEffDen_fcal[0][2] = (TH1D*) h_muonID_MSEffDen_fcal[0][1]->Clone ("h_muonID_MSEffDen_fcal_PbPb18_sign2");
    h_muonID_MSEffNum_fcal[1][2] = (TH1D*) h_muonID_MSEffNum_fcal[1][1]->Clone ("h_muonID_MSEffNum_fcal_PbPb15_sign2");
    h_muonID_MSEffDen_fcal[1][2] = (TH1D*) h_muonID_MSEffDen_fcal[1][1]->Clone ("h_muonID_MSEffDen_fcal_PbPb15_sign2");
    h_muonID_MSEffNum_fcal[0][2]->Add (h_muonID_MSEffNum_fcal[0][0], -1);
    h_muonID_MSEffDen_fcal[0][2]->Add (h_muonID_MSEffDen_fcal[0][0], -1);
    h_muonID_MSEffNum_fcal[1][2]->Add (h_muonID_MSEffNum_fcal[1][0], -1);
    h_muonID_MSEffDen_fcal[1][2]->Add (h_muonID_MSEffDen_fcal[1][0], -1);

    t_muonID_MSEff_fcal_PbPb18[2] = new TEfficiency (*(h_muonID_MSEffNum_fcal[0][2]), *(h_muonID_MSEffDen_fcal[0][2]));
    t_muonID_MSEff_fcal_PbPb15[2] = new TEfficiency (*(h_muonID_MSEffNum_fcal[1][2]), *(h_muonID_MSEffDen_fcal[1][2]));
  }



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

  TGAE* g_muonTrigEffPtpp = TEff2TGAE (t_muonTrigEff_pt_pp);
  TGAE* g_muonTrigEffPtPbPb18 = TEff2TGAE (t_muonTrigEff_pt_PbPb18);

  g_muonTrigEffPtpp->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  g_muonTrigEffPtPbPb18->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  g_muonTrigEffPtpp->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEffPtPbPb18->GetYaxis ()->SetTitle ("Muon trigger efficiency");

  g_muonTrigEffPtpp->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonTrigEffPtPbPb18->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonTrigEffPtpp->GetXaxis ()->SetRangeUser (5, 80);
  g_muonTrigEffPtPbPb18->GetXaxis ()->SetRangeUser (5, 80);

  g_muonTrigEffPtpp->SetMarkerStyle (kOpenSquare);
  g_muonTrigEffPtpp->SetMarkerColor (kAzure+2);
  g_muonTrigEffPtpp->SetLineColor (kAzure+2);
  g_muonTrigEffPtPbPb18->SetMarkerStyle (kOpenCircle);
  g_muonTrigEffPtPbPb18->SetMarkerColor (kRed+1);
  g_muonTrigEffPtPbPb18->SetLineColor (kRed+1);

  p_1->cd ();
  g_muonTrigEffPtpp->Draw ("AP");
  g_muonTrigEffPtPbPb18->Draw ("P");

  myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.22, 0.82, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.22, 0.765, kBlack, "HLT_mu14", 0.04);

  TGAE* g_muonTrigEffEtapp = TEff2TGAE (t_muonTrigEff_eta_pp);
  TGAE* g_muonTrigEffEtaPbPb18 = TEff2TGAE (t_muonTrigEff_eta_PbPb18);

  g_muonTrigEffEtapp->GetXaxis ()->SetTitle ("#eta");
  g_muonTrigEffEtaPbPb18->GetXaxis ()->SetTitle ("#eta");
  g_muonTrigEffEtapp->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEffEtaPbPb18->GetYaxis ()->SetTitle ("Muon trigger efficiency");

  g_muonTrigEffEtapp->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonTrigEffEtaPbPb18->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonTrigEffEtapp->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  g_muonTrigEffEtaPbPb18->GetXaxis ()->SetRangeUser (-2.5, 2.5);

  g_muonTrigEffEtapp->SetMarkerStyle (kOpenSquare);
  g_muonTrigEffEtapp->SetMarkerColor (kAzure+2);
  g_muonTrigEffEtapp->SetLineColor (kAzure+2);
  g_muonTrigEffEtaPbPb18->SetMarkerStyle (kOpenCircle);
  g_muonTrigEffEtaPbPb18->SetMarkerColor (kRed+1);
  g_muonTrigEffEtaPbPb18->SetLineColor (kRed+1);

  p_2->cd ();
  g_muonTrigEffEtapp->Draw ("AP");
  g_muonTrigEffEtaPbPb18->Draw ("P");

  myText (0.61, 0.89, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.042);
  myMarkerTextNoLine (0.22, 0.893, kAzure+2, kOpenSquare, "2017 #it{pp}, 258 pb^{-1}", 1.5, 0.042);
  myMarkerTextNoLine (0.22, 0.835, kRed+1, kOpenCircle, "2018 Pb+Pb, 1.4 nb^{-1}", 1.5, 0.042);

  c_1->SaveAs ("Plots/MuonTrigEffs/pp_PbPb_pt_eta.pdf");



  TGAE* g_muonMed_IDEffPtpp = TEff2TGAE (t_muonMed_IDEff_pt_pp[1]);
  TGAE* g_muonMed_IDEffPtPbPb18 = TEff2TGAE (t_muonMed_IDEff_pt_PbPb18[1]);

  g_muonMed_IDEffPtpp->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  g_muonMed_IDEffPtPbPb18->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  g_muonMed_IDEffPtpp->GetYaxis ()->SetTitle ("#varepsilon (Medium | ID)");
  g_muonMed_IDEffPtPbPb18->GetYaxis ()->SetTitle ("#varepsilon (Medium | ID)");

  g_muonMed_IDEffPtpp->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonMed_IDEffPtPbPb18->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonMed_IDEffPtpp->GetXaxis ()->SetRangeUser (5, 80);
  g_muonMed_IDEffPtPbPb18->GetXaxis ()->SetRangeUser (5, 80);

  g_muonMed_IDEffPtpp->SetMarkerStyle (kOpenSquare);
  g_muonMed_IDEffPtpp->SetMarkerColor (kAzure+2);
  g_muonMed_IDEffPtpp->SetLineColor (kAzure+2);
  g_muonMed_IDEffPtPbPb18->SetMarkerStyle (kOpenCircle);
  g_muonMed_IDEffPtPbPb18->SetMarkerColor (kRed+1);
  g_muonMed_IDEffPtPbPb18->SetLineColor (kRed+1);

  p_1->cd ();
  p_1->Clear ();
  g_muonMed_IDEffPtpp->Draw ("AP");
  g_muonMed_IDEffPtPbPb18->Draw ("P");

  myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.22, 0.82, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);

  TGAE* g_muonMed_IDEffEtapp = TEff2TGAE (t_muonMed_IDEff_eta_pp[1]);
  TGAE* g_muonMed_IDEffEtaPbPb18 = TEff2TGAE (t_muonMed_IDEff_eta_PbPb18[1]);

  g_muonMed_IDEffEtapp->GetXaxis ()->SetTitle ("#eta");
  g_muonMed_IDEffEtaPbPb18->GetXaxis ()->SetTitle ("#eta");
  g_muonMed_IDEffEtapp->GetYaxis ()->SetTitle ("#varepsilon (Medium | ID)");
  g_muonMed_IDEffEtaPbPb18->GetYaxis ()->SetTitle ("#varepsilon (Medium | ID)");

  g_muonMed_IDEffEtapp->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonMed_IDEffEtaPbPb18->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonMed_IDEffEtapp->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  g_muonMed_IDEffEtaPbPb18->GetXaxis ()->SetRangeUser (-2.5, 2.5);

  g_muonMed_IDEffEtapp->SetMarkerStyle (kOpenSquare);
  g_muonMed_IDEffEtapp->SetMarkerColor (kAzure+2);
  g_muonMed_IDEffEtapp->SetLineColor (kAzure+2);
  g_muonMed_IDEffEtaPbPb18->SetMarkerStyle (kOpenCircle);
  g_muonMed_IDEffEtaPbPb18->SetMarkerColor (kRed+1);
  g_muonMed_IDEffEtaPbPb18->SetLineColor (kRed+1);

  p_2->cd ();
  p_2->Clear ();
  g_muonMed_IDEffEtapp->Draw ("AP");
  g_muonMed_IDEffEtaPbPb18->Draw ("P");

  myText (0.61, 0.89, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.042);
  myMarkerTextNoLine (0.22, 0.893, kAzure+2, kOpenSquare, "2017 #it{pp}, 258 pb^{-1}", 1.5, 0.042);
  myMarkerTextNoLine (0.22, 0.835, kRed+1, kOpenCircle, "2018 Pb+Pb, 1.4 nb^{-1}", 1.5, 0.042);

  c_1->SaveAs ("Plots/MuonMed_IDEffs/pp_PbPb_pt_eta.pdf");



  TCanvas* c_2 = new TCanvas ("c_2", "", 800, 600);
  FormatTH2Canvas (c_2, true);

  TH1* h2_muonTrigEff_eta_phi_pp = TEff2TH1 (t_muonTrigEff_eta_phi_pp);
  h2_muonTrigEff_eta_phi_pp->GetYaxis ()->SetTitleOffset (1.1);
  h2_muonTrigEff_eta_phi_pp->GetZaxis ()->SetRangeUser (0, 1);
  h2_muonTrigEff_eta_phi_pp->GetZaxis ()->SetTitle ("Muon trigger efficiency");
  h2_muonTrigEff_eta_phi_pp->Draw ("colz");
  myText (0.55, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.55, 0.84, kBlack, "#it{pp}, 5.02 TeV", 0.045);
  myText (0.55, 0.78, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.045);
  myText (0.55, 0.72, kBlack, "HLT_mu14", 0.045);
  c_2->SaveAs ("Plots/MuonTrigEffs/pp_eta_phi.pdf");

  TH1* h2_muonTrigEff_eta_phi_PbPb18 = TEff2TH1 (t_muonTrigEff_eta_phi_PbPb18);
  h2_muonTrigEff_eta_phi_PbPb18->GetYaxis ()->SetTitleOffset (1.1);
  h2_muonTrigEff_eta_phi_PbPb18->GetZaxis ()->SetRangeUser (0, 1);
  h2_muonTrigEff_eta_phi_PbPb18->GetZaxis ()->SetTitle ("Muon trigger efficiency");
  h2_muonTrigEff_eta_phi_PbPb18->Draw ("colz");
  myText (0.55, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.55, 0.84, kBlack, "Pb+Pb, 5.02 TeV", 0.045);
  myText (0.55, 0.78, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.045);
  myText (0.55, 0.72, kBlack, "HLT_mu14", 0.045);
  c_2->SaveAs ("Plots/MuonTrigEffs/PbPb18_eta_phi.pdf");


  TCanvas* c_3 = new TCanvas ("c_3", "", 800, 600);
  FormatTH2Canvas (c_3, false);
  c_3->cd ();
  TGAE* g_muonTrigEff_fcal_PbPb18 = TEff2TGAE (t_muonTrigEff_fcal_PbPb18);
  g_muonTrigEff_fcal_PbPb18->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [TeV]");
  g_muonTrigEff_fcal_PbPb18->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEff_fcal_PbPb18->GetXaxis ()->SetRangeUser (-0.5, 5);
  g_muonTrigEff_fcal_PbPb18->GetYaxis ()->SetRangeUser (0.4, 1);
  g_muonTrigEff_fcal_PbPb18->SetMarkerStyle (kOpenCircle);
  g_muonTrigEff_fcal_PbPb18->Draw ("AP");

  TF1* f_muonTrigEff_fcal_PbPb18 = new TF1 ("f_muonTrigEff_fcal_PbPb18", "[0]+[1]*x", 0, 5);
  g_muonTrigEff_fcal_PbPb18->Fit (f_muonTrigEff_fcal_PbPb18, "Q0");
  f_muonTrigEff_fcal_PbPb18->SetLineStyle (2);
  f_muonTrigEff_fcal_PbPb18->SetLineColor (kGreen+2);
  f_muonTrigEff_fcal_PbPb18->SetLineWidth (2);
  f_muonTrigEff_fcal_PbPb18->Draw ("same");

  myText (0.2, 0.35, kBlack, Form ("Slope = %s TeV^{-1}", FormatMeasurement (f_muonTrigEff_fcal_PbPb18->GetParameter (1), f_muonTrigEff_fcal_PbPb18->GetParError (1), 2)), 0.04);
  myText (0.2, 0.30, kBlack, Form ("y-int. = %s", FormatMeasurement (f_muonTrigEff_fcal_PbPb18->GetParameter (0), f_muonTrigEff_fcal_PbPb18->GetParError (0), 2)), 0.04);
  myText (0.2, 0.25, kBlack, Form ("#chi^{2} / dof = %.2f / %i", f_muonTrigEff_fcal_PbPb18->GetChisquare (), f_muonTrigEff_fcal_PbPb18->GetNDF ()), 0.04);

  myText (0.62, 0.4, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.62, 0.35, kBlack, "2018 Pb+Pb, 1.4 nb^{-1}", 0.04);
  myText (0.62, 0.3, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.62, 0.25, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.04);
  myMarkerTextNoLine (0.65, 0.875, kBlack, kOpenCircle, "HLT_mu14", 1.25, 0.04);
  c_3->SaveAs ("Plots/MuonTrigEffs/PbPb_cent.pdf");


  TH1D* h_electronTrigEffNum_pt[3];
  TH1D* h_electronTrigEffDen_pt[3];
  TH1D* h_electronTrigEffNum_eta[3];
  TH1D* h_electronTrigEffDen_eta[3];
  TH2D* h2_electronTrigEffNum_pt_eta[3];
  TH2D* h2_electronTrigEffDen_pt_eta[3];

  TH1D* h_electronTrigEffNum_fcal[2];
  TH1D* h_electronTrigEffDen_fcal[2];

  h_electronTrigEffNum_pt[0] = (TH1D*) inFile->Get ("h_electronTrigEffNum_pt_pp");
  h_electronTrigEffDen_pt[0] = (TH1D*) inFile->Get ("h_electronTrigEffDen_pt_pp");
  h_electronTrigEffNum_pt[1] = (TH1D*) inFile->Get ("h_electronTrigEffNum_pt_PbPb18");
  h_electronTrigEffDen_pt[1] = (TH1D*) inFile->Get ("h_electronTrigEffDen_pt_PbPb18");
  h_electronTrigEffNum_pt[2] = (TH1D*) inFile->Get ("h_electronTrigEffNum_pt_PbPb15");
  h_electronTrigEffDen_pt[2] = (TH1D*) inFile->Get ("h_electronTrigEffDen_pt_PbPb15");
  h_electronTrigEffNum_eta[0] = (TH1D*) inFile->Get ("h_electronTrigEffNum_eta_pp");
  h_electronTrigEffDen_eta[0] = (TH1D*) inFile->Get ("h_electronTrigEffDen_eta_pp");
  h_electronTrigEffNum_eta[1] = (TH1D*) inFile->Get ("h_electronTrigEffNum_eta_PbPb18");
  h_electronTrigEffDen_eta[1] = (TH1D*) inFile->Get ("h_electronTrigEffDen_eta_PbPb18");
  h_electronTrigEffNum_eta[2] = (TH1D*) inFile->Get ("h_electronTrigEffNum_eta_PbPb15");
  h_electronTrigEffDen_eta[2] = (TH1D*) inFile->Get ("h_electronTrigEffDen_eta_PbPb15");
  h2_electronTrigEffNum_pt_eta[0] = (TH2D*) inFile->Get ("h2_electronTrigEffNum_pt_eta_pp");
  h2_electronTrigEffDen_pt_eta[0] = (TH2D*) inFile->Get ("h2_electronTrigEffDen_pt_eta_pp");
  h2_electronTrigEffNum_pt_eta[1] = (TH2D*) inFile->Get ("h2_electronTrigEffNum_pt_eta_PbPb18");
  h2_electronTrigEffDen_pt_eta[1] = (TH2D*) inFile->Get ("h2_electronTrigEffDen_pt_eta_PbPb18");
  h2_electronTrigEffNum_pt_eta[2] = (TH2D*) inFile->Get ("h2_electronTrigEffNum_pt_eta_PbPb15");
  h2_electronTrigEffDen_pt_eta[2] = (TH2D*) inFile->Get ("h2_electronTrigEffDen_pt_eta_PbPb15");

  h_electronTrigEffNum_fcal[0] = (TH1D*) inFile->Get ("h_electronTrigEffNum_fcal_PbPb18");
  h_electronTrigEffDen_fcal[0] = (TH1D*) inFile->Get ("h_electronTrigEffDen_fcal_PbPb18");
  h_electronTrigEffNum_fcal[1] = (TH1D*) inFile->Get ("h_electronTrigEffNum_fcal_PbPb15");
  h_electronTrigEffDen_fcal[1] = (TH1D*) inFile->Get ("h_electronTrigEffDen_fcal_PbPb15");

  t_electronTrigEff_pt_pp = new TEfficiency (*(h_electronTrigEffNum_pt[0]), *(h_electronTrigEffDen_pt[0]));
  t_electronTrigEff_pt_PbPb18 = new TEfficiency (*(h_electronTrigEffNum_pt[1]), *(h_electronTrigEffDen_pt[1]));
  t_electronTrigEff_pt_PbPb15 = new TEfficiency (*(h_electronTrigEffNum_pt[2]), *(h_electronTrigEffDen_pt[2]));
  t_electronTrigEff_eta_pp = new TEfficiency (*(h_electronTrigEffNum_eta[0]), *(h_electronTrigEffDen_eta[0]));
  t_electronTrigEff_eta_PbPb18 = new TEfficiency (*(h_electronTrigEffNum_eta[1]), *(h_electronTrigEffDen_eta[1]));
  t_electronTrigEff_eta_PbPb15 = new TEfficiency (*(h_electronTrigEffNum_eta[2]), *(h_electronTrigEffDen_eta[2]));
  t_electronTrigEff_pt_eta_pp = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[0]), *(h2_electronTrigEffDen_pt_eta[0]));
  t_electronTrigEff_pt_eta_PbPb18 = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[1]), *(h2_electronTrigEffDen_pt_eta[1]));
  t_electronTrigEff_pt_eta_PbPb15 = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[2]), *(h2_electronTrigEffDen_pt_eta[2]));

  t_electronTrigEff_fcal_PbPb18 = new TEfficiency (*(h_electronTrigEffNum_fcal[0]), *(h_electronTrigEffDen_fcal[0]));
  t_electronTrigEff_fcal_PbPb15 = new TEfficiency (*(h_electronTrigEffNum_fcal[1]), *(h_electronTrigEffDen_fcal[1]));


  TGAE* g_electronTrigEffPtpp = TEff2TGAE (t_electronTrigEff_pt_pp);
  TGAE* g_electronTrigEffPtPbPb18 = TEff2TGAE (t_electronTrigEff_pt_PbPb18);
  TGAE* g_electronTrigEffPtPbPb15 = TEff2TGAE (t_electronTrigEff_pt_PbPb15);

  g_electronTrigEffPtpp->GetXaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  g_electronTrigEffPtPbPb18->GetXaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  g_electronTrigEffPtpp->GetYaxis ()->SetTitle ("Electron trigger efficiency");
  g_electronTrigEffPtPbPb18->GetYaxis ()->SetTitle ("Electron trigger efficiency");

  g_electronTrigEffPtpp->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_electronTrigEffPtPbPb18->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_electronTrigEffPtpp->GetXaxis ()->SetRangeUser (5, 80);
  g_electronTrigEffPtPbPb18->GetXaxis ()->SetRangeUser (5, 80);

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

  myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.22, 0.82, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.22, 0.765, kBlack, "HLT_e15_lhloose(_ion)_L1EM12", 0.04);

  TGAE* g_electronTrigEffEtapp = TEff2TGAE (t_electronTrigEff_eta_pp);
  TGAE* g_electronTrigEffEtaPbPb18 = TEff2TGAE (t_electronTrigEff_eta_PbPb18);

  g_electronTrigEffEtapp->GetXaxis ()->SetTitle ("#eta");
  g_electronTrigEffEtaPbPb18->GetXaxis ()->SetTitle ("#eta");
  g_electronTrigEffEtapp->GetYaxis ()->SetTitle ("Electron trigger efficiency");
  g_electronTrigEffEtaPbPb18->GetYaxis ()->SetTitle ("Electron trigger efficiency");

  g_electronTrigEffEtapp->GetYaxis ()->SetRangeUser (0.0, 1.1);
  g_electronTrigEffEtaPbPb18->GetYaxis ()->SetRangeUser (0.0, 1.1);
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

  //myText (0.61, 0.89, kBlack, "20 < #it{p}_{T}^{#mu} < 80 GeV", 0.042);
  myMarkerTextNoLine (0.22, 0.28, kAzure+2, kOpenSquare, "2017 #it{pp}, 258 pb^{-1}", 1.5, 0.042);
  myMarkerTextNoLine (0.22, 0.225, kRed+1, kOpenCircle, "2018 Pb+Pb, 1.4 nb^{-1}", 1.5, 0.042);

  c_1->SaveAs ("Plots/ElectronTrigEffs/pp_PbPb18_pt_eta.pdf");


  FormatTH2Canvas (c_2, true);
  c_2->cd ();
  TH1* h2_electronTrigEff_pt_eta_pp = TEff2TH1 (t_electronTrigEff_pt_eta_pp);
  h2_electronTrigEff_pt_eta_pp->GetYaxis ()->SetTitleOffset (1.1);
  h2_electronTrigEff_pt_eta_pp->GetZaxis ()->SetRangeUser (0, 1);
  h2_electronTrigEff_pt_eta_pp->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronTrigEff_pt_eta_pp->GetZaxis ()->SetTitle ("Electron trigger efficiency");
  h2_electronTrigEff_pt_eta_pp->Draw ("colz");
  myText (0.55, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.55, 0.84, kBlack, "2017 #it{pp}, 5.02 TeV", 0.045);
  myText (0.55, 0.78, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.045);
  myText (0.55, 0.72, kBlack, "HLT_e15_lhloose_L1EM12", 0.045);
  c_2->SaveAs ("Plots/ElectronTrigEffs/pp_pt_eta.pdf");

  TH1* h2_electronTrigEff_pt_eta_PbPb18 = TEff2TH1 (t_electronTrigEff_pt_eta_PbPb18);
  h2_electronTrigEff_pt_eta_PbPb18->GetYaxis ()->SetTitleOffset (1.1);
  h2_electronTrigEff_pt_eta_PbPb18->GetZaxis ()->SetRangeUser (0, 1);
  h2_electronTrigEff_pt_eta_PbPb18->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronTrigEff_pt_eta_PbPb18->GetZaxis ()->SetTitle ("Electron trigger efficiency");
  h2_electronTrigEff_pt_eta_PbPb18->Draw ("colz");
  myText (0.55, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.55, 0.84, kBlack, "2018 Pb+Pb, 5.02 TeV", 0.045);
  myText (0.55, 0.78, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.045);
  myText (0.55, 0.72, kBlack, "HLT_e15_lhloose_ion_L1EM12", 0.045);
  c_2->SaveAs ("Plots/ElectronTrigEffs/PbPb18_pt_eta.pdf");


  FormatTH2Canvas (c_3, false);
  c_3->cd ();
  TGAE* g_electronTrigEff_fcal_PbPb18 = TEff2TGAE (t_electronTrigEff_fcal_PbPb18);
  g_electronTrigEff_fcal_PbPb18->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [TeV]");
  g_electronTrigEff_fcal_PbPb18->GetYaxis ()->SetTitle ("Electron trigger efficiency");
  g_electronTrigEff_fcal_PbPb18->GetXaxis ()->SetRangeUser (-0.5, 5);
  g_electronTrigEff_fcal_PbPb18->GetYaxis ()->SetRangeUser (0.4, 1);
  g_electronTrigEff_fcal_PbPb18->SetMarkerStyle (kOpenCircle);
  g_electronTrigEff_fcal_PbPb18->Draw ("AP");

  TF1* f_electronTrigEff_fcal_PbPb18 = new TF1 ("f_electronTrigEff_fcal_PbPb18", "[0]+[1]*x", 0, 5);
  f_electronTrigEff_fcal_PbPb18->SetParameter (0, 1);
  f_electronTrigEff_fcal_PbPb18->SetParameter (1, 0);
  g_electronTrigEff_fcal_PbPb18->Fit (f_electronTrigEff_fcal_PbPb18, "Q0");
  f_electronTrigEff_fcal_PbPb18->SetLineStyle (2);
  f_electronTrigEff_fcal_PbPb18->SetLineColor (kGreen+2);
  f_electronTrigEff_fcal_PbPb18->SetLineWidth (2);
  f_electronTrigEff_fcal_PbPb18->Draw ("same");

  myText (0.2, 0.35, kBlack, Form ("Slope = %s TeV^{-1}", FormatMeasurement (f_electronTrigEff_fcal_PbPb18->GetParameter (1), f_electronTrigEff_fcal_PbPb18->GetParError (1), 2)), 0.04);
  myText (0.2, 0.30, kBlack, Form ("y-int. = %s", FormatMeasurement (f_electronTrigEff_fcal_PbPb18->GetParameter (0), f_electronTrigEff_fcal_PbPb18->GetParError (0), 2)), 0.04);
  myText (0.2, 0.25, kBlack, Form ("#chi^{2} / dof = %.2f / %i", f_electronTrigEff_fcal_PbPb18->GetChisquare (), f_electronTrigEff_fcal_PbPb18->GetNDF ()), 0.04);

  myText (0.62, 0.4, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.62, 0.35, kBlack, "2018 Pb+Pb, 1.7 nb^{-1}", 0.04);
  myText (0.62, 0.3, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.62, 0.25, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.04);
  myMarkerTextNoLine (0.65, 0.875, kBlack, kOpenCircle, "HLT_e15_lhloose_ion_L1EM12", 1.25, 0.04);
  c_3->SaveAs ("Plots/ElectronTrigEffs/PbPb18_cent.pdf");



  TH2D* h2_zmumuTrigEffNum_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffNum_pt_y_pp");
  TH2D* h2_zmumuTrigEffDen_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffDen_pt_y_pp");
  TH2D* h2_zmumuTrigEffNum_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffNum_pt_y_PbPb18");
  TH2D* h2_zmumuTrigEffDen_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffDen_pt_y_PbPb18");
  TH2D* h2_zmumuTrigEffNum_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffNum_pt_y_PbPb15");
  TH2D* h2_zmumuTrigEffDen_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffDen_pt_y_PbPb15");

  TH2D* h2_zmumuTrigEff_pt_y_pp = (TH2D*) h2_zmumuTrigEffNum_pt_y_pp->Clone ("h2_zmumuTrigEffNum_pt_y_pp");
  h2_zmumuTrigEff_pt_y_pp->Divide (h2_zmumuTrigEffDen_pt_y_pp);
  TH2D* h2_zmumuTrigEff_pt_y_PbPb18 = (TH2D*) h2_zmumuTrigEffNum_pt_y_PbPb18->Clone ("h2_zmumuTrigEffNum_pt_y_PbPb18");
  h2_zmumuTrigEff_pt_y_PbPb18->Divide (h2_zmumuTrigEffDen_pt_y_PbPb18);
  TH2D* h2_zmumuTrigEff_pt_y_PbPb15 = (TH2D*) h2_zmumuTrigEffNum_pt_y_PbPb15->Clone ("h2_zmumuTrigEffNum_pt_y_PbPb15");
  h2_zmumuTrigEff_pt_y_PbPb15->Divide (h2_zmumuTrigEffDen_pt_y_PbPb15);

  TH2D* h2_zeeTrigEffNum_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zeeTrigEffNum_pt_y_pp");
  TH2D* h2_zeeTrigEffDen_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zeeTrigEffDen_pt_y_pp");
  TH2D* h2_zeeTrigEffNum_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zeeTrigEffNum_pt_y_PbPb18");
  TH2D* h2_zeeTrigEffDen_pt_y_PbPb18 = (TH2D*) (inFile)->Get ("h2_zeeTrigEffDen_pt_y_PbPb18");
  TH2D* h2_zeeTrigEffNum_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zeeTrigEffNum_pt_y_PbPb15");
  TH2D* h2_zeeTrigEffDen_pt_y_PbPb15 = (TH2D*) (inFile)->Get ("h2_zeeTrigEffDen_pt_y_PbPb15");

  TH2D* h2_zeeTrigEff_pt_y_pp = (TH2D*) h2_zeeTrigEffNum_pt_y_pp->Clone ("h2_zeeTrigEffNum_pt_y_pp");
  h2_zeeTrigEff_pt_y_pp->Divide (h2_zeeTrigEffDen_pt_y_pp);
  TH2D* h2_zeeTrigEff_pt_y_PbPb18 = (TH2D*) h2_zeeTrigEffNum_pt_y_PbPb18->Clone ("h2_zeeTrigEffNum_pt_y_PbPb18");
  h2_zeeTrigEff_pt_y_PbPb18->Divide (h2_zeeTrigEffDen_pt_y_PbPb18);
  TH2D* h2_zeeTrigEff_pt_y_PbPb15 = (TH2D*) h2_zeeTrigEffNum_pt_y_PbPb15->Clone ("h2_zeeTrigEffNum_pt_y_PbPb15");
  h2_zeeTrigEff_pt_y_PbPb15->Divide (h2_zeeTrigEffDen_pt_y_PbPb15);


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
  c_4->SaveAs ("Plots/ZmumuTrigEffs/pp_pt_y.pdf");

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
  c_4->SaveAs ("Plots/ZmumuTrigEffs/PbPb_pt_y.pdf");

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
  c_4->SaveAs ("Plots/ZeeTrigEffs/pp_pt_y.pdf");

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
  c_4->SaveAs ("Plots/ZeeTrigEffs/PbPb_pt_y.pdf");

}

#endif
