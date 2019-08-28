#ifndef __PlotTriggerEffs_C__
#define __PlotTriggerEffs_C__

#include "Params.h"
//#include "ZTrackUtilities.h"

#include <ArrayTemplates.h>

#include <TEfficiency.h>

#include <iostream>

typedef TGraphAsymmErrors TGAE;

TEfficiency* t_muonTrigEff_pt[4];
TEfficiency* t_muonTrigEff_eta[4][9];
TEfficiency* t_muonTrigEff_fcal;
TEfficiency* t_muonTrigEff_eta_phi_periph;
TEfficiency* t_muonTrigEff_eta_phi_cent;
TEfficiency* t_muonTrigEff_eta_phi_pp;
TEfficiency* t_muonTrigEff_eta_phi_PbPb;

TEfficiency* t_electronTrigEff_pt_pp;
TEfficiency* t_electronTrigEff_pt_PbPb;
TEfficiency* t_electronTrigEff_eta_pp;
TEfficiency* t_electronTrigEff_eta_PbPb;
TEfficiency* t_electronTrigEff_pt_eta_pp;
TEfficiency* t_electronTrigEff_pt_eta_PbPb;
TEfficiency* t_electronTrigEff_fcal;


bool isPbPb = true;

void PlotTriggerEffs () {

  SetupDirectories ("TriggerEffTagAndProbe/", "ZTrackAnalysis/");
  TFile* inFile = new TFile (Form ("%s/Nominal/outFile.root", rootPath.Data ()), "read");

  TH1D* h_muonTrigEffNum_pt[4];
  TH1D* h_muonTrigEffDen_pt[4];
  TH1D* h_muonTrigEffNum_eta[4][9];
  TH1D* h_muonTrigEffDen_eta[4][9];
  TH1D* h_muonTrigEffNum_fcal;
  TH1D* h_muonTrigEffDen_fcal;

  h_muonTrigEffNum_pt[0] = (TH1D*) inFile->Get ("h_muonTrigEffNum_pt_periph");
  h_muonTrigEffDen_pt[0] = (TH1D*) inFile->Get ("h_muonTrigEffDen_pt_periph");
  h_muonTrigEffNum_pt[1] = (TH1D*) inFile->Get ("h_muonTrigEffNum_pt_cent");
  h_muonTrigEffDen_pt[1] = (TH1D*) inFile->Get ("h_muonTrigEffDen_pt_cent");
  h_muonTrigEffNum_pt[2] = (TH1D*) inFile->Get ("h_muonTrigEffNum_pt_pp");
  h_muonTrigEffDen_pt[2] = (TH1D*) inFile->Get ("h_muonTrigEffDen_pt_pp");
  h_muonTrigEffNum_pt[3] = (TH1D*) inFile->Get ("h_muonTrigEffNum_pt_PbPb");
  h_muonTrigEffDen_pt[3] = (TH1D*) inFile->Get ("h_muonTrigEffDen_pt_PbPb");

  t_muonTrigEff_pt[0] = new TEfficiency (*(h_muonTrigEffNum_pt[0]), *(h_muonTrigEffDen_pt[0]));
  t_muonTrigEff_pt[1] = new TEfficiency (*(h_muonTrigEffNum_pt[1]), *(h_muonTrigEffDen_pt[1]));
  t_muonTrigEff_pt[2] = new TEfficiency (*(h_muonTrigEffNum_pt[2]), *(h_muonTrigEffDen_pt[2]));
  t_muonTrigEff_pt[3] = new TEfficiency (*(h_muonTrigEffNum_pt[3]), *(h_muonTrigEffDen_pt[3]));

  for (int iPhi = 0; iPhi < 9; iPhi++) {
    h_muonTrigEffNum_eta[0][iPhi] = (TH1D*) inFile->Get (Form ("h_muonTrigEffNum_eta_iPhi%i_periph", iPhi));
    h_muonTrigEffDen_eta[0][iPhi] = (TH1D*) inFile->Get (Form ("h_muonTrigEffDen_eta_iPhi%i_periph", iPhi));
    h_muonTrigEffNum_eta[1][iPhi] = (TH1D*) inFile->Get (Form ("h_muonTrigEffNum_eta_iPhi%i_cent", iPhi));
    h_muonTrigEffDen_eta[1][iPhi] = (TH1D*) inFile->Get (Form ("h_muonTrigEffDen_eta_iPhi%i_cent", iPhi));
    h_muonTrigEffNum_eta[2][iPhi] = (TH1D*) inFile->Get (Form ("h_muonTrigEffNum_eta_iPhi%i_pp", iPhi));
    h_muonTrigEffDen_eta[2][iPhi] = (TH1D*) inFile->Get (Form ("h_muonTrigEffDen_eta_iPhi%i_pp", iPhi));
    h_muonTrigEffNum_eta[3][iPhi] = (TH1D*) inFile->Get (Form ("h_muonTrigEffNum_eta_iPhi%i_PbPb", iPhi));
    h_muonTrigEffDen_eta[3][iPhi] = (TH1D*) inFile->Get (Form ("h_muonTrigEffDen_eta_iPhi%i_PbPb", iPhi));
    t_muonTrigEff_eta[0][iPhi] = new TEfficiency (*(h_muonTrigEffNum_eta[0][iPhi]), *(h_muonTrigEffDen_eta[0][iPhi]));
    t_muonTrigEff_eta[1][iPhi] = new TEfficiency (*(h_muonTrigEffNum_eta[1][iPhi]), *(h_muonTrigEffDen_eta[1][iPhi]));
    t_muonTrigEff_eta[2][iPhi] = new TEfficiency (*(h_muonTrigEffNum_eta[2][iPhi]), *(h_muonTrigEffDen_eta[2][iPhi]));
    t_muonTrigEff_eta[3][iPhi] = new TEfficiency (*(h_muonTrigEffNum_eta[3][iPhi]), *(h_muonTrigEffDen_eta[3][iPhi]));
  }

  h_muonTrigEffNum_fcal = (TH1D*) inFile->Get ("h_muonTrigEffNum_fcal");
  h_muonTrigEffDen_fcal = (TH1D*) inFile->Get ("h_muonTrigEffDen_fcal");
  t_muonTrigEff_fcal = new TEfficiency (*(h_muonTrigEffNum_fcal), *(h_muonTrigEffDen_fcal));


  TH2D* h2_muonTrigEffNum_eta_phi_periph = new TH2D ("h_muonTrigEffNum_eta_phi_periph", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  TH2D* h2_muonTrigEffDen_eta_phi_periph = new TH2D ("h_muonTrigEffDen_eta_phi_periph", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  TH2D* h2_muonTrigEffNum_eta_phi_cent = new TH2D ("h_muonTrigEffNum_eta_phi_cent", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  TH2D* h2_muonTrigEffDen_eta_phi_cent = new TH2D ("h_muonTrigEffDen_eta_phi_cent", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  TH2D* h2_muonTrigEffNum_eta_phi_pp = new TH2D ("h_muonTrigEffNum_eta_phi_pp", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  TH2D* h2_muonTrigEffDen_eta_phi_pp = new TH2D ("h_muonTrigEffDen_eta_phi_pp", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  TH2D* h2_muonTrigEffNum_eta_phi_PbPb = new TH2D ("h_muonTrigEffNum_eta_phi_PbPb", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  TH2D* h2_muonTrigEffDen_eta_phi_PbPb = new TH2D ("h_muonTrigEffDen_eta_phi_PbPb", ";#eta;#phi;Counts", 25, -2.5, 2.5, 8, 0, 2*pi);
  for (int iPhi = 0; iPhi < 8; iPhi++) {
    TH1D* num = h_muonTrigEffNum_eta[0][iPhi];
    TH1D* den = h_muonTrigEffDen_eta[0][iPhi];
    for (int ix = 1; ix <= num->GetNbinsX (); ix++) {
      h2_muonTrigEffNum_eta_phi_periph->SetBinContent (ix, iPhi+1, num->GetBinContent (ix));
      h2_muonTrigEffDen_eta_phi_periph->SetBinContent (ix, iPhi+1, den->GetBinContent (ix));
    }
    num = h_muonTrigEffNum_eta[1][iPhi];
    den = h_muonTrigEffDen_eta[1][iPhi];
    for (int ix = 1; ix <= num->GetNbinsX (); ix++) {
      h2_muonTrigEffNum_eta_phi_cent->SetBinContent (ix, iPhi+1, num->GetBinContent (ix));
      h2_muonTrigEffDen_eta_phi_cent->SetBinContent (ix, iPhi+1, den->GetBinContent (ix));
    }
    num = h_muonTrigEffNum_eta[2][iPhi];
    den = h_muonTrigEffDen_eta[2][iPhi];
    for (int ix = 1; ix <= num->GetNbinsX (); ix++) {
      h2_muonTrigEffNum_eta_phi_pp->SetBinContent (ix, iPhi+1, num->GetBinContent (ix));
      h2_muonTrigEffDen_eta_phi_pp->SetBinContent (ix, iPhi+1, den->GetBinContent (ix));
    }
    num = h_muonTrigEffNum_eta[3][iPhi];
    den = h_muonTrigEffDen_eta[3][iPhi];
    for (int ix = 1; ix <= num->GetNbinsX (); ix++) {
      h2_muonTrigEffNum_eta_phi_PbPb->SetBinContent (ix, iPhi+1, num->GetBinContent (ix));
      h2_muonTrigEffDen_eta_phi_PbPb->SetBinContent (ix, iPhi+1, den->GetBinContent (ix));
    }
  }

  t_muonTrigEff_eta_phi_periph = new TEfficiency (*(h2_muonTrigEffNum_eta_phi_periph), *(h2_muonTrigEffDen_eta_phi_periph));
  t_muonTrigEff_eta_phi_cent = new TEfficiency (*(h2_muonTrigEffNum_eta_phi_cent), *(h2_muonTrigEffDen_eta_phi_cent));
  t_muonTrigEff_eta_phi_pp = new TEfficiency (*(h2_muonTrigEffNum_eta_phi_pp), *(h2_muonTrigEffDen_eta_phi_pp));
  t_muonTrigEff_eta_phi_PbPb = new TEfficiency (*(h2_muonTrigEffNum_eta_phi_PbPb), *(h2_muonTrigEffDen_eta_phi_PbPb));


  //TH2D* h2_muonTrigEffNum = (TH2D*) inFile->Get ("h2_muonTrigEffNum");
  //TH2D* h2_muonTrigEffDen = (TH2D*) inFile->Get ("h2_muonTrigEffDen");
  //TH2D* h2_electronTrigEffNum = (TH2D*) inFile->Get ("h2_electronTrigEffNum");
  //TH2D* h2_electronTrigEffDen = (TH2D*) inFile->Get ("h2_electronTrigEffDen");

  //TEfficiency* t_muonTrigEffPt = new TEfficiency (*(h2_muonTrigEffNum->ProjectionY ()), *(h2_muonTrigEffDen->ProjectionY ()));
  //TEfficiency* t_electronTrigEffPt = new TEfficiency (*(h2_electronTrigEffNum->ProjectionY ()), *(h2_electronTrigEffDen->ProjectionY ()));

  //TEfficiency* t_muonTrigEffFCal = new TEfficiency (*(h2_muonTrigEffNum->ProjectionX ("_px", h2_muonTrigEffNum->GetXaxis ()->FindBin (20), h2_muonTrigEffNum->GetNbinsY ())), *(h2_muonTrigEffDen->ProjectionX ("_px", h2_muonTrigEffDen->GetXaxis ()->FindBin (20), h2_muonTrigEffDen->GetNbinsY ())));
  //TEfficiency* t_electronTrigEffFCal = new TEfficiency (*(h2_electronTrigEffNum->ProjectionX ("_px", h2_electronTrigEffNum->GetXaxis ()->FindBin (20), h2_electronTrigEffNum->GetNbinsY ())), *(h2_electronTrigEffDen->ProjectionX ("_px", h2_electronTrigEffNum->GetXaxis ()->FindBin (20), h2_electronTrigEffDen->GetNbinsY ())));


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

  //TGAE* g_muonTrigEffPtPeriph = TEff2TGAE (t_muonTrigEff_pt[0]);
  //TGAE* g_muonTrigEffPtCent = TEff2TGAE (t_muonTrigEff_pt[1]);
  TGAE* g_muonTrigEffPtpp = TEff2TGAE (t_muonTrigEff_pt[2]);
  TGAE* g_muonTrigEffPtPbPb = TEff2TGAE (t_muonTrigEff_pt[3]);

  //g_muonTrigEffPtPeriph->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  //g_muonTrigEffPtCent->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  g_muonTrigEffPtpp->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  g_muonTrigEffPtPbPb->GetXaxis ()->SetTitle ("#it{p}_{T}^{#mu} [GeV]");
  //g_muonTrigEffPtPeriph->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  //g_muonTrigEffPtCent->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEffPtpp->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEffPtPbPb->GetYaxis ()->SetTitle ("Muon trigger efficiency");

  //g_muonTrigEffPtPeriph->GetYaxis ()->SetRangeUser (0.3, 1.1);
  //g_muonTrigEffPtCent->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonTrigEffPtpp->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonTrigEffPtPbPb->GetYaxis ()->SetRangeUser (0.3, 1.1);
  //g_muonTrigEffPtPeriph->GetXaxis ()->SetRangeUser (5, 80);
  //g_muonTrigEffPtCent->GetXaxis ()->SetRangeUser (5, 80);
  g_muonTrigEffPtpp->GetXaxis ()->SetRangeUser (5, 80);
  g_muonTrigEffPtPbPb->GetXaxis ()->SetRangeUser (5, 80);

  //g_muonTrigEffPtPeriph->SetMarkerStyle (kOpenCircle);
  //g_muonTrigEffPtPeriph->SetMarkerColor (kAzure+2);
  //g_muonTrigEffPtPeriph->SetLineColor (kAzure+2);
  //g_muonTrigEffPtCent->SetMarkerStyle (kOpenCircle);
  //g_muonTrigEffPtCent->SetMarkerColor (kRed+1);
  //g_muonTrigEffPtCent->SetLineColor (kRed+1);
  g_muonTrigEffPtpp->SetMarkerStyle (kOpenSquare);
  g_muonTrigEffPtpp->SetMarkerColor (kAzure+2);
  g_muonTrigEffPtpp->SetLineColor (kAzure+2);
  g_muonTrigEffPtPbPb->SetMarkerStyle (kOpenCircle);
  g_muonTrigEffPtPbPb->SetMarkerColor (kRed+1);
  g_muonTrigEffPtPbPb->SetLineColor (kRed+1);

  p_1->cd ();
  //g_muonTrigEffPtPeriph->Draw ("AP");
  //g_muonTrigEffPtCent->Draw ("P");
  g_muonTrigEffPtpp->Draw ("AP");
  g_muonTrigEffPtPbPb->Draw ("P");

  myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.22, 0.82, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.22, 0.765, kBlack, "HLT_mu14", 0.04);

  //TH1* h_muonTrigEffPtPeriph = TEff2TH1 (t_muonTrigEff_pt[0]);
  //TH1* h_muonTrigEffPtCent = TEff2TH1 (t_muonTrigEff_pt[1]);

  //h_muonTrigEffPtCent->Divide (h_muonTrigEffPtPeriph);
  //g_muonTrigEffPtCent = make_graph (h_muonTrigEffPtCent);

  //p_3->cd ();
  //g_muonTrigEffPtCent->GetYaxis ()->SetTitle ("Cent. / Periph.");
  //g_muonTrigEffPtCent->GetYaxis ()->CenterTitle (true);
  //g_muonTrigEffPtCent->SetMarkerStyle (kOpenCircle);
  //g_muonTrigEffPtCent->GetXaxis ()->SetRangeUser (5, 80);
  //g_muonTrigEffPtCent->GetYaxis ()->SetRangeUser (0.5, 1.5);
  //g_muonTrigEffPtCent->Draw ("AP");

  //TLine* l = new TLine (5, 1, 80, 1);
  //l->SetLineColor (kPink-8);
  //l->SetLineStyle (2);
  //l->Draw ("same");

  //for (int iPhi = 0; iPhi < 8; iPhi++) {
  //TGAE* g_muonTrigEffEtaPeriph = TEff2TGAE (t_muonTrigEff_eta[0][8]);
  //TGAE* g_muonTrigEffEtaCent = TEff2TGAE (t_muonTrigEff_eta[1][8]);
  TGAE* g_muonTrigEffEtapp = TEff2TGAE (t_muonTrigEff_eta[2][8]);
  TGAE* g_muonTrigEffEtaPbPb = TEff2TGAE (t_muonTrigEff_eta[3][8]);

  //g_muonTrigEffEtaPeriph->GetXaxis ()->SetTitle ("#eta");
  //g_muonTrigEffEtaCent->GetXaxis ()->SetTitle ("#eta");
  g_muonTrigEffEtapp->GetXaxis ()->SetTitle ("#eta");
  g_muonTrigEffEtaPbPb->GetXaxis ()->SetTitle ("#eta");
  //g_muonTrigEffEtaPeriph->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  //g_muonTrigEffEtaCent->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEffEtapp->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEffEtaPbPb->GetYaxis ()->SetTitle ("Muon trigger efficiency");

  //g_muonTrigEffEtaPeriph->GetYaxis ()->SetRangeUser (0.3, 1.1);
  //g_muonTrigEffEtaCent->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonTrigEffEtapp->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_muonTrigEffEtaPbPb->GetYaxis ()->SetRangeUser (0.3, 1.1);
  //g_muonTrigEffEtaPeriph->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  //g_muonTrigEffEtaCent->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  g_muonTrigEffEtapp->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  g_muonTrigEffEtaPbPb->GetXaxis ()->SetRangeUser (-2.5, 2.5);

  //g_muonTrigEffEtaPeriph->SetMarkerStyle (kOpenCircle);
  //g_muonTrigEffEtaPeriph->SetMarkerColor (kAzure+2);
  //g_muonTrigEffEtaPeriph->SetLineColor (kAzure+2);
  //g_muonTrigEffEtaCent->SetMarkerStyle (kOpenCircle);
  //g_muonTrigEffEtaCent->SetMarkerColor (kRed+1);
  //g_muonTrigEffEtaCent->SetLineColor (kRed+1);
  g_muonTrigEffEtapp->SetMarkerStyle (kOpenSquare);
  g_muonTrigEffEtapp->SetMarkerColor (kAzure+2);
  g_muonTrigEffEtapp->SetLineColor (kAzure+2);
  g_muonTrigEffEtaPbPb->SetMarkerStyle (kOpenCircle);
  g_muonTrigEffEtaPbPb->SetMarkerColor (kRed+1);
  g_muonTrigEffEtaPbPb->SetLineColor (kRed+1);

  p_2->cd ();
  //g_muonTrigEffEtaPeriph->Draw ("AP");
  //g_muonTrigEffEtaCent->Draw ("P");
  g_muonTrigEffEtapp->Draw ("AP");
  g_muonTrigEffEtaPbPb->Draw ("P");

  myText (0.61, 0.89, kBlack, "20 < #it{p}_{T}^{#mu} < 80 GeV", 0.042);
  myMarkerTextNoLine (0.22, 0.893, kAzure+2, kOpenSquare, "2017 #it{pp}, 258 pb^{-1}", 1.5, 0.042);
  myMarkerTextNoLine (0.22, 0.835, kRed+1, kOpenCircle, "2018 Pb+Pb, 1.4 nb^{-1}", 1.5, 0.042);

  //TH1* h_muonTrigEffEtaPeriph = TEff2TH1 (t_muonTrigEff_eta[0][8]);
  //TH1* h_muonTrigEffEtaCent = TEff2TH1 (t_muonTrigEff_eta[1][8]);
  //h_muonTrigEffEtaCent->Divide (h_muonTrigEffEtaPeriph);
  //g_muonTrigEffEtaCent = make_graph (h_muonTrigEffEtaCent);

  //p_4->cd ();
  //g_muonTrigEffEtaCent->SetMarkerStyle (kOpenCircle);
  //g_muonTrigEffEtaCent->GetYaxis ()->SetTitle ("Cent. / Periph.");
  //g_muonTrigEffEtaCent->GetYaxis ()->CenterTitle (true);
  //g_muonTrigEffEtaCent->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  //g_muonTrigEffEtaCent->GetYaxis ()->SetRangeUser (0.5, 1.5);
  //g_muonTrigEffEtaCent->Draw ("AP");

  //l = new TLine (-2.5, 1, 2.5, 1);
  //l->SetLineColor (kPink-8);
  //l->SetLineStyle (2);
  //l->Draw ("same");
  //}
  c_1->SaveAs ("Plots/MuonTrigEffs/pp_PbPb_pt_eta.pdf");

  TCanvas* c_2 = new TCanvas ("c_2", "", 800, 600);
  FormatTH2Canvas (c_2, true);

  //t_muonTrigEff_eta_phi_periph->Draw ("colz");
  TH1* h2_muonTrigEff_eta_phi_cent = TEff2TH1 (t_muonTrigEff_eta_phi_cent);
  h2_muonTrigEff_eta_phi_cent->GetYaxis ()->SetTitleOffset (1.1);
  h2_muonTrigEff_eta_phi_cent->GetZaxis ()->SetRangeUser (0, 1);
  h2_muonTrigEff_eta_phi_cent->GetZaxis ()->SetTitle ("Muon trigger efficiency");
  h2_muonTrigEff_eta_phi_cent->Draw ("colz");
  myText (0.55, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.55, 0.84, kBlack, "Pb+Pb 0-10%", 0.045);
  myText (0.55, 0.78, kBlack, "20 < #it{p}_{T}^{#mu} < 80 GeV", 0.045);
  myText (0.55, 0.72, kBlack, "HLT_mu14", 0.045);
  c_2->SaveAs ("Plots/MuonTrigEffs/cent_eta_phi.pdf");

  TH1* h2_muonTrigEff_eta_phi_periph = TEff2TH1 (t_muonTrigEff_eta_phi_periph);
  h2_muonTrigEff_eta_phi_periph->GetYaxis ()->SetTitleOffset (1.1);
  h2_muonTrigEff_eta_phi_periph->GetZaxis ()->SetRangeUser (0, 1);
  h2_muonTrigEff_eta_phi_periph->GetZaxis ()->SetTitle ("Muon trigger efficiency");
  h2_muonTrigEff_eta_phi_periph->Draw ("colz");
  myText (0.55, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.55, 0.84, kBlack, "Pb+Pb 60-80%", 0.045);
  myText (0.55, 0.78, kBlack, "20 < #it{p}_{T}^{#mu} < 80 GeV", 0.045);
  myText (0.55, 0.72, kBlack, "HLT_mu14", 0.045);
  c_2->SaveAs ("Plots/MuonTrigEffs/periph_eta_phi.pdf");

  TH1* h2_muonTrigEff_eta_phi_pp = TEff2TH1 (t_muonTrigEff_eta_phi_pp);
  h2_muonTrigEff_eta_phi_pp->GetYaxis ()->SetTitleOffset (1.1);
  h2_muonTrigEff_eta_phi_pp->GetZaxis ()->SetRangeUser (0, 1);
  h2_muonTrigEff_eta_phi_pp->GetZaxis ()->SetTitle ("Muon trigger efficiency");
  h2_muonTrigEff_eta_phi_pp->Draw ("colz");
  myText (0.55, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.55, 0.84, kBlack, "#it{pp}, 5.02 TeV", 0.045);
  myText (0.55, 0.78, kBlack, "20 < #it{p}_{T}^{#mu} < 80 GeV", 0.045);
  myText (0.55, 0.72, kBlack, "HLT_mu14", 0.045);
  c_2->SaveAs ("Plots/MuonTrigEffs/pp_eta_phi.pdf");

  TH1* h2_muonTrigEff_eta_phi_PbPb = TEff2TH1 (t_muonTrigEff_eta_phi_PbPb);
  h2_muonTrigEff_eta_phi_PbPb->GetYaxis ()->SetTitleOffset (1.1);
  h2_muonTrigEff_eta_phi_PbPb->GetZaxis ()->SetRangeUser (0, 1);
  h2_muonTrigEff_eta_phi_PbPb->GetZaxis ()->SetTitle ("Muon trigger efficiency");
  h2_muonTrigEff_eta_phi_PbPb->Draw ("colz");
  myText (0.55, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.55, 0.84, kBlack, "Pb+Pb, 5.02 TeV", 0.045);
  myText (0.55, 0.78, kBlack, "20 < #it{p}_{T}^{#mu} < 80 GeV", 0.045);
  myText (0.55, 0.72, kBlack, "HLT_mu14", 0.045);
  c_2->SaveAs ("Plots/MuonTrigEffs/PbPb_eta_phi.pdf");


  TCanvas* c_3 = new TCanvas ("c_3", "", 800, 600);
  FormatTH2Canvas (c_3, false);
  c_3->cd ();
  TGAE* g_muonTrigEff_fcal = TEff2TGAE (t_muonTrigEff_fcal);
  g_muonTrigEff_fcal->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [TeV]");
  g_muonTrigEff_fcal->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEff_fcal->GetXaxis ()->SetRangeUser (-0.5, 5);
  g_muonTrigEff_fcal->GetYaxis ()->SetRangeUser (0.4, 1);
  g_muonTrigEff_fcal->SetMarkerStyle (kOpenCircle);
  g_muonTrigEff_fcal->Draw ("AP");

  TF1* f_muonTrigEff_fcal = new TF1 ("f_muonTrigEff_fcal", "[0]+[1]*x", 0, 5);
  g_muonTrigEff_fcal->Fit (f_muonTrigEff_fcal, "Q0");
  f_muonTrigEff_fcal->SetLineStyle (2);
  f_muonTrigEff_fcal->SetLineColor (kGreen+2);
  f_muonTrigEff_fcal->SetLineWidth (2);
  f_muonTrigEff_fcal->Draw ("same");

  myText (0.2, 0.35, kBlack, Form ("Slope = %s TeV^{-1}", FormatMeasurement (f_muonTrigEff_fcal->GetParameter (1), f_muonTrigEff_fcal->GetParError (1), 2)), 0.04);
  myText (0.2, 0.30, kBlack, Form ("y-int. = %s", FormatMeasurement (f_muonTrigEff_fcal->GetParameter (0), f_muonTrigEff_fcal->GetParError (0), 2)), 0.04);
  myText (0.2, 0.25, kBlack, Form ("#chi^{2} / dof = %.2f / %i", f_muonTrigEff_fcal->GetChisquare (), f_muonTrigEff_fcal->GetNDF ()), 0.04);

  myText (0.62, 0.4, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.62, 0.35, kBlack, "2018 Pb+Pb, 1.4 nb^{-1}", 0.04);
  myText (0.62, 0.3, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.62, 0.25, kBlack, "20 < #it{p}_{T}^{#mu} < 80 GeV", 0.04);
  myMarkerTextNoLine (0.65, 0.875, kBlack, kOpenCircle, "HLT_mu14", 1.25, 0.04);
  c_3->SaveAs ("Plots/MuonTrigEffs/PbPb_cent.pdf");


  TH1D* h_electronTrigEffNum_pt[2];
  TH1D* h_electronTrigEffDen_pt[2];
  TH1D* h_electronTrigEffNum_eta[2];
  TH1D* h_electronTrigEffDen_eta[2];
  TH2D* h2_electronTrigEffNum_pt_eta[2];
  TH2D* h2_electronTrigEffDen_pt_eta[2];

  TH1D* h_electronTrigEffNum_fcal;
  TH1D* h_electronTrigEffDen_fcal;

  h_electronTrigEffNum_pt[0] = (TH1D*) inFile->Get ("h_electronTrigEffNum_pt_pp");
  h_electronTrigEffDen_pt[0] = (TH1D*) inFile->Get ("h_electronTrigEffDen_pt_pp");
  h_electronTrigEffNum_pt[1] = (TH1D*) inFile->Get ("h_electronTrigEffNum_pt_PbPb");
  h_electronTrigEffDen_pt[1] = (TH1D*) inFile->Get ("h_electronTrigEffDen_pt_PbPb");
  h_electronTrigEffNum_eta[0] = (TH1D*) inFile->Get ("h_electronTrigEffNum_eta_pp");
  h_electronTrigEffDen_eta[0] = (TH1D*) inFile->Get ("h_electronTrigEffDen_eta_pp");
  h_electronTrigEffNum_eta[1] = (TH1D*) inFile->Get ("h_electronTrigEffNum_eta_PbPb");
  h_electronTrigEffDen_eta[1] = (TH1D*) inFile->Get ("h_electronTrigEffDen_eta_PbPb");
  h2_electronTrigEffNum_pt_eta[0] = (TH2D*) inFile->Get ("h2_electronTrigEffNum_pt_eta_pp");
  h2_electronTrigEffDen_pt_eta[0] = (TH2D*) inFile->Get ("h2_electronTrigEffDen_pt_eta_pp");
  h2_electronTrigEffNum_pt_eta[1] = (TH2D*) inFile->Get ("h2_electronTrigEffNum_pt_eta_PbPb");
  h2_electronTrigEffDen_pt_eta[1] = (TH2D*) inFile->Get ("h2_electronTrigEffDen_pt_eta_PbPb");

  h_electronTrigEffNum_fcal = (TH1D*) inFile->Get ("h_electronTrigEffNum_fcal");
  h_electronTrigEffDen_fcal = (TH1D*) inFile->Get ("h_electronTrigEffDen_fcal");

  t_electronTrigEff_pt_pp = new TEfficiency (*(h_electronTrigEffNum_pt[0]), *(h_electronTrigEffDen_pt[0]));
  t_electronTrigEff_pt_PbPb = new TEfficiency (*(h_electronTrigEffNum_pt[1]), *(h_electronTrigEffDen_pt[1]));
  t_electronTrigEff_eta_pp = new TEfficiency (*(h_electronTrigEffNum_eta[0]), *(h_electronTrigEffDen_eta[0]));
  t_electronTrigEff_eta_PbPb = new TEfficiency (*(h_electronTrigEffNum_eta[1]), *(h_electronTrigEffDen_eta[1]));
  t_electronTrigEff_pt_eta_pp = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[0]), *(h2_electronTrigEffDen_pt_eta[0]));
  t_electronTrigEff_pt_eta_PbPb = new TEfficiency (*(h2_electronTrigEffNum_pt_eta[1]), *(h2_electronTrigEffDen_pt_eta[1]));

  t_electronTrigEff_fcal = new TEfficiency (*(h_electronTrigEffNum_fcal), *(h_electronTrigEffDen_fcal));


  TGAE* g_electronTrigEffPtpp = TEff2TGAE (t_electronTrigEff_pt_pp);
  TGAE* g_electronTrigEffPtPbPb = TEff2TGAE (t_electronTrigEff_pt_PbPb);

  g_electronTrigEffPtpp->GetXaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  g_electronTrigEffPtPbPb->GetXaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  g_electronTrigEffPtpp->GetYaxis ()->SetTitle ("Electron trigger efficiency");
  g_electronTrigEffPtPbPb->GetYaxis ()->SetTitle ("Electron trigger efficiency");

  g_electronTrigEffPtpp->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_electronTrigEffPtPbPb->GetYaxis ()->SetRangeUser (0.3, 1.1);
  g_electronTrigEffPtpp->GetXaxis ()->SetRangeUser (5, 80);
  g_electronTrigEffPtPbPb->GetXaxis ()->SetRangeUser (5, 80);

  g_electronTrigEffPtpp->SetMarkerStyle (kOpenSquare);
  g_electronTrigEffPtpp->SetMarkerColor (kAzure+2);
  g_electronTrigEffPtpp->SetLineColor (kAzure+2);
  g_electronTrigEffPtPbPb->SetMarkerStyle (kOpenCircle);
  g_electronTrigEffPtPbPb->SetMarkerColor (kRed+1);
  g_electronTrigEffPtPbPb->SetLineColor (kRed+1);

  p_1->cd ();
  g_electronTrigEffPtpp->Draw ("AP");
  g_electronTrigEffPtPbPb->Draw ("P");

  myText (0.22, 0.89, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.22, 0.82, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.22, 0.765, kBlack, "HLT_e15_lhloose(_ion)_L1EM12", 0.04);

  TGAE* g_electronTrigEffEtapp = TEff2TGAE (t_electronTrigEff_eta_pp);
  TGAE* g_electronTrigEffEtaPbPb = TEff2TGAE (t_electronTrigEff_eta_PbPb);

  g_electronTrigEffEtapp->GetXaxis ()->SetTitle ("#eta");
  g_electronTrigEffEtaPbPb->GetXaxis ()->SetTitle ("#eta");
  g_electronTrigEffEtapp->GetYaxis ()->SetTitle ("Electron trigger efficiency");
  g_electronTrigEffEtaPbPb->GetYaxis ()->SetTitle ("Electron trigger efficiency");

  g_electronTrigEffEtapp->GetYaxis ()->SetRangeUser (0.0, 1.1);
  g_electronTrigEffEtaPbPb->GetYaxis ()->SetRangeUser (0.0, 1.1);
  g_electronTrigEffEtapp->GetXaxis ()->SetRangeUser (-2.5, 2.5);
  g_electronTrigEffEtaPbPb->GetXaxis ()->SetRangeUser (-2.5, 2.5);

  g_electronTrigEffEtapp->SetMarkerStyle (kOpenSquare);
  g_electronTrigEffEtapp->SetMarkerColor (kAzure+2);
  g_electronTrigEffEtapp->SetLineColor (kAzure+2);
  g_electronTrigEffEtaPbPb->SetMarkerStyle (kOpenCircle);
  g_electronTrigEffEtaPbPb->SetMarkerColor (kRed+1);
  g_electronTrigEffEtaPbPb->SetLineColor (kRed+1);

  p_2->cd ();
  p_2->Clear ();
  g_electronTrigEffEtapp->Draw ("AP");
  g_electronTrigEffEtaPbPb->Draw ("P");

  //myText (0.61, 0.89, kBlack, "20 < #it{p}_{T}^{#mu} < 80 GeV", 0.042);
  myMarkerTextNoLine (0.22, 0.28, kAzure+2, kOpenSquare, "2017 #it{pp}, 258 pb^{-1}", 1.5, 0.042);
  myMarkerTextNoLine (0.22, 0.225, kRed+1, kOpenCircle, "2018 Pb+Pb, 1.4 nb^{-1}", 1.5, 0.042);

  c_1->SaveAs ("Plots/ElectronTrigEffs/pp_PbPb_pt_eta.pdf");


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

  TH1* h2_electronTrigEff_pt_eta_PbPb = TEff2TH1 (t_electronTrigEff_pt_eta_PbPb);
  h2_electronTrigEff_pt_eta_PbPb->GetYaxis ()->SetTitleOffset (1.1);
  h2_electronTrigEff_pt_eta_PbPb->GetZaxis ()->SetRangeUser (0, 1);
  h2_electronTrigEff_pt_eta_PbPb->GetYaxis ()->SetTitle ("#it{p}_{T}^{e} [GeV]");
  h2_electronTrigEff_pt_eta_PbPb->GetZaxis ()->SetTitle ("Electron trigger efficiency");
  h2_electronTrigEff_pt_eta_PbPb->Draw ("colz");
  myText (0.55, 0.90, kBlack, "#bf{#it{ATLAS}} Internal", 0.05);
  myText (0.55, 0.84, kBlack, "2018 Pb+Pb, 5.02 TeV", 0.045);
  myText (0.55, 0.78, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.045);
  myText (0.55, 0.72, kBlack, "HLT_e15_lhloose_ion_L1EM12", 0.045);
  c_2->SaveAs ("Plots/ElectronTrigEffs/PbPb_pt_eta.pdf");


  FormatTH2Canvas (c_3, false);
  c_3->cd ();
  TGAE* g_electronTrigEff_fcal = TEff2TGAE (t_electronTrigEff_fcal);
  g_electronTrigEff_fcal->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [TeV]");
  g_electronTrigEff_fcal->GetYaxis ()->SetTitle ("Electron trigger efficiency");
  g_electronTrigEff_fcal->GetXaxis ()->SetRangeUser (-0.5, 5);
  g_electronTrigEff_fcal->GetYaxis ()->SetRangeUser (0.4, 1);
  g_electronTrigEff_fcal->SetMarkerStyle (kOpenCircle);
  g_electronTrigEff_fcal->Draw ("AP");

  TF1* f_electronTrigEff_fcal = new TF1 ("f_electronTrigEff_fcal", "[0]+[1]*x", 0, 5);
  f_electronTrigEff_fcal->SetParameter (0, 1);
  f_electronTrigEff_fcal->SetParameter (1, 0);
  g_electronTrigEff_fcal->Fit (f_electronTrigEff_fcal, "Q0");
  f_electronTrigEff_fcal->SetLineStyle (2);
  f_electronTrigEff_fcal->SetLineColor (kGreen+2);
  f_electronTrigEff_fcal->SetLineWidth (2);
  f_electronTrigEff_fcal->Draw ("same");

  myText (0.2, 0.35, kBlack, Form ("Slope = %s TeV^{-1}", FormatMeasurement (f_electronTrigEff_fcal->GetParameter (1), f_electronTrigEff_fcal->GetParError (1), 2)), 0.04);
  myText (0.2, 0.30, kBlack, Form ("y-int. = %s", FormatMeasurement (f_electronTrigEff_fcal->GetParameter (0), f_electronTrigEff_fcal->GetParError (0), 2)), 0.04);
  myText (0.2, 0.25, kBlack, Form ("#chi^{2} / dof = %.2f / %i", f_electronTrigEff_fcal->GetChisquare (), f_electronTrigEff_fcal->GetNDF ()), 0.04);

  myText (0.62, 0.4, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.62, 0.35, kBlack, "2018 Pb+Pb, 1.7 nb^{-1}", 0.04);
  myText (0.62, 0.3, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.62, 0.25, kBlack, "#it{p}_{T}^{e} > 20 GeV", 0.04);
  myMarkerTextNoLine (0.65, 0.875, kBlack, kOpenCircle, "HLT_e15_lhloose_ion_L1EM12", 1.25, 0.04);
  c_3->SaveAs ("Plots/ElectronTrigEffs/PbPb_cent.pdf");



  TFile* outFile = new TFile (Form ("%s/Nominal/triggerEfficiencies.root", rootPath.Data ()), "recreate");
  
  t_muonTrigEff_pt[0]->Write ("t_muonTrigEff_pt_periph");
  t_muonTrigEff_pt[1]->Write ("t_muonTrigEff_pt_cent");
  t_muonTrigEff_pt[2]->Write ("t_muonTrigEff_pt_pp");
  t_muonTrigEff_pt[3]->Write ("t_muonTrigEff_pt_PbPb");
  t_muonTrigEff_eta_phi_periph->Write ("t_muonTrigEff_eta_phi_periph");
  t_muonTrigEff_eta_phi_cent->Write ("t_muonTrigEff_eta_phi_cent");
  t_muonTrigEff_eta_phi_pp->Write ("t_muonTrigEff_eta_phi_pp");
  t_muonTrigEff_eta_phi_PbPb->Write ("t_muonTrigEff_eta_phi_PbPb");
  t_muonTrigEff_fcal->Write ("t_muonTrigEff_fcal");
 
  t_electronTrigEff_pt_pp->Write ("t_electronTrigEff_pt_pp");
  t_electronTrigEff_pt_PbPb->Write ("t_electronTrigEff_pt_PbPb");
  t_electronTrigEff_eta_pp->Write ("t_electronTrigEff_eta_pp");
  t_electronTrigEff_eta_PbPb->Write ("t_electronTrigEff_eta_PbPb");
  t_electronTrigEff_pt_eta_pp->Write ("t_electronTrigEff_pt_eta_pp");
  t_electronTrigEff_pt_eta_PbPb->Write ("t_electronTrigEff_pt_eta_PbPb");
  t_electronTrigEff_fcal->Write ("t_electronTrigEff_fcal");

  outFile->Close ();



  TH2D* h2_zmumuTrigEffNum_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffNum_pt_y_pp");
  TH2D* h2_zmumuTrigEffDen_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffDen_pt_y_pp");
  TH2D* h2_zmumuTrigEffNum_pt_y_PbPb = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffNum_pt_y_PbPb");
  TH2D* h2_zmumuTrigEffDen_pt_y_PbPb = (TH2D*) (inFile)->Get ("h2_zmumuTrigEffDen_pt_y_PbPb");

  TH2D* h2_zmumuTrigEff_pt_y_pp = (TH2D*) h2_zmumuTrigEffNum_pt_y_pp->Clone ("h2_zmumuTrigEffNum_pt_y_pp");
  h2_zmumuTrigEff_pt_y_pp->Divide (h2_zmumuTrigEffDen_pt_y_pp);
  TH2D* h2_zmumuTrigEff_pt_y_PbPb = (TH2D*) h2_zmumuTrigEffNum_pt_y_PbPb->Clone ("h2_zmumuTrigEffNum_pt_y_PbPb");
  h2_zmumuTrigEff_pt_y_PbPb->Divide (h2_zmumuTrigEffDen_pt_y_PbPb);

  TH2D* h2_zeeTrigEffNum_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zeeTrigEffNum_pt_y_pp");
  TH2D* h2_zeeTrigEffDen_pt_y_pp = (TH2D*) (inFile)->Get ("h2_zeeTrigEffDen_pt_y_pp");
  TH2D* h2_zeeTrigEffNum_pt_y_PbPb = (TH2D*) (inFile)->Get ("h2_zeeTrigEffNum_pt_y_PbPb");
  TH2D* h2_zeeTrigEffDen_pt_y_PbPb = (TH2D*) (inFile)->Get ("h2_zeeTrigEffDen_pt_y_PbPb");

  TH2D* h2_zeeTrigEff_pt_y_pp = (TH2D*) h2_zeeTrigEffNum_pt_y_pp->Clone ("h2_zeeTrigEffNum_pt_y_pp");
  h2_zeeTrigEff_pt_y_pp->Divide (h2_zeeTrigEffDen_pt_y_pp);
  TH2D* h2_zeeTrigEff_pt_y_PbPb = (TH2D*) h2_zeeTrigEffNum_pt_y_PbPb->Clone ("h2_zeeTrigEffNum_pt_y_PbPb");
  h2_zeeTrigEff_pt_y_PbPb->Divide (h2_zeeTrigEffDen_pt_y_PbPb);


  TCanvas* c_4 = new TCanvas ("c_4", "", 1000, 800);
  FormatTH2Canvas (c_4, false);
  c_4->SetLogy ();

  h2_zmumuTrigEff_pt_y_pp->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zmumuTrigEff_pt_y_pp->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zmumuTrigEff_pt_y_pp->GetZaxis ()->SetTitle ("Z#rightarrow#mu#mu trigger efficiency");
  h2_zmumuTrigEff_pt_y_pp->GetXaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_pp->GetYaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_pp->GetZaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_pp->GetZaxis ()->SetRangeUser (0.5, 1);
  h2_zmumuTrigEff_pt_y_pp->Draw ("lego2");
  c_4->SaveAs ("Plots/ZmumuTrigEffs/pp_pt_y.pdf");

  TCanvas* c_5 = new TCanvas ("c_5", "", 1000, 800);
  FormatTH2Canvas (c_5, false);
  c_5->SetLogy ();

  h2_zmumuTrigEff_pt_y_PbPb->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zmumuTrigEff_pt_y_PbPb->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zmumuTrigEff_pt_y_PbPb->GetZaxis ()->SetTitle ("Z#rightarrow#mu#mu trigger efficiency");
  h2_zmumuTrigEff_pt_y_PbPb->GetXaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_PbPb->GetYaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_PbPb->GetZaxis ()->SetLabelSize (0.04);
  h2_zmumuTrigEff_pt_y_PbPb->GetZaxis ()->SetRangeUser (0.5, 1);
  h2_zmumuTrigEff_pt_y_PbPb->Draw ("lego2");
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
  h2_zeeTrigEff_pt_y_pp->GetZaxis ()->SetRangeUser (0.5, 1);
  h2_zeeTrigEff_pt_y_pp->Draw ("lego2");
  c_4->SaveAs ("Plots/ZeeTrigEffs/pp_pt_y.pdf");

  TCanvas* c_7 = new TCanvas ("c_7", "", 1000, 800);
  FormatTH2Canvas (c_7, false);
  c_7->SetLogy ();

  h2_zeeTrigEff_pt_y_PbPb->GetXaxis ()->SetTitle ("y_{Z}");
  h2_zeeTrigEff_pt_y_PbPb->GetYaxis ()->SetTitle ("#it{p}_{T}^{Z} [GeV]");
  h2_zeeTrigEff_pt_y_PbPb->GetZaxis ()->SetTitle ("Z#rightarrow ee trigger efficiency");
  h2_zeeTrigEff_pt_y_PbPb->GetXaxis ()->SetLabelSize (0.04);
  h2_zeeTrigEff_pt_y_PbPb->GetYaxis ()->SetLabelSize (0.04);
  h2_zeeTrigEff_pt_y_PbPb->GetZaxis ()->SetLabelSize (0.04);
  h2_zeeTrigEff_pt_y_PbPb->GetZaxis ()->SetRangeUser (0.5, 1);
  h2_zeeTrigEff_pt_y_PbPb->Draw ("lego2");
  c_4->SaveAs ("Plots/ZeeTrigEffs/PbPb_pt_y.pdf");

}

#endif
