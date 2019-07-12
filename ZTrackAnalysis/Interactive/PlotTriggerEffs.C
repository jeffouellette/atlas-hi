#ifndef __PlotTriggerEffs_C__
#define __PlotTriggerEffs_C__

#include "Params.h"
#include "ZTrackUtilities.h"

#include <ArrayTemplates.h>

#include <TEfficiency.h>

#include <iostream>

typedef TGraphAsymmErrors TGAE;

void PlotTriggerEffs () {

  SetupDirectories ("TriggerEffTagAndProbe/", "ZTrackAnalysis/");
  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");

  TH2D* h2_muonTrigEffNum = (TH2D*) inFile->Get ("h2_muonTrigEffNum");
  TH2D* h2_muonTrigEffDen = (TH2D*) inFile->Get ("h2_muonTrigEffDen");
  TH2D* h2_electronTrigEffNum = (TH2D*) inFile->Get ("h2_electronTrigEffNum");
  TH2D* h2_electronTrigEffDen = (TH2D*) inFile->Get ("h2_electronTrigEffDen");

  TEfficiency* t_muonTrigEffPt = new TEfficiency (*(h2_muonTrigEffNum->ProjectionY ()), *(h2_muonTrigEffDen->ProjectionY ()));
  TEfficiency* t_electronTrigEffPt = new TEfficiency (*(h2_electronTrigEffNum->ProjectionY ()), *(h2_electronTrigEffDen->ProjectionY ()));

  TEfficiency* t_muonTrigEffFCal = new TEfficiency (*(h2_muonTrigEffNum->ProjectionX ("_px", h2_muonTrigEffNum->GetXaxis ()->FindBin (20), h2_muonTrigEffNum->GetNbinsY ())), *(h2_muonTrigEffDen->ProjectionX ("_px", h2_muonTrigEffDen->GetXaxis ()->FindBin (20), h2_muonTrigEffDen->GetNbinsY ())));
  TEfficiency* t_electronTrigEffFCal = new TEfficiency (*(h2_electronTrigEffNum->ProjectionX ("_px", h2_electronTrigEffNum->GetXaxis ()->FindBin (20), h2_electronTrigEffNum->GetNbinsY ())), *(h2_electronTrigEffDen->ProjectionX ("_px", h2_electronTrigEffNum->GetXaxis ()->FindBin (20), h2_electronTrigEffDen->GetNbinsY ())));

  TCanvas* c = new TCanvas ("c_muonTrigEffPt", "", 800, 600);
  TGAE* g_muonTrigEffPt = t_muonTrigEffPt->CreateGraph ();
  g_muonTrigEffPt->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
  g_muonTrigEffPt->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEffPt->GetYaxis ()->SetRangeUser (0, 1.08);
  g_muonTrigEffPt->Draw ("AP");

  myText (0.24, 0.9, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.24, 0.85, kBlack, "2018 Pb+Pb, 1.4 nb^{-1}", 0.04);
  myText (0.24, 0.8, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myMarkerText (0.65, 0.9, kBlack, kFullCircle, "HLT_mu14", 1.25, 0.04);

  c->SaveAs ("Plots/muonTrigEffPt.pdf");

  c = new TCanvas ("c_muonTrigEffFCal", "", 800, 600);
  TGAE* g_muonTrigEffFCal = t_muonTrigEffFCal->CreateGraph ();
  g_muonTrigEffFCal->GetXaxis ()->SetTitle ("#Sigma#it{E}_{T}^{FCal} [TeV]");
  g_muonTrigEffFCal->GetYaxis ()->SetTitle ("Muon trigger efficiency");
  g_muonTrigEffFCal->GetYaxis ()->SetRangeUser (0, 1.08);
  g_muonTrigEffFCal->Draw ("AP");

  myText (0.62, 0.4, kBlack, "#bf{#it{ATLAS}} Internal", 0.04);
  myText (0.62, 0.35, kBlack, "2018 Pb+Pb, 1.4 nb^{-1}", 0.04);
  myText (0.62, 0.3, kBlack, "#sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.62, 0.25, kBlack, "#it{p}_{T}^{#mu} > 20 GeV", 0.04);
  myMarkerText (0.65, 0.9, kBlack, kFullCircle, "HLT_mu14", 1.25, 0.04);

  c->SaveAs ("Plots/muonTrigEffFCal.pdf");

  c = new TCanvas ("c_electronTrigEffPt", "", 800, 600);
  t_electronTrigEffPt->Draw ();

  c = new TCanvas ("c_electronTrigEffFCal", "", 800, 600);
  t_electronTrigEffFCal->Draw ();
}

#endif
