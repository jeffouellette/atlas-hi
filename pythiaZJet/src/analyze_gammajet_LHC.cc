#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TLatex.h>

#include <iostream>
#include <sstream>
#include <string>

#include <math.h>

//#include "eps09_cxx/eps09.h"

#include <AtlasStyle.h>
#include <AtlasUtils.h>

#include <GlobalParams.h>
#include <Utilities.h>

using namespace atlashi;

string FormatCounts (int counts) {
  if (counts < 1000) return "";
  else if (1000 <= counts && counts < 10000) {
    string countsStr = FormatMeasurement (counts, 0, 1);
    countsStr = countsStr.substr(0, 1) + "k";
    return countsStr;
  }
  else if (10000 <= counts && counts < 100000) {
    string countsStr = FormatMeasurement (counts, 0, 2);
    countsStr = countsStr.substr(0, 2) + "k";
    return countsStr;
  }
  else if (100000 <= counts && counts < 1000000) {
    string countsStr = FormatMeasurement (counts, 0, 3);
    countsStr = countsStr.substr(0, 3) + "k";
    return countsStr;
  }
  else return "";
}

int main() {

  SetAtlasStyle();


  //TFile* f = new TFile("output/gammajet_RHIC_2400k.root", "READ");
  TFile* f = new TFile("output/gammajet_2400k.root", "READ");

  TTree* t = (TTree*) f->Get("tree");

  const int nPtBins = 200;
  const double* ptBins = logspace (2, 200, nPtBins);

  TH1D* h_gamma_pt = new TH1D ("h_gamma_pt", "", nPtBins, ptBins);
  h_gamma_pt->Sumw2 ();
  TH1D* h_jet_pt = new TH1D ("h_jet_pt", "", nPtBins, ptBins);
  h_jet_pt->Sumw2 ();

  TH1D* h_njet = new TH1D ("h_njet", "", 30, 0, 30);
  h_njet->Sumw2 ();

  TH2D* h_jet_etaphi = new TH2D ("h_jet_etaphi", "", 50, -5, 5, 80, -pi, pi);
  h_jet_etaphi->Sumw2 ();

  TH1D* h_xgammajet_5g7 = new TH1D ("h_xgammajet_5g7", "", 50, 0, 2.0);
  h_xgammajet_5g7->Sumw2 ();
  TH1D* h_xgammajet_7g10 = new TH1D ("h_xgammajet_7g10", "", 50, 0, 2.0);
  h_xgammajet_7g10->Sumw2 ();
  TH1D* h_xgammajet_10g16 = new TH1D ("h_xgammajet_10g16", "", 50, 0, 2.0);
  h_xgammajet_10g16->Sumw2 ();
  TH1D* h_xgammajet_16g25 = new TH1D ("h_xgammajet_16g25", "", 50, 0, 2.0);
  h_xgammajet_16g25->Sumw2 ();

  int code;
  int id1;
  int id2;
  float x1pdf;
  float x2pdf;
  float Q;
  bool isValence1;
  bool isValence2;

  int photon_n, jet_n;
  vector<float>* photon_pt = nullptr, *photon_eta = nullptr, *photon_phi = nullptr;
  vector<float>* jet_pt = nullptr, *jet_eta = nullptr, *jet_phi = nullptr, *jet_e = nullptr;

  t->SetBranchAddress ("code",  &code);
  t->SetBranchAddress ("id1",   &id1);
  t->SetBranchAddress ("id2",   &id2);
  t->SetBranchAddress ("x1pdf", &x1pdf);
  t->SetBranchAddress ("x2pdf", &x2pdf);
  t->SetBranchAddress ("Q",     &Q);

  t->SetBranchAddress ("photon_n",   &photon_n);
  t->SetBranchAddress ("photon_pt",  &photon_pt);
  t->SetBranchAddress ("photon_eta", &photon_eta);
  t->SetBranchAddress ("photon_phi", &photon_phi);

  t->SetBranchAddress ("jet_n",   &jet_n);
  t->SetBranchAddress ("jet_pt",  &jet_pt);
  t->SetBranchAddress ("jet_eta", &jet_eta);
  t->SetBranchAddress ("jet_phi", &jet_phi);
  t->SetBranchAddress ("jet_e",   &jet_e);

  const int nEvts = t->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {

    t->GetEntry (iEvt);

    h_njet->Fill (jet_n);

    for (int iP = 0; iP < photon_n; iP++) {

      const float ppt = photon_pt->at (iP);

      h_gamma_pt->Fill (ppt);

      for (int iJ = 0; iJ < jet_n; iJ++) {
        if (jet_pt->at (iJ) < 3) continue;

        if (DeltaPhi (photon_phi->at (iP), jet_phi->at (iJ)) < 7*pi/8) continue;

        h_jet_pt->Fill (jet_pt->at (iJ));
        h_jet_etaphi->Fill (jet_eta->at (iJ), jet_phi->at (iJ));

        if (5 < ppt && ppt < 7) {
          h_xgammajet_5g7->Fill (jet_pt->at (iJ) / ppt);
        } else if (7 < ppt && ppt < 10) {
          h_xgammajet_7g10->Fill (jet_pt->at (iJ) / ppt);
        } else if (10 < ppt && ppt < 16) {
          h_xgammajet_10g16->Fill (jet_pt->at (iJ) / ppt);
        } else if (16 < ppt && ppt < 25) {
          h_xgammajet_16g25->Fill (jet_pt->at (iJ) / ppt);
        }

      }

    }

  }

  
  Color_t colors[4] = {kBlue, kBlack, kRed, kGreen+2};

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots pT distributions for Z's
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_gamma_pt = new TCanvas("canvas_gamma_pt", "", 800, 600);
  canvas_gamma_pt->cd();

  gPad->SetLogx();
  gPad->SetLogy();

  //h_gamma_pt->Scale (1./h_gamma_pt->Integral (), "width");

  h_gamma_pt->SetLineColor (colors[0]);
  h_gamma_pt->SetMarkerColor (colors[0]);
  h_gamma_pt->GetXaxis ()->SetTitle ("#it{p}_{T}^{ #gamma} [GeV]");
  h_gamma_pt->GetYaxis ()->SetTitle ("1/N_{#gamma} dN/dx");
  h_gamma_pt->Draw ("hist");

  myText (0.26, 0.3, colors[0], "No parton", 0.06);
  myText (0.26, 0.22, colors[1], "Parton produced", 0.06);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots pT distributions for jets
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_jet_pt = new TCanvas("canvas_jet_pt", "", 800, 600);
  canvas_jet_pt->cd();

  gPad->SetLogx();
  gPad->SetLogy();

  const float sf = h_gamma_pt->Integral ();
  h_jet_pt->Scale (1./sf, "width");

  h_jet_pt->SetLineColor (colors[0]);
  h_jet_pt->SetMarkerColor (colors[0]);
  h_jet_pt->GetXaxis ()->SetTitle ("#it{p}_{T}^{ jet} [GeV]");
  h_jet_pt->GetYaxis ()->SetTitle ("1/N_{#gamma} dN/dx");
  h_jet_pt->Draw ("hist");

  myText (0.26, 0.3, colors[0], "No parton", 0.06);
  myText (0.26, 0.22, colors[1], "Parton produced", 0.06);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots pT distributions for jets
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_xgammajet = new TCanvas("canvas_xgammajet", "", 800, 600);
  canvas_xgammajet->cd();

  h_xgammajet_5g7->Scale (1. / (h_gamma_pt->Integral (h_gamma_pt->FindBin (5), h_gamma_pt->FindBin (7))), "width");
  h_xgammajet_7g10->Scale (1. / (h_gamma_pt->Integral (h_gamma_pt->FindBin (7), h_gamma_pt->FindBin (10))), "width");
  h_xgammajet_10g16->Scale (1. / (h_gamma_pt->Integral (h_gamma_pt->FindBin (10), h_gamma_pt->FindBin (16))), "width");
  h_xgammajet_16g25->Scale (1. / (h_gamma_pt->Integral (h_gamma_pt->FindBin (16), h_gamma_pt->FindBin (25))), "width");

  double max = 0, min = 1e30;
  if (h_xgammajet_5g7->GetMaximum () > max) max = h_xgammajet_5g7->GetMaximum ();
  if (h_xgammajet_7g10->GetMaximum () > max) max = h_xgammajet_7g10->GetMaximum ();
  if (h_xgammajet_10g16->GetMaximum () > max) max = h_xgammajet_10g16->GetMaximum ();
  if (h_xgammajet_16g25->GetMaximum () > max) max = h_xgammajet_16g25->GetMaximum ();
  if (h_xgammajet_5g7->GetMinimum (0) < min) min = h_xgammajet_5g7->GetMinimum (0);
  if (h_xgammajet_7g10->GetMinimum (0) < min) min = h_xgammajet_7g10->GetMinimum (0);
  if (h_xgammajet_10g16->GetMinimum (0) < min) min = h_xgammajet_10g16->GetMinimum (0);
  if (h_xgammajet_16g25->GetMinimum (0) < min) min = h_xgammajet_16g25->GetMinimum (0);
  max *= 1.3;
  min *= 0.7;

  h_xgammajet_5g7->GetYaxis ()->SetRangeUser (min, max);
  h_xgammajet_7g10->GetYaxis ()->SetRangeUser (min, max);
  h_xgammajet_10g16->GetYaxis ()->SetRangeUser (min, max);
  h_xgammajet_16g25->GetYaxis ()->SetRangeUser (min, max);

  h_xgammajet_5g7->SetLineWidth (2);
  h_xgammajet_7g10->SetLineWidth (2);
  h_xgammajet_10g16->SetLineWidth (2);
  h_xgammajet_16g25->SetLineWidth (2);

  h_xgammajet_5g7->SetLineColor (colors[0]);
  h_xgammajet_5g7->SetMarkerColor (colors[0]);
  h_xgammajet_7g10->SetLineColor (colors[1]);
  h_xgammajet_7g10->SetMarkerColor (colors[1]);
  h_xgammajet_10g16->SetLineColor (colors[2]);
  h_xgammajet_10g16->SetMarkerColor (colors[2]);
  h_xgammajet_16g25->SetLineColor (colors[3]);
  h_xgammajet_16g25->SetMarkerColor (colors[3]);

  h_xgammajet_5g7->GetXaxis ()->SetTitle ("#it{x}_{#gamma}^{ jet}");
  h_xgammajet_5g7->GetYaxis ()->SetTitle ("1/N_{#gamma} dN/dx");
  h_xgammajet_7g10->GetXaxis ()->SetTitle ("#it{x}_{#gamma}^{ jet}");
  h_xgammajet_7g10->GetYaxis ()->SetTitle ("1/N_{#gamma} dN/dx");
  h_xgammajet_10g16->GetXaxis ()->SetTitle ("#it{x}_{#gamma}^{ jet}");
  h_xgammajet_10g16->GetYaxis ()->SetTitle ("1/N_{#gamma} dN/dx");
  h_xgammajet_16g25->GetXaxis ()->SetTitle ("#it{x}_{#gamma}^{ jet}");
  h_xgammajet_16g25->GetYaxis ()->SetTitle ("1/N_{#gamma} dN/dx");

  h_xgammajet_5g7->Draw ("hist");
  h_xgammajet_7g10->Draw ("same hist");
  h_xgammajet_10g16->Draw ("same hist");
  h_xgammajet_16g25->Draw ("same hist");

  myText (0.66, 0.89, colors[0], "5 < #it{p}_{T}^{ #gamma} < 7 GeV", 0.05); 
  myText (0.66, 0.82, colors[1], "7 < #it{p}_{T}^{ #gamma} < 10 GeV", 0.05); 
  myText (0.66, 0.75, colors[2], "10 < #it{p}_{T}^{ #gamma} < 16 GeV", 0.05); 
  myText (0.66, 0.68, colors[3], "16 < #it{p}_{T}^{ #gamma} < 25 GeV", 0.05); 

  //outFile->Close();
  //if (outFile) delete outFile;
  
  return 0;
}

void analyze_gammajet () {
  main ();
} 
