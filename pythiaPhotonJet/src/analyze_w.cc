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


  TFile* f = new TFile("output/wjet_1500k.root", "READ");

  TTree* t = (TTree*) f->Get("tree");

  const int NYBINS = 3;
  double YLOBINS[] = {-2.37, -1.37, 1.56};
  double YHIBINS[] = {-1.56, 1.37, 2.37};

  const int NPTBINS = 17;
  double PTBINS[] = {25,35,45,55,65,75,85,105,125,150,175,200,250,300,350,400,470,550};
  //double alpha = pow( 10.0, 1 / 10.0 );

  //const int NPTBINS2 = 170;
  //double PTBINS2[ NPTBINS2+1 ];
  //double alpha2 = pow( 10.0, 1 / 100.0 );
  //for (int n = 0; n < NPTBINS2+1; n++) {
  //  PTBINS2[ n ] = 10 * pow( alpha2, n );
  //  std::cout << n << " " << PTBINS2[ n ] << std::endl;
  //}

  const double* xA_bins = logspace (2e-4, 2, 200);
  const double* Q_bins = linspace (1, 22500, 200);

  TH2D* hist_xA_Q_mu_minus[3];
  hist_xA_Q_mu_minus[0] = new TH2D("hist_xA_Q_mu_minus_y0", "", 200, xA_bins, 200, Q_bins);
  hist_xA_Q_mu_minus[0]->Sumw2();
  hist_xA_Q_mu_minus[1] = new TH2D("hist_xA_Q_mu_minus_y1", "", 200, xA_bins, 200, Q_bins);
  hist_xA_Q_mu_minus[1]->Sumw2();
  hist_xA_Q_mu_minus[2] = new TH2D("hist_xA_Q_mu_minus_y2", "", 200, xA_bins, 200, Q_bins);
  hist_xA_Q_mu_minus[2]->Sumw2();

  TGraph* graph_xA_Q_mu_minus[3];
  graph_xA_Q_mu_minus[0] = new TGraph();
  graph_xA_Q_mu_minus[1] = new TGraph();
  graph_xA_Q_mu_minus[2] = new TGraph();

  TH2D* hist_xA_Q_mu_plus[3];
  hist_xA_Q_mu_plus[0] = new TH2D("hist_xA_Q_mu_plus_y0", "", 200, xA_bins, 200, Q_bins);
  hist_xA_Q_mu_plus[0]->Sumw2();
  hist_xA_Q_mu_plus[1] = new TH2D("hist_xA_Q_mu_plus_y1", "", 200, xA_bins, 200, Q_bins);
  hist_xA_Q_mu_plus[1]->Sumw2();
  hist_xA_Q_mu_plus[2] = new TH2D("hist_xA_Q_mu_plus_y2", "", 200, xA_bins, 200, Q_bins);
  hist_xA_Q_mu_plus[2]->Sumw2();

  TGraph* graph_xA_Q_mu_plus[3];
  graph_xA_Q_mu_plus[0] = new TGraph();
  graph_xA_Q_mu_plus[1] = new TGraph();
  graph_xA_Q_mu_plus[2] = new TGraph();

  TH1D* hist_pT_mu_minus[3];
  hist_pT_mu_minus[0] = new TH1D("hist_pT_mu_minus_y0", "", NPTBINS, PTBINS);
  hist_pT_mu_minus[0]->Sumw2();
  hist_pT_mu_minus[1] = new TH1D("hist_pT_mu_minus_y1", "", NPTBINS, PTBINS);
  hist_pT_mu_minus[1]->Sumw2();
  hist_pT_mu_minus[2] = new TH1D("hist_pT_mu_minus_y2", "", NPTBINS, PTBINS);
  hist_pT_mu_minus[2]->Sumw2();

  TH1D* hist_pT_mu_plus[3];
  hist_pT_mu_plus[0] = new TH1D("hist_pT_mu_plus_y0", "", NPTBINS, PTBINS);
  hist_pT_mu_plus[0]->Sumw2();
  hist_pT_mu_plus[1] = new TH1D("hist_pT_mu_plus_y1", "", NPTBINS, PTBINS);
  hist_pT_mu_plus[1]->Sumw2();
  hist_pT_mu_plus[2] = new TH1D("hist_pT_mu_plus_y2", "", NPTBINS, PTBINS);
  hist_pT_mu_plus[2]->Sumw2();

  TH1D* hist_xA_mu_minus[3];
  hist_xA_mu_minus[0] = new TH1D("hist_xA_mu_minus_y0", "", 200, xA_bins);
  hist_xA_mu_minus[0]->Sumw2();
  hist_xA_mu_minus[1] = new TH1D("hist_xA_mu_minus_y1", "", 200, xA_bins);
  hist_xA_mu_minus[1]->Sumw2();
  hist_xA_mu_minus[2] = new TH1D("hist_xA_mu_minus_y2", "", 200, xA_bins);
  hist_xA_mu_minus[2]->Sumw2();

  TH1D* hist_xA_mu_plus[3];
  hist_xA_mu_plus[0] = new TH1D("hist_xA_mu_plus_y0", "", 200, xA_bins);
  hist_xA_mu_plus[0]->Sumw2();
  hist_xA_mu_plus[1] = new TH1D("hist_xA_mu_plus_y1", "", 200, xA_bins);
  hist_xA_mu_plus[1]->Sumw2();
  hist_xA_mu_plus[2] = new TH1D("hist_xA_mu_plus_y2", "", 200, xA_bins);
  hist_xA_mu_plus[2]->Sumw2();

  int id2;
  float x2pdf;
  float x1pdf;
  float Q;

  int muon_minus_n;
  float muon_minus_pt[20];
  float muon_minus_eta[20];
  float muon_minus_phi[20];

  int muon_plus_n;
  float muon_plus_pt[20];
  float muon_plus_eta[20];
  float muon_plus_phi[20];

  t->SetBranchAddress("id2",&id2);
  t->SetBranchAddress("x2pdf",&x2pdf);
  t->SetBranchAddress("x1pdf",&x1pdf);
  t->SetBranchAddress("Q",&Q);

  t->SetBranchAddress("mu_minus_n",&muon_minus_n);
  t->SetBranchAddress("mu_minus_pt",muon_minus_pt);
  t->SetBranchAddress("mu_minus_eta",muon_minus_eta);
  t->SetBranchAddress("mu_minus_phi",muon_minus_phi);

  t->SetBranchAddress("mu_plus_n",&muon_plus_n);
  t->SetBranchAddress("mu_plus_pt",muon_plus_pt);
  t->SetBranchAddress("mu_plus_eta",muon_plus_eta);
  t->SetBranchAddress("mu_plus_phi",muon_plus_phi);

  for (int e = 0; e < t->GetEntries(); e++) {

    t->GetEntry(e);

    for (int m = 0; m < muon_minus_n; m++) {
      // beam-1 is the nucleus
      // so boost is in the direction of beam-2...
      TLorentzVector muon;
      muon.SetPtEtaPhiM (muon_minus_pt[m], muon_minus_eta[m], muon_minus_phi[m], muon_mass);
      const float muon_y = muon.Rapidity() + 0.465; // boost to CM frame

      int y_bin = -1;
      if (muon_y > -2.37 && muon_y < -1.56) y_bin = 0;
      if (muon_y > -1.37 && muon_y < +1.37) y_bin = 1;
      if (muon_y > +1.56 && muon_y < +2.37) y_bin = 2;
      if (y_bin == -1) continue;

      const double factorOfTen = pow (10, y_bin - 1);

      graph_xA_Q_mu_minus[y_bin]->SetPoint(graph_xA_Q_mu_minus[y_bin]->GetN(), x2pdf, Q * factorOfTen);
      hist_xA_Q_mu_minus[y_bin]->Fill (x2pdf, Q);
      hist_pT_mu_minus[y_bin]->Fill (muon_minus_pt[m]);
      hist_xA_mu_minus[y_bin]->Fill (x2pdf);
    }

    for (int m = 0; m < muon_plus_n; m++) {
      // beam-1 is the nucleus
      // so boost is in the direction of beam-2...
      TLorentzVector muon;
      muon.SetPtEtaPhiM (muon_plus_pt[m], muon_plus_eta[m], muon_plus_phi[m], muon_mass);
      const float muon_y = muon.Rapidity() + 0.465; // boost to CM frame

      int y_bin = -1;
      if (muon_y > -2.37 && muon_y < -1.56) y_bin = 0;
      if (muon_y > -1.37 && muon_y < +1.37) y_bin = 1;
      if (muon_y > +1.56 && muon_y < +2.37) y_bin = 2;
      if (y_bin == -1) continue;

      const double factorOfTen = pow (10, y_bin - 1);

      graph_xA_Q_mu_plus[y_bin]->SetPoint(graph_xA_Q_mu_plus[y_bin]->GetN(), x2pdf, Q * factorOfTen);
      hist_xA_Q_mu_minus[y_bin]->Fill (x2pdf, Q);
      hist_pT_mu_plus[y_bin]->Fill (muon_plus_pt[m]);
      hist_xA_mu_plus[y_bin]->Fill (x2pdf);
    }

  }

  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots pT distributions for mu+ and mu-
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_mu_pt = new TCanvas("canvas_mu_pt", "", 800, 600);
  canvas_mu_pt->cd();

  gPad->SetLogx();
  gPad->SetLogy();

  Color_t colors[3] = {kBlue, kBlack, kRed};

  for (int i = 0; i < 3; i++) {
    hist_pT_mu_plus[i]->Scale(1., "width");
    hist_pT_mu_plus[i]->SetMarkerColor(colors[i]);
    hist_pT_mu_plus[i]->SetLineColor(colors[i]);
    hist_pT_mu_plus[i]->SetMarkerStyle(kFullCircle);
    hist_pT_mu_plus[i]->GetXaxis()->SetTitle ("#it{p}_{T}^{#mu} #left[GeV#right]");
    hist_pT_mu_plus[i]->GetYaxis()->SetTitle ("Counts / GeV");
    hist_pT_mu_plus[i]->GetYaxis()->SetRangeUser(1e-2, 1e3);
    if (i == 0) hist_pT_mu_plus[i]->Draw("e1");
    else hist_pT_mu_plus[i]->Draw("e1 same");
  }
  for (int i = 0; i < 3; i++) {
    hist_pT_mu_minus[i]->Scale(1., "width");
    hist_pT_mu_minus[i]->SetMarkerColor(colors[i]);
    hist_pT_mu_minus[i]->SetLineColor(colors[i]);
    hist_pT_mu_minus[i]->SetMarkerStyle(kOpenCircle);
    hist_pT_mu_minus[i]->GetXaxis()->SetTitle ("#it{p}_{T}^{#mu} #left[GeV#right]");
    hist_pT_mu_minus[i]->GetYaxis()->SetTitle ("Counts / GeV");
    hist_pT_mu_minus[i]->GetYaxis()->SetRangeUser(1e-2, 1e3);
    hist_pT_mu_minus[i]->Draw("e1 same");
  }
  for (int i = 0; i < 3; i++) {
    myMarkerText (0.65, 0.85-0.07*i, colors[i], kFullCircle, Form("%g < #it{y} < %g", YLOBINS[i], YHIBINS[i]), 1.25, 0.04);
  }

  canvas_mu_pt->SaveAs ("muonPtSpectrum.pdf");


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots xA 1-d distributions for mu+
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_xA_mu_plus = new TCanvas("canvas_xA_mu_plus", "", 800, 600);
  canvas_xA_mu_plus->cd();

  gPad->SetLogx();

  int nplus = 0;
  for (int i = 0; i < 3; i++) {
    hist_xA_mu_plus[i]->Rebin(2);
    const double integral = hist_xA_mu_plus[i]->Integral();
    nplus += (int)integral;
    hist_xA_mu_plus[i]->Scale (1./integral);
    hist_xA_mu_plus[i]->SetMarkerColor(colors[i]);
    hist_xA_mu_plus[i]->SetLineColor(colors[i]);
    hist_xA_mu_plus[i]->GetXaxis()->SetTitle ("#it{x}_{A}");
    hist_xA_mu_plus[i]->GetYaxis()->SetTitle ("Counts / Total"); 
    hist_xA_mu_plus[i]->GetYaxis()->SetRangeUser (0, 0.1);
    if (i == 0) hist_xA_mu_plus[i]->Draw("hist");
    else hist_xA_mu_plus[i]->Draw("same hist");
  }

  myText (0.2, 0.9, kBlack, "Pythia8 #it{pp} #it{W}^{#plus} #rightarrow #mu^{#plus} + #nu_{#mu}", 0.04);
  myText (0.2, 0.84, kBlack, "#it{p}_{T}^{#mu} > 25 GeV", 0.04);
  myText (0.2, 0.78, kBlack, Form("1.5M events, %s #mu^{#plus}'s", FormatCounts(nplus).c_str()), 0.04);

  for (int i = 0; i < 3; i++) {
    //myMarkerText (0.65, 0.78+0.06*i, kBlack, Form("%g < #it{y} < %g (#it{Q} #times %g)", YLOBINS[i], YHIBINS[i], pow(10, i-1)), 0.04);
    myMarkerText (0.6, 0.78+0.06*i, colors[i], kFullCircle, Form("%g < #it{y} < %g (#it{Q} #times %g)", YLOBINS[i], YHIBINS[i], pow(10, i-1)), 1.25, 0.04);
  }

  canvas_xA_mu_plus->SaveAs ("muPlus_xA.pdf");


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots xA 1-d distributions for mu+
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_xA_mu_minus = new TCanvas("canvas_xA_mu_minus", "", 800, 600);
  canvas_xA_mu_minus->cd();

  gPad->SetLogx();

  int nminus = 0;
  for (int i = 0; i < 3; i++) {
    hist_xA_mu_minus[i]->Rebin(2);
    const double integral = hist_xA_mu_minus[i]->Integral();
    nminus += (double)integral;
    hist_xA_mu_minus[i]->Scale (1./integral);
    hist_xA_mu_minus[i]->SetMarkerColor(colors[i]);
    hist_xA_mu_minus[i]->SetLineColor(colors[i]);
    hist_xA_mu_minus[i]->GetXaxis()->SetTitle ("#it{x}_{A}");
    hist_xA_mu_minus[i]->GetYaxis()->SetTitle ("Counts / Total"); 
    hist_xA_mu_minus[i]->GetYaxis()->SetRangeUser (0, 0.1);
    if (i == 0) hist_xA_mu_minus[i]->Draw("hist");
    else hist_xA_mu_minus[i]->Draw("same hist");
  }

  myText (0.2, 0.9, kBlack, "Pythia8 #it{pp} #it{W}^{#minus} #rightarrow #mu^{#minus} + #bar{#nu}_{#mu}", 0.04);
  myText (0.2, 0.84, kBlack, "#it{p}_{T}^{#mu} > 25 GeV", 0.04);
  myText (0.2, 0.78, kBlack, Form("1.5M events, %s #mu^{#minus}'s", FormatCounts(nminus).c_str()), 0.04);

  for (int i = 0; i < 3; i++) {
    //myText (0.65, 0.78+0.06*i, kBlack, Form("%g < #it{y} < %g (#it{Q} #times %g)", YLOBINS[i], YHIBINS[i], pow(10, i-1)), 0.04);
    myMarkerText (0.6, 0.78+0.06*i, colors[i], kFullCircle, Form("%g < #it{y} < %g (#it{Q} #times %g)", YLOBINS[i], YHIBINS[i], pow(10, i-1)), 1.25, 0.04);
  }

  canvas_xA_mu_minus->SaveAs ("muMinus_xA.pdf");


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots xA vs Q 2-d distributions for mu+
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_xA_Q_mu_plus = new TCanvas("canvas_xA_Q_mu_plus", "", 800, 600);
  canvas_xA_Q_mu_plus->cd();
  //FormatTH2Canvas (canvas_xA_Q_mu_plus, false);
  gPad->SetLogx ();
  gPad->SetLogy ();

  nplus = 0;
  for (int i = 0; i < 3; i++) {
    graph_xA_Q_mu_plus[i]->SetMarkerColor(colors[i]);
    graph_xA_Q_mu_plus[i]->SetMarkerStyle(kDot);
    graph_xA_Q_mu_plus[i]->GetXaxis()->SetTitle ("#it{x}_{A}");
    graph_xA_Q_mu_plus[i]->GetYaxis()->SetTitle ("#sqrt{#it{Q}^{2}} #left[GeV#right]"); 
    graph_xA_Q_mu_plus[i]->GetXaxis()->SetLimits (1e-4, 2);
    graph_xA_Q_mu_plus[i]->GetYaxis()->SetRangeUser (1, 22500);
    nplus += graph_xA_Q_mu_plus[i]->GetN();
  }
  graph_xA_Q_mu_plus[0]->Draw("ap");
  graph_xA_Q_mu_plus[1]->Draw("p"); 
  graph_xA_Q_mu_plus[2]->Draw("p");

  myText (0.2, 0.9, kBlack, "Pythia8 #it{pp} #it{W}^{#plus} #rightarrow #mu^{#plus} + #nu_{#mu}", 0.04);
  myText (0.2, 0.84, kBlack, "#it{p}_{T}^{#mu} > 25 GeV", 0.04);
  myText (0.2, 0.78, kBlack, Form("1.5M events, %s #mu^{#plus}'s", FormatCounts(nplus).c_str()), 0.04);

  for (int i = 0; i < 3; i++) {
    myText (0.55, 0.26+0.18*i, kBlack, Form("%g < #it{y} < %g (#it{Q} #times %g)", YLOBINS[i], YHIBINS[i], pow(10, i-1)), 0.04);
    //myMarkerText (0.65, 0.78+0.06*i, colors[i], kFullCircle, Form("%g < #it{y} < %g (#it{Q} #times %g)", YLOBINS[i], YHIBINS[i], pow(10, i-1)), 1.25, 0.04);
  }

  canvas_xA_Q_mu_plus->SaveAs ("muPlus_xA_Q.pdf");


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots xA vs Q 2-d distributions for mu-
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_xA_Q_mu_minus = new TCanvas("canvas_xA_Q_mu_minus", "", 800, 600);
  canvas_xA_Q_mu_minus->cd();

  gPad->SetLogx();
  gPad->SetLogy();

  nminus = 0;
  for (int i = 0; i < 3; i++) {
    graph_xA_Q_mu_minus[i]->SetMarkerColor(colors[i]);
    graph_xA_Q_mu_minus[i]->SetMarkerStyle(kDot);
    graph_xA_Q_mu_minus[i]->GetXaxis()->SetTitle ("#it{x}_{A}");
    graph_xA_Q_mu_minus[i]->GetYaxis()->SetTitle ("#sqrt{#it{Q}^{2}} #left[GeV#right]"); 
    graph_xA_Q_mu_minus[i]->GetXaxis()->SetLimits (1e-4, 2);
    graph_xA_Q_mu_minus[i]->GetYaxis()->SetRangeUser (1, 22500);
    nminus += graph_xA_Q_mu_minus[i]->GetN();
  }

  graph_xA_Q_mu_minus[0]->Draw("ap");
  graph_xA_Q_mu_minus[1]->Draw("p"); 
  graph_xA_Q_mu_minus[2]->Draw("p");

  myText (0.2, 0.9, kBlack, "Pythia8 #it{pp} #it{W}^{#minus} #rightarrow #mu^{#minus} + #bar{#nu}_{#mu}", 0.04);
  myText (0.2, 0.84, kBlack, "#it{p}_{T}^{#mu} > 25 GeV", 0.04);
  myText (0.2, 0.78, kBlack, Form("1.5M events, %s #mu^{#minus}'s", FormatCounts(nminus).c_str()), 0.04);
  for (int i = 0; i < 3; i++) {
    myText (0.55, 0.26+0.18*i, kBlack, Form("%g < #it{y} < %g (#it{Q} #times %g)", YLOBINS[i], YHIBINS[i], pow(10, i-1)), 0.04);
    //myMarkerText (0.65, 0.78+0.06*i, colors[i], kFullCircle, Form("%g < #it{y} < %g (#it{Q} #times %g)", YLOBINS[i], YHIBINS[i], pow(10, i-1)), 1.25, 0.04);
  }

  canvas_xA_Q_mu_minus->SaveAs ("muMinus_xA_Q.pdf");


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Save canvases, histograms, graphs to a root file
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TFile* outFile = new TFile ("muons.root", "RECREATE");
  outFile->cd();

  canvas_mu_pt->Write();
  canvas_xA_Q_mu_plus->Write();
  canvas_xA_Q_mu_minus->Write();

  for (int i = 0; i < 3; i++) {
    hist_pT_mu_plus[i]->Write();
    hist_pT_mu_minus[i]->Write();
    graph_xA_Q_mu_plus[i]->Write(Form("graph_xA_Q_mu_plus_y%i", i));
    graph_xA_Q_mu_minus[i]->Write(Form("graph_xA_Q_mu_minus_y%i", i));
  }


  outFile->Close();
  if (outFile) delete outFile;
  
  return 0;
}
