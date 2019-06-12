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

bool ProducesParton (const int code) {
  return (code == 221 || code == 231 || code == 232) && (code != 241 && code != 242 && code != 243 && code != 244);
}

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


  //TFile* f = new TFile("output/zjet_1_50k.root", "READ");
  TFile* f = new TFile("output/zjet_240k.root", "READ");

  TTree* t = (TTree*) f->Get("tree");

  const int nPtBins = 200;
  const double* ptBins = logspace (2, 200, nPtBins);

  TH1D* h_noparton_z_pt = new TH1D ("h_noparton_z_pt", "", nPtBins, ptBins);
  h_noparton_z_pt->Sumw2 ();
  TH1D* h_noparton_jet_pt = new TH1D ("h_noparton_jet_pt", "", nPtBins, ptBins);
  h_noparton_jet_pt->Sumw2 ();
  TH1D* h_parton_z_pt = new TH1D ("h_parton_z_pt", "", nPtBins, ptBins);
  h_parton_z_pt->Sumw2 ();
  TH1D* h_parton_jet_pt = new TH1D ("h_parton_jet_pt", "", nPtBins, ptBins);
  h_parton_jet_pt->Sumw2 ();

  TH1D* h_noparton_njet = new TH1D ("h_noparton_njet", "", 30, 0, 30);
  h_noparton_njet->Sumw2 ();
  TH1D* h_parton_njet = new TH1D ("h_parton_njet", "", 30, 0, 30);
  h_parton_njet->Sumw2 ();

  TH2D* h_noparton_jet_etaphi = new TH2D ("h_noparton_jet_etaphi", "", 50, -5, 5, 80, -pi, pi);
  h_noparton_jet_etaphi->Sumw2 ();
  TH2D* h_parton_jet_etaphi = new TH2D ("h_parton_jet_etaphi", "", 50, -5, 5, 80, -pi, pi);
  h_parton_jet_etaphi->Sumw2 ();

  TH1D* h_jet_25z = new TH1D ("h_jet_25z", "", 25, 0, 100);
  TH1D* h_jet_30z = new TH1D ("h_jet_30z", "", 25, 0, 100);
  TH1D* h_jet_40z = new TH1D ("h_jet_40z", "", 25, 0, 100);
  TH1D* h_jet_50z = new TH1D ("h_jet_50z", "", 25, 0, 100);
  h_jet_25z->Sumw2 ();
  h_jet_30z->Sumw2 ();
  h_jet_40z->Sumw2 ();
  h_jet_50z->Sumw2 ();

  TH1D* h_trk_25z = new TH1D ("h_trk_25z", "", 39, 2, 80);
  TH1D* h_trk_30z = new TH1D ("h_trk_30z", "", 39, 2, 80);
  TH1D* h_trk_40z = new TH1D ("h_trk_40z", "", 39, 2, 80);
  TH1D* h_trk_50z = new TH1D ("h_trk_50z", "", 39, 2, 80);
  h_trk_25z->Sumw2 ();
  h_trk_30z->Sumw2 ();
  h_trk_40z->Sumw2 ();
  h_trk_50z->Sumw2 ();

  TH1D* h_xzjet_5z10 = new TH1D ("h_xzjet_5z10", "", 50, 0, 2.0);
  h_xzjet_5z10->Sumw2 ();
  TH1D* h_xzjet_10z20 = new TH1D ("h_xzjet_10z20", "", 50, 0, 2.0);
  h_xzjet_10z20->Sumw2 ();
  TH1D* h_xzjet_20z40 = new TH1D ("h_xzjet_20z40", "", 50, 0, 2.0);
  h_xzjet_20z40->Sumw2 ();
  TH1D* h_xzjet_40z80 = new TH1D ("h_xzjet_40z80", "", 50, 0, 2.0);
  h_xzjet_40z80->Sumw2 ();

  int code;
  int id1;
  int id2;
  float x1pdf;
  float x2pdf;
  float Q;
  bool isValence1;
  bool isValence2;

  int z_n, jet_n, part_n;
  vector<float>* z_pt = nullptr, *z_eta = nullptr, *z_phi = nullptr, *z_m = nullptr;
  vector<float>* jet_pt = nullptr, *jet_eta = nullptr, *jet_phi = nullptr, *jet_e = nullptr;
  vector<float>* part_pt = nullptr, *part_eta = nullptr, *part_phi = nullptr;

  t->SetBranchAddress ("code",  &code);
  t->SetBranchAddress ("id1",   &id1);
  t->SetBranchAddress ("id2",   &id2);
  t->SetBranchAddress ("x1pdf", &x1pdf);
  t->SetBranchAddress ("x2pdf", &x2pdf);
  t->SetBranchAddress ("Q",     &Q);

  t->SetBranchAddress ("z_n",   &z_n);
  t->SetBranchAddress ("z_pt",  &z_pt);
  t->SetBranchAddress ("z_eta", &z_eta);
  t->SetBranchAddress ("z_phi", &z_phi);
  t->SetBranchAddress ("z_m",   &z_m);

  t->SetBranchAddress ("jet_r04_n",   &jet_n);
  t->SetBranchAddress ("jet_r04_pt",  &jet_pt);
  t->SetBranchAddress ("jet_r04_eta", &jet_eta);
  t->SetBranchAddress ("jet_r04_phi", &jet_phi);
  t->SetBranchAddress ("jet_r04_e",   &jet_e);

  t->SetBranchAddress ("part_n",    &part_n);
  t->SetBranchAddress ("part_pt",   &part_pt);
  t->SetBranchAddress ("part_eta",  &part_eta);
  t->SetBranchAddress ("part_phi",  &part_phi);

  const int nEvts = t->GetEntries ();
  for (int iEvt = 0; iEvt < nEvts; iEvt++) {

    t->GetEntry (iEvt);

    const bool producesParton = ProducesParton (code);

    if (!producesParton) {
      h_noparton_njet->Fill (jet_n);
    } else {
      h_parton_njet->Fill (jet_n);
    }

    for (int iZ = 0; iZ < z_n; iZ++) {

      const float zpt = z_pt->at (iZ);

      if (!producesParton) {
        h_noparton_z_pt->Fill (zpt);
      } else {
        h_parton_z_pt->Fill (zpt);
      }

      TLorentzVector highestPtJet;
      highestPtJet.SetPxPyPzE (0, 0, 0, 0);
      for (int iJ = 0; iJ < jet_n; iJ++) {
        if (jet_pt->at (iJ) < 3) continue;

        if (DeltaPhi (z_phi->at (iZ), jet_phi->at (iJ)) < 7*pi/8) continue;

        if (jet_pt->at (iJ) > highestPtJet.Pt ())
          highestPtJet.SetPtEtaPhiE (jet_pt->at (iJ), jet_eta->at (iJ), jet_phi->at (iJ), jet_e->at (iJ));

        if (!producesParton) {
          h_noparton_jet_pt->Fill (jet_pt->at (iJ));
          h_noparton_jet_etaphi->Fill (jet_eta->at (iJ), jet_phi->at (iJ));
        } else {
          h_parton_jet_pt->Fill (jet_pt->at (iJ));
          h_parton_jet_etaphi->Fill (jet_eta->at (iJ), jet_phi->at (iJ));
        }

        if (5 < zpt && zpt < 10) {
          h_xzjet_5z10->Fill (jet_pt->at (iJ) / zpt);
        } else if (10 < zpt && zpt < 20) {
          h_xzjet_10z20->Fill (jet_pt->at (iJ) / zpt);
        } else if (20 < zpt && zpt < 40) {
          h_xzjet_20z40->Fill (jet_pt->at (iJ) / zpt);
        } else if (40 < zpt && zpt < 80) {
          h_xzjet_40z80->Fill (jet_pt->at (iJ) / zpt);
        }

        if (25 < zpt)
          h_jet_25z->Fill (jet_pt->at (iJ));
        if (30 < zpt)
          h_jet_30z->Fill (jet_pt->at (iJ));
        if (40 < zpt)
          h_jet_40z->Fill (jet_pt->at (iJ));
        if (50 < zpt)
          h_jet_50z->Fill (jet_pt->at (iJ));

      }

      float max_pt = 2;
      for (int iP = 0; iP < part_pt->size (); iP++) {
        if (part_pt->at (iP) < 2) continue;

        if (DeltaPhi (z_phi->at (iZ), part_phi->at (iP)) < 7*pi/8) continue;

        if (DeltaR (part_eta->at (iP), highestPtJet.Eta (), part_phi->at (iP), highestPtJet.Phi ()) < 0.4)
          continue;

        if (part_pt->at (iP) > max_pt) {
          max_pt = part_pt->at (iP);
        }
      }

      if (25 < zpt)
        h_trk_25z->Fill (max_pt);
      if (30 < zpt)
        h_trk_30z->Fill (max_pt);
      if (40 < zpt)
        h_trk_40z->Fill (max_pt);
      if (50 < zpt)
        h_trk_50z->Fill (max_pt);

    }

  }

  
  Color_t colors[4] = {kBlue, kBlack, kRed, kGreen+2};

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots pT distributions for Z's
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_z_pt = new TCanvas("canvas_z_pt", "", 800, 600);
  canvas_z_pt->cd();

  gPad->SetLogx();
  gPad->SetLogy();

  //h_noparton_z_pt->Scale (1./h_noparton_z_pt->Integral (), "width");
  //h_parton_z_pt->Scale (1./h_parton_z_pt->Integral (), "width");

  h_noparton_z_pt->SetLineColor (colors[0]);
  h_parton_z_pt->SetLineColor (colors[1]);
  h_noparton_z_pt->SetMarkerColor (colors[0]);
  h_parton_z_pt->SetMarkerColor (colors[1]);
  h_noparton_z_pt->GetXaxis ()->SetTitle ("#it{p}_{T}^{ Z} [GeV]");
  h_parton_z_pt->GetXaxis ()->SetTitle ("#it{p}_{T}^{ Z} [GeV]");
  h_noparton_z_pt->GetYaxis ()->SetTitle ("1/N_{Z} dN/dx");
  h_parton_z_pt->GetYaxis ()->SetTitle ("1/N_{Z} dN/dx");
  h_noparton_z_pt->Draw ("hist");
  h_parton_z_pt->Draw ("same hist");

  myText (0.26, 0.3, colors[0], "No parton", 0.06);
  myText (0.26, 0.22, colors[1], "Parton produced", 0.06);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots pT distributions for jets
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_jet_pt = new TCanvas("canvas_jet_pt", "", 800, 600);
  canvas_jet_pt->cd();

  //gPad->SetLogx();
  gPad->SetLogy();

  //h_jet_25z->Scale (1./t->GetEntries ("z_pt>25"));
  //h_jet_30z->Scale (1./t->GetEntries ("z_pt>30"));
  //h_jet_40z->Scale (1./t->GetEntries ("z_pt>40"));
  //h_jet_50z->Scale (1./t->GetEntries ("z_pt>50"));

  h_jet_25z->SetLineColor (colors[0]);
  h_jet_25z->SetMarkerColor (colors[0]);
  h_jet_30z->SetLineColor (colors[1]);
  h_jet_30z->SetMarkerColor (colors[1]);
  h_jet_40z->SetLineColor (colors[2]);
  h_jet_40z->SetMarkerColor (colors[2]);
  h_jet_50z->SetLineColor (colors[3]);
  h_jet_50z->SetMarkerColor (colors[3]);
  h_jet_25z->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
  h_jet_30z->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
  h_jet_40z->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
  h_jet_50z->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
  h_jet_25z->GetYaxis ()->SetTitle ("N_{jet} / N_{Z}");
  h_jet_30z->GetYaxis ()->SetTitle ("N_{jet} / N_{Z}");
  h_jet_40z->GetYaxis ()->SetTitle ("N_{jet} / N_{Z}");
  h_jet_50z->GetYaxis ()->SetTitle ("N_{jet} / N_{Z}");
  h_jet_25z->Draw ("e1");
  h_jet_30z->Draw ("same e1");
  h_jet_40z->Draw ("same e1");
  h_jet_50z->Draw ("same e1");

  myText (0.7, 0.90, colors[0], "#it{p}_{T}^{Z} > 25 GeV (nom.)", 0.04);
  myText (0.7, 0.85, colors[1], "#it{p}_{T}^{Z} > 30 GeV", 0.04);
  myText (0.7, 0.80, colors[2], "#it{p}_{T}^{Z} > 40 GeV", 0.04);
  myText (0.7, 0.75, colors[3], "#it{p}_{T}^{Z} > 50 GeV", 0.04);

  //const float sf = h_parton_z_pt->Integral () + h_noparton_z_pt->Integral ();
  //h_noparton_jet_pt->Scale (1./sf, "width");
  //h_parton_jet_pt->Scale (1./sf, "width");

  //h_noparton_jet_pt->SetLineColor (colors[0]);
  //h_parton_jet_pt->SetLineColor (colors[1]);
  //h_noparton_jet_pt->SetMarkerColor (colors[0]);
  //h_parton_jet_pt->SetMarkerColor (colors[1]);
  //h_noparton_jet_pt->GetXaxis ()->SetTitle ("#it{p}_{T}^{ jet} [GeV]");
  //h_parton_jet_pt->GetXaxis ()->SetTitle ("#it{p}_{T}^{ jet} [GeV]");
  //h_noparton_jet_pt->GetYaxis ()->SetTitle ("1/N_{Z} dN/dx");
  //h_parton_jet_pt->GetYaxis ()->SetTitle ("1/N_{Z} dN/dx");
  //h_noparton_jet_pt->Draw ("hist");
  //h_parton_jet_pt->Draw ("same hist");

  //myText (0.26, 0.3, colors[0], "No parton", 0.06);
  //myText (0.26, 0.22, colors[1], "Parton produced", 0.06);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots pT distributions for charged hadrons
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_trk_pt = new TCanvas("canvas_trk_pt", "", 800, 600);
  canvas_trk_pt->cd();

  gPad->SetLogy();

  //h_trk_25z->Scale (1./t->GetEntries ("z_pt>25"));
  //h_trk_30z->Scale (1./t->GetEntries ("z_pt>30"));
  //h_trk_40z->Scale (1./t->GetEntries ("z_pt>40"));
  //h_trk_50z->Scale (1./t->GetEntries ("z_pt>50"));

  h_trk_25z->SetLineColor (colors[0]);
  h_trk_25z->SetMarkerColor (colors[0]);
  h_trk_30z->SetLineColor (colors[1]);
  h_trk_30z->SetMarkerColor (colors[1]);
  h_trk_40z->SetLineColor (colors[2]);
  h_trk_40z->SetMarkerColor (colors[2]);
  h_trk_50z->SetLineColor (colors[3]);
  h_trk_50z->SetMarkerColor (colors[3]);
  h_trk_25z->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
  h_trk_30z->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
  h_trk_40z->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
  h_trk_50z->GetXaxis ()->SetTitle ("#it{p}_{T} [GeV]");
  h_trk_25z->GetYaxis ()->SetTitle ("N_{trk} / N_{Z}");
  h_trk_30z->GetYaxis ()->SetTitle ("N_{trk} / N_{Z}");
  h_trk_40z->GetYaxis ()->SetTitle ("N_{trk} / N_{Z}");
  h_trk_50z->GetYaxis ()->SetTitle ("N_{trk} / N_{Z}");
  h_trk_25z->Draw ("e1");
  h_trk_30z->Draw ("same e1");
  h_trk_40z->Draw ("same e1");
  h_trk_50z->Draw ("same e1");

  myText (0.7, 0.90, colors[0], "#it{p}_{T}^{Z} > 25 GeV (nom.)", 0.04);
  myText (0.7, 0.85, colors[1], "#it{p}_{T}^{Z} > 30 GeV", 0.04);
  myText (0.7, 0.80, colors[2], "#it{p}_{T}^{Z} > 40 GeV", 0.04);
  myText (0.7, 0.75, colors[3], "#it{p}_{T}^{Z} > 50 GeV", 0.04);
  myText (0.4, 0.90, kBlack, "#Delta#phi > 7#pi/8", 0.04);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Plots pT distributions for jets
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* canvas_xzjet = new TCanvas("canvas_xzjet", "", 800, 600);
  canvas_xzjet->cd();

  h_xzjet_5z10->Scale (1. / (h_noparton_z_pt->Integral (h_noparton_z_pt->FindBin (5), h_noparton_z_pt->FindBin (10)) + h_parton_z_pt->Integral (h_parton_z_pt->FindBin (5), h_parton_z_pt->FindBin (10))), "width");
  h_xzjet_10z20->Scale (1. / (h_noparton_z_pt->Integral (h_noparton_z_pt->FindBin (10), h_noparton_z_pt->FindBin (20)) + h_parton_z_pt->Integral (h_parton_z_pt->FindBin (10), h_parton_z_pt->FindBin (20))), "width");
  h_xzjet_20z40->Scale (1. / (h_noparton_z_pt->Integral (h_noparton_z_pt->FindBin (20), h_noparton_z_pt->FindBin (40)) + h_parton_z_pt->Integral (h_parton_z_pt->FindBin (20), h_parton_z_pt->FindBin (40))), "width");
  h_xzjet_40z80->Scale (1. / (h_noparton_z_pt->Integral (h_noparton_z_pt->FindBin (40), h_noparton_z_pt->FindBin (80)) + h_parton_z_pt->Integral (h_parton_z_pt->FindBin (40), h_parton_z_pt->FindBin (80))), "width");

  double max = 0, min = 1e30;
  if (h_xzjet_5z10->GetMaximum () > max) max = h_xzjet_5z10->GetMaximum ();
  if (h_xzjet_10z20->GetMaximum () > max) max = h_xzjet_10z20->GetMaximum ();
  if (h_xzjet_20z40->GetMaximum () > max) max = h_xzjet_20z40->GetMaximum ();
  if (h_xzjet_40z80->GetMaximum () > max) max = h_xzjet_40z80->GetMaximum ();
  if (h_xzjet_5z10->GetMinimum (0) < min) min = h_xzjet_5z10->GetMinimum (0);
  if (h_xzjet_10z20->GetMinimum (0) < min) min = h_xzjet_10z20->GetMinimum (0);
  if (h_xzjet_20z40->GetMinimum (0) < min) min = h_xzjet_20z40->GetMinimum (0);
  if (h_xzjet_40z80->GetMinimum (0) < min) min = h_xzjet_40z80->GetMinimum (0);
  max *= 1.3;
  min *= 0.7;

  h_xzjet_5z10->GetYaxis ()->SetRangeUser (min, max);
  h_xzjet_10z20->GetYaxis ()->SetRangeUser (min, max);
  h_xzjet_20z40->GetYaxis ()->SetRangeUser (min, max);
  h_xzjet_40z80->GetYaxis ()->SetRangeUser (min, max);

  h_xzjet_5z10->SetLineWidth (2);
  h_xzjet_10z20->SetLineWidth (2);
  h_xzjet_20z40->SetLineWidth (2);
  h_xzjet_40z80->SetLineWidth (2);

  h_xzjet_5z10->SetLineColor (colors[0]);
  h_xzjet_5z10->SetMarkerColor (colors[0]);
  h_xzjet_10z20->SetLineColor (colors[1]);
  h_xzjet_10z20->SetMarkerColor (colors[1]);
  h_xzjet_20z40->SetLineColor (colors[2]);
  h_xzjet_20z40->SetMarkerColor (colors[2]);
  h_xzjet_40z80->SetLineColor (colors[3]);
  h_xzjet_40z80->SetMarkerColor (colors[3]);

  h_xzjet_5z10->GetXaxis ()->SetTitle ("#it{x}_{Z}^{ jet}");
  h_xzjet_5z10->GetYaxis ()->SetTitle ("1/N_{Z} dN/dx");
  h_xzjet_10z20->GetXaxis ()->SetTitle ("#it{x}_{Z}^{ jet}");
  h_xzjet_10z20->GetYaxis ()->SetTitle ("1/N_{Z} dN/dx");
  h_xzjet_20z40->GetXaxis ()->SetTitle ("#it{x}_{Z}^{ jet}");
  h_xzjet_20z40->GetYaxis ()->SetTitle ("1/N_{Z} dN/dx");
  h_xzjet_40z80->GetXaxis ()->SetTitle ("#it{x}_{Z}^{ jet}");
  h_xzjet_40z80->GetYaxis ()->SetTitle ("1/N_{Z} dN/dx");

  h_xzjet_5z10->Draw ("hist");
  h_xzjet_10z20->Draw ("same hist");
  h_xzjet_20z40->Draw ("same hist");
  h_xzjet_40z80->Draw ("same hist");

  myText (0.66, 0.89, colors[0], "5 < #it{p}_{T}^{ Z} < 10 GeV", 0.05); 
  myText (0.66, 0.82, colors[1], "10 < #it{p}_{T}^{ Z} < 20 GeV", 0.05); 
  myText (0.66, 0.75, colors[2], "20 < #it{p}_{T}^{ Z} < 40 GeV", 0.05); 
  myText (0.66, 0.68, colors[3], "40 < #it{p}_{T}^{ Z} < 80 GeV", 0.05); 


  //outFile->Close();
  //if (outFile) delete outFile;
  
  return 0;
}

void analyze_z_LHC () {
  main ();
} 
