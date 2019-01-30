#include "Params.h"

#include <GlobalParams.h>
#include <Utilities.h>
#include <ArrayTemplates.h>
#include <Jet.h>

#include <TLine.h>
#include <TH2D.h>
#include <TFile.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

#include <iostream>

using namespace atlashi;
using namespace std;

using namespace JetTrackAnalysis;

const int numRtrkBins = 50;
const double* rtrkBins = linspace (0, 2.0, numRtrkBins);
//const double* rtrkBins = linspace (0, 250, numRtrkBins);

const double pbins[17] = {20, 30, 40, 50, 60, 80, 100, 120, 140, 160, 180, 200, 240, 280, 320, 360, 400};
const int numpbins = sizeof (pbins) / sizeof (pbins[0]) - 1;

void RtrkHistMaker () {

  SetAtlasStyle ();

  SetupDirectories ("RtrkStudy/", "JetTrackAnalysis/");

  TH3D** rtrkHists = Get1DArray <TH3D*> (2); // iMC (0=signal, 1=overlay)
  for (int iMC = 0; iMC < 2; iMC++) {
    rtrkHists[iMC] = new TH3D (Form ("rtrkHist_%s", iMC == 0 ? "signal":"overlay"), "", numpbins, pbins, numetabins, etabins, numRtrkBins, rtrkBins);
    rtrkHists[iMC]->Sumw2 ();
  }

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");
  TTree* signalJets = NULL, *overlayJets = NULL;
  if (inFile) {
    signalJets = (TTree*)inFile->Get ("signalJets");
    overlayJets = (TTree*)inFile->Get ("overlayJets");
  }
  if (inFile == NULL || signalJets == NULL || overlayJets == NULL) {
    cout << "Error: In RtrkHistMaker.C: TTrees not obtained. Quitting." << endl;
    return;
  }

  Jet jet;
  float sum_trk_pt = 0, evtWeight = 0;
  signalJets->SetBranchAddress ("jet_pt",     &jet.Pt);
  signalJets->SetBranchAddress ("jet_eta",    &jet.Eta);
  signalJets->SetBranchAddress ("jet_phi",    &jet.Phi);
  signalJets->SetBranchAddress ("jet_e",      &jet.E);
  signalJets->SetBranchAddress ("sum_trk_pt", &sum_trk_pt);
  signalJets->SetBranchAddress ("evtWeight",  &evtWeight);

  overlayJets->SetBranchAddress ("jet_pt",     &jet.Pt);
  overlayJets->SetBranchAddress ("jet_eta",    &jet.Eta);
  overlayJets->SetBranchAddress ("jet_phi",    &jet.Phi);
  overlayJets->SetBranchAddress ("jet_e",      &jet.E);
  overlayJets->SetBranchAddress ("sum_trk_pt", &sum_trk_pt);
  overlayJets->SetBranchAddress ("evtWeight",  &evtWeight);

  const int nSignalJets = signalJets->GetEntries ();
  const int nOverlayJets = overlayJets->GetEntries ();

  for (int iJet = 0; iJet < nSignalJets; iJet++) {
    signalJets->GetEntry (iJet);
    rtrkHists[0]->Fill (jet.Pt, jet.Eta, sum_trk_pt / jet.Pt, evtWeight);
    //rtrkHists[0]->Fill (jet.Pt, jet.Eta, sum_trk_pt, evtWeight);
  }
  for (int iJet = 0; iJet < nOverlayJets; iJet++) {
    overlayJets->GetEntry (iJet);
    rtrkHists[1]->Fill (jet.Pt, jet.Eta, sum_trk_pt / jet.Pt, evtWeight);
    //rtrkHists[1]->Fill (jet.Pt, jet.Eta, sum_trk_pt, evtWeight);
  }

  if (signalJets) { delete signalJets; signalJets = NULL; }
  if (overlayJets) { delete overlayJets; overlayJets = NULL; }

  inFile->Close ();

  TFile* outFile = new TFile (Form ("%s/rtrkHists.root", rootPath.Data ()), "recreate");
  rtrkHists[0]->Write ();
  rtrkHists[1]->Write ();

  const double* rtrk_los = linspace (0.1, 0.1, numpbins);
  const double* rtrk_his = linspace (1.50, 1.50, numpbins);
  //const double* rtrk_los = linspace (1, 1, numpbins);
  //const double* rtrk_his = linspace (250, 250, numpbins);

  TH2D* signalProj2d = Project2D ("signalProj2d", rtrkHists[0], "x", "z", 1, numetabins);
  //signalProj2d->RebinY (20);
  TH1D* signalProf = GetProfileX ("signalProf", signalProj2d, numpbins, pbins, true, rtrk_los, rtrk_his);

  TH2D* overlayProj2d = Project2D ("overlayProj2d", rtrkHists[1], "x", "z", 1, numetabins);
  //overlayProj2d->RebinY (20);
  TH1D* overlayProf = GetProfileX ("overlayProf", overlayProj2d, numpbins, pbins, true, rtrk_los, rtrk_his);

  signalProf->GetYaxis ()->SetRangeUser (0, 1);
  overlayProf->GetYaxis ()->SetRangeUser (0, 1);
  //signalProf->GetYaxis ()->SetRangeUser (0, 150);
  //overlayProf->GetYaxis ()->SetRangeUser (0, 150);

  signalProf->GetXaxis ()->SetTitle ("Jet #it{p}_{T} #left[GeV#right]");
  signalProf->GetYaxis ()->SetTitle ("<R_{trk}> = #Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{Calo}");
  //signalProf->GetYaxis ()->SetTitle ("<#Sigma#it{p}_{T}^{trk}> #left[GeV#right]");
  overlayProf->GetXaxis ()->SetTitle ("Jet #it{p}_{T} #left[GeV#right]");
  overlayProf->GetYaxis ()->SetTitle ("<R_{trk}> = #Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{Calo}");
  //overlayProf->GetYaxis ()->SetTitle ("<#Sigma#it{p}_{T}^{trk}> #left[GeV#right]");

  TCanvas* c1 = new TCanvas ("c1", "c1", 800, 600);
  c1->cd ();

  signalProf->SetLineColor (kBlue);
  signalProf->SetMarkerColor (kBlue);

  overlayProf->SetLineColor (kRed);
  overlayProf->SetMarkerColor (kRed);

  signalProf->Draw ("e1 x0");
  overlayProf->Draw ("e1 x0 same");

  myText (0.6, 0.88, kBlue, "Pythia8", 0.04);
  myText (0.6, 0.81, kRed, "Pythia8 with Data Overlay", 0.04);

  c1->SaveAs (Form ("%s/RtrkComp.pdf", plotPath.Data ()));
  //c1->SaveAs (Form ("%s/SumPtTrkComp.pdf", plotPath.Data ()));

  TCanvas* c2 = new TCanvas ("c2", "c2", 800, 600);
  c2->cd ();

  signalProj2d->GetXaxis ()->SetTitle ("#it{p}_{T}^{calo} #left[GeV#right]");
  signalProj2d->GetYaxis ()->SetTitle ("R_{trk} = <#Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{calo}>");
  //signalProj2d->GetYaxis ()->SetTitle ("#Sigma#it{p}_{T}^{trk} #left[GeV#right]");

  signalProj2d->GetYaxis ()->SetRangeUser (0.1, 2);
  //signalProj2d->GetYaxis ()->SetRangeUser (5, 200);
  signalProj2d->Draw ("colz");
  myText (0.6, 0.88, kBlack, "Pythia8", 0.04);

  c2->SaveAs (Form ("%s/SignalRtrkDist.pdf", plotPath.Data ()));
  //c2->SaveAs (Form ("%s/SignalSumPtTrkDist.pdf", plotPath.Data ()));

  TCanvas* c3 = new TCanvas ("c3", "c3", 800, 600);
  c3->cd ();

  overlayProj2d->GetXaxis ()->SetTitle ("#it{p}_{T}^{calo} #left[GeV#right]");
  overlayProj2d->GetYaxis ()->SetTitle ("R_{trk} = <#Sigma#it{p}_{T}^{trk} / #it{p}_{T}^{calo}>");
  //overlayProj2d->GetYaxis ()->SetTitle ("#Sigma#it{p}_{T}^{trk} #left[GeV#right]");

  overlayProj2d->GetYaxis ()->SetRangeUser (0.1, 2);
  //overlayProj2d->GetYaxis ()->SetRangeUser (5, 200);
  overlayProj2d->Draw ("colz");
  myText (0.6, 0.88, kBlack, "Pythia8 with Data Overlay", 0.04);

  c3->SaveAs (Form ("%s/OverlayRtrkDist.pdf", plotPath.Data ()));
  //c3->SaveAs (Form ("%s/OverlaySumPtTrkDist.pdf", plotPath.Data ()));

}
