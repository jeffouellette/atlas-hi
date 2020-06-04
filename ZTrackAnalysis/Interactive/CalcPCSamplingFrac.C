#ifndef __CalcPCSamplingFrac_C__
#define __CalcPCSamplingFrac_C__

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include <AtlasUtils.h>

#include "Params.h"

using namespace std;

void CalcPCSamplingFrac () {

  TFile* infile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/MixingAnalysis/Nominal/fcals.root", "read");
  TTree* intree = (TTree*) infile->Get ("fcals");

  float fcal_et = 0;
  bool HLT_mb_sptrk_L1ZDC_A_C_VTE50 = false;
  bool HLT_noalg_pc_L1TE50_VTE600_0ETA49 = false;
  bool HLT_noalg_cc_L1TE600_0ETA49 = false;

  intree->SetBranchAddress ("fcal_et", &fcal_et);
  intree->SetBranchAddress ("HLT_mb_sptrk_L1ZDC_A_C_VTE50",      &HLT_mb_sptrk_L1ZDC_A_C_VTE50);
  intree->SetBranchAddress ("HLT_noalg_pc_L1TE50_VTE600.0ETA49", &HLT_noalg_pc_L1TE50_VTE600_0ETA49);
  intree->SetBranchAddress ("HLT_noalg_cc_L1TE600_0ETA49",       &HLT_noalg_cc_L1TE600_0ETA49);

  TFile* outfile = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/mbSamplingFracs.root", "recreate");
  TH1D* h_centrality_sampled_pc = new TH1D ("h_centrality_sampled_pc", ";Centrality [%];Events", 80, 0, 80);
  TH1D* h_centrality_sampled_cc = new TH1D ("h_centrality_sampled_cc", ";Centrality [%];Events", 80, 0, 80);
  TH1D* h_centrality_all_pc = new TH1D ("h_centrality_all_pc", ";Centrality [%];Events", 80, 0, 80);
  TH1D* h_centrality_all_cc = new TH1D ("h_centrality_all_cc", ";Centrality [%];Events", 80, 0, 80);

  long nEvt = intree->GetEntries ();
  for (long iEvt = 0; iEvt < nEvt; iEvt++) {
    if (nEvt > 100 && iEvt % (nEvt / 100) == 0)
      cout << iEvt / (nEvt / 100) << "\% done...\r" << flush;

    intree->GetEntry (iEvt);

    const int centrality = GetPercentileCentrality (fcal_et); 
    if (HLT_mb_sptrk_L1ZDC_A_C_VTE50 || HLT_noalg_pc_L1TE50_VTE600_0ETA49)
      h_centrality_sampled_pc->Fill (centrality);
    if (HLT_noalg_cc_L1TE600_0ETA49)
      h_centrality_sampled_cc->Fill (centrality);
  }

  infile->Close ();
  if (infile) delete infile;
  infile = nullptr;
  intree = nullptr;

  TFile* infile2 = new TFile ("/atlasgpfs01/usatlas/data/jeff/ZTrackAnalysis/rootFiles/EventSkims/Nominal/fcals.root", "read");
  TTree* intree2 = (TTree*) infile2->Get ("fcals");

  intree2->SetBranchAddress ("fcal_et", &fcal_et);
  intree2->SetBranchAddress ("HLT_mb_sptrk_L1ZDC_A_C_VTE50",      &HLT_mb_sptrk_L1ZDC_A_C_VTE50);
  intree2->SetBranchAddress ("HLT_noalg_pc_L1TE50_VTE600.0ETA49", &HLT_noalg_pc_L1TE50_VTE600_0ETA49);
  intree2->SetBranchAddress ("HLT_noalg_cc_L1TE600_0ETA49",       &HLT_noalg_cc_L1TE600_0ETA49);

  nEvt = intree2->GetEntries ();
  for (long iEvt = 0; iEvt < nEvt; iEvt++) {
    if (nEvt > 100 && iEvt % (nEvt / 100) == 0)
      cout << iEvt / (nEvt / 100) << "\% done...\r" << flush;

    intree2->GetEntry (iEvt);

    const int centrality = GetPercentileCentrality (fcal_et); 
    if (HLT_mb_sptrk_L1ZDC_A_C_VTE50 || HLT_noalg_pc_L1TE50_VTE600_0ETA49)
      h_centrality_all_pc->Fill (centrality);
    if (HLT_noalg_cc_L1TE600_0ETA49)
      h_centrality_all_cc->Fill (centrality);
  }

  infile2->Close ();
  if (infile2) delete infile2;
  infile2 = nullptr;
  intree2 = nullptr;

  outfile->cd ();


  TCanvas* c = new TCanvas ("c", "", 800, 800);
  gPad->SetLogy ();

  h_centrality_sampled_pc->SetLineColor (kAzure-1);
  h_centrality_sampled_cc->SetLineColor (kRed+1);
  h_centrality_all_pc->SetLineColor (kAzure-1);
  h_centrality_all_cc->SetLineColor (kRed+1);

  h_centrality_sampled_pc->SetLineStyle (1);
  h_centrality_sampled_cc->SetLineStyle (1);
  h_centrality_all_pc->SetLineStyle (2);
  h_centrality_all_cc->SetLineStyle (2);

  //const float max = 7.0 * fmax (h_centrality_sampled_pc->GetMaximum (), h_centrality_sampled_cc->GetMaximum ());
  //const float min = 10. * fmin (h_centrality_sampled_pc->GetMinimum (0), h_centrality_sampled_cc->GetMinimum (0));
  const float max = 5e8;
  const float min = 2e3;

  h_centrality_sampled_pc->GetYaxis ()->SetRangeUser (min, max);
  h_centrality_sampled_cc->GetYaxis ()->SetRangeUser (min, max);
  h_centrality_all_pc->GetYaxis ()->SetRangeUser (min, max);
  h_centrality_all_cc->GetYaxis ()->SetRangeUser (min, max);

  h_centrality_sampled_pc->Draw ("hist");
  h_centrality_sampled_cc->Draw ("hist same");
  h_centrality_all_pc->Draw ("hist same");
  h_centrality_all_cc->Draw ("hist same");

  double ccAvgErr, pcAvgErr;

  double ccAvg = (h_centrality_sampled_cc->IntegralAndError (h_centrality_sampled_cc->FindFixBin (30), h_centrality_sampled_cc->FindFixBin (51), ccAvgErr)/ 22.);
  ccAvgErr = ccAvgErr / 22.;
  double pcAvg = (h_centrality_sampled_pc->IntegralAndError (h_centrality_sampled_pc->FindFixBin (58), h_centrality_sampled_pc->FindFixBin (77), pcAvgErr)/ 20.);
  pcAvgErr = pcAvgErr / 20.;

  myText (0.22, 0.86, kBlack, "#bf{#it{ATLAS}} Internal", 0.055);
  myText (0.22, 0.80, kBlack, "Pb+Pb, #sqrt{s_{NN}} = 5.02 TeV", 0.04);
  myText (0.22, 0.74, kBlack, Form ("Eff. PC prescale: %.4f #pm %.4f", pcAvg / ccAvg, fabs (pcAvg/ccAvg) * sqrt (pow (ccAvgErr/ccAvg, 2) + pow (pcAvgErr/pcAvg, 2))), 0.04);

  myText (0.20, 0.36, kBlack, "All evts", 0.035);
  myText (0.32, 0.36, kBlack, "Sample", 0.035);
  myLineText (0.290, 0.315, kAzure-1, 2, "", 2.0, 0.054);
  myLineText (0.415, 0.315, kAzure-1, 1, "", 2.0, 0.054);
  myLineText (0.290, 0.255, kRed+1, 2, "", 2.0, 0.054);
  myLineText (0.415, 0.255, kRed+1, 1, "", 2.0, 0.054);
  myText (0.43, 0.30, kAzure-1, "PC stream", 0.04);
  myText (0.43, 0.24, kRed+1, "CC stream", 0.04);


  c->SaveAs ("/atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Interactive/Plots/MBCentralitySampling.pdf");
 

  h_centrality_sampled_pc->Write ();
  h_centrality_sampled_cc->Write ();
  h_centrality_all_pc->Write ();
  h_centrality_all_cc->Write ();
  outfile->Close ();
}

#endif
