#ifndef __PlotEventWeights_C__
#define __PlotEventWeights_C__

#include "Params.h"

using namespace atlashi;

void PlotEventWeights () {

  TFile* inFile = new TFile ("../rootFiles/DataAnalysis/Nominal/data18hi.root", "read");

  bool isEE = false;
  float event_weight = 0;
  float z_pt = 0;
  float fcal_et = 0;

  TH1D**** h_event_weight = Get3DArray <TH1D*> (2, numCentBins, nPtZBins);
  for (short iSpc : {0, 1}) {
    for (short iCent = 0; iCent < numCentBins; iCent++) {
      for (short iPtZ = 0; iPtZ < nPtZBins; iPtZ++) {
        h_event_weight[iSpc][iCent][iPtZ] = new TH1D (Form ("h_event_weight_iSpc%i_iCent%i_iPtZ%i", iSpc, iCent, iPtZ), "1 / #varepsilon_{trig} #times #varepsilon_{ID}", 80, 0.9, 2.5);
        h_event_weight[iSpc][iCent][iPtZ]->Sumw2 ();
      }
    }
  }

  TTree* inTree = (TTree*) inFile->Get ("PbPbZTrackTree");
  inTree->SetBranchAddress ("isEE",         &isEE);
  inTree->SetBranchAddress ("event_weight", &event_weight);
  inTree->SetBranchAddress ("z_pt",         &z_pt);
  inTree->SetBranchAddress ("fcal_et",      &fcal_et);

  for (int iEvt = 0; iEvt < inTree->GetEntries (); iEvt++) {
    inTree->GetEntry (iEvt);

    const short iCent = GetCentBin (fcal_et);
    if (iCent < 1 || iCent > numCentBins-1)
      continue;
    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;
    const short iSpc = (isEE ? 0 : 1);

    h_event_weight[iSpc][iCent][iPtZ]->Fill (event_weight);
  }

  inTree = (TTree*) inFile->Get ("ppZTrackTree");
  inTree->SetBranchAddress ("isEE",         &isEE);
  inTree->SetBranchAddress ("event_weight", &event_weight);
  inTree->SetBranchAddress ("z_pt",         &z_pt);

  for (int iEvt = 0; iEvt < inTree->GetEntries (); iEvt++) {
    inTree->GetEntry (iEvt);

    const short iCent = 0;
    const short iPtZ = GetPtZBin (z_pt);
    if (iPtZ < 0 || iPtZ > nPtZBins-1)
      continue;
    const short iSpc = (isEE ? 0 : 1);

    h_event_weight[iSpc][iCent][iPtZ]->Fill (event_weight);
  }



  TCanvas* c1 = new TCanvas ("c1", "", 800, 800);
  c1->Divide (2, 2);
  for (short iCent = 0; iCent < numCentBins; iCent++) {
    c1->cd (iCent+1);
    for (short iPtZ = 2; iPtZ < nPtZBins; iPtZ++) {
      for (short iSpc : {0, 1}) {
        TH1D* h = h_event_weight[iSpc][iCent][iPtZ];
        h->GetXaxis ()->SetTitle ("1 / #varepsilon_{trig} #times #varepsilon_{ID}");
        h->GetYaxis ()->SetTitle ("Counts");
        h->SetLineColor (colors[iPtZ-1]);
        h->SetMarkerColor (colors[iPtZ-1]);
        h->SetLineStyle (iSpc+1);
        h->SetLineWidth (1);

        h->Draw (iPtZ == 2 && iSpc == 0 ? "hist" : "hist same");
      }
    }
  }

  c1->cd (1);
  myText (0.23, 0.87, kBlack, "#bf{#it{ATLAS}} Internal", 0.075);
  myText (0.23, 0.80, kBlack, "#it{pp}, #sqrt{s_{NN}} = 5.02 TeV", 0.06); 
  myText (0.53, 0.65, kBlack, "15 < #it{p}_{T}^{Z} < 30 GeV", 0.05);
  myText (0.53, 0.58, kBlack, "30 < #it{p}_{T}^{Z} < 60 GeV", 0.05);
  myText (0.53, 0.51, kBlack, "#it{p}_{T}^{Z} > 60 GeV", 0.05);
  myText (0.35, 0.71, kBlack, "ee", 0.05);
  myText (0.43, 0.71, kBlack, "#mu#mu", 0.05);
  myLineText (0.40, 0.665, colors[1], 1, "", 2.0, 0.036);
  myLineText (0.49, 0.665, colors[1], 2, "", 2.0, 0.036);
  myLineText (0.40, 0.595, colors[2], 1, "", 2.0, 0.036);
  myLineText (0.49, 0.595, colors[2], 2, "", 2.0, 0.036);
  myLineText (0.40, 0.525, colors[3], 1, "", 2.0, 0.036);
  myLineText (0.49, 0.525, colors[3], 2, "", 2.0, 0.036);
  c1->cd (2); myText (0.23, 0.87, kBlack, "Pb+Pb, 30-80%", 0.06);
  c1->cd (3); myText (0.23, 0.87, kBlack, "Pb+Pb, 10-30%", 0.06);
  c1->cd (4); myText (0.23, 0.87, kBlack, "Pb+Pb, 0-10%", 0.06);
}

#endif
