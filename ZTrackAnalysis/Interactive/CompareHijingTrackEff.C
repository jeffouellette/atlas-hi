#ifndef __CompareHijingTrackEff_C__
#define __CompareHijingTrackEff_C__

#include "Params.h"

#include <Utilities.h>
#include <ArrayTemplates.h>

void CompareHijingTrackEff () {

  const float multBins[6] = {-0.5, 499.5, 699.5, 1499.5, 1799.5, 3399.5};
  const int _numMultBins = (sizeof (multBins) / sizeof (multBins[0])) - 1;
  const int numMultBins = _numMultBins - 2;

  double newPtTrkBins[6] = {1, 1.5, 2, 6, 10, 60};
  int newNPtTrkBins = (sizeof (newPtTrkBins) / sizeof (newPtTrkBins[0])) - 1;

  const TString dir = "/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/TrackingEfficiencies/Nominal/";

  TFile* pythiaFile = new TFile (dir + "trackingEfficienciesMult_18.root", "read");
  //TFile* pythiaFile = new TFile (dir + "PbPb_Hijing_15.root", "read");

  TH1D** h_pythia_trk_effs = Get1DArray <TH1D*> (_numMultBins);
  for (int iMult = 0; iMult < _numMultBins; iMult++) {
    TH1D* num;
    TH1D* den;
    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {

      if (iEta == 0) {
        num = (TH1D*) pythiaFile->Get (Form ("h_trk_eff_num_iMult%i_iEta%i", iMult, iEta));
        den = (TH1D*) pythiaFile->Get (Form ("h_trk_eff_den_iMult%i_iEta%i", iMult, iEta));
      }
      else {
        num->Add ((TH1D*) pythiaFile->Get (Form ("h_trk_eff_num_iMult%i_iEta%i", iMult, iEta)));
        den->Add ((TH1D*) pythiaFile->Get (Form ("h_trk_eff_den_iMult%i_iEta%i", iMult, iEta)));
      }

      //if (iMult > 0) {
      //  num->Rebin (2);
      //  den->Rebin (2);
      //}
    }

    RebinSomeBins (num, newNPtTrkBins, newPtTrkBins);
    RebinSomeBins (den, newNPtTrkBins, newPtTrkBins);

    h_pythia_trk_effs[iMult] = (TH1D*) num->Clone (Form ("h_pythia_trk_eff_iMult%i", iMult));
    h_pythia_trk_effs[iMult]->Divide (den);
  }


  TFile* hijingFile = new TFile (dir + "PbPb18_Zee_Hijing_18.root", "read");
  TH1D** h_hijing_trk_effs = Get1DArray <TH1D*> (_numMultBins);
  for (int iMult = 0; iMult < _numMultBins; iMult++) {
    TH1D* num;
    TH1D* den;
    for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {

      if (iEta == 0) {
        num = (TH1D*) hijingFile->Get (Form ("h_trk_eff_num_iMult%i_iEta%i", iMult, iEta));
        den = (TH1D*) hijingFile->Get (Form ("h_trk_eff_den_iMult%i_iEta%i", iMult, iEta));
      }
      else {
        num->Add ((TH1D*) hijingFile->Get (Form ("h_trk_eff_num_iMult%i_iEta%i", iMult, iEta)));
        den->Add ((TH1D*) hijingFile->Get (Form ("h_trk_eff_den_iMult%i_iEta%i", iMult, iEta)));
      }

      //if (iMult > 0) {
      //  num->Rebin (2);
      //  den->Rebin (2);
      //}
    }

    RebinSomeBins (num, newNPtTrkBins, newPtTrkBins);
    RebinSomeBins (den, newNPtTrkBins, newPtTrkBins);

    h_hijing_trk_effs[iMult] = (TH1D*) num->Clone (Form ("h_hijing_trk_eff_iMult%i", iMult));
    h_hijing_trk_effs[iMult]->Divide (den);
  }

  TCanvas* c = new TCanvas ("c", "", 1400, 600);
  c->Draw ();
  const double dPadY = 0.4;
  const int axisTextSize = 23;

  for (int _iMult = 0; _iMult < _numMultBins; _iMult++) {
    int iMult = _iMult;
    if (_iMult == 1)
      continue;
    else if (_iMult > 1)
      iMult--;
    if (_iMult == 3)
      continue;
    else if (_iMult > 3)
      iMult--;

    c->cd ();

    const char* uPadName = Form ("uPad_%i", iMult);
    const char* dPadName = Form ("dPad_%i", iMult);

    TPad* uPad = new TPad (uPadName, "", (1./(numMultBins))*(iMult), dPadY, (1./(numMultBins))*(iMult+1), 1);
    TPad* dPad = new TPad (dPadName, "", (1./(numMultBins))*(iMult), 0, (1./(numMultBins))*(iMult+1), dPadY);

    uPad->SetTopMargin (0.04);
    uPad->SetBottomMargin (0);
    uPad->SetLeftMargin (0.17);
    uPad->SetRightMargin (0.06);
    dPad->SetTopMargin (0);
    dPad->SetBottomMargin (0.30);
    dPad->SetLeftMargin (0.17);
    dPad->SetRightMargin (0.06);
    uPad->Draw ();
    dPad->Draw ();

    TH1D* h = nullptr;


    uPad->cd ();
    uPad->SetLogx ();
    //for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {
      h = h_pythia_trk_effs[_iMult];

      h->SetLineColor (kAzure+2);
      h->SetMarkerColor (kAzure+2);
      h->SetMarkerStyle (kOpenSquare);

      h->GetYaxis ()->SetRangeUser (0.3, 1.08);

      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      h->GetYaxis ()->SetTitle ("Reconstruction Eff.");

      h->GetXaxis ()->SetTitleFont (43);
      h->GetXaxis ()->SetTitleSize (22);
      h->GetXaxis ()->SetLabelFont (43);
      h->GetXaxis ()->SetLabelSize (18);

      h->GetYaxis ()->SetTitleFont (43);
      h->GetYaxis ()->SetTitleSize (22);
      h->GetYaxis ()->SetLabelFont (43);
      h->GetYaxis ()->SetLabelSize (18);

      h->GetXaxis ()->SetTitleOffset (1.6 * h->GetXaxis ()->GetTitleOffset ());
      h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

      h->GetXaxis ()->SetMoreLogLabels ();
      

      h = h_hijing_trk_effs[_iMult];

      h->SetLineColor (kRed+1);
      h->SetMarkerColor (kRed+1);
      h->SetMarkerStyle (kOpenCircle);

      h->GetYaxis ()->SetRangeUser (0.3, 1.08);

      h->GetXaxis ()->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      h->GetYaxis ()->SetTitle ("Reco. Eff.");

      h->GetXaxis ()->SetTitleFont (43);
      h->GetXaxis ()->SetTitleSize (22);
      h->GetXaxis ()->SetLabelFont (43);
      h->GetXaxis ()->SetLabelSize (18);

      h->GetYaxis ()->SetTitleFont (43);
      h->GetYaxis ()->SetTitleSize (22);
      h->GetYaxis ()->SetLabelFont (43);
      h->GetYaxis ()->SetLabelSize (18);

      h->GetXaxis ()->SetTitleOffset (1.6 * h->GetXaxis ()->GetTitleOffset ());
      h->GetYaxis ()->SetTitleOffset (1.2 * h->GetYaxis ()->GetTitleOffset ());

      h->GetXaxis ()->SetMoreLogLabels ();
     
 
      h_pythia_trk_effs[_iMult]->Draw ("e1");
      h_hijing_trk_effs[_iMult]->Draw ("e1 same");


      myText (0.22, 0.88, kBlack, Form ("%i #leq N_{ch}^{>0.5GeV} #leq %i", (int)(multBins[_iMult]+0.5), (int)(multBins[_iMult+1]+0.5)), 0.08);

      if (iMult == 0)
        myText (0.25, 0.10, kBlack, "Only #pi^{#pm}'s", 0.08);
 
      //if (iMult == 1)
      //  myMarkerTextNoLine (0.65, 0.28-0.05*iEta, colors[iEta], kOpenCircle, Form ("%g < |#eta| < %g", etaTrkBins[iEta], etaTrkBins[iEta+1]), 1.2, 0.05);

      if (iMult == 1) {
        //myMarkerTextNoLine (0.25, 0.13, kBlack, kOpenCircle, "2018 cond.", 1.2, 0.05);
        //myMarkerTextNoLine (0.25, 0.08, kBlack, kOpenSquare, "2015 cond.", 1.2, 0.05);
        myMarkerTextNoLine (0.25, 0.13, kRed+1, kOpenCircle, "Pythia Z#rightarrow ee with HIJING Overlay", 1.2, 0.05);
        myMarkerTextNoLine (0.25, 0.08, kAzure+2, kOpenSquare, "Pythia Z#rightarrow ll with Data Overlay", 1.2, 0.05);
      }
    //}


    dPad->cd ();
    dPad->SetLogx ();
    TH1D* rat;
    //for (int iEta = 0; iEta < numEtaTrkBins; iEta++) {

      rat = (TH1D*) h_pythia_trk_effs[_iMult]->Clone (Form ("h_ratio_trk_effs_iMult%i", _iMult));
      rat->Divide (h_hijing_trk_effs[_iMult]);
      
      rat->SetLineColor (kAzure+2);
      rat->SetMarkerColor (kAzure+2);

      rat->GetYaxis ()->SetRangeUser (0.9, 1.1);
      rat->GetXaxis ()->SetTitle ("#it{p}_{T}^{truth} [GeV]");
      //rat->GetYaxis ()->SetTitle ("2015 / 2018");
      rat->GetYaxis ()->SetTitle ("PYTHIA / HIJING");
      rat->GetYaxis ()->CenterTitle ();
      rat->GetXaxis ()->SetTitleOffset (1.4*rat->GetXaxis ()->GetTitleOffset ());

      rat->Draw ("e1");

    //}

    TLine* l = new TLine (rat->GetBinLowEdge (1), 1, rat->GetBinLowEdge (rat->GetNbinsX ()) + rat->GetBinWidth (rat->GetNbinsX ()), 1);
    l->SetLineColor (kPink-8);
    l->SetLineStyle (2);
    l->SetLineWidth (2);
    l->Draw ("same");
  }

  c->SaveAs ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/Plots/ZTrackAnalysis/TrackingEfficienciesComp.pdf");




  TFile* f_mcEW = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/MCAnalysis/eventWeightsFile.root", "read");
  TH1D* h_PbPbNchDist_mc = (TH1D*) f_mcEW->Get ("h_PbPbNchDist_mc");

  TFile* f_hijingEW = new TFile ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/rootFiles/MinbiasAnalysis/PbPb_Hijing_18_eventWeights.root", "read");
  TH1D* h_PbPbNchDist_hijing = (TH1D*) f_hijingEW->Get ("h_PbPbNchDist_hijing");

  TCanvas* c2 = new TCanvas ("c2", "", 1200, 800);
  c2->SetLogy ();
  h_PbPbNchDist_mc->SetLineColor (kAzure+2);
  h_PbPbNchDist_mc->SetMarkerColor (kAzure+2);
  h_PbPbNchDist_mc->SetMarkerStyle (kOpenCircle);

  h_PbPbNchDist_mc->GetXaxis ()->SetTitle ("N_{ch}");
  h_PbPbNchDist_mc->GetYaxis ()->SetTitle ("Normalized entries");
  h_PbPbNchDist_mc->GetYaxis ()->SetRangeUser (1e-4, 0.2);
  
  h_PbPbNchDist_hijing->SetLineColor (kRed+1);
  h_PbPbNchDist_hijing->SetMarkerColor (kRed+1);
  h_PbPbNchDist_hijing->SetMarkerStyle (kOpenSquare);

  h_PbPbNchDist_hijing->GetXaxis ()->SetTitle ("N_{ch}");
  h_PbPbNchDist_hijing->GetYaxis ()->SetTitle ("Normalized entries");
  h_PbPbNchDist_hijing->GetYaxis ()->SetRangeUser (1e-4, 0.2);

  h_PbPbNchDist_mc->Draw ("e1");
  h_PbPbNchDist_hijing->Draw ("e1 same");

  TLine* l1 = new TLine (500, 1e-4, 500, 0.2);
  l1->SetLineStyle (2);
  l1->SetLineWidth (2);
  l1->Draw ("same");
  TLine* l2 = new TLine (700, 1e-4, 700, 0.2);
  l2->SetLineStyle (2);
  l2->SetLineWidth (2);
  l2->Draw ("same");
  TLine* l3 = new TLine (1500, 1e-4, 1500, 0.2);
  l3->SetLineStyle (2);
  l3->SetLineWidth (2);
  l3->Draw ("same");
  TLine* l4 = new TLine (1800, 1e-4, 1800, 0.2);
  l4->SetLineStyle (2);
  l4->SetLineWidth (2);
  l4->Draw ("same");
  TLine* l5 = new TLine (3400, 1e-4, 3400, 0.2);
  l5->SetLineStyle (2);
  l5->SetLineWidth (2);
  l5->Draw ("same");

  myMarkerText (0.55, 0.88, kAzure+2, kOpenCircle, "Pythia Z#rightarrow ll with Data Overlay", 1.25, 0.04);
  myMarkerText (0.55, 0.81, kRed+1, kOpenSquare, "Pythia Z#rightarrow ee with HIJING Overlay", 1.25, 0.04);

  c2->SaveAs ("/Users/jeffouellette/Research/atlas-hi/ZTrackAnalysis/Plots/ZTrackAnalysis/NchComp.pdf");
}

#endif
