#include "JetTrackEtaPhiHist.h"
#include "Params.h"

#include <GlobalParams.h>
#include <Utils.h>

#include <TLine.h>
#include <TH2D.h>
#include <TFile.h>

#include <AtlasStyle.h>
#include <AtlasUtils.h>

using namespace atlashi;

namespace JetTrackAnalysis {

void JetTrackEtaPhiHist () {

  SetAtlasStyle ();

  SetupDirectories ("JetTrackEtaPhi/", "JetTrackAnalysis/");

  TFile* inFile = new TFile (Form ("%s/outFile.root", rootPath.Data ()), "read");
  //TFile* inFile = new TFile (Form ("%s/outFile_mc.root", rootPath.Data ()), "read");

  //vector<int> runNumbers (0);
  //for (short i = 0; i < sizeof (full_run_list)/sizeof (full_run_list[0]); i++) runNumbers.push_back (full_run_list[i]);

  TH2D* JetTrackDeltaEtaDeltaPhi[numCentBins][2] = {};
  TH2D* JetPtSpectrum[numCentBins] = {};
  TH2D* JetCounts[numCentBins] = {};
  TH2D* TrackPtSpectrum[numCentBins] = {};
  TH2D* TrackCounts[numCentBins] = {};

  TH1D* TrackSelectionHist[numCentBins] = {};
  TH1D* JetSelectionHist[numCentBins] = {};

  TH2D* DijetDeltaEtaDeltaPhi[numCentBins] = {};
  TH2D* DijetEta1Eta2[numCentBins] = {};
  TH1D* DijetAj[numCentBins] = {};

  double ajmax = 0;
  double ajmin = 1e100;
  for (short iCent = 0; iCent < numCentBins; iCent++) {
   for (short iCut = 0; iCut < 2; iCut++) {
    JetTrackDeltaEtaDeltaPhi[iCent][iCut] = (TH2D*)inFile->Get (Form ("JetTrackDeltaEtaDeltaPhi_cent%i_cut%i", iCent, iCut));
    //JetTrackDeltaEtaDeltaPhi[iCent][iCut]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));
   }
   JetPtSpectrum[iCent] = (TH2D*)inFile->Get (Form ("JetPtSpectrum_cent%i", iCent));
   JetPtSpectrum[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

   JetCounts[iCent] = (TH2D*)inFile->Get (Form ("JetCounts_cent%i", iCent)); // divide by centrality bin width
   JetCounts[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

   TrackPtSpectrum[iCent] = (TH2D*)inFile->Get (Form ("TrackPtSpectrum_cent%i", iCent));
   TrackPtSpectrum[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

   TrackCounts[iCent] = (TH2D*)inFile->Get (Form ("TrackCounts_cent%i", iCent));
   TrackCounts[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

   TrackSelectionHist[iCent] = (TH1D*)inFile->Get (Form ("TrackSelectionHist_cent%i", iCent));
   JetSelectionHist[iCent] = (TH1D*)inFile->Get (Form ("JetSelectionHist_cent%i", iCent));

   DijetDeltaEtaDeltaPhi[iCent] = (TH2D*)inFile->Get (Form ("DijetDeltaEtaDeltaPhi_cent%i", iCent));
   DijetDeltaEtaDeltaPhi[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

   DijetEta1Eta2[iCent] = (TH2D*)inFile->Get (Form ("DijetEta1Eta2_cent%i", iCent));
   DijetEta1Eta2[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));

   DijetAj[iCent] = (TH1D*)inFile->Get (Form ("DijetAj_cent%i", iCent));
   DijetAj[iCent]->Scale (1. / (float)(centBins[iCent+1]-centBins[iCent]));
   if (DijetAj[iCent]->GetMaximum () > ajmax)
    ajmax = DijetAj[iCent]->GetMaximum ();
   if (DijetAj[iCent]->GetMinimum () > ajmin)
    ajmin = DijetAj[iCent]->GetMinimum ();
  }


  const Color_t colors[9] = {kBlack, kRed, kBlue, kMagenta, 8, kCyan+1, kOrange, kViolet, kGray};

  //TCanvas* c = new TCanvas ("c", "", 900, 800);
  //gStyle->SetPalette (kRainBow);
  ////gPad->SetLogz();
  //JetTrackDeltaEtaDeltaPhi[0]->RebinX (2);
  //JetTrackDeltaEtaDeltaPhi[0]->RebinY (2);
  //JetTrackDeltaEtaDeltaPhi[0]->Draw ("lego2");
  //c->SaveAs (Form ("%s/JetTrackEtaPhiCorr.pdf", plotPath.Data ()));

  TCanvas* DijetDeltaEtaDeltaPhiCanvas = new TCanvas ("DijetDeltaEtaDeltaPhiCanvas", "", 800, 600);
  FormatTH2Canvas (DijetDeltaEtaDeltaPhiCanvas);

  TCanvas* DijetDeltaPhiCanvas = new TCanvas ("DijetDeltaPhiCanvas", "", 800, 600);

  TCanvas* DijetEta1Eta2Canvas = new TCanvas ("DijetEta1Eta2Canvas", "", 800, 600);
  FormatTH2Canvas (DijetEta1Eta2Canvas);

  TCanvas* DijetAjCanvas = new TCanvas ("DijetAjCanvas", "", 800, 600);
  for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
   DijetDeltaEtaDeltaPhiCanvas->cd ();
   gPad->SetLogz ();

   DijetDeltaEtaDeltaPhi[iCent]->GetYaxis ()->SetTitleOffset (1);
   DijetDeltaEtaDeltaPhi[iCent]->GetZaxis ()->SetTitle ("Relative counts");
   DijetDeltaEtaDeltaPhi[iCent]->Draw ("colz");

   myText (0.7, 0.88, kBlack, Form ("%g-%g%%", centBins[iCent], centBins[iCent+1]), 0.045);
   DijetDeltaEtaDeltaPhiCanvas->SaveAs (Form ("%s/DijetDeltaEtaDeltaPhi_cent%i.pdf", plotPath.Data (), iCent));

   DijetDeltaPhiCanvas->cd ();
   gPad->SetLogy ();
   TH1D* th = DijetDeltaEtaDeltaPhi[iCent]->ProjectionY ();
   th->GetYaxis ()->SetRangeUser (1e6, 1e9);
   th->SetLineColor (colors[iCent]);
   th->SetMarkerColor (colors[iCent]);
   th->GetYaxis ()->SetTitle ("Relative counts");
   th->GetYaxis ()->SetTitleOffset (1);
   th->GetXaxis ()->SetTitleOffset (1);
  
   if (iCent == numCentBins-1) {
    th->DrawCopy ("e1");
    TLine* line1 = new TLine (3*pi/4, 1e6, 3*pi/4, 1e9);
    line1->SetLineStyle (7);
    line1->Draw ();
    TLine* line2 = new TLine (5*pi/4, 1e6, 5*pi/4, 1e9);
    line2->SetLineStyle (7);
    line2->Draw ();
   }
   else
    th->DrawCopy ("same e1");
   if (th) { delete th; th = NULL; }
   myMarkerText (0.25, 0.9 - iCent*0.05, colors[iCent], kFullCircle, Form ("%g-%g%%", centBins[iCent], centBins[iCent+1]), 1.25, 0.045/gPad->GetHNDC ());

   DijetEta1Eta2Canvas->cd ();
   gPad->SetLogz ();
   DijetEta1Eta2[iCent]->GetXaxis ()->SetTitle ("#eta^{lead}_{det}");
   DijetEta1Eta2[iCent]->GetXaxis ()->SetTitleOffset (1);
   DijetEta1Eta2[iCent]->GetYaxis ()->SetTitle ("#eta^{sublead}_{det}");
   DijetEta1Eta2[iCent]->GetYaxis ()->SetTitleOffset (1);
   DijetEta1Eta2[iCent]->Draw ("colz");
   
   myText (0.7, 0.88, kBlack, Form ("%g-%g%%", centBins[iCent], centBins[iCent+1]), 0.045);
   DijetEta1Eta2Canvas->SaveAs (Form ("%s/DijetEta1Eta2corr_cent%i.pdf", plotPath.Data (), iCent));

   DijetAjCanvas->cd ();
   gPad->SetLogy ();
   DijetAj[iCent]->GetYaxis ()->SetRangeUser (0.5*ajmin, 2*ajmax);
   DijetAj[iCent]->SetMarkerColor (colors[iCent]);
   DijetAj[iCent]->SetLineColor (colors[iCent]);
   if (iCent == numCentBins-1)
    DijetAj[iCent]->Draw ("e1");
   else
    DijetAj[iCent]->Draw ("same e1");
  }
  DijetDeltaPhiCanvas->SaveAs (Form ("%s/DijetDeltaPhi.pdf", plotPath.Data ()));
  DijetAjCanvas->SaveAs (Form ("%s/DijetAjCorr.pdf", plotPath.Data ()));


  TCanvas* JetTrackPhiCanvas = new TCanvas ("JetTrackPhiCanvas", "", 1200, 600);
  TPad* leftPad = new TPad ("leftPad", "", 0, 0, 0.50, 1);
  TPad* rightPad = new TPad ("rightPad", "", 0.50, 0, 1, 1);
  leftPad->SetLeftMargin (0.20);
  leftPad->SetRightMargin (0.02);
  rightPad->SetLeftMargin (0.20);
  rightPad->SetRightMargin (0.02);
  JetTrackPhiCanvas->Draw ();
  leftPad->Draw ();
  rightPad->Draw ();

  for (short iCut = 0; iCut < 2; iCut++) {
   for (short iEta = 0; iEta < 2; iEta++) {
    const int eta_lo = (iEta == 0 ? 0 : 2);
    const int eta_hi = (iEta == 0 ? 2 : 6);

    if (iEta == 0)
     leftPad->cd ();
    else
     rightPad->cd ();
    gPad->SetLogy (false);

    TH1D* hArr [numCentBins] = {};

    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
     TH1D* thisHist;
     thisHist = JetTrackDeltaEtaDeltaPhi[iCent][iCut]->ProjectionY (Form ("JetTrackDeltaEtaDeltaPhi_py%i_cent%i", iEta, iCent), eta_lo*10 + 1, eta_hi*10);
     float scale = JetSelectionHist[iCent]->GetBinContent (11);

     thisHist->Scale (0.1);

     thisHist->Rebin (4);
     thisHist->Scale (0.25);
     //if (iEta >= 1) {
     // thisHist->Rebin (2); // rebin again by 2
     // thisHist->Scale (0.5);
     //}

     thisHist->Scale (1. / scale);

     //thisHist->GetYaxis ()->SetRangeUser (0.9*ymin, 1.1*ymax);
     const float min = thisHist->GetMinimum ();
     for (int ix = 1; ix <= thisHist->GetNbinsX (); ix++) {
      thisHist->SetBinContent (ix, thisHist->GetBinContent (ix) - min); 
     }
     //thisHist->Scale (1. / thisHist->Integral ());

     //float max = 0;
     //for (int ix = thisHist->GetNbinsX () / 2; ix < thisHist->GetNbinsX (); ix++) {
     // if (max < thisHist->GetBinContent (ix))
     //  max = thisHist->GetBinContent (ix);
     //}
     //if (max > 0)
     // thisHist->Scale (1. / max);

     //float scale = thisHist->Integral (thisHist->GetNbinsX () / 2, thisHist->GetNbinsX ());
     //if (scale > 0)
     // thisHist->Scale (1. / scale);

     //thisHist->Scale (1. / scale);

     hArr[iCent] = thisHist;
    }

    double max = 0;
    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
     if (hArr[iCent]->GetMaximum () > max)
      max = hArr[iCent]->GetMaximum ();
    }

    for (short iCent = numCentBins-1; iCent >= 0; iCent--) {
     TH1D* thisHist = hArr[iCent];

     thisHist->GetYaxis ()->SetRangeUser (0, 1.1* max);
     thisHist->GetYaxis ()->SetTitle ("1/N_{jet} dN_{trk}/d#Delta#phi,  ZYAM");
     thisHist->GetXaxis ()->SetTitle ("#phi_{det}^{Jet} - #phi_{det}^{Trk}");
     thisHist->GetXaxis ()->SetTitleOffset (1);
     thisHist->GetYaxis ()->SetTitleOffset (1.6);
     thisHist->SetLineColor (colors[iCent]);
     thisHist->SetMarkerColor (colors[iCent]);
     if (iCent == numCentBins-1)
      thisHist->DrawCopy ("e1 x0");
     else
      thisHist->DrawCopy ("same e1 x0");
     if (thisHist) { delete thisHist; hArr[iCent] = NULL; }

     if (iEta == 0) {
      if (iCut == 0) {
       myMarkerText (0.74, 0.72 - iCent*0.05, colors[iCent], kFullCircle, Form ("%g-%g%%", centBins[iCent], centBins[iCent+1]), 1.25, 0.045/gPad->GetHNDC ());
      }
      else {
       myMarkerText (0.74, 0.72 - iCent*0.05, colors[iCent], kFullCircle, Form ("%g-%g%%", centBins[iCent], centBins[iCent+1]), 1.25, 0.045/gPad->GetHNDC ());
      }
     }
    }

    if (iCut == 0) {
     if (iEta == 0) {
      myText (0.65, 0.89, kBlack, Form ("%g < #left|#Delta#eta^{Jet}_{Trk}#right| < %g", (float)eta_lo, (float)eta_hi), 0.045/gPad->GetHNDC ()); 
      myText (0.65, 0.81, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/gPad->GetHNDC ());
     }
     else {
      myText (0.25, 0.89, kBlack, Form ("%g < #left|#Delta#eta^{Jet}_{Trk}#right| < %g", (float)eta_lo, (float)eta_hi), 0.045/gPad->GetHNDC ()); 
      myText (0.25, 0.81, kBlack, "#it{p}_{T}^{lead} > 30 GeV", 0.045/gPad->GetHNDC ());
     }
    }
    else {
     if (iEta == 0) {
      myText (0.65, 0.89, kBlack, Form ("%g < #left|#Delta#eta^{Jet}_{Trk}#right| < %g", (float)eta_lo, (float)eta_hi), 0.045/gPad->GetHNDC ()); 
      myText (0.65, 0.81, kBlack, "#bf{#it{ATLAS}} Internal", 0.045/gPad->GetHNDC ());
     }
     else {
      myText (0.25, 0.89, kBlack, Form ("%g < #left|#Delta#eta^{Jet}_{Trk}#right| < %g", (float)eta_lo, (float)eta_hi), 0.045/gPad->GetHNDC ()); 
      myText (0.25, 0.81, kBlack, "#left|#Delta#eta_{sublead.}^{Trk}#right| > 2", 0.045/gPad->GetHNDC ());
      myText (0.25, 0.73, kBlack, "#it{p}_{T}^{lead} > 30 GeV", 0.045/gPad->GetHNDC ());
      myText (0.25, 0.65, kBlack, "#it{p}_{T}^{sublead} > 30 GeV", 0.045/gPad->GetHNDC ());
     }
    }
   }
   JetTrackPhiCanvas->SaveAs (Form ("%s/JetTrackPhiCorr_cut%i.pdf", plotPath.Data (), iCut));
  }


  TCanvas* JetPtCanvas = new TCanvas ("JetPtCanvas", "", 1800, 800);
  JetPtCanvas->Divide (4, 2, -1, -1);
  JetPtCanvas->cd ();
  JetPtCanvas->Draw ();
  for (short iEta = 0; iEta < numetabins; iEta++) {
   JetPtCanvas->cd (iEta+1);
   gPad->SetLogy (true);
   for (short iCent = 0; iCent < numCentBins; iCent++) {

    TH1D* thisHist = JetPtSpectrum[iCent]->ProjectionX (Form ("_py%i", iEta), iEta, iEta);
    //thisHist->Rebin (2);

    thisHist->GetXaxis ()->SetTitle ("#it{p}_{T}^{jet} [GeV]");
    thisHist->GetYaxis ()->SetTitle ("d^{2}N / d#it{p}_{T} d(Cent) [GeV^{-1}]");
    thisHist->GetXaxis ()->SetTitleOffset (1);
    thisHist->GetYaxis ()->SetTitleOffset (1);
    thisHist->SetLineColor (colors[iCent]);
    thisHist->SetMarkerColor (colors[iCent]);
    if (iCent == 0)
     thisHist->DrawCopy ("e1 x0");
    else
     thisHist->DrawCopy ("same e1 x0");
    if (thisHist) delete thisHist;

    if (iEta == 1)
     myMarkerText (0.1, 0.45 - iCent*0.06, colors[iCent], kFullCircle, Form ("%g%%-%g%%", centBins[iCent], centBins[iCent+1]), 1.25, 0.03/gPad->GetHNDC ());
   }
   if (iEta == 0)
    myText (0.25, 0.94, kBlack, "#bf{#it{ATLAS}} Internal", 0.03/gPad->GetHNDC ());
   myText (0.6, 0.92, kBlack, Form ("%g < #eta^{Jet}_{det} < %g", etabins[iEta], etabins[iEta+1]), 0.03/gPad->GetHNDC ());
  }
  JetPtCanvas->SaveAs (Form ("%s/JetPtSpectra.pdf", plotPath.Data ()));


  TCanvas* TrackPtCanvas = new TCanvas ("TrackPtCanvas", "", 2250, 800);
  TrackPtCanvas->Divide (5, 2, -1, -1);
  TrackPtCanvas->cd ();
  TrackPtCanvas->Draw ();
  for (short iEta = 0; iEta < 10; iEta++) {
   TrackPtCanvas->cd (iEta+1);
   gPad->SetLogx (true);
   gPad->SetLogy (true);
   for (short iCent = 0; iCent < numCentBins; iCent++) {

    TH1D* thisHist = TrackPtSpectrum[iCent]->ProjectionX (Form ("_py%i", iEta), iEta+1, iEta+1);
    thisHist->GetYaxis ()->SetTitle ("N^{cent}_{evt} / N^{total}_{evt} d^{2}N / d#it{p}_{T} d#eta #left[GeV^{-1}#right]");
    thisHist->GetXaxis ()->SetTitleOffset (1);
    thisHist->GetYaxis ()->SetTitleOffset (1);
    thisHist->SetLineColor (colors[iCent]);
    thisHist->SetMarkerColor (colors[iCent]);
    if (iCent == 0)
     thisHist->DrawCopy ("e1 x0");
    else
     thisHist->DrawCopy ("same e1 x0");
    if (thisHist) delete thisHist;

    if (iEta == 1)
     myMarkerText (0.1, 0.45 - iCent*0.06, colors[iCent], kFullCircle, Form ("%g%%-%g%%", centBins[iCent], centBins[iCent+1]), 1.25, 0.03/gPad->GetHNDC ());
   }
   if (iEta == 0)
    myText (0.25, 0.94, kBlack, "#bf{#it{ATLAS}} Internal", 0.03/gPad->GetHNDC ());
   myText (0.6, 0.92, kBlack, Form ("%g < #eta^{Trk}_{det} < %g", -2.5 + 0.5*iEta, -2 + 0.5*iEta), 0.03/gPad->GetHNDC ());
  }
  TrackPtCanvas->SaveAs (Form ("%s/TrackPtSpectra.pdf", plotPath.Data ()));

}

}
